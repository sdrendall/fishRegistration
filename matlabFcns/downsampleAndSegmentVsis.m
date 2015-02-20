function downsampleAndSegmentVsis(vsiPath, outputPath, varargin)
    %% downsampleAndSegmentVsis(vsiPath, outputPath, [options])
    %
    % Loads the image at vsiPath, locates the brain section using active contour segmentation, 
    %  and downsamples the image such that the pysical scale of each pixel in the resulting image
    %  is the same as the pixel scale of images in the allen brain atlas.  This is intended as a 
    %  preprocessing step to format input images for registration to the allen brain atlas.
    %
    % The downsampled image is saved at outputPath
    %
    % Additional options can be specified following positional arguments to obtain more desireable results:
    %  'flip' - the image will be flipped vertically once it is downsampled
    %  'flop' - the image will be flipped horizontally once it is downsampled
    %  'registrationPlane' planeNo - planeNo is an integer specifying the plane of the vsi image to be downsampled and saved.
    %    This should be the plane number of the image that is to be registered to the allen brain atlas.  Ideally, this image
    %    should depict some sort of nuclear stain, such as dapi, hoechst, or nissl
    %  'segmentationPlane' planeNo - The plane (channel) to be used to segment the brain section in the image from the background.
    %    Ideally this plane contains an image with high tissue autofluorescence (green), and relatively sparse staining'
    
    args = parseVarargin(nargin-2, varargin);

    %% Find the brain section in the segmentation image
    segIm = loadPlaneFromVsi(vsiPath, args.segmentationPlane);
    segIm = downsampleToAtlasScale(segIm);
    segIm = mat2gray(segIm);
    section = findBrainSection_AC(segIm);

    %% Load and Process the registration image
    % Downsample, segment brain section, and flip, flop
    regIm = loadPlaneFromVsi(vsiPath, args.registrationPlane);
    regIm = downsampleToAtlasScale(regIm);
    regIm(~section) = 0;
    regIm = permuteImage(regIm, args);

    %% Pad the resulting image with zeros so it matches the size of the ABA
    regIm = padToAtlasSize(regIm);
    
    %% Convert to 8-bit then save 
    regIm = convertToUint8(regIm);
    saveImage(regIm, outputPath);
    
function args = parseVarargin(argc, argv)
    % Set defaults
    args = struct( ...
        'flip', false, ...
        'flop', false, ...
        'registrationPlane', 1, ...
        'segmentationPlane', 2);

    % Iterate though the input arguments until none are left
    i = 1;
    while i <= argc;
        switch lower(argv{i})
        case 'flip'
            args.flip = true;
            i = i + 1;
        case 'flop'
            args.flop = true;
            i = i + 1;
        case 'registrationplane'
            args.registrationPlane = argv{i + 1};
            i = i + 2;
        case 'segmentationplane'
            args.segmentationPlane = argv{i + 1};
            i = i + 2;
        otherwise
            i = i + 1;
        end
    end

function im = loadPlaneFromVsi(vsiPath, plane)
    disp(['Loading plane ', num2str(plane), ' from ', vsiPath, '.....'])
    r = bfGetReader(vsiPath);
    im = bfGetPlane(r, plane);

function im = downsampleToAtlasScale(im)
    disp('Downsizing image.....')
    vsiPixelSize = .64497; % in um
    atlasPixelSize = 25; % in um
    atlasScale = vsiPixelSize/atlasPixelSize;
    im = imresize(im, atlasScale);

function im = permuteImage(im, args)
    if args.flip && args.flop
        im = flipflop(im);
    elseif args.flip
        im = flip(im);
    elseif args.flop
        im = flip(im, 2);
    end

function im = flipflop(im)
    % flips the image vertically and horizontally
    im = im(end:-1:1, end:-1:1);

function im = padToAtlasSize(im)
    disp('Padding to atlas size.....')
    atlasSize = [320, 456];
    marginSize = ceil((atlasSize - size(im))./2);
    marginSize(marginSize < 0) = 0; % Remove 0s if im is larger than atlasSize in either dimension
    im = padarray(im, marginSize);

function [m, s] = getBackgroundStats(im)
    disp('Computing background noise statistics.....')

    % Use the corners
    q{1} = im(1:500, 1:500);
    q{2} = im(1:500, end-500:end);
    q{3} = im(end-500:end, 1:500);
    q{4} = im(end-500:end, end-500:end);

    m = zeros([1 4]);
    for i = 1:4
        m(i) = mean(q{i}(:));
    end

    m = median(m);
    s = sqrt(m)*2;  % Photon arrival is poisson, each photon is worth 2 in these images


function im = convertToUint8(im)
    disp('Converting to uint8.....')
    im = mat2gray(im); % Normalize, set all values 0 - 1
    im = uint8(im*255); % Spread over 255, convert to uint8

function im = uint16ToUint8(im)
    disp('Converting uint16 to uint8.....')
    im = double(im)./(2^16 - 1); % Set all values 0 - 1
    im = uint8(im*255); % Spread over 255, convert to uint8

function showMinMax(im)
    disp(['min: ', num2str(min(im(:)))])
    disp(['max: ', num2str(max(im(:)))])

function saveImage(im, outputPath)
    % Calls the appropriate save function based on the file extension of the given output path
    [~, ~, ext] = fileparts(outputPath);
    if strcmpi(ext, '.mhd')
        formatAndSaveMhd(im, outputPath)
    else
        imwrite(im, outputPath)
    end

function formatAndSaveMhd(data, outputPath)
    % Generates an ImageType object from the given image data, then saves it to outputPath in the mhd format
    disp('Generating image metadata.....')
    img = generateMhdImage(data);
    disp(['Saving Downsampled Image to ', outputPath])
    write_mhd(outputPath, img,'ElementType', 'uint8', 'NDims', 2);

function img = generateMhdImage(data)
    % Generates an ImageType object from the given image data
    s = size(data);
    origin = [0 0];
    spacing = [25 25];
    orientation =[1 0; 0 1];

    img = ImageType(s, origin, spacing, orientation);
    img.data = data;

%% Unused Funtions %%
% These are functions that were used in previous versions of this script

function im = convertZerosToNoise(im, m, s)
    z = im(im == 0);
    z = m + s .* randn(size(z));
    im(im == 0) = z;

function outputPath = generateoutputPath(jsonPath, vsiPath)
    disp('Generating output path.....')
    dataDir = fileparts(jsonPath);
    [~, baseVsiName] = fileparts(vsiPath);
    outputPath = fullfile(dataDir, [baseVsiName, '_downsampled.png']);

function relPath = generateRelativePath(outputPath)
    pathCoord = strfind(outputPath, '.registrationData');
    relPath = outputPath(pathCoord:end);