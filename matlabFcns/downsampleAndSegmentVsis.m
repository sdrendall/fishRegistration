function downsampleAndSegmentVsis(jsonPath, varargin)
    %% downsampleAndSegmentVsis(jsonPath, options)
    %

    % TODO, help text!!!
    
    args = parseVarargin(nargin, varargin);

    %% Load image data from the JSON file
    imageData = loadjson(jsonPath);

    for i = 1:length(imageData)
        %% Load data
        currentImData = imageData{i};

        %% Find the brain section in the segmentation image
        segIm = loadPlaneFromVsi(currentImData.vsiPath, args.segmentationPlane);
        segIm = downsampleToAtlasScale(segIm);
        segIm = mat2gray(segIm);
        section = findBrainSection_AC(segIm);

        %% Load and Process the registration image
        % Downsample, segment brain section, and flip, flop
        regIm = loadPlaneFromVsi(currentImData.vsiPath, args.registrationPlane);
        regIm = downsampleToAtlasScale(regIm);
        regIm(~section) = 0;
        if args.flip && args.flop
            regIm = flipflop(regIm);
        elseif args.flip
            regIm = flip(regIm);
        elseif args.flop
            regIm = flop(regIm);
        end

        %% Pad the resulting image with zeros so it matches the size of the ABA
        regIm = padToAtlasSize(regIm);

        %% Save as 8-bit mhd
        %regIm = convertToUint8(regIm);
        %outputPath = generateoutputPath(jsonPath, s.vsiPath);
        formatAndSaveMhd(regIm, outputPath);
        imwrite(regIm, outputPath)

        %% Update JSON file
        disp('Updating JSON file.....')

        % Paths stored in metadata.json should be relative to the experimentPath
        relativePath = generateRelativePath(outputPath);
        s.downsampledImagePath = relativePath;
        imageData{i} = s;
        savejson('', imageData, jsonPath);
    end

function args = parseVarargin(argc, argv)
    % Set defaults
    args = struct( ...
        'flip', false, ...
        'flop', false, ...
        'registrationPlane', 1, ...
        'segmentationPlane', 2)

    % Iterate though the input arguments until none are left
    i = 1;
    while i <= argc - 1;
        switch lower(argv{i})
        case 'flip'
            args.flip = true;
            i = i + 1;
        case 'flop'
            args.flop = true;
            i = i + 1;
        case 'registrationplane'
            args.registrationPlane = str2num(argv{i + 1});
            i = i + 2;
        case 'segmentationplane'
            args.segmentationPlane = str2num(argv{i + 1});
            i = i + 2;
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

function im = flipflop(im)
    % flips the image vertically and horizontally
    im = im(end:-1:1, end:-1:1);

function im = flip(im)
    % flips the image vertically
    im = im(end:-1:1, :);

function im = flop(im)
    % flips the image horizontally
    im = im(:, end:-1:1);

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

function formatAndSaveMhd(data, outputPath)
    disp('Generating Image Metadata.....')
    img = generateMhdImage(data);
    disp(['Saving Downsampled Image to ', outputPath])
    write_mhd(outputPath, img,'ElementType', 'uint8', 'NDims', 2);

function img = generateMhdImage(data)
    s = size(data);
    origin = [0 0];
    spacing = [25 25];
    orientation =[1 0; 0 1];

    img = ImageType(s, origin, spacing, orientation);
    img.data = data;