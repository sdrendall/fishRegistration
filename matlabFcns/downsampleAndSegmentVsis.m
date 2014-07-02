function downsampleAndSegmentVsis(jsonPath)
    %% downsampleAndSegmentVsis(jsonPath)
    %
    % Parses JSON file located at jsonPath
    %
    % File should consist of a cell array of structures, each with the field vsiPath
    %   which denotes the .vsi file to be preprocessed
    %
    % Preprocessing of each file specified includes:
    %   1. The Dapi (or nissl) channel is loaded from each vsi
    %   2. Brain Section is determined using otsu thresholding
    %   3. Non-brain pixels are set to 0
    %   4. Image is downsized such that each pixel represents 25um x 25um of physical space
    %   5. Image is padded with zeros to be at least the same size as coronal images in the Allen Brain Atals (456 x 320)
    %   6. Image is saved as an 8-bit png file in the directory containing the JSON file at jsonPath
    %
    %   The JSON file is updated after each iteration to include the path to the downsampled image

    %% Load image data from the JSON file
    imageData = loadjson(jsonPath);

    for i = 1:length(imageData)
        %% Load data
        s = imageData{i};
        im = loadDapiFromVsi(s.vsiPath);

        %% Segment Section
        disp('Detecting brain section.....')
        brain = findBrainSection(im);
        im(~brain) = 0;

        %% Downsample, flip, flop, and pad
        im = downsampleToAtlasScale(im);
        im = flipflop(im);
        im = padToAtlasSize(im);

        %% Save as 8-bit png
        pngPath = generatePngPath(jsonPath, s.vsiPath);
        im = convertToUint8(im);
        imwrite(im, pngPath);

        %% Update JSON file
        disp('Updating JSON file.....')
        s.downsampledImagePath = pngPath;
        imageData{i} = s;
        savejson('', imageData, jsonPath);
    end


function im = loadDapiFromVsi(vsiPath)
    disp(['Loading ', vsiPath, '.....'])
    r = bfGetReader(vsiPath);
    im = bfGetPlane(r, 1);
    im = double(im);

function im = downsampleToAtlasScale(im)
    disp('Downsizing image.....')
    vsiPixelSize = .64497; % in um
    atlasPixelSize = 25; % in um
    atlasScale = vsiPixelSize/atlasPixelSize;
    im = imresize(im, atlasScale);

function im = flipflop(im)
    im = im(end:-1:1, end:-1:1);

function im = padToAtlasSize(im)
    disp('Padding to atlas size.....')
    atlasSize = [320, 456];
    marginSize = ceil((atlasSize - size(im))./2);
    marginSize(marginSize < 0) = 0; % Remove 0s if im is larger than atlasSize in either dimension
    im = padarray(im, marginSize);

function pngPath = generatePngPath(jsonPath, vsiPath)
    disp('Generating output path.....')
    jsonDir = fileparts(jsonPath);
    [~, baseVsiName] = fileparts(vsiPath);
    pngPath = fullfile(jsonDir, [baseVsiName, '.png']);

function im = convertToUint8(im)
    disp('Converting to uint8.....')
    im = mat2gray(im); % Normalize, set all values 0 - 1
    im = uint8(im*255); % Spread over 255, convert to uint8