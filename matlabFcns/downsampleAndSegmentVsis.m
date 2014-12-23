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
    %   2. Image is normalized so that it fills the full range of the output bitdepth
    %   3. Image is downsized such that each pixel represents 25um x 25um of physical space
    %   4. Image is padded with zeros to be at least the same size as coronal images in the Allen Brain Atals (456 x 320)
    %   5. Zeros are replaced with poisson noise
    %   6. Image is saved as an 8-bit png file in the directory containing the JSON file at jsonPath
    %
    %   The JSON file is updated after each iteration to include the path to the downsampled image

    %% Load image data from the JSON file
    imageData = loadjson(jsonPath);

    for i = 1:length(imageData)
        %% Load data
        s = imageData{i};
        im = loadDapiFromVsi(s.vsiPath);
        showMinMax(im)

        %% Calculate background noise statistics
        [avg, sDev] = getBackgroundStats(im);  

        %% Downsample, flip, flop, and pad
        im = downsampleToAtlasScale(im);
        showMinMax(im)
        im = flipflop(im);
        showMinMax(im)
        im = padToAtlasSize(im);
        showMinMax(im)
        im = convertZerosToNoise(im, avg, sDev);
        showMinMax(im)

        %% Save as 8-bit png
        im = convertToUint8(im);
        showMinMax(im)
        pngPath = generatePngPath(jsonPath, s.vsiPath);
        imwrite(im, pngPath);

        %% Update JSON file
        disp('Updating JSON file.....')

        % Paths stored in metadata.json should be relative to the experimentPath
        relativePath = generateRelativePath(pngPath);
        s.downsampledImagePath = relativePath;
        imageData{i} = s;
        savejson('', imageData, jsonPath);
    end


function im = loadDapiFromVsi(vsiPath)
    disp(['Loading ', vsiPath, '.....'])
    r = bfGetReader(vsiPath);
    im = bfGetPlane(r, 1);

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

function pngPath = generatePngPath(jsonPath, vsiPath)
    disp('Generating output path.....')
    dataDir = fileparts(jsonPath);
    [~, baseVsiName] = fileparts(vsiPath);
    pngPath = fullfile(dataDir, [baseVsiName, '_downsampled.png']);

function relPath = generateRelativePath(pngPath)
    pathCoord = strfind(pngPath, '.registrationData');
    relPath = pngPath(pathCoord:end);

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