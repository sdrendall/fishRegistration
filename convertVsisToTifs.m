function convertVsisToTifs(startingImagePath, startingSavePath)

disp('Starting......')

if ~exist('startingImagePath', 'var')
    startingImagePath = '/hms/scratch1/sr235/ccValidation_03-18-14';
end

if ~exist('startingSavePath', 'var')
    startingSavePath = '/hms/scratch1/sr235/ccValidation_03-18-14_hemisphere_tifs';
end

% Find VSIs 
disp('Searching for Vsis.......')
vsiFiles = findVsis(startingImagePath, startingSavePath);

disp('Calculating Output Image Dimensions........')
maxDims = findLargestImageSize(vsiFiles);


for i = 1:length(vsiFiles)
    % check for existing file
    [~, nameNoExt] = fileparts(vsiFiles(i).name);
    writeName = [nameNoExt, '.tif'];
    destPath = fullfile(vsiFiles(i).dataPath, writeName);
    disp(['Checking for ', writeName])
    if ~exist(destPath, 'file')
        disp(['Converting ', vsiFiles(i).name, 'to tiff ', writeName])
        tic
        convertToTif(vsiFiles(i), maxDims)
        toc
    else
        disp(['Found ', writeName, '.  Moving to next image'])
    end
end

disp('Process Complete!')

function dims = findLargestImageSize(vsiFiles)
    dims = [0,0];
    for i = 1:length(vsiFiles)
        reader = bfGetReader(vsiFiles(i).path);
        h = reader.getSizeY;
        w = reader.getSizeX;
        if h > dims(1)
            dims(1) = h;
        end
        if w > dims(2)
            dims(2) = w;
        end
    end


function convertToTif(vsiFile, targetDimensions)
    % Load VSI
    vsi = bfopen(vsiFile.path);

    % Construct RGB image
    marginSize = targetDimensions - size(vsi{1,1}{2,1});
    rgb = zeros([targetDimensions, 3], 'uint8');
    rgb(:,:,2) = toUint8(padarray(mat2gray(vsi{1,1}{2,1}), marginSize, 'post'));
    rgb(:,:,3) = toUint8(padarray(mat2gray(vsi{1,1}{1,1}), marginSize, 'post'));

    % Save to tif
    [~, nameNoExt] = fileparts(vsiFile.name);
    writeName = [nameNoExt, '.tif'];
    bfsave(rgb, fullfile(vsiFile.dataPath, writeName))

    clear rgb vsi


function im = toUint8(im)
    % Converts im to uint8, normalizes to 0-255
    im = im - min(im(:));
    im = im ./ max(im(:));
    im = uint8(im .* 255);


function vsiFiles = findVsis(locationPath, savePath)
    vsiFiles = [];

    % Search locationPath
    locationPathContents = dirNoDot(locationPath);
    if length(locationPathContents) == 0
        return
    end
    
    % Create vsiData objects for files with .vsi extensions
    for i = 1:length(locationPathContents)
        [~, ~, ext] = fileparts(locationPathContents(i).name);
        if any(strcmpi(ext, '.vsi'))
            if ~exist('vsiFiles', 'var') || isempty(vsiFiles)
                vsiFiles = vsiData(locationPathContents(i), locationPath, savePath);
            else
                vsiFiles(end + 1) = vsiData(locationPathContents(i), locationPath, savePath);
            end
        end
    end


    % Recursively search directories
    directories = locationPathContents([locationPathContents(:).isdir]);
    for i = 1:length(directories)
        % Construct paths to vsi files and to data locations
        vsiPath = fullfile(locationPath, directories(i).name);
        dumpPath = fullfile(savePath, directories(i).name);
        if isempty(vsiFiles)
            vsiFiles = findVsis(vsiPath, dumpPath);
        else
            incomingVsiFiles = findVsis(vsiPath, dumpPath);
            if ~isempty(incomingVsiFiles)
                vsiFiles = [vsiFiles, incomingVsiFiles];
            end
        end
    end
