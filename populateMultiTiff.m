function populateMultiTiff(imagePaths, baseOutput, channels)
%% reads images from imagePaths and writes them to a single tiff file at baseOutput
% For now, it just writes the green and blue channels
%
% Intended for use with tiff files

if ~exist('channels', 'var')
    channels = [2 3];
end

outputPath = cell(size(channels));

% Calculate Target Dimensions
[maxX, maxY] = getLargestDimensions(imagePaths);

for i = 1:length(channels)
    % read in the first image
    im = imread(imagePaths{1}, channels(i));
    im = removeZeros(padToMax(im, maxX, maxY));
    outputPath{i} = formatOutputPath(baseOutput, channels(i));
    imwrite(im, outputPath{i});
end

if length(imagePaths) < 2
    return
end

parfor chan = 1:length(channels)
    currChannel = (channels(chan));
    for i = 2:length(imagePaths)
        currIm = imread(imagePaths{i}, currChannel);
        currIm = removeZeros(padToMax(currIm, maxX, maxY));
        imwrite(currIm, outputPath{chan}, 'WriteMode', 'append')
    end
end



function im = padToMax(im, maxX, maxY)
    marginSize = [maxY, maxX] - size(im);
    im = padarray(im, marginSize, 'post');

function im = removeZeros(im)
        nonZero = im(:) ~= 0;
        im(~nonZero) = min(im(nonZero));
        im = mat2gray(im)

function [maxX, maxY] = getLargestDimensions(paths)
    for i = 1:length(paths)
        info(i) = imfinfo(paths{i})
    end
    maxX = max([info(:).Width]);
    maxY = max([info(:).Height]);

function outputPath = formatOutputPath(basePath, channel)
    [baseDir, baseName] = fileparts(basePath);
    switch channel
        case 1
            c = 'red';
        case 2
            c = 'green';
        case 3
            c = 'blue';
        otherwise
            c = '';
    end
    outputPath = fullfile(baseDir, [baseName, '_', c, '.tif']);