function outputPaths = cropRemoveZerosAndResize(imagePaths, cropCoordinatesPath)
%% Crops images at imagePaths, removes zeros, then resize

    % Load 'crop'
    disp('Loading crop coordinates.....')
    load(cropCoordinatesPath)

    disp('Loading readers.....')
    for i = 1:length(imagePaths)
        disp([num2str(i), '.....'])
        readers(i) = bfGetReader(imagePaths{i});
    end

    disp('calculating output dimensions.....')
    [maxX, maxY] = calculateMaxImageSize(readers, crop);
    clear readers

    for i = 1:length(imagePaths)
        disp(['processing image ', num2str(i), '.....'])
        tic
        r = bfGetReader(imagePaths{i});
        im = zeros([maxY, maxX, 3]);
        disp('loading channel 2.....')
        im(:,:,2) = cropAndProcess(bfGetPlane(r, 2), crop(i), maxX, maxY);
        disp('loading channel 3.....')
        im(:,:,3) = cropAndProcess(bfGetPlane(r, 1), crop(i), maxX, maxY);
        [baseDir, baseName] = fileparts(imagePaths{i});
        outputPaths{i} = fullfile(baseDir, [baseName, '_preProc.tif']);
        imwrite(im, outputPaths{i})
        disp(['processed ', baseName, ' in ', num2str(toc), ' seconds'])
    end


function im = cropAndProcess(im, crop, maxX, maxY)
    disp('cropping.....')
    im = cropROI(im, crop);
    disp('clearing zeros.....')
    im = removeZeros(im);
    disp('padding.....')
    im = padToMax(im, maxX, maxY);

function im = removeZeros(im)
        nonZero = im(:) ~= 0;
        im(~nonZero) = min(im(nonZero));
        im = mat2gray(im);

function im = padToMax(im, maxX, maxY)
    marginSize = [maxY, maxX] - size(im);
    im = padarray(im, marginSize, 'post', 'replicate');

function [maxX, maxY] = calculateMaxImageSize(r, crop_all)
    for i = 1:length(r)
        [sizeX(i), sizeY(i)] = getCropSize(r(i), crop_all(i));
    end
    maxX = max(sizeX);
    maxY = max(sizeY);

function im = cropROI(im, crop)
    % Find start points
    [nr, nc] = size(im);
    minX = round(min(crop.x) * nc);
    minY = round(min(crop.y) * nr);
    maxX = round(max(crop.x) * nc);
    maxY = round(max(crop.y) * nr);
    im = im(minY:maxY, minX:maxX);

function [sizeX, sizeY] = getCropSize(r, crop)
    minX = round(min(crop.x) * r.getSizeX);
    minY = round(min(crop.y) * r.getSizeY);
    maxX = round(max(crop.x) * r.getSizeX);
    maxY = round(max(crop.y) * r.getSizeY);
    sizeX = maxX - minX + 1;
    sizeY = maxY - minY + 1;