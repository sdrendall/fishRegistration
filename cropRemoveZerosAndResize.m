function outputPaths = cropRemoveZerosAndResize(imagePaths, cropCoordinatesPath)
%% Crops images at imagePaths, removes zeros, then resize

    % Load 'crop'
    disp('Loading crop coordinates.....')
    load(cropCoordinatesPath)

    for i = 1:length(imagePaths)
        readers(i) = bfGetReader(imagePaths{i});
    end

    disp('calculating output dimensions.....')
    [maxX, maxY] = calculateMaxImageSize(readers, crop)
    clear readers

    for i = 1:length(imagePaths)
        disp(['processing image ', num2str(i), '.....'])
        tic
        r = bfGetReader(imagePaths{i});
        im = zeros([maxY, maxX, 3]);
        im(:,:,2) = cropAndProcess(bfGetPlane(r, 2), r, crop, maxX, maxY);
        im(:,:,3) = cropAndProcess(bfGetPlane(r, 1), r, crop, maxX, maxY);
        [baseDir, baseName] = fileparts(imagePaths{i});
        outputPaths{i} = fullfile(baseDir, [baseName, '_preProc.tif']);
        imwrite(im, outputPaths{i})
        disp(['processed ', baseName, ' in ', num2str(toc), ' seconds'])
    end


function im = cropAndProcess(im, r, crop, maxX, maxY)
    disp('cropping.....')
    im = cropROI(im, r, crop);
    disp('padding.....')
    im = padToMax(im, maxX, maxY);
    disp('clearing zeros.....')
    im = removeZeros(im);

function im = removeZeros(im)
        nonZero = im(:) ~= 0;
        im(~nonZero) = min(im(nonZero));
        im = mat2gray(im)

function im = padToMax(im, maxX, maxY)
    marginSize = [maxY, maxX] - size(im);
    im = padarray(im, marginSize, 'post');

function [maxX, maxY] = calculateMaxImageSize(r, crop)
    for i = 1:length(r)
        [sizeX(i), sizeY(i)] = getCropSize(r(i), crop(i));
    end
    maxX = max(sizeX);
    maxY = max(sizeY);

function im = cropROI(im, r, crop)
    [sizeX, sizeY] = getCropSize(r, crop);
    % Find start points
    [nr, nc] = size(im);
    minX = round(min(crop.x) * nc);
    minY = round(min(crop.y) * nr);
    im = im(minY:minY+sizeY-1, minX:minX+sizeX-1);

function [sizeX, sizeY] = getCropSize(r, crop)
    xScaleCoEff = (max(crop.x) - min(crop.x));
    yScaleCoEff = (max(crop.y) - min(crop.y));
    sizeX = round(r.getSizeX * xScaleCoEff);
    sizeY = round(r.getSizeY * yScaleCoEff);
