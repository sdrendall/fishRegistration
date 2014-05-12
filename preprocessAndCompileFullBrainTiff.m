function preprocessAndCompileFullBrainTiff(r)

    thumb = cell(size(r));
    % load thumbs
    disp('loading thumbnails....')
    for i = 1:length(r)
        disp(['  ', num2str(i), '....'])
        thumb{i} = imresize(bfGetPlane(r(i), 1), .0625);
    end

    % Get Crop coords
    disp('specify save location for crop coordinates')
    cropSave = uigetfile();

    disp('please specify crop coordinates')
    figure
    for i = 1:length(thumb)
        imshow(thumb{i})
        [x, y] = ginput(4);
        [ty, tx] = size(thumb{i});
        % Normalize x,y to thumbnail size
        crop(i).x = x ./ tx;
        crop(i).y = y ./ ty;
    end

    disp('saving coordinates.....')
    save(cropSave, 'crop')

    % Create individual, preprocessed images
    disp('specify output for cropped images.....')
    outputDir = uigetdir();
    if ~exist('outputDir', 'dir')
        mkdir(outputDir)
    end

    disp('cropping and saving images.....')
    parfor i = 1:length(r)
        r(i).setSeries(0)
        xScaleCoEff = (max(crop(i).x) - min(crop(i).x));
        yScaleCoEff = (max(crop(i).y) - min(crop(i).y));
        sizeX = round(r(i).getSizeX * xScaleCoEff);
        sizeY = round(r(i).getSizeY * yScaleCoEff);
        im = zeros([sizeY, sizeX, 3]);
        im(:,:,2) = cropROI(bfGetPlane(r(i), 2), crop(i).x, crop(i).y);
        im(:,:,3) = cropROI(bfGetPlane(r(i), 1), crop(i).x, crop(i).y);
        outputPath = fullfile(outputPath, ['preProc_', num2str(i), '.tif']);
        imwrite(im, outputPath)
    end

    disp('compiling multitiff....')
    % Save to MultiTiff
    f = dir(fullfile(outputDir, '*.tif'));
    multiOut = fullfile(outputDir, 'fullBrain');
    populateMultiTiff({f(:).name}, multiOut);


function im = cropROI(im, x, y)
    % Determine Ultimate Image Size
    xScaleCoEff = (max(crop(i).x) - min(crop(i).x));
    yScaleCoEff = (max(crop(i).y) - min(crop(i).y));
    sizeX = round(r(i).getSizeX * xScaleCoEff);
    sizeY = round(r(i).getSizeY * yScaleCoEff);
    % Find start points
    [nr, nc] = size(im);
    minX = round(min(x) * nc);
    minY = round(min(y) * nr);
    im = im(minY:minY+sizeY-1, minX:minX+sizeX-1);