function generateTiffThumbnails(imagePath)

    % Recurse for cells
    if iscell(imagePath)
        for i = 1:length(imagePath)
            generateTiffThumbnails(imagePath{i});
        end
    end

    % load and downsample images
    r = bfGetReader(imagePath);
    thumb = imresize(bfGetPlane(r, 1), .0625);

    % save
    [baseDir, baseName] = fileparts(imagePath);
    outputPath = fullfile(baseDir, [baseName, '_downsampled.tif']);
    imwrite(thumb, outputPath)