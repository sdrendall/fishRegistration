function getCropCoordinates(imPaths, outputPath)
%% Get's user input crop coordinates and saves them to a .mat file

    for i = 1:length(imPaths)
        thumb{i} = imread(imPaths{i});
    end

    figure
    for i = 1:length(thumb)
        imshow(thumb{i})
        [x, y] = ginput;
        [ty, tx]= size(thumb{i});
        x = x ./ tx;
        y = y ./ ty;
        crop(i).x = assertBoundries(x)
        crop(i).y = assertBoundries(y)
    end

    save(outputPath, 'crop')