function getCropCoordinates(imPaths, outputPath)
%% Get's user input crop coordinates and saves them to a .mat file

    for i = 1:length(imPaths)
        thumb{i} = imread(imPaths{i});
        thumb{i} = mat2gray(thumb{i}(:,:, 3));
    end

    figure
    for i = 1:length(thumb)
        imshow(thumb{i})
        [x, y] = ginput
        [ty, tx]= size(thumb{i})
        x = x ./ tx;
        y = y ./ ty;
        crop(i).x = assertBoundries(x);
        crop(i).y = assertBoundries(y);
        crop(i).x
        crop(i).y
    end

    save(outputPath, 'crop')

function c = assertBoundries(c)
    c(c > 1) = 1;
    c(c <= 0) = [];