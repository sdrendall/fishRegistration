function mask = findBrainSection(im)
    %% mask = findBrainSection(im)
    % 
    % Returns a binary mask of the brain section in im

    %% Downsample im
    scaleFactor = .005;
    sIm = imresize(im, scaleFactor);

    %% Segment Brain section
    thresh = graythresh(sIm(:));
    mask = im2bw(sIm, thresh);
    %debugShow(mask, 'raw mask')

    %% Clean up the mask
    mask = imfill(mask, 'holes');
    %debugShow(mask, 'filled')

    mask = imerode(mask, strel('disk', 5));
    %debugShow(mask, 'after erosion')

    %% Get largest object
    mask = getLargest(mask);
    %debugShow(mask, 'largest');

    %% Dialate
    mask = imdilate(mask, strel('disk', 6));
    %debugShow(mask, 'dilated')


    %% Resize
    mask = imresize(mask, size(im));

    %% Hull
    mask = bwconvhull(mask, 'Union');
    %debugShow(mask, 'final')


    %% Close
    %mask = imdilate(mask, strel('disk', 300));
    %%debugShow(mask, 'final')


function mask = getLargest(mask)
    %% label regions
    labMask = bwlabel(mask);

    %% Find largest region
    p = regionprops(labMask, 'Area');
    areas = [p(:).Area];
    largest = find(areas == max(areas));

    %% Set all other regions to zero
    mask(labMask ~= largest) = 0;


function debugShow(im, titleStr)
    figure, imshow(im, [])
    if exist('titleStr', 'var')
        title(titleStr)
    end