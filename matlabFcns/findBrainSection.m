function mask = findBrainSection(im)
    %% mask = findBrainSection(im)
    % 
    % Returns a binary mask of the brain section in im
    im = mat2gray(im);

    %% Downsample im
    disp('Downsampling and low pass filtering image.....')
    scaleFactor = .02;
    sIm = imresize(im, scaleFactor);

    %% Segment Brain section
    disp('Binarizing image.....')
    thresh = graythresh(sIm(sIm ~= 0));
    mask = im2bw(sIm, thresh);
    %debugShow(mask, 'raw mask')

    %% Close
    disp('Cleaning binary image.....')
    mask = imclose(mask, strel('disk', 6));
    %debugShow(mask, 'closed')

    %% Clean up the mask
    mask = imfill(mask, 'holes');
    %debugShow(mask, 'filled')

    %% Erode
    mask = imerode(mask, strel('disk', 20));
    %debugShow(mask, 'after erosion')

    %% Get largest object
    disp('Choosing brain segment.....')
    mask = getLargest(mask);
    %debugShow(mask, 'largest');

    %% Dialate
    mask = imdilate(mask, strel('disk', 20));
    %debugShow(mask, 'dilated')

    %% Resize
    disp('Resizing binary image.....')
    mask = imresize(mask, size(im));
    %debugShow(mask, 'final')
    

function mask = getLargest(mask)
    %% label regions
    labMask = bwlabel(mask);

    %% Find largest region
    p = regionprops(labMask, 'Area');
    areas = [p(:).Area];
    largest = find(areas == max(areas));

    %% Set all other regions to zero
    if ~isempty(largest)
        mask(labMask ~= largest) = 0;
    end

function im = hackyBackgroundSubtract(im)
    [nr, nc] = size(im);
    fineSize = round(sqrt((nr*nc)/200)/(2*pi)) % radius
    courseSize = round(sqrt((nr*nc)/6)) % square edge
    fine = imfilter(im, fspecial('disk', fineSize));
    course = imfilter(im, fspecial('average', courseSize));
    im = fine - course;

function debugShow(im, titleStr)
    figure, imshow(im, [])
    if exist('titleStr', 'var')
        title(titleStr)
    end