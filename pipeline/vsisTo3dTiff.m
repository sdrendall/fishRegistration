function vsisTo3dTiff(vsiPath, tiffPath, dicomPath)
    %% vsisTo3dTiff(vsiPath)
    %
    % Sequentially loads each .vsi file at vsiPath and stacks the contained
    % images to a 3D .tiff file.
    %
    % Each image is segmented to find the brain section, and all pixels not
    % in the identified section are set to 0.
    %
    % The x and y dimensions of the output image are all set to the largest
    % found within the specified files.
    dicomScale = .25;

    %% Find vsi files -- using bfReaders
    readers = findVsis(vsiPath);

    %% Get max X and Y coordinates
    [maxX, maxY] = getMaxXandY(readers);
    dicomX = round(maxX*dicomScale);
    dicomY = round(maxY*dicomScale);

    %% Compile 3D tiff
    % Initialize output container
    dicomImage = zeros([dicomY, dicomX, length(readers)]);
    tiffWriter = TiffWriter(tiffPath);
    for i = 1:length(readers)
        %% Load dapi channel
        im = bfGetPlane(readers(i), 1);
        im = convertToUint8(im);

        %% Find brain section and clear surrounding area
        section = findBrainSection(im);
        im(~section) = 0;

        %% Pad image to reach max dimensions
        im = padToMax(im, maxX, maxY);

        %% Append to 3d Image
        tiffWriter.append(im)
        dicomImage(:,:,i) = imresize(im, [dicomY, dicomX]);
    end

    %% Write dicom image
    dicomwrite(dicomImage, dicomPath)


function im = convertToUint8(im)
    im = mat2gray(im); % Normalize, set all values 0 - 1
    im = uint8(im*255); % Spread over 255, convert to uint8

function readers = findVsis(vsiPath)
    %% Returns bfReader objects for each vsi at vsiPath
    f = dir(vsiPath);
    vsiDir = fileparts(vsiPath); 
    for i = 1:length(f)
        readers(i) = bfGetReader(fullfile(vsiDir, f(i).name));
    end

function [maxX, maxY] = getMaxXandY(readers)
    %% Returns the maximum X and Y dimensions found in images referenced
    % by readers
    xVals = zeros(size(readers));
    yVals = xVals;
    for i = 1:length(readers)
        xVals(i) = readers(i).getSizeX;
        yVals(i) = readers(i).getSizeY;
    end
    maxX = max(xVals);
    maxY = max(yVals);

function im = padToMax(im, maxX, maxY)
    %% Pads im with zeros to the given maxX and maxY dimensiosn
    marginSize = [maxY, maxX] - size(im);
    im = padarray(im, marginSize, 'post');