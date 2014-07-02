function assignBregmaCoordinates(experimentPath)

    %% Load JSON data
    jsonPath = fullfile(experimentPath, 'data.json');
    data = loadjson(jsonPath);

    for i = 1:length(data)
        %% Load each downsampled image
        currImData = data{i};
        im = imread(currImData.downsampledImagePath);

        %% Get Bregma Coordinate - disqualify images?
        currImData.bregmaCoord = getBregmaCoordinate(im);

        %% Update json file
        data{i} = currImData;
        savejson('', data, jsonPath);
    end


function coord = getBregmaCoordinate(im)
    %% Display image
    figure, imshow(im, [])
    title('Image Values Adjusted for Easier Visibility')

    %% Get Coordinate
    coord = queryUserForCoordinate();

function coord = queryUserForCoordinate()
    response = inputdlg('Please Enter Bregma Coord (mm)');
    coord = str2double(response);
    if isnan(coord)
        warndlg('Please Input a Number!', 'Unacceptable Input')
        coord = queryUserForCoordinate();
    end