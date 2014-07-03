function assignBregmaCoordinates(experimentPath)

    %% Load JSON data
    jsonPath = fullfile(experimentPath, '.metadata.json');
    data = loadjson(jsonPath);

    for i = 1:length(data)
        %% Load each downsampled image
        currImData = data{i};
        im = imread(currImData.downsampledImagePath);

        %% Get Bregma Coordinate or disqualify images
        response = getUserInput(im);
        currImData.bregmaCoord = str2double(response);
        currImData.exclude = isexclude(response);

        %% Update json file
        data{i} = currImData;
        savejson('', data, jsonPath);
    end


function rsp = getUserInput(im)
    %% Display image and get user response
    figure, imshow(im, [])
    title('Image Values Adjusted for Easier Visibility')

    rsp = inputdlg('Please Enter Bregma Coord (mm) or "exclude"');

    %% Ensure that response is valid
    if isnan(coord) && ~isexclude(rsp)
        warndlg('Please Input a Number or "exclude"!', 'Unacceptable Input')
        rsp = getUserInput(im);
    end

function bool = isexclude(str)
    bool = strcmpi(rsp, 'exclude');