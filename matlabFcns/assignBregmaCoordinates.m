function assignBregmaCoordinates(experimentPath)

    %% Load JSON data
    jsonPath = fullfile(experimentPath, '.registrationData', 'metadata.json');
    data = loadjson(jsonPath);

    for i = 1:length(data)
        %% Load each downsampled image
        currImData = data{i};
        im = imread(currImData.downsampledImagePath);

        [~, name] = fileparts(currImData.downsampledImagePath);

        %% Get Bregma Coordinate or disqualify images
        response = getUserInput(im, name);
        currImData.bregmaCoord = str2double(response);
        currImData.atlasCoord = bregmaToAtlas(currImData.bregmaCoord)
        currImData.exclude = isexclude(response);

        %% Update json file
        data{i} = currImData;
        savejson('', data, jsonPath);
    end


function rsp = getUserInput(im, name)
    %% Display image and get user response
    imshow(im, [])
    %title('Image Values Adjusted for Easier Visibility')
    title(name)


    rsp = inputdlg('Please Enter Bregma Coord (mm) or "exclude"');
    
    % if cancel is pressed
    if isempty(rsp)
        b = questdlg('Exit Bregma Coordinate Assignment?', 'Exit?', 'Yes');
        switch b
        case 'Yes'
            exit
        otherwise
            rsp = getUserInput(im, name);
        end
    end

    %% Ensure that response is valid, empty cell if cancel is pressed
    if isnan(str2double(rsp)) && ~isexclude(rsp)
        disp('Please Input a Number or "exclude"!')
        rsp = getUserInput(im, name);
    end

function bool = isexclude(str)
    bool = strcmpi(str, 'exclude');

function atlasCoord = bregmaToAtlas(bCoord)
    bregmaInAtlas = 5525; % in um
    bCoord = mm2um(bCoord);
    atlasCoord = bregmaInAtlas - bCoord;

function um = mm2um(mm)
    um = mm*1000;