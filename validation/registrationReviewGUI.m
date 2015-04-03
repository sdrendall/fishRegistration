function varargout = registrationReviewGUI(varargin)
    % REGISTRATIONREVIEWGUI MATLAB code for registrationReviewGUI.fig
    %      REGISTRATIONREVIEWGUI, by itself, creates a new REGISTRATIONREVIEWGUI or raises the existing
    %      singleton*.
    %
    %      H = REGISTRATIONREVIEWGUI returns the handle to a new REGISTRATIONREVIEWGUI or the handle to
    %      the existing singleton*.
    %
    %      REGISTRATIONREVIEWGUI('CALLBACK',hObject,eventData,handles,...) calls the local
    %      function named CALLBACK in REGISTRATIONREVIEWGUI.M with the given input arguments.
    %
    %      REGISTRATIONREVIEWGUI('Property','Value',...) creates a new REGISTRATIONREVIEWGUI or raises the
    %      existing singleton*.  Starting from the left, property value pairs are
    %      applied to the GUI before registrationReviewGUI_OpeningFcn gets called.  An
    %      unrecognized property name or invalid value makes property application
    %      stop.  All inputs are passed to registrationReviewGUI_OpeningFcn via varargin.
    %
    %      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
    %      instance to run (singleton)".
    %
    % See also: GUIDE, GUIDATA, GUIHANDLES
    
    % Edit the above text to modify the response to help registrationReviewGUI
    
    % Last Modified by GUIDE v2.5 24-Mar-2015 18:38:31
    
    % Begin initialization code - DO NOT EDIT
    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
                       'gui_Singleton',  gui_Singleton, ...
                       'gui_OpeningFcn', @registrationReviewGUI_OpeningFcn, ...
                       'gui_OutputFcn',  @registrationReviewGUI_OutputFcn, ...
                       'gui_LayoutFcn',  [] , ...
                       'gui_Callback',   []);
    if nargin && ischar(varargin{1})
        gui_State.gui_Callback = str2func(varargin{1});
    end
    
    if nargout
        [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
    else
        gui_mainfcn(gui_State, varargin{:});
    end
    % End initialization code - DO NOT EDIT
    
    
% --- Executes just before registrationReviewGUI is made visible.
function registrationReviewGUI_OpeningFcn(hObject, eventdata, handles, varargin)
    % This function has no output args, see OutputFcn.
    % varargin   command line arguments to registrationReviewGUI (see VARARGIN)
    json.startup()

    % Choose default command line output for registrationReviewGUI
    handles.output = hObject;

    % Default values
    handles.imageToDisplay = 'slice';

    % Update handles structure
    guidata(hObject, handles);
    
    % UIWAIT makes registrationReviewGUI wait for user response (see UIRESUME)
    % uiwait(handles.figure1);
    
% --- Outputs from this function are returned to the command line.
function varargout = registrationReviewGUI_OutputFcn(hObject, eventdata, handles) 
    % varargout  cell array for returning output args (see VARARGOUT);
    
    % Get default command line output from handles structure
    varargout{1} = handles.output;


% --- Executes on button press in loadData.
function loadData_Callback(hObject, eventdata, handles)
    % Loads data from metadata.json, loads the first image set, 
    %  and refreshes the gui
    try
       handles = loadJsonData(handles);
    catch err
        disp('loadData failed')
        return
    end
    handles = setCurrentImageSet(handles, 1);
    handles = loadCurrentImageSet(handles);
    refreshGUI(handles)
    guidata(hObject, handles)


% --- Executes on button press in nextIm.
function nextIm_Callback(hObject, eventdata, handles)
    handles.jsonData{handles.currentSetKey}.reviewed = true;
    updateJson(handles);
    handles = loadNextImage(handles);
    refreshGUI(handles)
    guidata(hObject, handles)


% --- Executes on button press in prevIm.
function prevIm_Callback(hObject, eventdata, handles)
    handles.jsonData{handles.currentSetKey}.reviewed = true;
    updateJson(handles);
    handles = loadPreviousImage(handles);
    refreshGUI(handles)
    guidata(hObject, handles)


% --- Executes when selected object is changed in sliceRating.
function sliceRating_SelectionChangeFcn(hObject, eventdata, handles)
    % Update the metadata with the value of the selected radiobutton
    currentSet = getCurrentSet(handles);
    switch get(eventdata.NewValue, 'Tag')
        case 'fullyUsable'
            currentSet.sliceUsable = 'yes';
        case 'partiallyUsable'
            currentSet.sliceUsable = 'partially';
        case 'notUsable'
            currentSet.sliceUsable = 'no';
    end       
    handles = updateCurrentSet(handles, currentSet);
    guidata(hObject, handles)


function registrationQualityScore_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of registrationQualityScore as text
%        str2double(get(hObject,'String')) returns contents of registrationQualityScore as a double
        currentSet = getCurrentSet(handles);
        currentSet.registrationQualityScore = str2double(get(hObject, 'String'));
        handles = updateCurrentSet(handles, currentSet);
        guidata(hObject, handles)    


% --- Executes during object creation, after setting all properties.
function registrationQualityScore_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


% --- Executes on button press in showReference.
function showReference_Callback(hObject, eventdata, handles)
    handles = setImageToDisplay(handles, 'reference');
    refreshDisplay(handles)
    guidata(hObject, handles)


% --- Executes on button press in showSlice.
function showSlice_Callback(hObject, eventdata, handles)
    handles = setImageToDisplay(handles, 'slice');
    refreshDisplay(handles)
    guidata(hObject, handles)


% --- Executes on button press in showLabels.
function showLabels_Callback(hObject, eventdata, handles)
    handles = setImageToDisplay(handles, 'labels');
    refreshDisplay(handles)
    guidata(hObject, handles)


% --- Executes on button press in displayExclusionMask_checkbox.
function displayExclusionMask_checkbox_Callback(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of displayExclusionMask_checkbox
    refreshDisplay(handles)


% --- Executes on button press in showOverlay.
function showOverlay_Callback(hObject, eventdata, handles)
    handles = setImageToDisplay(handles, 'overlay');
    refreshDisplay(handles)
    guidata(hObject, handles)


% --- Executes on button press in overlaySlice.
function overlaySlice_Callback(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of overlaySlice
    handles = generateOverlay(handles);
    if strcmp(handles.imageToDisplay, 'overlay')
        refreshDisplay(handles)
    end
    guidata(hObject, handles)


% --- Executes on button press in overlayReference.
function overlayReference_Callback(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of overlayReference
    handles = generateOverlay(handles);
    if strcmp(handles.imageToDisplay, 'overlay')
        refreshDisplay(handles)
    end
    guidata(hObject, handles)


function commentsBox_Callback(hObject, eventdata, handles)
% hObject    handle to commentsBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of commentsBox as text
%        str2double(get(hObject,'String')) returns contents of commentsBox as a double
    currentSet = getCurrentSet(handles);
    currentSet.qcComments = get(hObject, 'String');
    handles = updateCurrentSet(handles, currentSet);
    guidata(hObject, handles)    


% --- Executes during object creation, after setting all properties.
function commentsBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to commentsBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


% --- Executes on button press in excludeLeft_checkbox.
function excludeLeft_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to excludeLeft_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of excludeLeft_checkbox
    currentSet = getCurrentSet(handles);
    currentSet.excludeLeftHemisphere = get(hObject, 'Value');
    handles = updateCurrentSet(handles, currentSet);
    guidata(hObject, handles)


% --- Executes on button press in excludeRight_checkbox.
function excludeRight_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to excludeRight_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of excludeRight_checkbox
    currentSet = handles.jsonData{handles.currentSetKey};
    currentSet.excludeRightHemisphere = get(hObject, 'Value');
    handles = updateCurrentSet(handles, currentSet);
    guidata(hObject, handles)


% --- Executes on button press in excludeRegions_button.
function excludeRegions_button_Callback(hObject, eventdata, handles)
% hObject    handle to excludeRegions_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    axes(handles.mainDisplay)
    button = 1;
    % Bleh.. breaks when button is empty (enter is pressed)
    while true
        [x, y, button] = ginput(1);
        if isempty(button)
            break
        else
            id = getIdFromSelection(handles, x, y);
            hem = getHemisphereFromSelection(handles, x, y);
            handles = updateExclusionIdsAndMask(handles, id, hem, button);
        end
    end
    guidata(hObject, handles)


% --- Executes when selected object is changed in filterReviewed_buttonGroup.
function filterReviewed_buttonGroup_SelectionChangedFcn(hObject, eventdata, handles)
    % There are two modes for storing the current state of the program.
    % When reviewed and unreviewed sections are being shown, an index to the current set number is stored
    % When only unreviewed sections are being shown, a list of unreviewed keys is kept, when a section is reviewed
    %  it's key is popped from the list.  When the list is exhausted, a message is displayed
    %
    % The current key always remains the same.

    keys = getCleanJsonDataKeys(handles);
    currentKeyInSet = keys == handles.currentSetKey;

    switch get(hObject, 'Tag')
        case 'showReviewed_radio'
            handles.currentKeyIndex = find(currentKeyInSet, 1);
        case 'skipReviewed_radio'
            keys(currentKeyInSet) = [];
    end

    handles.jsonKeys = keys;
    handles.numberOfKeys = length(keys);
    guidata(hObject, handles)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STATE CHANGING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function handles = loadNextImage(handles)
    switch get(get(handles.filterReviewed_buttonGroup, 'SelectedObject'), 'Tag')
    case 'showReviewed_radio'
        handles = incrementCurrentKeyIndex(handles);
    case 'skipReviewed_radio'
        handles = popCurrentSetKey(handles);
    end
    handles = loadCurrentImageSet(handles);
    refreshGUI(handles)


function handles = loadPreviousImage(handles)
    switch get(get(handles.filterReviewed_buttonGroup, 'SelectedObject'), 'Tag')
    case 'showReviewed_radio'
        handles = decrementCurrentKeyIndex(handles);
    case 'skipReviewed_radio'
        handles = popCurrentSetKey(handles);
    end
    handles = loadCurrentImageSet(handles);
    refreshGUI(handles)


function handles = incrementCurrentKeyIndex(handles)
    if handles.currentKeyIndex == handles.numberOfKeys
        nextSetKeyIndex = 1;
    else
        nextSetKeyIndex = handles.currentKeyIndex + 1;
    end
    handles = setCurrentImageSet(handles, nextSetKeyIndex);


function handles = decrementCurrentKeyIndex(handles)
    if handles.currentKeyIndex == 1
        nextSetKeyIndex = handles.numberOfKeys;
    else 
        nextSetKeyIndex = handles.currentKeyIndex - 1;
    end
    handles = setCurrentImageSet(handles, nextSetKeyIndex);


function handles = popCurrentSetKey(handles)
    if ~isempty(handles.jsonKeys)
        handles.currentSetKey = handles.jsonKeys(end);
        handles.jsonKeys(end) = [];
    else
        msgbox('All images have been reviewed')
    end
    handles.numberOfKeys = length(handles.jsonKeys);


function handles = setImageToDisplay(handles, value)
    handles.imageToDisplay = value;


function handles = setCurrentImageSet(handles, index)
    handles.currentKeyIndex = index;
    key = handles.jsonKeys(index);
    handles.currentSetKey = key;


function currentSet = getCurrentSet(handles)
    currentSet = handles.jsonData{handles.currentSetKey};


function handles = updateCurrentSet(handles, set)
    handles.jsonData{handles.currentSetKey} = set;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DATA LOADING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% These functions load data from the current set (loading isn't really the best word, should have said get)
%  They aren't responsible for changing the state of the program, or the appearance of the GUI


function handles = loadCurrentImageSet(handles)
    % Loads all of the data associated with the current image set to the handles object
    % There are a large number of 'load' functions to handle setting defaults in each case
    handles = loadImages(handles);
    handles = loadScores(handles);
    handles = loadComments(handles);
    handles = loadExclusionSettings(handles);


function handles = loadScores(handles)
    % Loads the Registration Quality Scores associated with an image
    %  sets the score to 5 if no field is found
    currentSet = getCurrentSet(handles);
    if ~isfield(currentSet, 'sliceUsable')
        currentSet.sliceUsable = 'no';
    end
    if ~isfield(currentSet, 'registrationQualityScore')
        currentSet.registrationQualityScore = 5;
    end
    handles = updateCurrentSet(handles, currentSet);


function handles = loadComments(handles)
    % Loads the comments field associated with the current image set
    currentSet = getCurrentSet(handles);
    if ~isfield(currentSet, 'qcComments')
        return
    else
        set(handles.commentsBox, 'String', currentSet.qcComments)
    end
    handles = updateCurrentSet(handles, currentSet);


function handles = loadExclusionSettings(handles)
    % Loads the exclusion settings associated with an image set
    currentSet = getCurrentSet(handles);
    if ~isfield(currentSet, 'excludeLeftHemisphere')
        currentSet.excludeLeftHemisphere = false;
    end
    if ~isfield(currentSet, 'excludeRightHemisphere')
        currentSet.excludeRightHemisphere = false;
    end
    if ~isfield(currentSet, 'regionIdsToExclude')
        currentSet.regionIdsToExclude = [];
    end
    handles = updateCurrentSet(handles, currentSet);
    handles = generateExclusionMask(handles);


function handles = loadImages(handles)
    % Loads the set of images specified in one entry of the json metadata file
    % Also regenerates the overlay image stored in handles.  Ensures that all images in handles are current
    currentSet = getCurrentSet(handles);
    handles = loadReferenceImage(handles, currentSet);
    handles = loadSliceImage(handles, currentSet);
    handles = loadReferenceLabels(handles, currentSet);
    handles = generateOverlay(handles);


function handles = loadReferenceImage(handles, setData)
    imPath = setData.registeredAtlasReferenceImagePath;
    handles.refIm = loadRegistrationImage(handles, imPath);


function handles = loadSliceImage(handles, setData)
    imPath = setData.downsampledImagePath;
    handles.sliceIm = loadRegistrationImage(handles, imPath);
    handles.hemiIm = true(size(handles.sliceIm));  % TODO: Temporary - for debugging


function handles = loadReferenceLabels(handles, setData)
    imPath = setData.registeredAtlasLabelsPath;
    handles.labelIm = loadRegistrationImage(handles, imPath);


function im = loadRegistrationImage(handles, imPath)
    % image paths specified in metadata.json are relative to the 'experimentPath' used in the registration pipeline
    % since the folder containing these images may be out of context, it is assumed that the images to be evaluated
    % are in the same folder as the .json file the metadata is being read from.  This function attempts to load the 
    % image with the filename specified in imPath, located in the directory where the json file was read from.
    imName = getFilename(imPath);
    try
        im = readImage(fullfile(handles.basePath, imName));
    catch err
        warning(['Could not read image at ', fullfile(handles.basePath, imName), ' perhaps it doesnt exist.']);
        im = 0;
    end


function im = readImage(filepath)
    % Uses the appropriate load function based on the filetype found at the filepath
    disp(['loading ', filepath])
    [~, ~, ext] = fileparts(filepath);
    if strcmpi(ext, '.mhd')
        imageObj = read_mhd(filepath);
        im = imageObj.data;
    else
        im = imread(filepath);
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% JSON FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function handles = loadJsonData(handles)
    % Load data from the json file in the specified location
    try
        [handles.basePath, handles.jsonPath] = getJsonPath();
    catch err
        disp(err.message)
        % rethrow(err)
    end

    disp(['Loading data from ', handles.jsonPath])
    jsonData = json.read(handles.jsonPath);

    % The new JSON API I am using will try to import a list of json objects as a structure array.
    %  This is not ideal for a number of reasons, first and formost being that this code expects a cell
    %  array of structures.  Here, I convert the jsonData array to a cell if it is not one.
    if ~iscell(jsonData)
        structArray = jsonData;
        jsonData = cell(size(jsonData));
        for i = 1:length(jsonData)
            jsonData(i) = {structArray(i)};
        end
    end

    % Set the reviewed field to false for data sets that don't have it
    for i = 1:length(jsonData)
        if ~isfield(jsonData{i}, 'reviewed')
            jsonData{i}.reviewed = false;
        end
    end
     
    handles.jsonData = jsonData;
    handles.jsonKeys = getCleanJsonDataKeys(handles);
    handles.numberOfKeys = length(handles.jsonKeys);


function updateJson(handles)
    json.write(handles.jsonData, handles.jsonPath);


function [basePath, jsonPath] = getJsonPath()
    % Query user for the experiment path
    basePath = uigetdir();
    % Cancel loading if the user closes the dialog box
    if ~ischar(basePath) % uigetdir returns 0 if the action is cancelled.  
        throw(MException('LoadActionTerminated', 'LoadData was terminated by the user'))
    end
    jsonPath = searchForJson(basePath);


function jsonPath = searchForJson(basePath)
    % Finds json files in the directory at basePath
    jsonFiles = dir(fullfile(basePath, '*.json')); % I'll make this recursive if things get more complicated
    jsonFile = chooseJsonFile(jsonFiles);          % If multiple files are returned, one must be chosen.  This will also handle warnings
    jsonPath = fullfile(basePath, jsonFile.name);


function jsonFile = chooseJsonFile(jsonFiles)
    % Find the json file named 'metadata.json', if that doesnt exist, choose another one
    % Currently not safe for jsonFiles coming from multiple directories
    if length(jsonFiles) > 1
        warning('Multiple .json files found in the spcified directory')
        isMetadataDotJson = strcmpi({jsonFiles(:).name}, 'metadata.json');
        if any(isMetadataDotJson)
            jsonFile = jsonFiles(isMetadataDotJson);
        else
            warning('No file named metadata.json found in the spcified directory')
            jsonFile = jsonFiles(1);
        end
    elseif isempty(jsonFiles)
        warning('No .json file found in the specified directory')
        jsonFile = jsonFiles;
    else
        jsonFile = jsonFiles;
    end


function filename = getFilename(filePath)
    % Returns the filename, with extension, of a file located at filePath
    [~, name, extension] = fileparts(filePath);
    filename = [name, extension];


function keys = getCleanJsonDataKeys(handles)
    % Generates an array of indicies corresponding to the metadata for usable image sets
    % Eliminating these sets would cause them to be removed when the json file is updated
    % Using the currently stored data.

    % Images that have been marked to be excluded, or unsuccessful registration attempts should be ignored
    inputData = handles.jsonData;
    keys = [];
    for i = 1:length(inputData)
        if includeDataSet(inputData{i}, handles)
            keys(end + 1) = i;
        end
    end


function include_bool = includeDataSet(dataSet, handles)
    try
        include_bool = ~dataSet.exclude && dataSet.registrationSuccessful;

        if strcmp(get(get(handles.filterReviewed_buttonGroup, 'SelectedObject'), 'Tag'), 'skipReviewed_radio')
            include_bool = include_bool && ~dataSet.reviewed;
        end

    catch err
        warning('Error encountered when loading ', dataSet.vsiPath, ' so it will not be loaded.')
        include_bool = false;
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% DISPLAY FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function handles = generateOverlay(handles)
    % Checks the status of the overlay checkboxes, and generates an overlay image
    % Does not refresh the display
    overlay = zeros([size(handles.refIm), 3]);
    if get(handles.overlaySlice, 'value')
        overlay(:,:,2) = mat2gray(handles.sliceIm);
    end
    if get(handles.overlayReference, 'value')
        overlay(:,:,1) = mat2gray(handles.refIm);
    end
    handles.overlayIm = overlay;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% REFRESH FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% These functions update the appearance of the gui

function refreshGUI(handles)
    % Convenience function to refresh everything
    refreshDisplay(handles)
    refreshScores(handles)
    refreshExclusionSettings(handles)


function refreshDisplay(handles)
    % refreshes the currently displayed image set
    axes(handles.mainDisplay)

    switch handles.imageToDisplay
        case 'slice'
            imshow(handles.sliceIm, [])
        case 'reference'
            imshow(handles.refIm, [])
        case 'labels'
            imagesc(handles.labelIm)
            colormap(lines)
            axis off
            axis equal
        case 'overlay'
            imshow(handles.overlayIm)
    end

    if get(handles.displayExclusionMask_checkbox, 'Value') && isfield(handles, 'exclusionMask')
        alphamask(handles.exclusionMask, [1 0 0], 0.5);
    end


function refreshScores(handles)
    currentSet = getCurrentSet(handles);
    set(handles.registrationQualityScore,'String', num2str(currentSet.registrationQualityScore))
    switch currentSet.sliceUsable
        case 'yes'
            set(handles.fullyUsable, 'Value', 1)
        case 'partially'
            set(handles.partiallyUsable, 'Value', 1)
        case 'no'
            set(handles.notUsable, 'Value', 1)
    end


function refreshExclusionSettings(handles)
    currentSet = getCurrentSet(handles);
    set(handles.excludeLeft_checkbox, 'Value', currentSet.excludeLeftHemisphere);
    set(handles.excludeRight_checkbox, 'Value', currentSet.excludeRightHemisphere);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% REGION EXCLUSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function id = getIdFromSelection(handles, x, y);
    x = round(x);
    y = round(y);
    maskSize = size(handles.exclusionMask);
    [y, x] = constrainToSize(maskSize, y, x);
    id = handles.labelIm(y, x);


function hem = getHemisphereFromSelection(handles, x, y)
    x = round(x);
    y = round(y);
    maskSize = size(handles.exclusionMask);
    [y, x] = constrainToSize(maskSize, y, x);
    hem = handles.hemiIm(y, x);


function [irow, icol] = constrainToSize(s, irow, icol);
% Given the size of a 2d array, a row index, and a column index, ensures that the 
% Row and column indicies are within the bounds of the given array such that they can
% be used to index that array
    if irow < 1
        irow = 1;
    elseif irow > s(1)
        irow = s(1);
    end
    if icol < 1
        icol = 1;
    elseif icol > s(2)
        icol = s(2);
    end


function handles = updateExclusionIdsAndMask(handles, id, hem, button)
    currentSet = getCurrentSet(handles);
    % button == 1 indicates a left click, this adds a region to the mask (excludes a region)
    if button == 1 && isempty(currentSet.regionIdsToExclude) || ~ismember([id, hem], currentSet.regionIdsToExclude, 'rows')
        currentSet.regionIdsToExclude(end + 1, :) = [id, hem];
        handles.exclusionMask(handles.labelIm == id & handles.hemiIm == hem) = 1;
    % button == 2 indicates a right click, removing a region from the exclusion list
    elseif button == 3 && ismember([id, hem], currentSet.regionIdsToExclude, 'rows')
        entryToRemove = currentSet.regionIdsToExclude(:, 1) == id & currentSet.regionIdsToExclude(:, 2) == hem;
        currentSet.regionIdsToExclude(entryToRemove, :) = [];
        handles.exclusionMask(handles.labelIm == id & handles.hemiIm == hem) = 0;
    end
    handles = updateCurrentSet(handles, currentSet);
    refreshDisplay(handles)
        

function handles = generateExclusionMask(handles)
    % Generates a mask corresponding to the regions listed in regionIdsToExclude
    currentSet = getCurrentSet(handles);
    mask = false(size(handles.sliceIm));
    if ~isempty(currentSet.regionIdsToExclude)
        for i = 1:size(currentSet.regionIdsToExclude, 1)
            structId = currentSet.regionIdsToExclude(i, 1);
            hemiId = currentSet.regionIdsToExclude(i, 2);
            mask = mask | (handles.labelIm == structId & handles.hemiIm == hemiId);
        end
    end
    handles.exclusionMask = mask;
