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
    
    % Last Modified by GUIDE v2.5 29-Jan-2015 16:53:29
    
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

% --- Executes on button press in nextIm.
function nextIm_Callback(hObject, eventdata, handles)
    updateJson(handles);
    if handles.currentKeyIndex == handles.numberOfKeys
        nextSetKeyIndex = 1;
    else
        nextSetKeyIndex = handles.currentKeyIndex + 1;
    end
    handles = setCurrentImageSet(hObject, handles, nextSetKeyIndex);
    handles = loadCurrentImageSet(hObject, handles);
    handles = loadScores(hObject, handles);
    refreshDisplay(handles)
    refreshScores(handles)

% --- Executes on button press in prevIm.
function prevIm_Callback(hObject, eventdata, handles)
    updateJson(handles);
    if handles.currentKeyIndex == 1
        nextSetKeyIndex = handles.numberOfKeys;
    else 
        nextSetKeyIndex = handles.currentKeyIndex - 1;
    end
    handles = setCurrentImageSet(hObject, handles, nextSetKeyIndex);
    handles = loadCurrentImageSet(hObject, handles);
    handles = loadScores(hObject, handles);
    refreshDisplay(handles)
    refreshScores(handles)

% --- Executes on button press in showReference.
function showReference_Callback(hObject, eventdata, handles)
    handles = setImageToDisplay(hObject, handles, 'reference');
    refreshDisplay(handles)

% --- Executes on button press in showSlice.
function showSlice_Callback(hObject, eventdata, handles)
    handles = setImageToDisplay(hObject, handles, 'slice');
    refreshDisplay(handles)

% --- Executes on button press in showLabels.
function showLabels_Callback(hObject, eventdata, handles)
    handles = setImageToDisplay(hObject, handles, 'labels');
    refreshDisplay(handles)

% --- Executes when selected object is changed in sliceRating.
function sliceRating_SelectionChangeFcn(hObject, eventdata, handles)
    % Update the metadata with the value of the selected radiobutton
    currentSet = handles.jsonData{handles.currentSetNo};
    switch get(eventdata.NewValue, 'Tag')
        case 'fullyUsable'
            currentSet.sliceUsable = 'yes';
        case 'partiallyUsable'
            currentSet.sliceUsable = 'partially';
        case 'notUsable'
            currentSet.sliceUsable = 'no';
    end       
    handles.jsonData{handles.currentSetNo} = currentSet;
    guidata(hObject, handles)

function registrationQualityScore_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of registrationQualityScore as text
%        str2double(get(hObject,'String')) returns contents of registrationQualityScore as a double
        currentSet = handles.jsonData{handles.currentSetNo};
        currentSet.registrationQualityScore = str2double(get(hObject, 'String'));
        handles.jsonData{handles.currentSetNo} = currentSet;
        guidata(hObject, handles)    

% --- Executes during object creation, after setting all properties.
function registrationQualityScore_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

% --- Executes on button press in loadData.
function loadData_Callback(hObject, eventdata, handles)
    try
       handles = loadJsonData(hObject, handles);
    catch err
        disp('loadData failed')
        return
    end
    handles = setCurrentImageSet(hObject, handles, 1);
    handles = loadCurrentImageSet(hObject, handles);
    handles = loadScores(hObject, handles);
    refreshDisplay(handles)
    refreshScores(handles)

% --- Executes on button press in overlaySlice.
function overlaySlice_Callback(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of overlaySlice
    handles = regenerateOverlay(hObject, handles);
    if strcmp(handles.imageToDisplay, 'overlay')
        refreshDisplay(handles)
    end

% --- Executes on button press in overlayReference.
function overlayReference_Callback(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of overlayReference
    handles = regenerateOverlay(hObject, handles);
    if strcmp(handles.imageToDisplay, 'overlay')
        refreshDisplay(handles)
    end

% --- Executes on button press in showOverlay.
function showOverlay_Callback(hObject, eventdata, handles)
    handles = setImageToDisplay(hObject, handles, 'overlay');
    refreshDisplay(handles)

function handles = setImageToDisplay(hObject, handles, value)
    handles.imageToDisplay = value;
    guidata(hObject, handles)

function handles = setCurrentImageSet(hObject, handles, index)
    handles.currentKeyIndex = index;
    key = handles.jsonKeys(index);
    handles.currentSetNo = key;
    guidata(hObject, handles)

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

function refreshScores(handles)
    currentSet = handles.jsonData{handles.currentSetNo};
    set(handles.registrationQualityScore,'String', num2str(currentSet.registrationQualityScore))
    switch currentSet.sliceUsable
        case 'yes'
            set(handles.fullyUsable, 'Value', 1)
        case 'partially'
            set(handles.partiallyUsable, 'Value', 1)
        case 'no'
            set(handles.notUsable, 'Value', 1)
    end

function handles = regenerateOverlay(hObject, handles)
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
    guidata(hObject, handles)

function handles = loadImageSet(hObject, handles, setData)
    % Loads the set of images specified in one entry of the json metadata file
    % Also regenerates the overlay image stored in handles.  Ensures that all images in handles are current
    handles = loadReferenceImage(hObject, handles, setData);
    handles = loadSliceImage(hObject, handles, setData);
    handles = loadReferenceLabels(hObject, handles, setData);
    handles = regenerateOverlay(hObject, handles);

function handles = loadCurrentImageSet(hObject, handles)
    currentSet = handles.jsonData{handles.currentSetNo};
    handles = loadImageSet(hObject, handles, currentSet);

function handles = loadReferenceImage(hObject, handles, setData)
    imPath = setData.registeredAtlasReferenceImagePath;
    handles.refIm = getRegistrationImage(handles, imPath);
    guidata(hObject, handles);

function handles = loadSliceImage(hObject, handles, setData)
    imPath = setData.downsampledImagePath;
    handles.sliceIm = getRegistrationImage(handles, imPath);
    guidata(hObject, handles);

function handles = loadReferenceLabels(hObject, handles, setData)
    imPath = setData.registeredAtlasLabelsPath;
    handles.labelIm = getRegistrationImage(handles, imPath);
    guidata(hObject, handles);

function im = getRegistrationImage(handles, imPath)
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

function handles = loadScores(hObject, handles)
    currentSet = handles.jsonData{handles.currentSetNo};
    if ~isfield(currentSet, 'sliceUsable')
        currentSet.sliceUsable = 'no';
    end
    if ~isfield(currentSet, 'registrationQualityScore')
        currentSet.registrationQualityScore = 5;
    end
    handles.jsonData{handles.currentSetNo} = currentSet;
    guidata(hObject, handles)
    refreshScores(handles)

function updateJson(handles)
    savejson('', handles.jsonData, handles.jsonPath);

function handles = loadJsonData(hObject, handles)
    % Load data from the json file in the specified location
    try
        [handles.basePath, handles.jsonPath] = getJsonPath();
    catch err
        disp(err.message)
        rethrow(err)
    end
    disp(['Loading data from ', handles.jsonPath])
    handles.jsonData = loadjson(handles.jsonPath);
    handles.jsonKeys = getCleanJsonDataKeys(handles.jsonData);
    handles.numberOfKeys = length(handles.jsonKeys);

    % Update the handles object
    guidata(hObject, handles)

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
    jsonFile = chooseJsonFile(jsonFiles);              % If multiple files are returned, one must be chosen.  This will also handle warnings
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

function keys = getCleanJsonDataKeys(inputData)
    % Generates an array of indicies corresponding to the metadata for usable image sets
    % Eliminating these sets would cause them to be removed when the json file is updated
    % Using the currently stored data.
    keys = [];
    for i = 1:length(inputData)
        if ~inputData{i}.exclude && inputData{i}.registrationSuccessful
            keys(end + 1) = i;
        end
    end

function filename = getFilename(filePath)
    % Returns the filename, with extension, of a file located at filePath
    [~, name, extension] = fileparts(filePath);
    filename = [name, extension];
