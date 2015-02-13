function cropRegionsUsingAtlas(experimentPath, ids, varargin)

	%% Parse input
	args = parseVarargin(varargin);

	%% Load metadata
	if ~exist('experimentPath', 'var')
		disp('Experiment Path Not Specified. Please Choose Experiment Path')
		experimentPath = uigetdir;
	end
	
	cd(experimentPath)
	jsonPath = generateJsonPath(experimentPath);
	regData = loadRegistrationData(jsonPath);

	%% Crop Each Image
	for i = 1:length(regData)
		%% Load Image Data
		currentSet = regData{i};
		annotations = readImage(fullfile(experimentPath, currentSet.registeredAtlasLabelsPath));
	
		%% Generate Region Mask
		mask = getRegionsById(annotations, TeA_ids);
		%% Only continue if there are actually pixels to crop
		if any(mask(:))
			section = loadvsi(fullfile(experimentPath, currentSet.vsiPath));
			
			fullSizeMask = rescaleMask(mask, section);
			teaRegions = regionprops(fullSizeMask, 'BoundingBox', 'Image');
		
			croppedRegions = cropRegions(section, teaRegions);
		
			outputBase = generateBaseOutputPath(currentSet.vsiPath);
			saveCroppedRegions(croppedRegions, outputBase);
		else
			disp(['Region not found in ', currentSet.vsiPath, '.  Advancing to next image.'])
		end
	end