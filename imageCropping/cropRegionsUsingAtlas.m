function cropRegionsUsingAtlas(ids, varargin)
	% cropRegionsUsingAtlas(ids, [options])
	%
	% 'slice order' - 1 x n array containing the permutation of the indicies 1:n 
	%  		correpsonding to the order in which planes in the vsi files should be loaded
	%       defaults to [3 2 1]
	% 'split channels' - boolean, output channels individually, or as an rgb image
	%					default false
	% 'experiment path' - string

	%% Parse input
	args = parseVarargin(varargin, nargin);

	%% Load metadata
	if isempty(args.exp_path)
		disp('Experiment Path Not Specified. Please Choose Experiment Path')
		args.exp_path = uigetdir;
	end
	
	cd(args.exp_path)
	jsonPath = generateJsonPath(args.exp_path);
	regData = loadRegistrationData(jsonPath);

	%% Crop Each Image
	for i = 1:length(regData)
		currentSet = regData{i};
		try
			%% Load Image Data
			annotations = readImage(fullfile(args.exp_path, currentSet.registeredAtlasLabelsPath));
		
			%% Generate Region Mask
			mask = getRegionsById(annotations, ids);
			%% Only continue if there are actually pixels to crop
			if any(mask(:))
				planes = unpackvsi(fullfile(args.exp_path, currentSet.vsiPath), args);
				
				fullSizeMask = rescaleMask(mask, planes{1}); % Resize the mask to the size of a plane
				regionProperties = regionprops(fullSizeMask, 'BoundingBox', 'Image');
			
				croppedRegions = cropRegionsFromPlanes(planes, regionProperties);
			
				outputBase = generateBaseOutputPath(currentSet.vsiPath);
				saveCroppedRegions(croppedRegions, outputBase);
			else
				disp(['Region not found in ', currentSet.vsiPath, '.  Advancing to next image.'])
			end
		catch err
			disp(['Warning!, cropping failed for ', currentSet.vsiPath]);
		end
	end


function args = parseVarargin(argv, argc)
	args = struct( ...
		'slice_order', [3, 2, 1], ...
		'split_channels', false,  ...
		'exp_path', '');
	if argc < 2
		return
	end

	i = 1;
	while i <= argc - 1
		switch lower(argv{i})
		case 'slice order'
			args.slice_order = argv{i + 1};
			i = i + 2;
		case 'split channels'
			args.split_channels = true;
			i = i + 1;
		case 'experiment path'
			args.exp_path = argv{i + 1};
			i = i + 2;
		end
	end


function mask = rescaleMask(mask, section) 
	disp('Resizing Binary Mask.....')
	[nr, nc, ~] = size(section);
	sectionSize = [nr, nc];
	atlasScale = .64497/25;
	maskSize = size(mask);

	% Determine image dimension ratio
	startPoint = ceil((size(mask) - sectionSize.*atlasScale)./2);
	endPoint = ceil(startPoint + sectionSize.*atlasScale);
	startPoint(startPoint < 1) = 1;
	endPoint(endPoint > maskSize) = maskSize(endPoint > maskSize);

	% Rescale
	mask = mask(startPoint(1):endPoint(1), startPoint(2):endPoint(2));
	mask = imresize(mask, [nr, nc], 'nearest');

	% Increasing the size of a binary image leaves 'pixelated' edges.
	% To remedy this, the morphological opening operation can be applied to smooth sharp edges
	% The size of the structuring element used in the opening operation is determined using the
	% Ratio between the number of pixels in the first dimension of the image
	pixelRatio = nr/(endPoint(1) - startPoint(1));
	mask = imopen(mask, strel('disk', round(pixelRatio)));


function jsonPath = generateJsonPath(experimentPath)
	disp('Generating Json Path.....')
	jsonPath = fullfile(experimentPath, '/.registrationData/', 'metadata.json');


function saveCroppedRegions(croppedRegions, basePath)
	disp('Saving Cropped Regions.....')
	for i = 1:length(croppedRegions)
		savePath = [basePath, '_croppedTeA_', num2str(i), '.tiff'];
		disp(['Saving ', savePath])
		imwrite(croppedRegions{i}, savePath)
	end


function basePath = generateBaseOutputPath(vsiPath)
	disp('Generating Output Paths.....')
	[baseDir, baseName] = fileparts(vsiPath);
	basePath = fullfile(baseDir, baseName);


function croppedRegions = cropRegionsFromPlanes(planes, regions)
	croppedRegions = {};
	for i = 1:length(planes)
		croppedRegions = [croppedRegions; cropRegions(planes{i}, regions)];
	end


function croppedRegions = cropRegions(im, regions)
	disp('Cropping Image.....')
	croppedRegions = cell([length(regions), 1]);
	for i = 1:length(regions)
		bb = regions(i).BoundingBox;
		croppedRegion = im(bb(2):bb(2)+bb(4)-1, bb(1):bb(1) + bb(3)-1, :);
		% Set non-region pixels to 0
		for j = 1:size(croppedRegion, 3)
			croppedRegion(:,:,j) = croppedRegion(:,:,j).*regions(i).Image;
		end
		croppedRegions{i} = croppedRegion;
	end


function mask = getRegionsById(annotations, ids)
	disp('Identifying Region Pixels.....')
	mask = false(size(annotations));
	for i = 1:length(ids)
		mask = mask | (annotations == ids(i));
	end


function regData = loadRegistrationData(jsonPath)
	% Read the data from the specified json file then filter out the useless images
	disp('Loading Metadata.....')
	regJson = loadjson(jsonPath);
	disp('Disqualifying Image Sets.....')
	regData = {};
	for i = 1:length(regJson)
		set = regJson{i};
		if ~set.exclude && set.registrationSuccessful
			regData{end + 1} = set;
		end
	end


function im = readImage(filepath)
    % Uses the appropriate load function based on the filetype found at the filepath
    [~, ~, ext] = fileparts(filepath);
    if strcmpi(ext, '.mhd')
        imageObj = read_mhd(filepath);
        im = imageObj.data';
    else
        im = imread(filepath);
    end


function planes = unpackvsi(filepath, args)
	disp(['Loading ', filepath, '.....'])
	r = bfGetReader(filepath);
	np = r.getImageCount();
	if args.split_channels
		planes = cell([np, 1]);
		for i = 1:np
			planes{i} = mat2gray(bfGetPlane(r, i));
		end
	else
		planes = {loadRgbVsi(r, args)};
	end

function im = loadRgbVsi(reader, args)
	nr = reader.getSizeY;
	nc = reader.getSizeX;
	np = reader.getImageCount;
	im = zeros([nr, nc, 3]);

	for i = 1:np
		im(:,:,args.slice_order(i)) = mat2gray(bfGetPlane(reader, i));
	end	