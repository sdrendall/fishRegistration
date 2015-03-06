function cropRegionsUsingAtlas(vsiPath, labelsPath, ids, varargin)
	% cropRegionsUsingAtlas(vsiPath, labelsPath, ids, [options])
	%
	% 'slice order' - 1 x n array containing the permutation of the indicies 1:n 
	%  		correpsonding to the order in which planes in the vsi files should be loaded
	%       defaults to [3 2 1]
	% 'split channels' - boolean, output channels individually, or as an rgb image
	%					default false
	% 'experiment path' - string


	%% Parse input
	args = parseVarargin(varargin, nargin - 3);

	%% Crop Each Image
	try
		%% Load Image Data
		annotations = readImage(fullfile(args.exp_path, labelsPath));

		%% Generate Region Mask
		mask = getRegionsById(annotations, ids);
		%% Only continue if there are actually pixels to crop
		if any(mask(:))
			planes = unpackvsi(fullfile(args.exp_path, vsiPath), args);

			fullSizeMask = rescaleMask(mask, planes{1}); % Resize the mask to the size of a plane
			regionProperties = regionprops(fullSizeMask, 'BoundingBox', 'Image');

			croppedRegions = cropRegionsFromPlanes(planes, regionProperties);

			outputBase = generateBaseOutputPath(vsiPath, args);
			saveCroppedRegions(croppedRegions, outputBase, args);
		else
			disp(['Region not found in ', vsiPath, '.  Advancing to next image.'])
		end
	catch err
		disp(['Warning!, cropping failed for ', vsiPath]);
		rethrow(err)
	end


function args = parseVarargin(argv, argc)
	args = struct( ...
		'slice_order', [3, 2, 1], ...
		'split_channels', false,  ...
		'exp_path', '', ...
		'out_path', '', ...
		'region_name', '');

	i = 1;
	while i <= argc
		switch lower(argv{i})
		case 'slice order'
			args.slice_order = argv{i + 1};
			i = i + 2;
		case 'split channels'
			args.split_channels = argv{i + 1};
			i = i + 2;
		case 'experiment path'
			args.exp_path = argv{i + 1};
			i = i + 2;
		case 'output path'
		    args.out_path = argv{i + 1};
		    i = i + 2;
		case 'region name'
		    args.region_name = argv{i + 1};
		    i = i + 2;
		otherwise
		    disp('Argument unrecognized')
		    disp(argv{i})
		    i = i + 1;
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
	mask = imopen(mask, strel('disk', 2*round(pixelRatio)));


function jsonPath = generateJsonPath(experimentPath)
	disp('Generating Json Path.....')
	jsonPath = fullfile(experimentPath, '/.registrationData/', 'metadata.json');


function saveCroppedRegions(croppedRegions, basePath, args)
	disp('Saving Cropped Regions.....')
	ensureSaveLocation(basePath);
	for i = 1:length(croppedRegions)
		savePath = [basePath, '_cropped_', args.region_name, '_', num2str(i), '.tiff'];
		disp(['Saving ', savePath])
		imwrite(croppedRegions{i}, savePath)
	end


function ensureSaveLocation(savePath)
    saveDir = fileparts(savePath);
    [~, ~] = mkdir(saveDir)  % Nullifying the 2nd output prevents Warning messages from being printed if the dir exists


function basePath = generateBaseOutputPath(vsiPath, args)
    % Creates the base output path based on the vsiPath, and the specified output location.
    % The directory structure, relative to the experiment path, is used to maintain some level of organization
    % args.out_path is relative to args.exp_path.
    % I'm not sure what will happen if args.exp_path is behind a symbolic link.
	disp('Generating Output Paths.....')
	[baseDir, baseName] = fileparts(vsiPath);
	basePath = fullfile(args.exp_path, args.out_path, baseDir, baseName);


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
        im = imageObj.data;
    else
        im = imread(filepath);
    end


function planes = unpackvsi(filepath, args)
	disp(['Loading ', filepath, '.....'])
	r = bfGetReader(filepath);
	np = r.getImageCount();
	if args.split_channels
	    disp('Loading Channels seperately')
		planes = cell([np, 1]);
		for i = 1:np
			planes{i} = mat2gray(bfGetPlane(r, i));
		end
	else
	    disp('Loading as RGB image')
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