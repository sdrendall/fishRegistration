function cropRegionsUsingAtlas(vsiPath, labelsPath, hemispherePath, ids, varargin)
	% cropRegionsUsingAtlas(vsiPath, labelsPath, hemispherePath, ids, [options])
	%
	% 'region name' - string, The name of the region being cropped. Used in the filename of output images
	% 'slice order' - 1 x n array containing the permutation of the indicies 1:n 
	%  		correpsonding to the order in which planes in the vsi files should be loaded
	%       defaults to [3 2 1]
	% 'split channels' - boolean, output channels individually, or as an rgb image
	%					default false
	% 'experiment path' - string
	% 'hemisphere' - 'left', 'right', or 'both', Designates the hemisphere images should
	%               be cropped from.  Defaults to 'both'
	% 'exclusions' - N x 2 array, Designates N regions that should be excluded from crop operations.
	%               First column corresponds to the region id, 2nd column corresponds to the hemisphere id
	% 'flip' - Corrects vertically flipped .vsi files
	% 'flop' - Corrects horizontally flipped .vsi files


	%% Parse input
	args = parseVarargin(varargin, nargin - 4);

	%% Crop Each Image
	try
		%% Load Image Data
		annotations = readImage(fullfile(args.exp_path, labelsPath));
		hemisphereLabels = readImage(fullfile(args.exp_path, hemispherePath));

		%% Generate Region Mask
		mask = generateRegionMask(annotations, hemisphereLabels, ids, args);

		%% Only continue if there are actually pixels to crop
		if any(mask(:))
			planes = unpackvsi(fullfile(args.exp_path, vsiPath), args);
			planes = permutePlanes(planes, args);

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
		'region_name', '', ...
		'hemisphere', 'both', ...
		'exclusions', [], ...
		'flip', false, ...
		'flop', false);

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
		case 'hemisphere'
		    args.hemisphere = lower(argv{i + 1});
		    i = i + 2;
		case 'exclusions'
		    args.exclusions = argv{i + 1};
		    i = i + 2;
		case 'flip'
			args.flip = true;
			i = i + 1;
		case 'flop'
			args.flop = true;
			i = i + 1;
		otherwise
		    disp('Argument unrecognized')
		    disp(argv{i})
		    i = i + 1;
		end
	end


function planes = permutePlanes(planes, args)
	for i = 1:length(planes)
		planes{i} = permuteImage(planes{i}, args);
	end


function im = permuteImage(im, args)
    if args.flip && args.flop
        im = flip(im);
        im = flip(im, 2);
    elseif args.flip
        im = flip(im);
    elseif args.flop
        im = flip(im, 2);
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
    [~, ~] = mkdir(saveDir);  % Nullifying the 2nd output prevents Warning messages from being printed if the dir exists


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
		% Normalize values in the cropped region
		croppedRegion(regions(i).Image) = mat2gray(croppedRegion(regions(i).Image));
		croppedRegions{i} = croppedRegion;
	end


function mask = generateRegionMask(annotations, hemisphereLabels, ids, args)
	annoMask = generateAnnotationsMask(annotations, ids);
	hemisphereMask = generateHemisphereMask(hemisphereLabels, args.hemisphere);
	exclusionMask = generateExclusionMask(annotations, hemisphereLabels, args.exclusions);

    mask = annoMask & hemisphereMask & ~exclusionMask;


function mask = generateAnnotationsMask(annotations, ids)
    % Generates a boolean mask of the registered annotations that cooresponds to ids in the list of ids to be cropped
	disp('Identifying Region Pixels.....')
	mask = false(size(annotations));
	for i = 1:length(ids)
		mask = mask | (annotations == ids(i));
	end


function mask = generateHemisphereMask(hemiLabels, hemisToInclude)
    % Generates a boolean mask from the registered hemisphere lables corresponding to the hemispheres listed in
    %   hemisToInclude
    switch hemisToInclude
    % Left is 2, right is 1
    case 'left'
        mask = hemiLabels == 2;
    case 'right'
        mask = hemiLabels == 1;
    case 'both'
        mask = hemiLabels ~= 0; % The only values in hemi labels are 0, 1 and 2
    end


function mask = generateExclusionMask(annotations, hemiLabels, exclusionList)
    % Generates a boolean mask, where true corresponds to pixels that should be excluded from cropping operations.
    mask = false(size(annotations));
    for i = 1:size(exclusionList, 1)
        mask = mask | (annotations == exclusionList(i, 1) & hemiLabels == exclusionList(i, 2));
    end


function im = readImage(filepath)
    % Loads an image using the appropriate load function based on the filetype inferred from filepath
    [~, ~, ext] = fileparts(filepath);
    if strcmpi(ext, '.mhd')
        imageObj = read_mhd(filepath);
        im = imageObj.data;
    else
        im = imread(filepath);
    end


function planes = unpackvsi(filepath, args)
    % Loads each plane from a vsi image into a cell array
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
    % Loads a vsi image as an RGB image.  Maps vsi planes to RGB channels based on the slice order argument
	nr = reader.getSizeY;
	nc = reader.getSizeX;
	np = reader.getImageCount;
	im = zeros([nr, nc, 3]);

	for i = 1:np
		im(:,:,args.slice_order(i)) = double(bfGetPlane(reader, i));
	end

	im(:,:,3) = mat2gray(im(:,:,3));
	im(:,:,2) = im(:,:,2)*20/2^16;