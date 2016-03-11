function p_usable = display_registration_quality_distribution(metadata_path)

	data = json.read(metadata_path);

	data = filter_data(data);

	% Calculate the probablily that a given region in an image is usable, for each image
	p_usable = zeros(size(data));
	fully_usable = ones(size(data));
	region_indicies = zeros(size(data));

	for i = 1:length(data)
		curr_set = data{i};
		region_indicies(i) = curr_set.atlasIndex;

		if strcmpi(curr_set.sliceUsable, 'no')
			p_usable(i) = 0;
		else
			p_usable(i) = get_region_usable_probability(curr_set);
		end
	end

	figure
	stem(region_indicies, fully_usable, 'r', 'linewidth', 1.5)
	hold on
	stem(region_indicies, p_usable, 'color', [.0, .5, .0], 'linewidth', 1.5)
	title('Probability that a Region is Usable in a Given Image')
	xlabel('(<-- Front of Brain)     Atlas Index     (Back of Brain -->)')
	ylabel('P(usable)')

function filtered_data = filter_data(input_data)
	filtered_data = cell(size(input_data));
	f_ind = 1;

	for i = 1:length(input_data)
		% The json library used here may load json object lists as struct arrays or cell arrays.  
		% We only want cells
		if iscell(input_data)
			curr_data = input_data{i};
		else
			curr_data = input_data(i);
		end

		% only include data sets that aren't excluded
		if ~curr_data.exclude
			filtered_data{f_ind} = curr_data;
			f_ind = f_ind + 1;
		end
	end

	% Clear extra cells
	if f_ind <= i
		% f_ind is always 1 higher than the last index written to
		filtered_data(f_ind : end) = [];
		disp([num2str( f_ind - i ), ' images were excluded.'])
	end


function p_usable = get_region_usable_probability(set)
	[~, labels_name, labels_extention] = fileparts(set.registeredAtlasLabelsPath);
	labels_path = [labels_name, labels_extention];
	labels = read_image(labels_path);
	unique_regions = length(unique(labels(:))) * 2;
	number_of_exclusions = double(size(set.regionIdsToExclude, 1));
	p_usable = 1 - number_of_exclusions/unique_regions;


function im = read_image(filepath)
    % Loads an image using the appropriate load function based on the filetype inferred from filepath
    [~, ~, ext] = fileparts(filepath);
    if strcmpi(ext, '.mhd')
        imageObj = read_mhd(filepath);
        im = imageObj.data;
    else
        im = imread(filepath);
    end