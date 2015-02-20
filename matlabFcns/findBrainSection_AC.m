function section = findBrainSection_AC(im, mask)
	%% function section = findBrainSection_AC(im, [mask])
	% 
	% Attempts to locate the brain section in 2D image 'im'
	%  using active contour segmentation 
	%
	% A binary mask, 'mask' can be specified to use as a starting point for the 
	%  segmentation process.  Otherwise, the boundry of the image will be used.

	if ~exist('mask', 'var')
		mask = true(size(im));
	end

	% Segment with active contour
	disp('Determining brain section.....')
	section = activecontour(im, mask, 2000);

	% Clean up the resulting image
	disp('Cleaning up brain section.....')
	% Morphological Processing
	section = imopen(section, strel('disk', 5));
	section = imdilate(section, strel('disk', 3));

	% Remove objects on the edges of the mask, excluding the largest segment (which is probably the brain)
	props = regionprops(section);
	% Only attempt to remove border sections if more than one region is present
	if length(props) > 1
		% Determine which region is the largest (hencforth the 'brain')
		largestRegionLabel = find([props(:).Area] == max([props(:).Area]), 1, 'first');
		% Compute a labeled image, the intensity of the pixels in this image correspond to the index in props storing information about that region
		labeledImage = bwlabel(section);
		% Create an image only containing the brain
		brain = section;
		brain(labeledImage ~= largestRegionLabel) = 0;
		% Union the brain mask with a mask not containing any regions along the edge
		section = brain | imclearborder(section);
	end

