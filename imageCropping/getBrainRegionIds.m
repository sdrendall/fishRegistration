function ids = getBrainRegionIds(name, jsonPath)
	% ids = getBrainRegionIds(structureName, structureDataPath)
	%
	% Returns the ids corresponding to the named brain structure, and all of its substructures in the allen brain atlas
	%
	% structureName is the name of the brain structure in the allen brain atlas (case sensitive)
	% structureDataPath is the path to the json file containing the allen brain atlas structure data

	%  which can be found here: TODO (aba api call)

	if ~exist('jsonPath', 'var')
		jsonPath = '/home/sam/Dropbox/grayLab/allenReferenceAtlas_mouseCoronal/structureData.json';
	end

	% load json file
	data = loadjson(jsonPath);
	
	% find the stucture corresponding to the given name
	[structure, structureFound] = getStructureByName(name, data);

	% get the ids of the named stucture and all children
	if structureFound
		ids = getIdsFromStructure(structure);
	else
		error(['Could not find structure with name: ', name])
	end


function [structure, structureFound] = getStructureByName(structureName, structure)
	% Checks if the given structure has the given name, then checks that structure's children if it doesn't
	if strcmp(structure.name, structureName)
		structureFound = true;
	else
		if hasChildren(structure)
			[structure, structureFound] = checkChildrenForStructureName(structureName, structure.children);
		else
			structureFound = false;
		end
	end


function [child, childFound] = checkChildrenForStructureName(name, children)
	% Iterates through a given cell array of children, and checks to see if they have the given name
	for i = 1:length(children)
		% check each child and its children
		[child, childFound] = getStructureByName(name, children{i});
		if childFound
			break
		end
	end


function ids = getIdsFromStructure(structure)
	% Returns a list of id numbers corresponding to the given structure, and all of it's children (and children's children....)
	if hasChildren(structure)
		ids = [structure.id, getChildIds(structure.children)];
	else
		ids = structure.id;
	end


function childIds = getChildIds(children)
	% Returns a list of all of the ids identifying a child structure and it's children
	childIds = [];
	for i = 1:length(children)
		childIds = [childIds, getIdsFromStructure(children{i})];
	end


function answer = hasChildren(structure)
	answer = size(structure.children, 1) > 0;