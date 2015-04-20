import json
from itertools import imap, chain


class StructureFinder:
    """
    I find the structure in the allen brain atlas corresponding to the given trait
    """

    def __init__(self, structure_data_path):
        self.structureData = self._load_structure_data(structure_data_path)

    def get_ids_by_acronym(self, acronym):
        return self._get_ids_by_attribute('acronym', acronym)

    def get_ids_by_structure_name(self, name):
        return self._get_ids_by_attribute('name', name)

    def get_ids_by_id(self, id_no):
        return self._get_ids_by_attribute('id', id_no)

    def _get_ids_by_attribute(self, attribute, value):
        """
        'Template' function for obtaining the set of structure ids that identify a structure by an attribute
        :param attribute: Attribute to identify structure by
        :param value: Value of the attribute
        :return: List of ids corresponding to the ids of the identified structure and all of it's progeny
        """
        structure = self.search_structure_data_for_attribute(attribute, value)
        if structure is not None:
            return self.get_ids_from_structure(structure)
        else:
            raise StructureNotFoundError(attribute, value)

    def search_structure_data_for_attribute(self, attribute, value):
        """
        Searches the structure data stored in self.structureData for a structure who's attribute == value
        :param attribute: The attribute to be checked
        :param value: The desired value of that attribute
        :return: The structure contained in structureData who's attribute == value
        """
        return self._check_structure_for_attribute(self.structureData, attribute, value)

    def _check_structure_for_attribute(self, structure, attribute, value):
        """
        Checks to see if a structure has the desired value.  If it doesn't, checks that structure's children
        :param structure: The root structure to be searched
        :param attribute: The attribute to be checked
        :param value: The desired value of the attribute
        :return: The structure who's attribute has the desired value.  None if that structure can't be found
        """
        if structure[attribute] == value:
            return structure
        else:
            return self._search_children_for_attribute(structure['children'], attribute, value)

    def _search_children_for_attribute(self, children, attribute, value):
        """
        Checks if each the specified attribute of each child in children has the desired value
        :param children: The structures to be checked
        :param attribute: The attribute to be checked
        :param value: The desired value
        :return: The structure of the child containing the desired attribute, None if that child isn't found
        """
        for child in children:
            structure = self._check_structure_for_attribute(child, attribute, value)
            if structure is not None:
                return structure
        return None

    def get_ids_from_structure(self, structure):
        """
        :param structure: A structure from the allen brain atlas structure data
        :return: A list of ids corresponding to the given structure and all of its descendants
        """
        return [id_no for id_no in self._generate_structure_ids(structure)]

    def _generate_structure_ids(self, structure):
        """
        A generator function that yields the id numbers of a given structure from the allen brain atlas,
            and all of its descendants
        :param structure: A structure from the allen brain atlas structure data
        """
        yield structure['id']
        for id_no in chain(*imap(self._generate_structure_ids, structure['children'])):
            yield id_no

    @staticmethod
    def get_structure_property_generator_function(structure_property):

        def structure_property_generator_function(structure):
            yield structure[structure_property]
            for property_value in chain(*imap(structure_property_generator_function, structure['children'])):
                yield property_value

        return structure_property_generator_function

    @staticmethod
    def _load_structure_data(path):
        f = open(path, 'r')
        data = json.load(f)
        f.close()
        return data


class StructureNotFoundError(Exception):
    """
    Raised when a structure cannot be found in the allen brain atlas
    """

    def __init__(self, identifier_type, identifier):
        self.msg = 'Structure with %s %s could not be found.' % (identifier_type, identifier)


def main():
    structure_data_path = '/home/sam/Dropbox/grayLab/allenReferenceAtlas_mouseCoronal/structureData.json'
    finder = StructureFinder(structure_data_path)

    print "Finding ids by acronym....."
    tea_ids = finder.get_ids_by_acronym('TEa')
    print tea_ids

    print "Finding ids by name....."
    iso_ids = finder.get_ids_by_structure_name('Isocortex')
    print iso_ids

    if all((id in iso_ids for id in tea_ids)):
        print 'success!'
    else:
        print 'failure :('

if __name__ == '__main__':
    main()