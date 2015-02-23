import os
import jsonTools
from itertools import imap


class CoordinateImporter(jsonTools.MetadataHandler):
    """
    Subclass of the MetadataHandler that knows how to interpret the various forms that allen brain atlas coordinates may
     be given as.

    Slices marked as 'damaged' or without an input coordinate are marked to be excluded
    """
    def __init__(self, experiment_path=os.getcwd()):
        jsonTools.MetadataHandler.__init__(self, experiment_path)
        # Initialize mode
        self.mode = 'index'
        self.update_entry_with_coord = self.update_entry_with_atlas_index
        # Load metadata
        self.load_metadata()
        # Exclude each entry by default
        for entry in self.metadata:
            entry['exclude'] = True

    def import_csv_data(self, csv_lines):
        # Map entries, and update Json File
        formatted_data = self.format_csv_data(csv_lines)
        self.add_input_to_metadata(formatted_data)

    def format_csv_data(self, csv_lines):
        """
        :param csv_lines - an iterable that returns lines in csv format, where the first item is a tag associated with
            a .vsi file in the experiment, and the rest of the items describe the rostral-caudal position of the slice
            in the brain, or the word 'damaged'
        :returns A list of InputData objects, each with the tag set to the tag found in the csv file, and appropriate
            values for exclude and coord depending on the contents of the corresponding csv line
        """
        # Some iterator stuff because I was bored
        iparsed_lines = imap(self.parse_csv_line, csv_lines)
        exclude_item = lambda inp: inp.lower() == 'damaged' or inp.lower() == 'exclude'
        return [InputData(tag=tag, exclude=any(imap(exclude_item, items)), coord=mean(self.iclean_input(items)))
                for tag, items in iparsed_lines]

    @staticmethod
    def parse_csv_line(line):
        line = remove_whitespace(line)
        items = line.split(',')
        tag = items.pop(0)
        return tag, items

    @staticmethod
    def iclean_input(items):
        for val in items:
            if is_number(val):
                yield float(val)

    def add_input_to_metadata(self, input_data):
        # map tags to metadata
        for item in input_data:
            entry = self.find_corresponding_metadata_entry(item.tag)
            if entry is not None and item.data['coord'] is not None:
                entry['exclude'] = item.data['exclude']
                self.update_entry_with_coord(entry, item.data['coord'])

        # warn about dqs
        self.warn_about_dqs()
        # update metadata.json
        self.update_metadata()

    def find_corresponding_metadata_entry(self, tag):
        for entry in self.metadata:
            if tag.lower() in entry['vsiPath'].lower():
                return entry

        # If an entry can't be found, print a warning and return none
        print "WARNING: No json entry found for " + tag
        return None

    def warn_about_dqs(self):
        for entry in self.metadata:
            if entry['exclude']:
                print "WARNING: " + entry['vsiPath'] + " will be excluded because it was not assigned an atlas coordinate"

    def set_mode_bregma(self):
        self.mode = 'bregma'
        self.update_entry_with_coord = self.update_entry_with_bregma_coord

    def set_mode_atlas(self):
        self.mode = 'atlas'
        self.update_entry_with_coord = self.update_entry_with_atlas_coord

    def set_mode_index(self):
        self.mode = 'index'
        self.update_entry_with_coord = self.update_entry_with_atlas_index

    @staticmethod
    def update_entry_with_bregma_coord(entry, coord):
        entry['bregmaCoord'] = float(coord)
        entry['atlasCoord'] = bregma_to_atlas(entry['bregmaCoord'])
        entry['atlasIndex'] = atlas_to_index(entry['atlasCoord'])

    @staticmethod
    def update_entry_with_atlas_coord(entry, coord):
        entry['atlasCoord'] = float(coord)
        entry['atlasIndex'] = atlas_to_index(entry['atlasCoord'])
        entry['bregmaCoord'] = atlas_to_bregma(entry['atlasCoord'])

    @staticmethod
    def update_entry_with_atlas_index(entry, coord):
        entry['atlasIndex'] = float(round(coord))
        entry['atlasCoord'] = index_to_atlas(entry['atlasIndex'])
        entry['bregmaCoord'] = atlas_to_bregma(entry['atlasCoord'])


class InputData:
    """
    I hold data to add to one entry of the metadata.json file.
    My tag identifies the metadata entry I should be added to.
    """

    # There's probably a better way to make these defaults but I can't figure it out.
    #  The empty dict will never be accessed, so it doesn't matter that it is mutable
    def __init__(self, tag, iterable={}, **kwargs):
        self.tag = tag
        self.data = dict(iterable, **kwargs)

## ----- CONVERSION FUNCTIONS ----- ##
bregmaInAtlas = 5525  # in um


def atlas_to_index(aCoord):
    ind = round(aCoord/25)
    return verify_aba_index(ind)


def atlas_to_bregma(aCoord):
    bCoord = bregmaInAtlas - aCoord
    return um2mm(bCoord)


def bregma_to_atlas(bCoord):
    bCoord = mm2um(bCoord)
    return bregmaInAtlas - bCoord


def index_to_atlas(aInd):
    return aInd*25


def verify_aba_index(ind):
    if ind < 0:
        print "WARNING: Given atlas coordinate corresponds to a negative index!"
        print "Defaulting to 0 (first slice)..."
        ind = 0
    elif ind > 528:
        print "WARNING: Given atlas coordinate corresponds to an index outside of the Allen Brain Atlas space!"
        print "Defaulting to 528 (last slice)..."
        ind = 528
    return ind


def mm2um(mm):
    return mm*1000


def um2mm(um):
    return um/1000


def mean(iterable):
    vals = map(float, iterable)
    n = len(vals)
    if n == 0:
        return None
    else:
        return sum(vals)/n


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


### ---- MISC FUNCTIONS ---- ###

def remove_whitespace(line):
    return ''.join(line.split())