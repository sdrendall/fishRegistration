import os, json
from sys import argv

## ----- CONVERSION FUNCTIONS ----- ##
bregmaInAtlas = 5525 # in um


def atlas2Index(aCoord):
    ind = round(aCoord/25)
    return verifyAbaIndex(ind)


def atlas2Bregma(aCoord):
    bCoord = bregmaInAtlas - aCoord
    return um2mm(bCoord)


def bregma2Atlas(bCoord):
    bCoord = mm2um(bCoord)
    return bregmaInAtlas - bCoord


def index2Atlas(aInd):
    return aInd*25


def verifyAbaIndex(ind):
    if ind < 0:
        print "WARNING: Given atlas coordinate corresponds to a negative index!"
        print "Defaulting to 0 (first slice)..."
        ind = 0
    elif ind > 528:
        print "WARNING: Given atlas coordinate corresponds to an index outside of the Allen Brain Atlas space!"
        print "Defaulting to 528 (last slice)..."
        ind = 528
    return ind


def getAvgCoord(self, c1, c2):
    cAv = mean(c1, c2)
    return round(cAv)


def mm2um(mm):
    return mm*1000


def um2mm(um):
    return um/1000


def mean(a, b):
    a = float(a)
    b = float(b)
    return (a + b)/2


class ArgumentParser:

    def __init__(self, argv=argv):
        self.argv = argv
        self.checkArgs()
        self.setExperimentPath()

    def checkArgs(self):
        if len(self.argv) > 3:
            print "Too many input arguments!"
            self.printProperUsage()
            raise Exception('Inproper Input')
        elif len(self.argv) == 3 and not os.path.isdir(os.path.abspath(self.argv[2])):
            print "Specified experimentPath does not exist!"
            self.printProperUsage()
            raise Exception('Inproper Input')
        elif len(self.argv) == 2 and not os.path.exists(os.path.abspath(self.argv[1])):
            print "Specified csvPath does not exist!"
            self.printProperUsage()
            raise Exception('Inproper Input')
        elif len(self.argv) == 1:
            print "Not enough input arguments!"
            self.printProperUsage()
            raise Exception('Inproper Input')

    def printProperUsage(self):
        print "Proper Usage:"
        print "%s csvPath [experimentPath]" % self.argv[0]

    def getExperimentPath(self):
        if len(self.argv) == 3:
            expPath = os.path.abspath(self.argv[1])
        else:
            expPath = os.path.abspath('.')
        return expPath

    def setExperimentPath(self, expPath=None):
        if expPath is None: expPath = self.getExperimentPath() # Won't let me use self.method() for default values
        self.checkExperimentPath(expPath)
        print 'Experiment Path: ' + expPath
        self.expPath = expPath

    def getCsvPath(self):
        return os.path.abspath(self.argv[1])

    def generateJsonPath(self, expPath=None):
        if expPath is None: expPath = self.expPath
        return os.path.join(expPath, '.registrationData', 'metadata.json')

    def checkExperimentPath(self, expPath):
        jsonPath = self.generateJsonPath(expPath)
        if not os.path.isfile(jsonPath):
            print "%s not found!" % jsonPath
            raise Exception('Invalid Experiment Path:' + expPath)


class CoordinateImporter:

    def __init__(self, csvPath, jsonPath):
        self.jsonPath = jsonPath
        self.setMode_index()
        self.csv = open(csvPath)
        self.dictList = self.getDictList()

    def importFromCsv(self, csvHeaderLines=0):
        # Remove Header Lines from the input csv
        self.popCsvHeaders(csvHeaderLines)

        # By default, exclude entries
        for entry in self.dictList: entry['exclude'] = 1

        # Map entries, and update Json File
        self.mapCsvToJson()
        self.warnAboutDqs()
        self.updateJson()
    
    def parseLine(self, line):
        line = removeWhitespace(line)
        items = line.split(',')
        tag = items.pop(0)
        return tag, items

    def cleanInput(self, items):
        return [float(val) for val in items if val.isdigit()]

    def mapCsvToJson(self):
        for line in self.csv:
            # Construct vsi filename from tag
            tag, inp = self.parseLine(line)
            # Search json for corresponding vsi entry
            entry = self.findCorrespondingJsonEntry(tag)
            self.updateEntry(entry, inp, tag)

    def findCorrespondingJsonEntry(self, tag):
        for entry in self.dictList:
            if tag in entry['vsiPath']:
                return entry
        # Return None if no entry matches
        return None

    def updateEntry(self, entry, inp, tag):
        # Disqualify damaged sections
        isDamaged = lambda inp: inp.lower() == 'damaged'
        for item in inp: 
            if isDamaged(item):
                return

        # The input must be cleaned AFTER it is searched for 'damaged'.
        # Cleaning removes all non-numeric strings, and converts all remaining strings to floats
        inp = self.cleanInput(inp)
        # entry is None if no matching entry json entry is found
        if entry is not None:
            coord = mean(inp)
            self.updateEntryWithCoord(entry, coord)
            entry['exclude'] = 0
        else:
            # If an entry is not found in the json file, warn the user
            print "WARNING: No json entry found for " + tag
            
    def warnAboutDqs(self):
        for entry in self.dictList:
            if entry['exclude'] == 1:
                print "WARNING: " + entry['vsiPath'] + " will be excluded because it was not assigned an atlas coordinate"
    
    def popCsvHeaders(self, numLines):
        for i in xrange(numLines): self.csv.readline()
    
    def getDictList(self):
        jsonFile = open(self.jsonPath)
        dictList = json.load(jsonFile)
        jsonFile.close()
        return dictList
    
    def updateJson(self):
        jsonFile = open(self.jsonPath, 'w')
        json.dump(self.dictList, jsonFile, sort_keys=True, indent=4)
        jsonFile.close()

    def setMode_bregma(self):
        self.mode = 'bregma'
        self.updateEntryWithCoord = self.updateEntryWithBregmaCoord

    def setMode_atlas(self):
        self.mode = 'atlas'
        self.updateEntryWithCoord = self.updateEntryWithAtlasCoord

    def setMode_index(self):
        self.mode = 'index'
        self.updateEntryWithCoord = self.updateEntryWithAtlasIndex

    def updateEntryWithBregmaCoord(self, entry, coord):
        entry['bregmaCoord'] = float(coord)
        entry['atlasCoord'] = bregma2Atlas(entry['bregmaCoord'])
        entry['atlasIndex'] = atlas2Index(entry['atlasCoord'])

    def updateEntryWithAtlasCoord(self, entry, coord):
        entry['atlasCoord'] = float(coord)
        entry['atlasIndex'] = atlas2Index(entry['atlasCoord'])
        entry['bregmaCoord'] = atlas2Bregma(entry['atlasCoord'])
    
    def updateEntryWithAtlasIndex(self, entry, coord):
        entry['atlasIndex'] = float(round(coord))
        entry['atlasCoord'] = index2Atlas(entry['atlasIndex'])
        entry['bregmaCoord'] = atlas2Bregma(entry['atlasCoord'])


### ---- MISC FUNCTIONS ---- ###

def removeWhitespace(line):
    return ''.join(line.split())

def mean(l):
    return sum(l)/float(len(l))