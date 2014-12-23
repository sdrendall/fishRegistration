#! /usr/bin/python

import os, json, sys
from sys import argv

def checkArgs():
    if len(argv) > 3:
        print "Too many input arguments!"
        printProperUsage()
        return False
    elif len(argv) == 3 and not os.path.isdir(os.path.abspath(argv[2])):
        print "Specified experimentPath does not exist!"
        printProperUsage()
        return False
    elif len(argv) == 2 and not os.path.exists(os.path.abspath(argv[1])):
        print "Specified csvPath does not exist!"
        printProperUsage()
        return False    
    elif len(argv) == 1:
        print "Not enough input arguments!"
        printProperUsage()
        return False
    return True

def printProperUsage():
    print "Proper Usage:"
    print "%s csvPath [experimentPath]" % argv[0]

def getExperimentPath():
    if len(argv) == 3:
        return os.path.abspath(argv[1])
    else:
        return os.path.abspath('.')

def generateJsonPath(expPath):
    return os.path.join(expPath, '.registrationData' ,'metadata.json')


def isInvalidExperimentPath(expPath):
    jsonPath = generateJsonPath(expPath)
    if not os.path.isfile(jsonPath):
        print "%s not found!" % jsonPath
        return True
    else:
        return False

def loadJson(jsonPath):
    jsonFile = open(jsonPath)
    dictList = json.load(jsonFile)
    jsonFile.close()
    return dictList

def updateJson(jsonPath, dictList):
    jsonFile = open(jsonPath, 'w')
    json.dump(dictList, jsonFile, sort_keys=True, indent=4)
    jsonFile.close()

def getCsvPath():
    return os.path.abspath(argv[1])

def findCorrespondingJsonEntry(tag, dictList):
    for entry in dictList:
        if tag in entry['vsiPath']:
            return entry
    # Return None if no entry matches
    return None

def updateBregmaCoord(entry, coord):
    entry['bregmaCoord'] = float(coord)
    entry['atlasCoord'] = bregmaToAtlas(coord)

def bregmaToAtlas(bCoord):
    bregmaInAtlas = 5525 # in um
    bCoord = mm2um(float(bCoord))
    return bregmaInAtlas - bCoord

def mm2um(mm):
    return mm*1000

def warnAboutDqs(dictList):
    for entry in dictList:
        if entry['exclude'] == 1:
            print "WARNING: " + entry['vsiPath'] + " will be excluded because it was not assigned an atlas coordinate"

def mapCsvToJson(csv, dictList):
    for line in csv:
        # Construct vsi filename from tag
        tag, coord = line.split(',')
        # Search json for corresponding vsi entry
        entry = findCorrespondingJsonEntry(tag, dictList)

        # findCorrespondingJsonEntry will return None if no entry is found
        if entry is not None:
            updateBregmaCoord(entry, coord)
            entry['exclude'] = 0
        else:
            # If an entry is not found in the json file, warn the user
            print "WARNING: No json entry found for " + tag

def main():
    # Check args
    argsValid = checkArgs()
    if not argsValid:
        return
    # Load experiment path
    expPath = getExperimentPath()
    if isInvalidExperimentPath(expPath):
        print "Invalid Experiment Path:"
        print expPath
        return

    print "expPath: ", expPath

    # Load JSON
    jsonPath = generateJsonPath(expPath)
    dictList = loadJson(jsonPath)

    # Load csv
    csvPath = getCsvPath()
    csv = open(csvPath)

    # Pop First Line
    csv.readline()

    # Json entries will be excluded by default, and included if a bregma coordinate is found
    for entry in dictList: entry['exclude'] = 1
        
    # Map csv entries to json entries
    mapCsvToJson(csv, dictList)

    # Warn about excluded entries
    warnAboutDqs(dictList)

    # Update json
    updateJson(jsonPath, dictList)


if __name__ == '__main__':
    main()