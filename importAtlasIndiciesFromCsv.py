#! /usr/bin/python

from pythonMods import importTools
from sys import argv

argParser = importTools.ArgumentParser(argv)
# Get input paths
jsonPath = argParser.generateJsonPath()
csvPath = argParser.getCsvPath()

importer = importTools.CoordinateImporter(csvPath=csvPath, jsonPath=jsonPath)
importer.setMode_index()
importer.importFromCsv(csvHeaderLines=2)