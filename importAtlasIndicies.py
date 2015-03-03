#!/usr/bin/python

from pythonMods import importTools
import argparse
import os


def generate_parser():
    parser = argparse.ArgumentParser(
        description='Imports allen index data from a specified location.  '
                    'Currently accepts .csv with rows beginning with a distinct '
                    'identification tag specifying the .vsi file the index information corresponds to, '
                    'and at least one column containing index information'
    )
    parser.add_argument('csvPath', help='The csv file to import data from.')
    parser.add_argument('-e', '--experimentPath', default=os.getcwd(),
                        help='The root directory containing the data to be processed.'
                             '  Defaults to the current directory.')

    import_modes = parser.add_mutually_exclusive_group(required=True)
    import_modes.add_argument('-i', '--useIndex', default=False, action='store_true',
                              help='Imported values are allen brain atlas slice indicies.')
    import_modes.add_argument('-c', '--useAtlasCoordinates', default=False, action='store_true',
                              help='Imported values correspond to physical rostral caudal coordinates in the aba.')
    import_modes.add_argument('-b', '--useBregmaCoordinates', default=False, action='store_true',
                              help='Imported values correspond to the paxinos bregma coordinates')
    return parser


def set_importer_mode(importer, args):
    if args.useIndex:
        importer.set_mode_index()
    elif args.useAtlasCoordinates:
        importer.set_mode_atlas()
    elif args.useBregmaCoordinates:
        importer.set_mode_bregma()


def main():
    parser = generate_parser()
    args = parser.parse_args()

    csv_file = open(args.csvPath, 'r')

    importer = importTools.CoordinateImporter(args.experimentPath)
    set_importer_mode(importer, args)
    importer.import_csv_data(csv_file)


if __name__ == '__main__':
    main()