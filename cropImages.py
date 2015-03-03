#! /usr/bin/python
__author__ = 'Sam Rendall'

import argparse
import os
import subprocess
from itertools import imap
from pythonMods import outputProcessing


def generate_parser():
    parser = argparse.ArgumentParser(
        description='A tool to crop brain regions from images registered to the allen brain atlas')
    parser.add_argument('-e', '--experimentPath', default=os.getcwd(),
                        help='Specifies the root directory containing the data to be processed.  '
                             'Defaults to the current directory.')
    parser.add_argument('-o', '--outputPath', help='The location to save the output to.', default=None)
    parser.add_argument('-r', '--regions', nargs='+',
                        help="A list of names or acronyms denoting the regions to be cropped. "
                             "Hint: Region names containing spaces should be wrapped in apostrophes 'like' 'this'")
    parser.add_argument('-s', '--splitChannels', action='store_const', const=1, default=0,
                        help='Save cropped images from each channel in the .vsi image as individual images')
    parser.add_argument('-i', '--sliceOrder', nargs=3, default=[1, 2, 3], type=int,
                        help='If channels are not split, this specifies the order that channels in the vsi Image '
                             'should be saved to the output rgb image.')
    parser.add_argument('-i', '--structureDataPath',
                        default=os.path.abspath('~/code/fishRegistration/structureData.json'),
                        help='The path to the .json file containing the allen brain atlas structure data. '
                             'Defaults to ~/code/fishRegistration/structureData.json')
    parser.add_argument('-b', '--useBatch', action='store_true', default=False,
                        help='Submit processes to an lsf cluster. Warning: Mistakes will take longer to catch, '
                             'prototype locally before using this option.')

    group = parser.add_mutually_exclusive_group(required=True,
                                                description='The type of identifier used to specify brain structures.')
    group.add_argument('-a', '--useAcronyms', action='store_true', default=False,
                       help='Identify brain regions using acronyms')
    group.add_argument('-n', '--useNames', action='store_true', default=False,
                       help='Identify brain regions using names')
    return parser


def generate_arg_string(ids, args):
    arg_string = "matlab -nosplash -nodesktop -r \"cropRegionsUsingAtlas([%s], 'slice order', {sliceOrder}, " \
                 "'split_channels', {splitChannels}, 'experiment path', {experimentPath}, " \
                 "'output path', {outputPath}, ); exit\""
    arg_string = arg_string   # Cast the list of structure ids to strings, add to argstr
    arg_string = arg_string.format(**vars(args))  # Unpack the args and insert them into the arg_string

    if args.useBatch:
        return 'bsub -q short -W 0:30 -R "rusage[mem=4000]" ' + arg_string
    else:
        return arg_string


def open_crop_process(arg_string, args):
    subprocess.call(arg_string, cwd=args.experimentPath, env=os.environ, shell=True)


def format_output_path(args):
    """
    Ensures that the output path in the args namespace is always relative to the experiment path
    This might be better handled by the argparse parser, but I'm not sure how.
    I'm also not sure how this method will handle symbolic links
    :param args: the namespace returned by the the parse_args method of an ArgumentParser
    :return: args, with the outputPath attribute made relative to the experimentPath attribute
    """
    if args.outputPath is None:
        args.outputPath = '.'
    else:
        args.outputPath = os.path.relpath(args.outputPath, args.experimentPath)

    return args


def main():
    parser = generate_parser()
    args = parser.parse_args()
    args = format_output_path(args)

    finder = outputProcessing.StructureFinder(args.structureDataPath)
    if args.useAcronyms:
        id_list_gen = (finder.get_ids_by_acronym(acronym) for acronym in args.regions)
    elif args.useNames:
        id_list_gen = (finder.get_ids_by_structure_name(name) for name in args.regions)
    else:
        print 'Input format must be specified'
        raise Exception()

    for id_list in id_list_gen:
        arg_string = generate_arg_string(id_list, args)
        open_crop_process(arg_string, args)