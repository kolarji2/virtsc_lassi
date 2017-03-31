#!/usr/bin/python2
# -*- coding: utf-8 -*-

import os
import sys
import time
import logging
import json
import analyzator
import settings
import io
import screen
import argparse

__name__ = 'virtsc_lassi'
__author__ = 'jiri1kolar'
__email__ = "jiri1kolar@gmail.com"

#logger setting
#create file handler
logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s [%(levelname)s] %(module)s - %(message)s',
        datefmt='%H:%M:%S')
log = logging.getLogger(__name__)

# create file handler which logs even debug messages
#fh=logging.FileHandler(__name__+time.clock().__str__()+'.log')
#fh.setLevel(logging.DEBUG)
## create console handler

# add the handlers to the logger
#log.addHandler(fh)
#log.addHandler(ch)

def main():
    root_directory = os.path.dirname(os.path.realpath(__file__))
    output_dest = 'results'
    #root_directory = '/home/jiri/mff_vs/virtsc_lassi'
    #for config_file in sys.argv[1:]:
    if args.config_file != "" and args.config_file is not None:
        log.info('Loading config file: %s', args.config_file)
        with open(root_directory + '/' + args.config_file) as input_stream:
            configuration = json.load(input_stream)
            if ('lw' in configuration):
                settings.local_weight_function = configuration['lw']
            if ('gw' in configuration):
                settings.global_weight_function = configuration['gw']
            if ('morgan' in configuration):
                settings.morgan_radius = configuration['morgan']
    if args.selection_file != "" and args.selection_file is not None:
        with open(root_directory + '/' + args.selection_file) as input_stream:
            configuration = json.load(input_stream)
            io.prepare_selection(configuration, output_dest, root_directory)
    if args.directory is not None:
        recognize_option={'recognize_sets': True,
                          'recognize_collection': True}
        screen_option={'method':args.method,
                       'padel-path':'padel',
                       'fragments':io.parseFragmentTypes(args.fragment_types),
                       'descriptors': args.descriptors_generator}
        screen_info=args.method
        if args.method!='descriptors':
            screen_info+=' '+args.fragment_types
        screen.recursive_screen(args.directory,args.output_directory,recognize_option,screen_option,screen_info)

if __name__ == 'virtsc_lassi':
    parser = argparse.ArgumentParser()
    parser.add_argument("-s","--selection-file",
                        metavar='FILE',
                        default="",
                        type=str,
                        help="Prepare selection from json file and save it in standard format to output directory.")
    parser.add_argument("-c", "--config-file",
                        metavar='FILE',
                        type=str,
                        help="Json configuration file")
    parser.add_argument("-d", "--directory",
                        metavar='FILE',
                        type=str,
                        help="Directory where program should look for molecule sets.")
    parser.add_argument("-od", "--output-directory",
                        metavar='FILE',
                        type=str,
                        default="",
                        help="Output directory for results, by default results are saved in same directory as data")
    parser.add_argument("-f", "--fragment-types",
                        help="Optional comma separated list of fragment types to extract (eg. \"tt.3,ecfp.2\")",
                        default="ecfp.2")
    parser.add_argument("-m", "--method",
                        help="Method of analyzing molecules. Allowed values are [fragments|descriptors|features]",
                        choices=['fragments', 'descriptors','features','pharm2d'],
                        default="fragments")
    parser.add_argument("-g", "--descriptors-generator",
                        help="Generator to be used to obtain fragments descriptors. Allowed values are [rdkit|padel]",
                        choices=['rdkit', 'padel'],
                        default="rdkit")
    args = parser.parse_args()
    main()
