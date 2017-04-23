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
settings.init_loggging()
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
    weight_fn_info = ""
    #root_directory = '/home/jiri/mff_vs/virtsc_lassi'
    #for config_file in sys.argv[1:]:
    if args.config_file != "" and args.config_file is not None:
        log.info('Loading config file: %s', args.config_file)
        with open(args.config_file) as input_stream:
            configuration = json.load(input_stream)
            if ('lw' in configuration):
                settings.local_weight_function = configuration['lw']
                weight_fn_info+=" lw: "+configuration['lw']
            if ('gw' in configuration):
                settings.global_weight_function = configuration['gw']
                weight_fn_info += " gw: " + configuration['gw']
    if args.selection_file != "" and args.selection_file is not None:
        with open(root_directory + '/' + args.selection_file) as input_stream:
            configuration = json.load(input_stream)
            io.prepare_selection(configuration, output_dest, root_directory)
    if args.selection_muv != "" and args.selection_muv is not None:
        opt={'nligands-to-train':args.nligands_to_train,
             'nrandom-set':10}
        io.prepare_selection_MUV(args.selection_muv,output_dest,opt)
    if args.directory is not None and args.directory!="":
        recognize_option={'recognize_sets': True,
                          'recognize_collection': True}
        screen_option={'method': args.method,
                       'padel-path': 'padel',
                       'fragments': io.parseFragmentTypes(args.fragment_types),
                       'descriptors': args.descriptors_generator,
                       'cos_th': args.threshold_similarity,
                       'lw': settings.local_weight_function,
                       'gw': settings.global_weight_function}
        screen_info=args.method
        if args.method!='descriptors':
            screen_info+=' '+args.fragment_types
        if args.method=='features':
            screen_info+='th_sim: '+args.threshold_similarity.__str__()
        screen_info+=weight_fn_info
        screen.recursive_screen(args.directory,args.output_directory,recognize_option,screen_option,screen_info)

if __name__ == 'virtsc_lassi':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-s","--selection-file",
                        metavar='FILE',
                        default="",
                        type=str,
                        help="Prepare selection from json file and save it in to output directory.")
    parser.add_argument("--selection-muv",
                        metavar='FILE',
                        default="",
                        type=str,
                        help="Prepare selection from MUV directory and save it to output directory.")
    parser.add_argument("--nligands-to-train",
                        help="Number of ligands choosen from MUV set as known-ligands",
                        type=int,
                        default=20)
    parser.add_argument("-c", "--config-file",
                        metavar='FILE',
                        type=str,
                        help="Json configuration file")
    parser.add_argument("-d", "--directory",
                        metavar='FILE',
                        type=str,
                        default="",
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
                        help="Method of analyzing molecules. Allowed values are:\n"
                             "\t[fragments] - each molecule is prepresented by fragment fingerprint\n"
                             "\t[descriptors] - each molecule is represented by descriptors\n"
                             "\t[features] - each molecule is split to fragments and each fragment is represented by descriptors",
                        choices=['fragments', 'descriptors','features'],
                        default="fragments")
    parser.add_argument("-g", "--descriptors-generator",
                        help="Generator to be used to obtain fragments descriptors.",
                        choices=['rdkit', 'padel'],
                        default="rdkit")
    parser.add_argument("-ths", "--threshold-similarity",
                        help="Threshold for fraction of compared molecules, for features only",
                        type=float,
                        default=-1.0)
    args = parser.parse_args()
    main()
