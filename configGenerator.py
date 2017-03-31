#!/usr/bin/python2
# -*- coding: utf-8 -*-

import os
import sys

__name__ = 'virtsc_lassi.configGenerator'
__author__ = 'jiri1kolar'
root_directory = os.path.dirname(os.path.realpath(__file__))

'''
with open(root_directory + '/collection_paths.data', 'r') as input:
    lw=['log','binary','tf']
    gw=['normal','entropy','idf','binary']
    for line in input:
        collection_path=line.split('\n')[0]
        dirs=collection_path.split('/')
        results_path='results/'+dirs[2]+'/'+dirs[3]+'/'+dirs[3]
        config_path='json/'+dirs[2]+'/'+dirs[3]
        if not os.path.exists(config_path):
            os.makedirs(config_path)
        for lwi in lw:
            for gwj in gw:
                with open(config_path+'/config_'+dirs[3]+lwi+gwj+'.json','w') as output:
                    output.write('{"collection_path":"'+collection_path+'",')
                    output.write('"files": {"data": "data.smi", "ligands": "ligands.smi","known-ligands": "known-ligands.smi", "known-decoys": "known-decoys.smi"},"output":"')
                    output.write((results_path+lwi+gwj+'","lw" : "{0:s}","gw" : "{1:s}","morgan" : {2:d}}}').format(lwi,gwj,3))
'''


