#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
from rdkit import Chem
import biochem_tools

__name__ = 'virtsc_lassi.features'
__author__ = 'jiri1kolar'
#load logger
log=logging.getLogger(__name__)




def fragments_extraction(file_path,output_path,fragment_types):
    extraction_options = {
        'kekule': False,
        'isomeric': False,
        'fragments': fragment_types
    }
    file_type = 'smi'
    if file_path.endswith('.sdf'):
        file_type = 'sdf'
    frag_json = output_path
    if output_path == "":
        frag_json = file_path[:-4] + '.frags.json'
    biochem_tools.extract_fragments([file_path], file_type, frag_json, extraction_options)
    return frag_json


def descriptors_extraction(frag_json, output_path, descriptors_generator, padel_path):
    desc_csv=output_path
    if output_path=="":
        desc_csv=frag_json[:-5] + ".csv"
    if descriptors_generator == "rdkit":
        biochem_tools.rdkit_compute_descriptors(frag_json, desc_csv, True)
    elif descriptors_generator == "padel":
        biochem_tools.padel_compute_descriptors(frag_json, desc_csv, True, padel_path)
    return desc_csv

def get_fragments_from_molecule(molecule,fragments_option):
    extraction_options = {
        'kekule': False,
        'isomeric': False,
        'fragments': fragments_option
    }
    return biochem_tools.extract_fragments_from_molecule(molecule, fragments_option, extraction_options)

def get_molecule_from_frag_smiles(str_smiles):
    sanitize_operation = Chem.SanitizeFlags.SANITIZE_ALL ^ \
                         Chem.SanitizeFlags.SANITIZE_KEKULIZE
    mol=Chem.MolFromSmiles(str_smiles, sanitize=False)
    Chem.SanitizeMol(mol, sanitizeOps=sanitize_operation)
    return mol

def generate_feature_desc(file_path,output_path,screen_option):
    frag_json=fragments_extraction(file_path,"",screen_option['fragments'])
    desc_csv=descriptors_extraction(frag_json,output_path,screen_option['descriptors'],screen_option['padel-path'])
    return desc_csv