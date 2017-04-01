# !/usr/bin/python2
# -*- coding: utf-8 -*-

from rdkit.Chem.AtomPairs import Pairs
from rdkit.Chem.AtomPairs import Torsions
from rdkit.Chem import AllChem
from rdkit import Chem
from rdkit.DataStructs import IntSparseIntVect as isvec
from rdkit.DataStructs import UIntSparseIntVect as uisvec
from rdkit.DataStructs import LongSparseIntVect as lsvec
from rdkit.Chem.Pharm2D import Gobbi_Pharm2D,Generate


import io
import features
import os
import subprocess

#local
import logging
import numpy as np
from biochem_tools.rdkit_descriptors import _names,_functions

__name__ = 'virtsc_lassi.descriptors'
__author__ = 'jiri1kolar'
#load logger
log=logging.getLogger(__name__)
th=1e-9



'''
Descriptor section
'''
def get_desc(molecule,screen_option):
    if screen_option['method']=='fragments':
        return get_fragments(molecule,screen_option['fragments'])
    elif screen_option['method']=='descriptors':
        return get_descriptors(molecule,screen_option['descriptors'])
    elif screen_option['method'] == 'pharm2d':
        return get_pharmacofores(molecule)
    elif screen_option['method'] == 'features':
        return get_descriptors(molecule,screen_option['descriptors'])
    else:
        log.error("Can not generate descriptors, unknown option!")
    return None

def get_descriptors(molecule,type):
    if type=='rdkit':
        return get_RDkit_Desc(molecule)
    elif type=='padel':
        return _desc_lib.get_desc(molecule)
    return None


def get_RDkit_Desc(molecule):
    desc = {}
    for name, fn in zip(_names, _functions):
        val = fn(molecule)
        if np.isnan(val) or abs(val) < th:
            continue
        desc[name] = val
    return desc

def get_fragments(molecule,types):
    desc = {}
    for item in types:
        if item['name'] == 'tt':
            desc.update(get_Torsion(molecule))
        elif item['name'] == 'ap':
            desc.update(get_AtomPairs(molecule))
        elif item['name'] == 'ecfp':
            desc.update(get_Morgan(molecule,item['size']))
    return desc


def get_AtomPairs(molecule):
    '''
    Compute atom pair descriptors(APs)
    :param molecule formula in smiles format:
    :return: sparse int array of frequency of APs
    '''
    ap=Pairs.GetAtomPairFingerprint(molecule)
    return isvec.GetNonzeroElements(ap)

def get_Torsion(molecule):
    '''
    Compute atom pair descriptors(APs)
    :param molecule formula in smiles format:
    :return: sparse int array of frequency of APs
    '''
    tt=Torsions.GetTopologicalTorsionFingerprint(molecule)
    return lsvec.GetNonzeroElements(tt)

def get_Morgan(molecule,radius):
    '''
    Compute morgan descriptors
    :param molecule formula RDKit.Mol:
    :return: dictionary of nonzero descriptors
    '''
    mrgn=AllChem.GetMorganFingerprint(molecule,radius)
    dict=uisvec.GetNonzeroElements(mrgn)
    mol={}
    for key, value in dict.iteritems():
        name=key.__str__()+'r'+radius.__str__()+'d'
        mol[name]=value
    return mol

'''
Support for 2D Pharmacophore fingerprints

'''
def get_pharmacofores(molecule):
    fp=Generate.Gen2DFingerprint(molecule, Gobbi_Pharm2D.factory)
    mol={}
    for onbit in list(fp.GetOnBits()):
        mol[onbit]=1.0
    return None


"""
Support for using Padel-Descriptors

store all descriptors in dictionary with smiles as keys
"""
_desc_lib=None
_library_rdy=False

def init_library(ligands_path,data_path,screen_option,recognize_option):
    global _desc_lib
    global _library_rdy
    _library_rdy=True
    padel_path=screen_option['padel-path']
    _desc_lib = desc_library(padel_path)
    if screen_option['method']=='features':
        init_feature_library(ligands_path,data_path,screen_option,recognize_option)
    else:
        init_molecule_library(ligands_path,data_path,screen_option,recognize_option)

def init_feature_library(ligands_path,data_path,screen_option,recognize_option):
    global _desc_lib
    if recognize_option['recognize_sets']:
        ligands_path_padel = ligands_path[:-21] + 'ligands.frags.csv'
        data_path_padel = data_path[:-12] + 'data.frags.csv'
    else:
        ligands_path_padel = ligands_path[:-4] + '.frags.csv'
        data_path_padel = data_path[:-4] + '.frags.csv'
    if not os.path.isfile(ligands_path_padel):
        print(ligands_path_padel)
        ligands_path_padel=features.generate_feature_desc(ligands_path,ligands_path_padel, screen_option)
        print(ligands_path_padel)
    if not os.path.isfile(data_path_padel):
        data_path_padel=features.generate_feature_desc(data_path,data_path_padel, screen_option)
    _desc_lib.load_padel_desc_file(ligands_path_padel)
    _desc_lib.load_padel_desc_file(data_path_padel)

def init_molecule_library(ligands_path,data_path,screen_option,recognize_option):
    global _desc_lib
    if recognize_option['recognize_sets']:
        ligands_path_padel = ligands_path[:-21] + 'ligands.padel.csv'
        data_path_padel = data_path[:-12] + 'data.padel.csv'
    else:
        ligands_path_padel = ligands_path[:-4] + '.padel.csv'
        data_path_padel = data_path[:-4] + '.padel.csv'
    if not os.path.isfile(ligands_path_padel):
        _desc_lib.generate_padel_desc(ligands_path, ligands_path_padel)
    if not os.path.isfile(data_path_padel):
        _desc_lib.generate_padel_desc(data_path, data_path_padel)
    _desc_lib.load_padel_desc_file(ligands_path_padel)
    _desc_lib.load_padel_desc_file(data_path_padel)


class desc_library:
    def __init__(self,padel_path):
        self.padel_path=padel_path
        self.library={}

    def generate_padel_desc(self,file_path,padel_output):
        mol_set = io.load_Molecules(file_path)
        padel_input = file_path[:-4] + '.tmp.smi'
        io.save_smi_for_padel(padel_input, mol_set)
        log.info('Executing PaDEL ...')
        thread = subprocess.Popen(
            ['java', '-jar', #'-Xmx1024m',
             self.padel_path + '/PaDEL-Descriptor.jar',
             # '-maxruntime', '5000',
             '-threads', '8',
             '-2d',
             '-dir', padel_input,
             '-file', padel_output])
        thread.wait()
        log.info('Executing PaDEL ... done')
        os.remove(padel_input)

    def load_padel_desc_file(self,file_path):
        log.info('Loading PaDEL descriptors for: %s',file_path)
        self.library.update(io.load_padel_desc_file(file_path))

    def load_desc_file(self,file_path):
        log.info('Loading descriptors to library for: %s', file_path)
        self.library.update(io.load_padel_desc_file(file_path))

    def get_desc(self,molecule):
        smile = Chem.MolToSmiles(molecule)
        if smile in self.library:
            return self.library[smile]
        log.error("Molecule not found!")
        return None
    def get_desc_from_smiles(self,smiles):
        if smiles in self.library:
            return self.library[smiles]
        return None