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
    else:
        log.error("Can not generate descriptors, unknown option!")
    return None

def get_descriptors(molecule,type):
    if type=='rdkit':
        return get_RDkit_Desc(molecule)
    elif type=='padel':
        return get_Padel_Desc(molecule)
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

def init_library(padel_path):
    global _desc_lib
    global _library_rdy
    _library_rdy=True
    _desc_lib=desc_library(padel_path)

def padel_generate_desc(file_path,padel_output):
    global desc_lib
    global lib_rdy
    if not lib_rdy:
        log.error("Padel not rdy, you have to set path to padel and initialize it by calling init_padel.")


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