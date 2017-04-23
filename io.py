#!/usr/bin/python2
# -*- coding: utf-8 -*-

import os
import json
import csv
import random
import fnmatch
import logging
import numpy as np
from rdkit import Chem

__name__ = 'virtsc_lassi.io'
__author__ = 'jiri1kolar'

#load logger
log=logging.getLogger(__name__)

#osetrit cteni a zapis proti spatnych vstupum


def create_parent_directory(path):
    """Create directory if it does not exists.

    :param path:
    :return:
    """
    dir_name = os.path.dirname(path)
    if not os.path.exists(dir_name) and not dir_name == "":
        os.makedirs(dir_name)

def find_files_recursively(directory, pattern):
    matches = []
    for root, dirnames, filenames in os.walk(directory):
        for filename in fnmatch.filter(filenames, pattern):
            matches.append(os.path.join(root, filename))
    return matches

def to_float(x):
    try:
        a = float(x)
        if np.isinf(a): a = float('nan')
    except ValueError:
        return float('nan')
    else:
        return a

def save_smi_file(file_path,moleculesSet):
    create_parent_directory(file_path)
    with open(file_path,'w') as stream:
        for mol in moleculesSet:
            stream.write(Chem.MolToSmiles(mol))
            stream.write('\n')


def load_smi_file(file_path):
    """Load smi from smi file and return them

    :param dataset_directory:
    :param file_name:
    :return: moleculeSet:RDKit.Mol[]
    """
    moleculeSet=[]
    log.info('Loading (SMI): %s',file_path)
    with open(file_path) as input_stream:
        for line in input_stream:
            line = line.strip()
            molecule = Chem.MolFromSmiles(line)
            if molecule is None:
                log.error('Invalid molecule detected. (Skipped)')
                continue
            # Molecules created from SMILES does not have any name,
            # so we use the SMILES as a name.
            molecule.SetProp('_Name', line)
            moleculeSet.append(molecule)
    return moleculeSet

def load_sdf_file(file_path):
    """Load molecules from sdc file and return them

    :param dataset_directory:
    :param file_name:
    :return: moleculeSet:RDKit.Mol[]
    """
    moleculeSet = []
    log.info('Loading (SDF): %s',file_path)
    for molecule in Chem.SDMolSupplier(file_path):
        if molecule is None:
            log.error('Invalid molecule detected. (Skipped)')
            continue
        moleculeSet.append(molecule)
    return moleculeSet

def load_Molecules(file_path):

    s=file_path.split('.')
    suffix=s[len(s)-1]
    moleculeSet=[]
    if (suffix=='sdf'):
        moleculeSet=load_sdf_file(file_path)
    elif (suffix=='smi'):
        moleculeSet=load_smi_file(file_path)
    else:
        log.error('Unknown molecule file format')
    return moleculeSet

def next_smile_molecule(mol_stream):
    line=mol_stream.readline()
    mol=line.split('\n')[0]
    return Chem.MolFromSmiles(mol)

def prepare_selection(configuration,output_dest,root_directory):
    progressNum = 0.0
    log.info("Preparing collection...")
    for i in range(0, len(configuration)):
        log.info("Creating selection from: %s", configuration[i]['selection'])
        collection_path = os.path.join(root_directory,'/instances',configuration[i]['collection_path'])
        selection_path = os.path.join(root_directory,'instances',configuration[i]['selection'])
        output_path = os.path.join(root_directory,output_dest,configuration[i]['selection'])
        selection_files = os.listdir(selection_path)
        decoys_path_tmp=""
        ligands_path_tmp=""
        decoys_set=[]
        ligands_set=[]
        for file in selection_files:
            output_path_spec = output_path + '/' + file.split('.')[0]
            with open(selection_path + '/' + file) as spec_cfg_stream:
                spec_cfg = json.load(spec_cfg_stream)
                # select data
                decoys_path = collection_path + '/' + spec_cfg['files'][0] + '.sdf'
                ligands_path = collection_path + '/' + spec_cfg['files'][1] + '.sdf'
                data_test = spec_cfg['data']['test']
                data_train_ligands = spec_cfg['data']['train']['ligands']
                # data_train_decoys=spec_cfg['data']['train']['decoys']
                decoys_path = decoys_path.encode('utf-8')
                ligands_path = ligands_path.encode('utf-8')
                if len(ligands_set)==0 or ligands_path_tmp!=ligands_path:
                    ligands_set = load_Molecules(ligands_path)
                if len(decoys_set)==0 or decoys_path_tmp!=decoys_path:
                    decoys_set = load_Molecules(decoys_path)
                decoys_path_tmp=decoys_path
                ligands_path_tmp=ligands_path
                known_ligands_set = load_MoleculesFromSelectionSDF(data_train_ligands, ligands_set, decoys_set,
                                                                      'Known-ligands')
                moleculeSet = load_MoleculesFromSelectionSDF(data_test, ligands_set, decoys_set, 'Data-set')
                # Create smi files for later processing
                known_ligands_path = output_path_spec + '/known-ligands.smi'
                ligands_path = output_path_spec + '/ligands.smi'
                data_path = output_path_spec + '/data.smi'
                save_smi_file(known_ligands_path, known_ligands_set)
                save_smi_file(ligands_path, ligands_set)
                save_smi_file(data_path, moleculeSet)

def prepare_selection_MUV(muv_path,output_dest,opt):
    '''
    prepare selection in form: MUV/???/known-ligands.smi'
                              MUV/???/ligands.smi
                              MUV/???/data.smi
    :param muv_path:
    :param output_dest:
    :return:
    '''
    logging.info("Preparing collection from MUV data sets...")
    nknown_actives = opt['nligands-to-train']
    nrandom_sets = opt['nrandom-set']
    logging.info("\tchoosing %d random ligands as known-ligands",nknown_actives)
    logging.info("\tgenerating %d random sets", nrandom_sets)
    muv_files= sorted(os.listdir(muv_path))
    muv_data_sets={}
    output_path=os.path.join(output_dest,'muv')

    random.seed(33)
    for file in muv_files:
        array=file.split('_')
        muv_id=array[-2]
        type=array[-1].split('.')[0]
        if muv_id not in muv_data_sets:
            muv_data_sets[muv_id]={}
        muv_data_sets[muv_id][type] = file
    for muv_id,data_set in muv_data_sets.iteritems():
        #load data
        print('\tProcessing set {0:s}...'.format(muv_id))
        ligands_path=os.path.join(muv_path,data_set['actives'])
        decoys_path=os.path.join(muv_path,data_set['decoys'])
        ligands_set=load_molecules_from_muv_data_file(ligands_path)
        decoys_set=load_molecules_from_muv_data_file(decoys_path)
        for i in range(nrandom_sets):
            zeros='00'
            output_path_spec=os.path.join(output_path,muv_id,zeros[(i+1)//10:]+str(i+1))
            known_ligands_path = os.path.join(output_path_spec,'known-ligands.smi')
            ligands_path = os.path.join(output_path_spec,'ligands.smi')
            data_path = os.path.join(output_path_spec,'data.smi')
            known_ligands_set=[]
            unknown_ligands_set=list(ligands_set)
            #randomly choose known ligands
            for j in range(nknown_actives):
                index=int(np.floor(random.random()*len(unknown_ligands_set)))
                known_ligands_set.append(unknown_ligands_set[index])
                unknown_ligands_set.pop(index)
            save_smi_file(ligands_path, ligands_set)
            save_smi_file(known_ligands_path,known_ligands_set)
            save_smi_file(data_path,decoys_set+unknown_ligands_set)

def load_molecules_from_muv_data_file(ligands_set):
    '''
    file is expected to be in format
    # PUBCHEM_COMPOUND_CID	ID	SMILES
    22162	MUV_859_A_1	Cl.CN(C)CCOC(=O)C(c1ccccc1)C1(O)CCCC1
    ...
    :param ligands_set:
    :return:
    '''
    molecules=[]
    with open(ligands_set,'r') as stream:
        lines=stream.readlines()
        for line in lines:
            smiles=line.split()[-1].strip()
            if (smiles=='SMILES'):
                continue
            mol=Chem.MolFromSmiles(smiles)
            if mol is not None:
                molecules.append(mol)
    return molecules

def load_MoleculesFromSelectionSDF(selection,ligandsSet,decoysSet,info):
    '''

    :param selection:
    :param ligandsSet:
    :param decoysSet:
    :param info:
    :return:selectedMol:RDKit.Mol[]
    '''
    selectedMol=[]
    selection_dic={}
    for target in selection:
        selection_dic[target['name']]=1
    for mol in ligandsSet:
        if not (type(mol) is Chem.rdchem.Mol):
            log.error("Cant load molecule from ligands, molecule is skipped")
            continue
        if mol.GetProp('_Name').__str__() in selection_dic:
            selectedMol.append(mol)

    for mol in decoysSet:
        if not (type(mol) is Chem.rdchem.Mol):
            log.error("Cant load molecule from decoys, molecule is skipped")
            continue
        if mol.GetProp('_Name').__str__() in selection_dic:
            selectedMol.append(mol)
    return selectedMol

def parseFragmentTypes(fragment_types):
    parsed_types = []
    for item in fragment_types.split(','):
        item_split = item.split('.')
        if not len(item_split) == 2:
            log.error('Invalid fragment type: %s', item)
            log.info('  Expected format {TYPE}.{SIZE}')
            return
        parsed_types.append({
            'name': item_split[0],
            'size': int(item_split[1])
        })
    return parsed_types

def load_padel_desc_file(file_path):
    th = 1e-9
    desc_names = []
    molecules = {}
    for row in csv.reader(open(file_path, "r")):
        if not row: continue
        if len(desc_names) == 0:
            for name in row[1:]:
                desc_names.append(name)
        else:
            if row[0] not in molecules:
                mol_desc={}
                for name,value in zip(desc_names,row[1:]):
                    num=to_float(value)
                    if np.isnan(num) or abs(num)<=th:
                        continue
                    mol_desc[name]=num
                molecules[row[0]]=mol_desc
    return molecules

def save_smi_for_padel(file_path,moleculesSet):
    create_parent_directory(file_path)
    with open(file_path,'w') as stream:
        for mol in moleculesSet:
            smiles=Chem.MolToSmiles(mol)
            stream.write(smiles)
            stream.write('\t')
            stream.write(smiles)
            stream.write('\n')

