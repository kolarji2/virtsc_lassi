# !/usr/bin/python2
# -*- coding: utf-8 -*-

import time
import logging
import numpy as np
from rdkit import Chem
from rdkit.ML.Scoring import Scoring


#internal modules
import similarity_an
import dmm_builder
import io
import os.path
import descriptors
import features

__name__ = 'virtsc_lassi.analyzator'
__author__ = 'jiri1kolar'
#load logger
log=logging.getLogger(__name__)

def compare(input_path_set,output_path,recognize_option,screen_option):
    known_ligands_path = input_path_set[0]
    data_path = input_path_set[1]
    ligands_path = input_path_set[2]

    '''
        TO DO
            BETTER implementation of PaDEL and features
    '''

    #prepare PaDel descriptors
    if screen_option['descriptors']=='padel':
        descriptors.init_library(ligands_path,data_path,screen_option,recognize_option)

    '''
    if screen_option['method'] == 'features':
        descriptors.init_library(screen_option['padel-path'])
        frag_json = features.fragments_extraction(data_path, screen_option['fragments'])
        data_frags_csv = features.descriptors_extraction(frag_json, screen_option['descriptors'],
                                                             screen_option['padel-path'])
        descriptors.desc_library_load_desc_file(data_frags_csv)
    '''
    #reset similariti lib
    similarity_an._similarity_library={}
    #generate Descriptors-Molecules-Matrix
    descNames,dmm=dmm_builder.build_matrix(known_ligands_path,screen_option)
    #compute weighted Descriptors-Molecules-Matrix and local and global weight coefficients
    [wdmm,weight]=similarity_an.compute_weighted_dmm(dmm)
    log.info('Number of analysed descriptors: %i', len(descNames))
    #compute Singular Value decomposition
    k = -1  # use all singular values
    U,S,VT=dmm_builder.svd(wdmm,k)
    #load tested data set and ligands for identification of correct results
    ligands_set=io.load_Molecules(ligands_path)
    io.create_parent_directory(output_path)
    with open(output_path + '/result.csv', 'w') as output_stream:
        compareData(U,S,VT,weight,data_path,descNames,output_stream,screen_option)
    #with open(output_path + '/ligands.csv', 'w') as output_stream:
    #    compareData(U,S,VT,ligands_path,output_stream)
    with open(output_path + '/result.csv', 'r') as input_stream:
        with open(output_path + '/result_sorted.csv', 'w') as output_stream:
            sortResults(input_stream,output_stream)
    log.info('Computing AUC and EF...')
    with open(output_path + '/result_sorted.csv', 'r') as input_stream:
        scores,col=make_scores(input_stream,ligands_set)
        with open(output_path + '/result_stat_data.csv', 'w') as output_stream:
            for row in scores:
                output_stream.write(str(row[0])+', '+str(row[1])+', '+str(row[2]+ '\n'))
        auc=Scoring.CalcAUC(scores,col)
        log.info('\tCurrent AUC: %.3f', auc)
        fraction=[0.5/100,1./100,1.5/100,2./100]
        ef=Scoring.CalcEnrichment(scores,col,fraction)
        with open(output_path + '/result_stat.csv', 'w') as output_stream:
            output_stream.write(str(auc)+'\n')
            i=0
            for result in ef:
                output_stream.write(str(fraction[i]) + ', ' + str(result)+'\n')
                i=i+1
    return auc

def compareData(U,S,VT,weight,data_path,descNames,output_stream,screen_option):
    '''
    Compute similarity of molecules in data_path to known actives
    :param U,S,VT,file_path:
    :return: resultSimilarity
    '''
    log.info('Computing similarity:')
    np.set_printoptions(threshold=np.nan)
    invS=np.linalg.inv(S)
    V=VT.transpose()
    #imax=min(900,len(moleculeSet))
    imax=file_len(data_path)

    with open(data_path, 'r') as mol_stream:
        log.info('Molecules to process: %i',imax)
        stopwatch=time.clock()
        startwatch=time.clock()
        percent=0.1
        for i in range(imax):
            if (i>percent*imax-2):
                log.info('\tProcessed: %i/%i Remaining time: %f min',i+1,imax,(imax-i-1)*(time.clock()-startwatch)/(i+1)/60)
                stopwatch=time.clock()
                percent+=0.3
            molecule=io.next_smile_molecule(mol_stream)
            cos_max=similarity_an.compute_similarity(U,invS,V,molecule,weight,descNames,screen_option)
            if cos_max == None:
                log.error("Wrong molecule")
                continue
                #result[molecule]=cos_max
            output_stream.write(cos_max.__str__()+', '+Chem.MolToSmiles(molecule)+'\n')

def sortResults(input_stream, output_stream):
    input=input_stream.read().split('\n')
    results={}
    input.pop(len(input)-1)
    for line in input:
        s=line.split(', ')
        key=s[1]
        value=s[0]
        results[key]=value
    sortedKeys=sorted(results, key=results.__getitem__, reverse=True)
    for key in sortedKeys:
            output_stream.write(results[key].__str__()+', '+ key +'\n')

def make_scores(input_stream,ligands_set):
    '''
    Compute similarity of molecules in data_path to known actives
    :param input_stream
    :param ligands_set
    :return: score:List[n,3]
    '''
    scores=[]
    ligands_smiles={}
    #generate ligands smiles
    for mol in ligands_set:
        ligands_smiles[Chem.MolToSmiles(mol)]=1
    input=input_stream.read().split('\n')
    for line in input:
        s=line.split(', ')
        if (len(s)==2):
            if s[1] in ligands_smiles:
                scores.append([s[0],True,s[1]])
            else:
                scores.append([s[0],False,s[1]])
    return scores,1

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1
