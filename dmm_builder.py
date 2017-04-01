# !/usr/bin/python2
# -*- coding: utf-8 -*-

import logging
import numpy as np
import descriptors
import io
import features
from scipy import linalg
import similarity_an

__name__ = 'virtsc_lassi.dmm_builder'
__author__ = 'jiri1kolar'
#load logger
log=logging.getLogger(__name__)
'''

Descriptors-Molecules/features-Matrix builder

'''


def build_matrix(known_ligands_path,screen_option):
    if screen_option['method']=='features':
        return build_feature_matrix(known_ligands_path,screen_option)
    else:
        known_ligands_set = io.load_Molecules(known_ligands_path)
        return build_molecule_matrix(known_ligands_set,screen_option)

def build_feature_matrix(file_path,screen_option):
    '''
    Make a matrix with features and their descriptors
    :param moleculeSet:
    :param screen_option:
    :return:
    '''

    log.info('Building Descriptor-Feature matrix...')
    desc_csv=features.generate_feature_desc(file_path,"",screen_option)
    mol_lib=io.load_padel_desc_file(desc_csv)
    #descriptors.desc_library_load_desc_file(desc_csv)
    descNames = {}
    dictList=[]
    # first scan for used descriptors
    j = 0
    for smiles,dict in mol_lib.iteritems():
        for key in dict.keys():
            if key not in descNames:
                descNames[key] = j
                j+=1
        dictList.append(dict)
    dfm=make_matrix(dictList,descNames)
    '''
    [wdfm, weight] = similarity_an.compute_weighted_dmm(dfm)
    U, S, VT = svd(wdfm, 0)
    V=VT.transpose()
    cos_th=0.6
    class_v=[[V[0,:]]]
    class_dfm=[[dfm[:,0]]]
    for i in range(V.shape[0]):
        vi=V[i,:]
        new_class = True
        for j in range(len(class_v)):
            add_to_class=True
            for vj in class_v[j]:
                curr_cos = np.dot(vi, vj) / np.linalg.norm(vi) / np.linalg.norm(vj)
                if curr_cos<cos_th:
                    add_to_class=False
            if add_to_class:
                class_v[j].append(vi)
                class_dfm[j].append(dfm[:,i])
                new_class=False
        if new_class:
            class_v.append([vi])
            class_dfm.append([dfm[:, i]])

    class_avg_dfm=[]
    for cls in class_dfm:
        avg=np.zeros(len(descNames))
        for vi in cls:
            avg+=vi
        avg=avg/len(cls)
        class_avg_dfm.append(avg)
    class_avg_dfm=np.matrix(class_avg_dfm).transpose()
    '''
    return descNames,dfm

def build_molecule_matrix(moleculeSet,screen_option):
    '''
    Make a matrix from input smiles file
    Use Atom Pair fingerprint descriptors
    :param file_path:
    :return: dictionary of used desc and its indexes,matrix(data)
    '''
    # load known ligands from file
    log.info('Building Descriptor-Molecule matrix...')
    numDesc=0
    descNames={}
    dictList=[]
    #first scan for used descriptors
    j=0
    for molecule in moleculeSet:
        dict=descriptors.get_desc(molecule,screen_option)
        for key in dict.keys():
            if key not in descNames:
                descNames[key] = j
                j+=1
        dictList.append(dict)
    j=0
    perm = np.random.permutation(len(descNames))
    for key in descNames.keys():
        descNames[key] = perm[j]
        j = j + 1

    dmm=make_matrix(dictList,descNames)
    #return sparse.csc_matrix(np.matrix(A).transpose())
    return descNames,dmm

def make_matrix(dictList,descNames):
    numDesc = len(descNames)
    numMol = len(dictList)
    # allocate matrix and compute desc
    dmm = np.zeros((numDesc, numMol))
    j = 0
    for dict in dictList:
        for key, value in dict.iteritems():
            i = descNames[key]
            dmm[i, j] = value
        j = j + 1
    return dmm

def svd(mat,k):
    '''
    Make SVD decomposition
    :param in_matrix:
    :return: U, S_diagonal, VT; that A=U.s.VT
    '''
    log.info('Making singular value decomposition')
    u,s,vt=np.linalg.svd(mat,full_matrices=False)
    #u,s,vt=linalg.svd(mat,full_matrices=False,check_finite=True)
    th=1e-9
    i=0
    nonzero_s=[]
    for si in s:
        if (abs(si)>th):
            nonzero_s.append(si)
            i+=1
    nonzero_u=u[:,:i]
    nonzero_vt=vt[:i,:]
    nonzero_s=np.diag(nonzero_s)
    return nonzero_u,nonzero_s,nonzero_vt