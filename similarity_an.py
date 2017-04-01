# !/usr/bin/python2
# -*- coding: utf-8 -*-
import logging
import math
import settings
import numpy as np
from rdkit import Chem

import features
import descriptors

__name__ = 'virtsc_lassi.similarity_an'
__author__ = 'jiri1kolar'
#load logger
log=logging.getLogger(__name__)
_similarity_library={}


def compute_similarity(U,invS,V,molecule,weight,dw,descNames,screen_option):
    if screen_option['method']=='features':
        return compute_similarity_features(U,invS,V,molecule,weight,dw,descNames,screen_option)
    else:
        return compute_similarity_molecules(U,invS,V,molecule,weight,dw,descNames,screen_option)

def compute_similarity_molecules(U,invS,V,molecule,weight,dw,descNames,screen_option):
    z = np.zeros(len(descNames))
    dict = descriptors.get_desc(molecule, screen_option)
    if dict == None:
        return None
    return compute_similarity_vector(U,invS,V,dict,weight,dw,descNames)

def compute_similarity_features(U,invS,V,molecule,weight,dw,descNames,screen_option):

    frags=features.get_fragments_from_molecule(molecule,screen_option['fragments'])
    smiles_list=[]
    for item in frags:
        if item['smiles'] not in smiles_list:
            smiles_list.append(item['smiles'])
    cos_max=0
    if len(smiles_list)>0:
        num_frag = (int)(np.ceil(screen_option['avg_frac']*len(frags)))
        cos_list = []
        for smiles in smiles_list:
            if smiles not in _similarity_library:
                frag = features.get_molecule_from_frag_smiles(str(smiles))
                dict=descriptors.get_desc(frag,screen_option)
                if dict==None:
                    log.error("Can not compute descriptors")
                curr_cos=compute_similarity_vector(U,invS,V,dict,weight,dw,descNames)
                _similarity_library[smiles]=curr_cos
            else:
                curr_cos=_similarity_library[smiles]
            cos_list.append(curr_cos)
        sorted_cos_list=np.sort(cos_list)
        cos_selected=sorted_cos_list[-num_frag:]
        return np.average(cos_selected)
    return -1


def compute_similarity_vector(U,invS,V,dict,weight,dw,descNames):
    z = np.zeros(len(descNames))
    for key, value in dict.iteritems():
        if (key in descNames):
            di = descNames[key]
            z[di] = get_lw(value,dw[di]) * weight[di]
    cos_max = -1
    x = np.dot(z.transpose(), np.dot(U, invS))
    for j in range(V.shape[0]):
        base = V[j, :].transpose()
        curr_cos = np.dot(x, base) / np.linalg.norm(x) / np.linalg.norm(base)
        if (cos_max < curr_cos):
            cos_max = curr_cos
    return cos_max


def compute_weighted_dmm(dmm):
    wdmm=dmm
    weight=[]
    dw=[]
    n=dmm.shape[1]
    for i in range(dmm.shape[0]):
        total_freq_i=0
        min=np.nanmin(dmm[i,:])
        max=np.nanmax(dmm[i,:])
        for j in range(dmm.shape[1]):
            if (abs(dmm[i,j])>0):
                total_freq_i+=abs(dmm[i,j]);
        #weighting functions
        weight.append(get_gw(dmm[i,:],total_freq_i,max-min))
        dw.append(min)

    for i in range(wdmm.shape[0]):
        for j in range(wdmm.shape[1]):
            wdmm[i,j]=weight[i]*get_lw(dmm[i, j],dw[i])
    return [wdmm,weight,dw]


def get_lw(tfij,dwi):
    #return value of local weight function according to global settings
    if settings.local_weight_function=='log':
        return get_logw(tfij)
    elif settings.local_weight_function == 'binary':
        return get_binaryw(tfij)
    elif settings.global_weight_function == 'max':
        return tfij-dwi
    else: # 'tf=default'
        return tfij

def get_gw(dmmi,gfi,max_freq_i):
    #return value of global weight function according to global settings
    if settings.global_weight_function=='entropy':
        return get_entropyw(dmmi,gfi)
    elif settings.global_weight_function == 'idf':
        return get_idfyw(dmmi,gfi)
    elif settings.global_weight_function == 'normal':
        return get_normalw(dmmi,gfi)
    elif settings.global_weight_function == 'max':
        return 1.0/max_freq_i
    else: # 'binary'=default
        return 1


'''
    Local weight functions
'''

def get_logw(tfij):
    return math.log(1.0+tfij)
def get_binaryw(tfij):
    if (abs(tfij)<settings.threshold):
        return 0.0
    return 1.0

'''
    Global weight functions
'''

def get_entropyw(dmmi,gfi):
    gi=1.0
    for j in range(len(dmmi)):
        if (dmmi[j] > settings.threshold):
            pij = 1.0*dmmi[j] / gfi
            gi = gi + pij * math.log(pij,2) / math.log(len(dmmi),2)
    return gi

def get_idfyw(dmmi,tfij):
    #return inverse document frequency
    dfi=0
    for j in range(len(dmmi)):
        if (dmmi[j] > settings.threshold):
            dfi+=1
    return math.log(len(dmmi)/(1.0+dfi),2)

def get_normalw(dmmi,tfij):
    #return normal value of row vector
    norm=np.linalg.norm(dmmi)
    if (norm<settings.threshold):
        return 0
    return 1.0/norm