'''
Process virtual screening on tested data sets
'''
import logging
import os
import json
#local
import io
import fnmatch
import analyzator

__name__ = 'virtsc_lassi.screen'
__author__ = 'jiri1kolar'
__email__ = "jiri1kolar@gmail.com"

#load logger
log=logging.getLogger(__name__)

def screen(directory,format,output_path,recognize_option,screen_option):
    '''
    process screening in directory where are located files:
    known-ligands.smi/sdf - training set of ligands
    ligands.smi/sdf - all ligands for later evaluation
    data.smi/sdf - mixed ligands and other molecules
    :param directory:
    :return:
    '''

    #file check, if files exist, and if it it smiles or sdf format
    if directory[-1]=='/':
        directory=directory[:-1]
    known_ligands_path = directory + '/known-ligands.'+format
    ligands_path = directory + '/ligands.'+format
    data_path = directory + '/data.'+format
    sets_path = [known_ligands_path, data_path, ligands_path]
    auc = analyzator.compare(sets_path, output_path,recognize_option,screen_option)
    return auc

def recursive_screen(root_data_dir,root_result_dir,recognize_option,screen_option,screen_info):
    '''
    search root_data_dir for all subdirectories that involves data files and process screening on them,
    save results to directory with same name as found, but in root_result_dir
    if root_result_dir is empty string, results are saved in the same dir as data was found
    Make an average statisics from screening and save it to root_result_dir
    :param root_data_dir:
    :param root_result_dir:
    :param recognize_option:
    :return:
    '''
    files=io.find_files_recursively(root_data_dir,'data.s[md][if]')

    avr_result_path=root_result_dir
    if len(avr_result_path)>0 and avr_result_path[-1]!='/':
        avr_result_path += '/'
    setList = {}
    if recognize_option['recognize_sets']:
        '''
        Assumes that dir tree is as folows:
        */set_name/dir1/data
        */set_name/dir2/data
        */set_name/dir3/data
        */set_name/any other dir/data
        ...
        '''
        for file_path in files:
            set_dir=file_path.split('/')[-3]
            if set_dir not in setList:
                setList[set_dir]=[]
            setList[set_dir].append(file_path)
    else:
        setList['default'] = files

    set_names=setList.keys()
    set_names.sort()
    for set_name in set_names:
        print(set_name)
        set=setList[set_name]
        collection_name = ""
        if recognize_option['recognize_collection']:
            if len(set) > 0:
                file_path = set[0].split('/')
                if len(file_path) >= 5:
                    collection_name = file_path[-5]
        avg_auc = 0.0
        vrc_auc = 0.0
        tot=len(set)
        for file in set:
            format=file[-3:]
            data_dir=file[:-8]
            if len(root_result_dir) > 0:
                result_dir=data_dir.replace(root_data_dir,'')
                if len(result_dir)>0 and result_dir[0]=='/':
                    result_dir=result_dir[1:]
                if root_result_dir[-1] == '/':
                    root_result_dir = root_result_dir[:-1]
                result_dir=root_result_dir+'/'+result_dir
            else:
                result_dir=data_dir
            auc=screen(data_dir,format,result_dir,recognize_option,screen_option)
            avg_auc+=auc
            vrc_auc+=auc*auc
        avg_auc=avg_auc/tot
        vrc_auc = vrc_auc / tot - (avg_auc * avg_auc)
        log.info("Average AUC: %.3f +- %.3f", avg_auc, vrc_auc)
        with open(avr_result_path +collection_name+screen_option['lw']+screen_option['gw']+'resultingAUC.csv', 'a') as output_stream:
             output_stream.write("{0:s}, {1:s}, {2:f}, {3:f}, {4:d}, {5:s} \n".format(collection_name, set_name, avg_auc, vrc_auc, tot,screen_info))
