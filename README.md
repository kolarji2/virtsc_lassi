# virtsc_lassi
Program for virtual screening using latent semantic indexing

For more information about usage:
python2 main.py -h

usage: main.py [-h] [-s FILE] [--selection-muv FILE]
               [--nligands-to-train NLIGANDS_TO_TRAIN] [-c FILE] [-d FILE]
               [-od FILE] [-f FRAGMENT_TYPES]
               [-m {fragments,descriptors,features}] [-g {rdkit,padel}]
               [-ths THRESHOLD_SIMILARITY]

optional arguments:
  -h, --help            show this help message and exit
  -s FILE, --selection-file FILE
                        Prepare selection from json file and save it in to output directory.
  --selection-muv FILE  Prepare selection from MUV directory and save it to output directory.
  --nligands-to-train NLIGANDS_TO_TRAIN
                        Number of ligands choosen from MUV set as known-ligands
  -c FILE, --config-file FILE
                        Json configuration file
  -d FILE, --directory FILE
                        Directory where program should look for molecule sets.
  -od FILE, --output-directory FILE
                        Output directory for results, by default results are saved in same directory as data
  -f FRAGMENT_TYPES, --fragment-types FRAGMENT_TYPES
                        Optional comma separated list of fragment types to extract (eg. "tt.3,ecfp.2"), special value ap.0 for atom-pairs, or tt.0 for topological torsion
  -m {fragments,descriptors,features}, --method {fragments,descriptors,features}
                        Method of analyzing molecules. Allowed values are:
                        	[fragments] - each molecule is prepresented by fragment fingerprint
                        	[descriptors] - each molecule is represented by descriptors
                        	[features] - each molecule is split to fragments and each fragment is represented by descriptors
  -g {rdkit,padel}, --descriptors-generator {rdkit,padel}
                        Generator to be used to obtain fragments descriptors.

Sample usage:
python2 main.py -f ecfp.2 -d results/8.0-8.5/ -c default.json -m features
