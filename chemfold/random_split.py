#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""random_split.py: Assigns a random fold number to each input.

Assigns a random fold number to each input for the ML model to compare the generated distribution of data and labels over the folds and the resulting performance with other methods like LFC, LSH or a scaffold-based fold splitting. This script is part of the folding single partner study of WP1 from the MELLODDY project. It is highly inspired by the script 'split_by_scaffold.py by Ansgar Schuffenhauer.
"""
#ATTENTION! WARNING! This is a beta version. Use on your own risk!


import argparse
import numpy as np
import pandas as pd
import json
import hashlib
import random
import hmac
import os, sys, logging

__author__ = "Lina Humbeck, Boehringer Ingelheim Pharma GmbH & Co.KG"
__contributors__ = "Ansgar Schuffenhauer, Novartis"

shandler = logging.StreamHandler(stream=sys.stdout)
formatter = logging.Formatter('%(asctime)s -%(name)s - %(levelname)s - %(message)s')
shandler.setFormatter(formatter)
shandler.setLevel(logging.INFO)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
logger.addHandler(shandler)

def hashed_fold_scaffold(scaff, secret, nfolds = 5):
    scaff = str(scaff).encode("ASCII")
    h = sha256([scaff], secret)
    random.seed(h, version=2)
    return random.randint(0, nfolds - 1)

def sha256(inputs, secret):
    """
    Encryption function using python's pre-installed packages HMAC and hashlib.
    We are using SHA256 (it is equal in security to SHA512).
    :param inputs: input strings
    :param secret: given pharma partner key
    :return:
    """
    m = hmac.new(secret, b'', hashlib.sha256)
    for i in inputs:
        m.update(i)
    return m.digest()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='compute fold vectors based on Scaffolds')
    parser.add_argument('--infolder', help="Input folder, working directory of melloddy_tuner, expected to conatin results and results_tmp subfolder ", type=str, required=True)
    parser.add_argument('--out', help="output folder", type=str, required=True)
    parser.add_argument('--params_file',help='path to parameters.json file',required=True)

    args = parser.parse_args()
    logger.info(args)

    with open(args.params_file) as cnf_f:
        params = json.load(cnf_f)

    #read standardized structures
    fname = os.path.join(args.infolder,'results_tmp','standardization','T2_standardized.csv')
    t2_standardized = pd.read_csv(fname)
    logger.info('read in standardized structure file: {}'.format(fname))

    t2_unique = t2_standardized[['canonical_smiles']].drop_duplicates()
    t2_scaff = t2_standardized.merge(t2_unique[['canonical_smiles']], on='canonical_smiles')
    logger.info('extracted unique compounds:{0}'.format(len(t2_scaff)))

    #now need to merg the T5 table, to get scaffold assignments to unique descriptors IDs
    t5 = pd.read_csv(os.path.join(args.infolder,'results_tmp','descriptors','mapping_table_T5.csv'))
    t2_joined = t2_scaff.merge(t5, on='input_compound_id')
    t6_scaff = t2_joined.groupby('descriptor_vector_id')['canonical_smiles'].min()

    #now we need to join to final file with continuous descriptor vector Ids
    t2_t11 = pd.read_csv(os.path.join(args.infolder,'results','T11.csv'))
    t2_t11_joined = t2_t11.join(t6_scaff,on='descriptor_vector_id')

    ## creating output path if it does not exist
    os.makedirs(args.out, exist_ok=True)

    #now we need to generate the hashes for the scaffolds
    key = params['key']['key'].encode("ASCII")
    nfolds = params['lsh']['nfolds'] 
    t2_t11_joined['fold'] = t2_t11_joined['canonical_smiles'].apply(lambda x : hashed_fold_scaffold(x, secret = key, nfolds = nfolds))
    t2_t11_joined.to_csv(os.path.join(args.out,'T2_T11_random_folds.csv'))
    logger.info('written out csv file with fold info')

    #now create numpy arrays for folds
    random_fold = np.zeros(len(t2_t11_joined['cont_descriptor_vector_id']))
    random_fold[t2_t11_joined['cont_descriptor_vector_id']] = t2_t11_joined['fold']

    np.save(os.path.join(args.out,'random_folds.npy'),random_fold)
    logger.info('written out numpy arrays with folds')

