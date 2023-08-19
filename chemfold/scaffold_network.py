import argparse
import numpy as np
import pandas as pd
import pickle

import json

import hashlib
import random
import hmac

from rdkit import Chem
from rdkit.Chem.Scaffolds import rdScaffoldNetwork, MurckoScaffold
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import PandasTools

import os, sys, logging

def has_unusual_ringsize(mol):
    return len([len(x) for x in mol.GetRingInfo().AtomRings() if len(x)>6 or len(x)<5 ])>0
    
def has_macrocycle(mol):
    return len([len(x) for x in mol.GetRingInfo().AtomRings() if len(x)>9])>0

def murcko_scaff_smiles(mol_smiles):
    mol = Chem.MolFromSmiles(mol_smiles)
    if mol is not None:
        return Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(mol))
    else:
        return None

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

    shandler = logging.StreamHandler(stream=sys.stdout)
    formatter = logging.Formatter('%(asctime)s -%(name)s - %(levelname)s - %(message)s')
    shandler.setFormatter(formatter)
    shandler.setLevel(logging.INFO)
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    logger.addHandler(shandler)


    parser = argparse.ArgumentParser(description='compute fold vectors based on Scaffolds')
    parser.add_argument('--infolder', help="Input folder, working directory of melloddy_tuner, expected to conatin results and results_tmp subfolder ", type=str, required=True)
    parser.add_argument('--out', help="output folder", type=str, required=True)
    parser.add_argument('--params_file',help='path to parameters.json file',required=True)
    parser.add_argument('--sn_flattenIsotopes',help = 'controls flattenIsotopes parameter of scaffold network',type=bool,default=True)
    parser.add_argument('--sn_includeGenericBondScaffolds',help = 'controls includeGenericBondScaffolds parameter of scaffold network',type=bool,default=False)
    parser.add_argument('--sn_includeGenericScaffolds',help = 'controls includeGenericScaffolds parameter of scaffold network',type=bool,default=False)
    parser.add_argument('--sn_includeScaffoldsWithAttachments',help = 'controls includeScaffoldsWithAttachments parameter of scaffold network',type=bool,default=False)
    parser.add_argument('--sn_includeScaffoldsWithoutAttachments',help = 'controls includeScaffoldsWithoutAttachments parameter of scaffold network',type=bool,default=True)
    parser.add_argument('--sn_pruneBeforeFragmenting',help = 'controls pruneBeforeFragmenting parameter of scaffold network',type=bool,default=True)
    parser.add_argument('--nrings', help="preferred number of rings for scaffold, defines the cut-level of scaffold network", type=int, default=3)


    args = parser.parse_args()
    logger.info(args)

    with open(args.params_file) as cnf_f:
        params = json.load(cnf_f)

    snparams = rdScaffoldNetwork.ScaffoldNetworkParams()
    snparams.flattenIsotopes=args.sn_flattenIsotopes
    snparams.includeGenericBondScaffolds=args.sn_includeGenericBondScaffolds
    snparams.includeGenericScaffolds=args.sn_includeGenericScaffolds
    snparams.includeScaffoldsWithAttachments = args.sn_includeScaffoldsWithAttachments
    snparams.includeScaffoldsWithoutAttachments = args.sn_includeScaffoldsWithoutAttachments
    snparams.pruneBeforeFragmenting = args.sn_pruneBeforeFragmenting



    #read standardized structures
    fname = os.path.join(args.infolder,'results_tmp','standardization','T2_standardized.csv')
    t2_standardized = pd.read_csv(fname)
    logger.info('read in standardized structure file: {}'.format(fname))

    #calculate Murcko Scaffold
    t2_standardized['murcko_scaff_smiles'] = t2_standardized['canonical_smiles'].apply(murcko_scaff_smiles)
    logger.info('created murcko scaff smiles column')

    #create unique Murcko Scaffolds dataframe
    t2_murcko_unique = t2_standardized[['murcko_scaff_smiles']].drop_duplicates()
    logger.info('extracted unique murcko scaffolds:{0}'.format(len(t2_murcko_unique)))



    #create scaffold network for each unique Murcko scaffold
    sn_dict = {}
    failed_list = []
    for row in t2_murcko_unique.itertuples():
        logger.debug('processing\t|{}|'.format(row.murcko_scaff_smiles))
        if row.murcko_scaff_smiles is not None and row.murcko_scaff_smiles != '' and pd.notnull( row.murcko_scaff_smiles):
            mol = Chem.MolFromSmiles(row.murcko_scaff_smiles)
            if mol is not None:
                try:
                    sn_dict[row.murcko_scaff_smiles] = rdScaffoldNetwork.CreateScaffoldNetwork([mol],snparams)
                except:
                    failed_list.append(row.murcko_scaff_smiles)
            else:
                failed_list.append(row.murcko_scaff_smiles)
        else:
            failed_list.append(row.murcko_scaff_smiles)       
    logger.info('generated scaffold network')

    ## creating output path if it does not exist
    os.makedirs(args.out, exist_ok=True)

    #save the scaffold netowrks
    with open(os.path.join(args.out,'T2_sn_dict.pkl'),'wb') as pklf:
        pickle.dump(sn_dict,pklf)
    logger.info('written pickle file for scaffold networks')

    #with open(os.path.join(args.out,'T2_sn_dict.pkl'),'rb') as pklf:
    #    sn_dict = pickle.load(pklf)
    #logger.info('read pickle file for scaffold networks')

    #write out failed list
    logger.info('encountered {} murcko scaffolds failing in scaffold network'.format(len(failed_list)))
    fname = os.path.join(args.out,'failed_sn_murcko_scaffolds.txt')
    pd.Series(failed_list).to_csv(fname,sep='\t')

    #make a unique data frame of scaffold network nodes, and clacukate properties
    all_nodes = set([])
    for sn in sn_dict.values():
        all_nodes |= set([str(n) for n in sn.nodes])
    logger.info('extracted nodes')    
    logger.info('number of nodes {}'.format(len(all_nodes)))

    node_df = pd.DataFrame({'node_smiles':list(all_nodes)})
    PandasTools.AddMoleculeColumnToFrame(node_df,'node_smiles','mol',includeFingerprints=False)
    node_df['num_rings'] = node_df['mol'].apply(Chem.rdMolDescriptors.CalcNumRings)
    node_df['num_rings_delta'] = (node_df['num_rings'] -args.nrings).abs() 
    node_df['num_rbonds'] = node_df['mol'].apply(Chem.rdMolDescriptors.CalcNumRotatableBonds)
    node_df['num_hrings'] = node_df['mol'].apply(Chem.rdMolDescriptors.CalcNumHeterocycles)
    node_df['num_arings'] = node_df['mol'].apply(Chem.rdMolDescriptors.CalcNumAromaticRings)
    node_df['num_bridge'] =  node_df['mol'].apply(Chem.rdMolDescriptors.CalcNumBridgeheadAtoms)
    node_df['num_spiro'] =  node_df['mol'].apply(Chem.rdMolDescriptors.CalcNumSpiroAtoms)
    node_df['has_macrocyle'] = node_df['mol'].apply(has_macrocycle)
    node_df['has_unusual_ring_size'] = node_df['mol'].apply(has_unusual_ringsize)
    logger.info('calculated node properties')   

    #the follwing lines define the precedence of nodes, we aim at ~3 ring scaffolds
    priority_cols =['num_rings_delta','has_macrocyle','num_rbonds','num_bridge','num_spiro','has_unusual_ring_size','num_hrings', 'num_arings','node_smiles']
    priority_asc = [True,False,True,False,False,False,False,True,True]
    node_df.sort_values(priority_cols, ascending = priority_asc, inplace =True)
    node_df['priority'] = np.arange(len(node_df))
    node_df.set_index('node_smiles',inplace=True)
    logger.info('sorted nodes by priority')

    with open(os.path.join(args.out,'T2_sn_nodes.pkl'),'wb') as pklf:
        pickle.dump(node_df,pklf)
    logger.info('written pickle file for node table')
    #with open(args.out,'T2_sn_nodes.pkl'),'rb') as pklf:
    #    node_df = pickle.load(pklf)
    #logger.info('read pickle file for node table')

    #assign for eac Murcko scaffold the preferred sn scaffolds
    for row in t2_murcko_unique.itertuples():
        if row.murcko_scaff_smiles in sn_dict:
            nodelist = [str(n) for n in sn_dict[row.murcko_scaff_smiles].nodes]
            node_data = node_df.loc[nodelist]
            t2_murcko_unique.loc[row.Index,'sn_scaff_smiles'] = node_data['priority'].idxmin()
    logger.info('assigned preferred scaffolds')

    with open(os.path.join(args.out,'T2_sn_murcko_unique.pkl'),'wb') as pklf:
        pickle.dump(t2_murcko_unique,pklf)
    logger.info('written pickle file for unique Murcko scaffolds')

    #join back to original scaffolds
    t2_scaff = t2_standardized.merge(t2_murcko_unique[['murcko_scaff_smiles','sn_scaff_smiles']],on='murcko_scaff_smiles')
    t2_scaff.to_csv(os.path.join(args.out,'T2_sn_scaff.csv'))
    logger.info('written scaffold assignment')

    #now need to merg the T5 table, to get scaffold assignments to unique descriptors IDs
    t5 = pd.read_csv(os.path.join(args.infolder,'results_tmp','descriptors','mapping_table_T5.csv'))
    t2_joined = t2_scaff.merge(t5,on='input_compound_id')

    # in rare cases of collisions 8diffeent scaffolds, but ientical descriptors we brea the ties by alpabetic precedence
    t6_scaff = t2_joined.groupby('descriptor_vector_id')[['murcko_scaff_smiles','sn_scaff_smiles']].min()

    #identify collisions: different scaffold, same descriptor vector
    nunique_murcko = t2_joined.groupby('descriptor_vector_id')['murcko_scaff_smiles'].nunique()
    collisions_murcko = t2_joined[t2_joined['descriptor_vector_id'].isin(nunique_murcko[nunique_murcko > 1].index)].sort_values('descriptor_vector_id')
    collisions_murcko.to_csv(os.path.join(args.out,'collisions_murcko.txt'),sep='\t')
    logger.info('Colliding Murcko scaffolds observed on {0} records and written out'.format(len(collisions_murcko)) )

    nunique_sn = t2_joined.groupby('descriptor_vector_id')['sn_scaff_smiles'].nunique()
    collisions_sn = t2_joined[t2_joined['descriptor_vector_id'].isin(nunique_sn[nunique_sn > 1].index)].sort_values('descriptor_vector_id')
    collisions_sn.to_csv(os.path.join(args.out,'collisions_sn.txt'),sep='\t')
    logger.info('Colliding SN scaffolds observed on {0} records and written out'.format(len(collisions_sn)) )

    #now we need to join to final file with continuous descriptor vector Ids
    t2_t11 = pd.read_csv(os.path.join(args.infolder,'results','T11.csv'))
    t2_t11_joined = t2_t11.join(t6_scaff,on='descriptor_vector_id')

    #now we need to generate the hashes for the scaffolds
    key = params['key']['key'].encode("ASCII")
    nfolds = params['lsh']['nfolds'] 
    t2_t11_joined['sn_scaff_fold'] = t2_t11_joined['sn_scaff_smiles'].apply(lambda x : hashed_fold_scaffold(x, secret = key, nfolds=nfolds))
    t2_t11_joined['murcko_scaff_fold'] = t2_t11_joined['murcko_scaff_smiles'].apply(lambda x : hashed_fold_scaffold(x, secret = key, nfolds = nfolds))
    t2_t11_joined.to_csv(os.path.join(args.out,'T2_T11_scaff_folds.csv'))
    logger.info('written out csv file with fold info')

    #now create numpy arrays for folds
    fold_sn = np.zeros(len(t2_t11_joined['cont_descriptor_vector_id']))
    fold_sn[t2_t11_joined['cont_descriptor_vector_id']] = t2_t11_joined['sn_scaff_fold']
    fold_murcko = np.zeros(len(t2_t11_joined['cont_descriptor_vector_id']))
    fold_murcko[t2_t11_joined['cont_descriptor_vector_id']] = t2_t11_joined['murcko_scaff_fold']

    np.save(os.path.join(args.out,'murcko_scaff_folds.npy'),fold_murcko)
    np.save(os.path.join(args.out,'sn_scaff_folds.npy'),fold_sn)
    logger.info('written out numpy arrays with folds')

