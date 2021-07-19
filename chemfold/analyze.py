import argparse
import logging
import sys
import analysis as cfa

#logging
shandler = logging.StreamHandler(stream=sys.stdout)
formatter = logging.Formatter('%(asctime)s -%(name)s - %(levelname)s - %(message)s')
shandler.setFormatter(formatter)
shandler.setLevel(logging.INFO)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
logger.addHandler(shandler)


parser = argparse.ArgumentParser(description="Suite of different folding methods for chemical structures suitable for privacy-preserving federated learning and analysis scripts.")
#parser.add_argument("--x", help="Descriptor file (matrix market or numpy)", type=str, required=True)
parser.add_argument("--inp", help="Folder with input files, e.g., descriptor file (matrix market or numpy, distances, y matrix, sparse chem model folder + json summary,Matrix of ECFP features [compounds x ecfp], .mtx file)", type=str, required=True)
parser.add_argument("--analysis", help="which analysis to run options are: [all, performance, imbalance, similarity]", choices=['all','performance','imbalance','similarity'], default='all')
parser.add_argument("--out", help="Output directory for the clusters file (.npy)", type=str, required=True)

#performance
parser.add_argument("--psc", help="Path to sparsechem models sparsechem/models", type=str, required=False)
parser.add_argument("--baseline_prefix", help="Baseline prefix", type=str, default='random')

#similarity
parser.add_argument("--maxsample", help="Maximal number of compound pairs to sample", type=int, default=10000000)
parser.add_argument("--batchsize", help="Number of compound pairs to precessed in one batch", type=int, default=50000)
parser.add_argument("--numbins", help="Number of bins of the histogram", type=int, default=10)
parser.add_argument("--precision", help="Maximal tolerated change of intra fold ratio per similarity bin between to batches to reach convergence",type=float,default=0.02)
parser.add_argument("--minpop", help="Minimal population of intra-fold samples per similarity bin to required to reach convergence",type=int, default=3)
parser.add_argument("--rseed",help="random seed to use", type=int, default = 123456)

args = parser.parse_args()
logger.info(args)


#call of analysis functions

if args.analysis:
	if args.analysis == 'imbalance' or args.analysis == 'all':
	#balance
		cfa.balance(args.inp,args.out)


	if args.analysis == 'performance' or args.analysis == 'all':
		#performance
		#before a model has to be trained with the baseline and the alternative fold splitting using SparseChem
		# We recognize the conditions by the prefix given as option during the training job submission. This prefix is found in the model json files. All  fold models under the same conditions are expected to have the same prefix. We now declare which condition is the baseline.

		#TODO: add an example call of sparsechem train
        #cd ./examples/chembl
		#python train.py --prefix [specify folding method here]


		#Analysis of the effect of the fold splits on model performance, using random folding as the baseline. It assumes you have created the models for 5 validation folds and reads in the json files created by sparsechem. Each folding scheme should use a diffrent sparsechem --prefix to distinguish them.
		cfa.performance(args.psc, args.baseline_prefix, args.out)


	if args.analysis == 'similarity' or args.analysis == 'all':
		#similarity
		#This analyzes the fraction of randomly chosen compound pairs per similariuy bins that come from the same fold. Ideally the fraction on intra-fold pairs for high similarity bins should be as high as possible.
		cfa.similarity(args.inp, args.out,args.maxsample,args.batchsize,args.numbins,args.precision,args.minpop,args.rseed)
