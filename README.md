# ChemFold

ChemFold provides several methods for computing train-validation-test splits, designed for both ordinary ML and federated ML tasks involving small molecules.
Following methods are included:
* Random split
* Sphere exclusion clustering based split
* Locality sensitive hashing (LSH) based split
* Scaffold trees

## Installation of ChemFold

ChemFold can be installed by cloning its repository. First we need two dependencies rdkit (at least 2021.03.3) and cython:
```bash
conda install -c conda-forge rdkit ## at least 2021.03.3
pip install cython
pip install -e .
```

## Example of running scaffold-network
```bash
python -m chemfold.scaffold_network --infolder <root folder of melloddy_tuner_output> \
                                    --out <folder to write output to> \
                                    --params_file <location of parameters.json>
```
For the `parameters.json` must contain values for `["key"]["key"]` and `["lsh"]["nfolds"]` entries.
This method also works in federated setting and guarantees consistent folding results across multiple parties.

## Example of running random split
Here is an example executing random split:
```bash
python -m chemfold.random_split --infolder <root folder of melloddy_tuner_output> \
                                --out <folder to write output to> \
                                --params_file <location of parameters.json>
```
This method also works in federated setting and guarantees consistent folding results across multiple parties.

## Example of running sphere exclusion clustering
```bash
python -m chemfold.sphere_exclusion --x chembl_23mini_x.npy \
                                    --dists 0.2 0.4 0.6 \
                                    --out chembl_23mini_clusters.npy
```
where `chembl_23mini_x.npy` is the sparse matrix (either scipy.sparse saved in .npy or matrix market, ending with `.mtx`) of descriptor values such as ECFPs. This command does sequence of sphere exclusion clusterings for distances 0.2, 0.4 and 0.6.

Note this implementation of sphere exclusion can be used only for non-federated setting and does **not** create consistent results if used by different parties.

## Example of getting Locality sensitive hashing
MELLODDY TUNER (see next section) directly offers a LSH-based fold splitting. The respective fold-file can be found in `files_4_ml/T11_fold_vector.npy`.

## Example of running analyses
After calculating the fold using one of the above functions analyses can be run. The performance difference (needs SparseChem (https://github.com/melloddy/SparseChem) runs before), label and data imbalance as well as similarity of chemical substances within one fold can be analyzed.
```bash
python -m chemfold.analyze --inp <ChemFold folder including the output files from the folding methods, e.g., sn_scaff_folds.npy> \
                           --out <folder to write output to> \
                           --psc <location of SparseChem folder inclusing the models folder>
                           --baseline_prefix <prefix used in SparseChem that corresponds to the baseline folding methods>
                           --analysis <which analyses to run, choose from 'all', 'performance', 'imbalance' and 'similarity'>
```


## Inputs from MELLODDY-TUNER
ChemFold uses as input processed molecular data, for which can be obtained from [MELLODDY-TUNER](https://github.com/melloddy/MELLODDY-TUNER).

Here we outline basic steps to get MELLODDY-TUNER installed and data preparation executed, for in-depth guide you can check MELLODDY-TUNER's manual.

### Clone MELLODDY-TUNER
First, clone the git repository from the MELLODDY gitlab repository:

```bash
git clone -b release/1.0 https://github.com/melloddy/MELLODDY-TUNER.git
```


### Create enviroment
Create your own enviroment from the given yml file with:

```bash
cd MELLODDY-TUNER
conda env create -f melloddy_pipeline_env.yml
```

Activate the environment:
```bash
conda activate melloddy_pipeline
```
 
### Install MELLODDY-TUNER
You have to install the melloddy-tuner package with pip:

```
pip install -e .
```


## 0. Prepare Input Files for MELLODDY-TUNER
See specifications in MELLODDY-TUNER [README.md](https://github.com/melloddy/MELLODDY-TUNER/blob/master/README.md) under section Prepare Input Files.


## 1. Run Data Prepration Script
```bash
python bin/prepare_4_melloddy.py \
--structure_file {path/to/your/structure_file_T2.csv} \
--activity_file {/path/to/your/activity_data_file_T4.csv} \
--weight_table {/path/to/your/weight_table_T3.csv} \
--config_file {/path/to/the/distributed/parameters.json} \
--output_dir {path/to/the/output_directory} \
--run_name {name of your current run} \
--number_cpu {number of CPUs to use} \
--ref_hash {path/to/the/provided/ref_hash.json} \
```

Example for ChEMBL:
```bash
python bin/prepare_4_melloddy.py \
--structure_file {path/to/your/chembl25_T2.csv} \
--activity_file {/path/to/your/chembl25_T4.csv} \
--weight_table {/path/to/your/chembl25_T3.csv} \
--config_file {/tests/structure_preparation_test/example_parameters.json} \
--output_dir {path/to/the/output_directory} \
--run_name {name of your current run} \
--number_cpu {number of CPUs to use} \
--ref_hash {/tests/structure_preparation_test/ref_hash.json}



