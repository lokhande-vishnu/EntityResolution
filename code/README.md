# Code

## Overview
- `F-MWSP/`: Code for solving minimum weight set packing problem with column generation and flexible dual optimal inequalities. Code in this directory is available in matlab. A python version is in preparation. F-MWSP is sometimes referred as ZEUS in the scripts.
- `dedupe/`: A python library which provides APIs for blocking and scoring in the entity resolution pipeline. We provide a customized dedupe with blocking and scoring functions separated.
- `patent_sample/`, `csv_example/`, `affiliations/`, `settlements/`, `music20k/`: The datasets on which experiments were conducted. 

## Installations
The code requires a customized dedupe. Please install it in the following way. Refer to the [dedupe](https://github.com/dedupeio/dedupe) page for further assistance.

It is recommended to install dedupe in a [virtual environment](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html).

```bash
conda create --name myenv python=3.6
source activate myenv
cd dedupe
pip install "numpy>=1.9"
pip install -r requirements.txt
cython src/*.pyx
pip install -e .
```

The F-MWSP code requires IBM CPLEX. Please refer to this [documentation](https://www.ibm.com/support/knowledgecenter/en/SSSA5P_12.9.0/ilog.odms.studio.help/Optimization_Studio/topics/COS_installing.html) for it's installation.


## QuickStart
Here we provide an example for the `patent_sample` dataset. Follow the same procedure for other datasets.

1. Perform blocking and scoring on _patent_sample_.
Adjust the experiment name within the _main.py_ file.
```bash
cd patent_sample
python main.py
cd ..
```
The code takes about a minute to execute.
This step creates `F_sample.mat` which contains a sparse graph with nodes representing the observations and the edges bearing a similarity measure between the nodes.

2. Run F-MWSP algorithm.
The code for this step is currently available in matlab.
Adjust the dataset name within the _main.m_ file appropriately.
```bash
cd F-MWSP/examples/entity_resolution/exec/
matlab --nodesktop --nosplash --nodisplay -r "main; exit"
cd ../../../../
```
This step generates `H_sample.mat` file in the dataset directory and contains a list of tuples of the form (observationID, clusterID).

3. Run the evaluation script.
```bash
cd patent_sample
python evaluation.py
cd ..
```
In this step we match the true cluster id's provided in the dataset with the cluster id's obtained by running F-MWSP algorithm. We compare the performance against hierarchical clustering baseline.

## Notes
The _affiliations_ and _music20k_ dataset take very long to run. Please use _main_low_memory.m_ file for these datasets.

Every dataset directory contains a _results_ folder. This contains the outputs of each step of Quickstart.