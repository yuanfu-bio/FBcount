# FBcount
Tools for Feature Barcode Counting. The Lite version of FBranger contributed by collaborator Fengrui Sui.


## Install

To access the latest version of the code, clone and install the repository:

```bash
git clone git@github.com:yuanfu-bio/FBcount.git
```
For easiest use, [create a conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-with-commands) and install requirements:

```bash
conda create -n <env_name> python=3.11
conda install -c conda-forge jq
conda activate <env_name>
cd FBcount
pip install -r ./requirements.txt
```