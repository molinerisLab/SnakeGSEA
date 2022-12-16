# SnakeGSEA

A Snakemake pipeline to perform fastGSEA after `bit_rnaseq` differential gene expression analysis.

SnakeGSEA can also be applyed on a pre-ranked list of genes produced from other analysis workflows (create the `contrast.rnk` files and list all the contrast names in the `config.yaml` file).

## Usage

Clone the repository
```
git clone git@github.com:molinerisLab/SnakeGSEA.git
```

Create the SnakeGSEA_Env environment with all dependencies
```
conda env create -f local/env/environment.yml -n SnakeGSEA_Env
```

Set the proper parameters in `config.yaml` file and then run
```
snakemake -p -j 4 all
```
