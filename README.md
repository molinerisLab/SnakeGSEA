# SnakeGSEA

A Snakemake pipeline to perform fGSEA after differential gene expression analysis.

It can be applyed on a pre-ranked list of genes too (you just need to create the `contrast.rnk` files and list all the contrast names in the `config.yaml` file).

## Usage

1. Set the proper parameters in `config.yaml` file
2. Run `snakemake -p -j 4 all`
