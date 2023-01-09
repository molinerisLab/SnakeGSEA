# SnakeGSEA

A Snakemake pipeline to perform fastGSEA after `bit_rnaseq` differential gene expression analysis.

SnakeGSEA uses gene set collections provided by MSigDB for [Human](http://www.gsea-msigdb.org/gsea/msigdb/collections.jsp) and [Mouse](https://www.gsea-msigdb.org/gsea/msigdb/mouse/collections.jsp).


## Usage

Clone the repository
```
git clone git@github.com:molinerisLab/SnakeGSEA.git SnakeGSEA
```

Move to SnakeGSEA directory and create the SnakeGSEA_Env environment with all dependencies
```
conda env create -f local/env/environment.yaml -n SnakeGSEA_Env
```

Move to `dataset/v1` directory and set the proper parameters in `config.yaml` file. Then run
```
snakemake -p -j N all
```

## Output
SnakeGSEA pipeline priduces the following files:

- multi_GSEA.link.header_added.xlsx
```
1       contrast          
2       msigdb_type
3       pathway
4       pval
5       padj
6       ES
7       NES
8       nMoreExtreme
9       size
10      shinySea_link
11      pathway_link
```

- multi_GSEA.html

Browser version of the multi_GSEA.link.header_added.xlsx table followed by the barplot visualization of top 5 and bottom 5 gene sets ordered according to NES value. Darker bars mean a significant association of the gene set, lighter bars are related to non significant associations. A bar plot is produced for each contrast and MSigDB collection combination.


## SnakeGSEA from other pipelines

SnakeGSEA can be applied on a generic pre-ranked list of genes produced from your custom analysis workflows.
Create the `contrast.rnk` files and list all the contrast names in the `config.yaml` file.
