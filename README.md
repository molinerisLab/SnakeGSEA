# SnakeGSEA

A Snakemake pipeline to perform fastGSEA after `bit_rnaseq` differential gene expression analysis.

## Installation

Clone the repository
```
git clone git@github.com:molinerisLab/SnakeGSEA.git SnakeGSEA
```

### Install dependencies 

Move to SnakeGSEA directory obtained with `git clone` and create the SnakeGSEA_Env environment with all dependencies
```
conda env create -f local/env/environment.yml -n SnakeGSEA_Env
```

### Download MSigDB data

SnakeGSEA uses gene set collections provided by MSigDB. Here the links for the download of the gene sets for [Human](http://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/2022.1.Hs/msigdb_v2022.1.Hs_files_to_download_locally.zip) and [Mouse](http://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/2022.1.Mm/msigdb_v2022.1.Mm_files_to_download_locally.zip)

## Usage

Move to `dataset/v1` directory and set the proper parameters in `config.yaml` file.
The parameters are documented in the file itself.

Pay attention to `MSIGDB_DIR` the variable contain the path that contain the downloaded MSigDB files.

To run the pipeline: 
```
snakemake -p -j N_CORES all
```
Put `N_CORES` to the number of CPUs you waunt to use.

## Output
SnakeGSEA pipeline produces the following files:

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

## Publish to a web server

```
snakemake -p -j 1 publish
```

After the publishg as been completed it is possible to interactively build GSEA plots using the [shinySea](https://github.com/molinerisLab/shinySea) app for each geneset.
In the html and xlsx ouput there are link that point to this app.


## SnakeGSEA from other pipelines

SnakeGSEA can be applied on a generic pre-ranked list of genes produced from your custom analysis workflows.

Create the tab separated `contrast.rnk` files:

```
1) gene_id        gene symbol (e.g., SOX1)
2) ranking_value  the numeric value considered for gene ordering (e.g., logFC)
```

in the current directory and then run 
```
snakemake -p -j N_CORES all
``` 

You do not have to modify the `config.yaml`. Eventual rnk files from bit_rna_seq pipelines will be ignored.

## Enrichment Score Weight

As reported in the [GSEA manual](https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideTEXT.htm#_GSEAPreranked_Page), the "classic" weight should be used:

> The GSEA PNAS 2005 paper introduced a method where a running sum statistic is incremented by the absolute value of the ranking metric when a gene belongs to the set. This method has proven to be efficient and facilitates intuitive interpretation of ranking metrics that reflect correlation of gene expression with phenotype. In the case of GSEAPreranked, you should make sure that this weighted scoring scheme applies to your choice of ranking statistic – i.e. the magnitude of the ranking metric is biologically meaningful. When in doubt, we recommend using a more conservative scoring approach by setting Enrichment statistic = ‘classic’.


