#!/usr/bin/env Rscript

library(fgsea)
library(optparse)
library(ggplot2)

# Functions ----
read.gmt <- function (file) 
{
  if (!grepl("\\.gmt$", file)[1]) {
    stop("Pathway information must be a .gmt file! Exit")
  }
  geneSetDB = readLines(file)
  geneSetDB = strsplit(geneSetDB, "\t")
  names(geneSetDB) = sapply(geneSetDB, "[", 1)
  geneSetDB = lapply(geneSetDB, "[", -1:-2)
  geneSetDB = lapply(geneSetDB, function(x) {
    x[which(x != "")]
  })
  return(geneSetDB)
}

# Options ----
option_list <- list( 
  make_option(c("-g", "--gmtfile"), action="store", default=NULL,
              help="Download the .gmt file"),
  make_option(c("-r", "--rnkfile"), action="store", default=NULL, 
              help="Download the .rnk file"),
  make_option(c("-n", "--nperm"), type="integer", default=1000, 
              help="Number of permutations to compute p-value",
              metavar="number"),
  make_option(c("-b", "--biggest"), type="integer", default=500, 
              help="Maximum genes in gene set",
              metavar="number"),
  make_option(c("-s", "--smallest"), type="integer", default=15, 
              help="Minimum genes in gene set",
              metavar="number"),
  make_option(c("-d", "--directory"), action ="store", default=".",
              help="Directory containing fgsea results, enrichment plots and fgsea Leading Edges Details"),
  make_option(c("-l", "--leadingedge"), action ="store", default=".",
              help="Subdirectory containing fgsea Leading Edges Details"),
  make_option(c("-f", "--filetype"), action ="store", default="csv",
              help="Determines file type for fgsea results (csv by default) and fgsea Leading Edges Details (tsv by default)"),
  make_option("--gseaParam", type="integer", default=0,
              help="GSEA parameter value, all gene-level stats are raised to the power of gseaParam before calculation of ES."),
  make_option(c("-j","--cores"), type="integer", default=1,
              help="Parallelization.")
)
opt <- parse_args(OptionParser(option_list = option_list))


# Script ----
# Check gseaParam
if(!opt$gseaParam %in% c(0,1))
  stop("gseaParam must be 0 or 1! Exit")

# Download pathway file
gmt <- read.gmt(paste0(opt$gmtfile))

# Download rank file
rnk <- read.table(opt$rnkfile, header = F, col.names = c("featureID", "value"))
rnk_vector <- unlist(rnk$value)
names(rnk_vector) <- rnk$featureID

# Create fGSEA an Leading_Edge directory
if(!dir.exists(opt$directory)) dir.create(opt$directory, recursive = T, showWarnings = FALSE)
if(!dir.exists(opt$leadingedge)) dir.create(opt$leadingedge, recursive = T, showWarnings = FALSE)

# Define output filename
outfile <- paste0(opt$rnkfile, ".", opt$gmtfile, ".fGSEA")

# Run fGSEA
# Add condition: if(simple or multilevel)
fgseaRes <- fgsea::fgseaSimple(pathways = gmt,
	stats    = rnk_vector,
	nperm    = opt$nperm,
	minSize  = opt$smallest,
	maxSize  = opt$biggest,
	gseaParam = opt$gseaParam, #https://www.biostars.org/p/263074/
	nproc = opt$cores,
)

# Check fGSEA output
if(!is.null(fgseaRes)){
  message(" -- fgseaSimple done")
} else {
  stop(" -- Error in fgseaSimple! Exit")
}

# Save results
if (opt$filetype == "excel") {
  outfile_xlsx <- paste0(opt$directory, "/", outfile,".xlsx")
  message(" -- saving to: ", outfile_xlsx)
  write.xlsx(subset(fgseaRes, select = -c(leadingEdge)), outfile_xlsx, col.names = T, row.names = F)
} else {
  outfile_tab <- paste0(opt$directory, "/", outfile,".tab")
  message(" -- saving to: ", outfile_tab)
  write.table(subset(fgseaRes, select = -c(leadingEdge)), outfile_tab, col.names = T, quote = F, row.names = F, sep="\t")
}

