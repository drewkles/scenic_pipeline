# Andrew Elkins
# 2019-07-16
# scenic_pipeline.R
# Running the scenic pipeline

# Make sure to do:
# module load R/3.5.0
# module python/3.7.2
# Rscript scenic_pipeline.R expression_matrix_rds cell_info_rds db_path --nCores 10 --title Scenic scBrain 


# Requirements
suppressPackageStartupMessages(library(Seurat))
library(argparse)
library(AUCell)
library(SCENIC)
library(reticulate)
# Required for running on hoffman
reticulate::use_python("/u/local/apps/python/3.7.2/bin/python3", required=TRUE)

options(stringsAsFactors=FALSE)


getArgs <- function(){
  parser <- ArgumentParser(description='Rscript to run SCENIC on human single cell data')
  parser$add_argument('expression_mat_rds', help=paste0('Path to expression matrix rds file (genes by cells)',
                                  ' Genes as rownames and Cells as colnames') )
  parser$add_argument('cell_info_rds', 
                      help=paste0('Path to cell info rds (clusters with cells as rownames)'))
  parser$add_argument('db_path', 
                      help=paste0('Put path to cis Target pathway,(e.g. ~/home/scenci/cisTarget_databases)'))
  parser$add_argument('--nCores', 
                      help=paste0('Number of cores to use'), type="integer", default=1)
  parser$add_argument('--title', help=paste0('add project title name'))
  
  parser$add_argument('--min_counts_gene', help=paste0('Filter by the total number of reads per gene.',
                                                       'By default it keeps only the genes with at least',
                                                       ' 6 UMI counts across all samples '))
  parser$add_argument('--min_samples', help=paste0('Filter by the number of cells in which the gene is',
                                                   'detected** (e.g. >0 UMI, or >1 log2(CPM)).',
                                                   'By default genes that are detected in at least 1%',
                                                   ' of the cells are kept'))
  parser$parse_args()
}

createDirs <- function(){
  if(!dir.exists("output")){
    dir.create("output")
  }
  if(!dir.exists("int")){
    dir.create("int")
  }
}

main <- function(args){
  # Run through main routine
  if(is.null(args)){
    stop('Need to pass arguments!')
  }
  message("Initializing...")
  createDirs()
  expr.mat <- readRDS(file=args$expression_mat_rds)
  cell.info <- readRDS(file=argscell_info_rds)
  
  # Catch and formating
  if(ncol(expr.mat) != nrow(cell.info)){
    stop('Number of cells in two files does not match!')
  }
  if(grepl("\\.",colnames(expr.mat)[1])){
    colnames(expr.mat) <- sapply(colnames(mat.expr), function(x) gsub("\\.", "-", x))
  }
  if(grepl("\\.",rownames(cell.info)[1] )){
    rownames(cell.info) <- sapply(rownames(cell.info), function(x) gsub("\\.", "-", x) )
  }
  if(all.equal(colnames(expr.mat),rownames(cell.info) )){
    stop('cell names from expression matrix and cell info do not match!')
  }
  
  colnames(cell.info)[1] <- "cluster"
  cell.info$nGene <- apply(expr.mat, 2, function(c) sum(c!=0))
  cell.info$nUMI <- apply(expr.mat, 2, sum)
  saveRDS(cell.info, file="int/cellInfo.Rds") # save cell info for later
  
  data(defaultDbNames)
  dbs <- defaultDbNames[["hgnc"]]
  s.title <- ifelse(!is.null(args$title), "Scenic Analysis", args$title)

  scenic.options <- initializeScenic(org="hgnc", dbDir=args$db_path,
                                     datasetTitle=s.title, nCores=nCores)
  scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
  saveRDS(scenicOptions, file="int/scenicOptions.Rds")
  
  message("  Saving scenic options")
  # min_counts_gene min_samples
  if(!is.null(args$min_counts_gene)){
    min.counts.gene <- args$min_count_gene
  }else{
    min.counts.gene <- 3*0.1*ncol(expr.mat)
  }
  if(!is.null(args$min_samples)){
    min.samples <- ncol(expr.mat) * 0.1
  } else{
    min.samples <- args$min_samples
  }
  
  message("...Filtering genes")
  genes.kept <- geneFiltering(expr.mat, scenicOptions=scenicOptions,
                             minCountsPerGene=3*0.1*ncol(expr.mat),
                             minSamples=ncol(expr.mat) * 0.1)
  saveRDS(genes.kept, file="int/1.1_genesKept.Rds")
  message("   Total kept genes: ", genes.kept)
  
  message("   Running correlation matrix")
  expr.mat.filtered <- expr.mat[genesKept, ]
  corrMat <- cor(t(expr.mat.filtered), method="spearman")
  saveRDS(corrMat, file=getIntName(scenicOptions, "corrMat"))
  
  message("   exporting for GRNBoost!")
  exportsForArboreto(expr.mat.filtered, scenicOptions=scenicOptions)
  message("\nGRNBoost")
  py_run_file("../scenic_pipeline/grn_boost.py")
  
  grnboost.output <- importArboreto("int/1.2_grnoutput.txt")
  colnames(grnboost.output) <- c("TF", "Target", "weight")
  saveRDS(grnboost.output, file="int/1.4_GENIE3_linkList.Rds")
  message("\nCreating coexpression newtork modules")
  expr.mat.log <- log2(small.nucseq.exp+1)
  runSCENIC_1_coexNetwork2modules(scenicOptions)
  message("\nCreating regulons")
  runSCENIC_2_createRegulons(scenicOptions)
  message("\nScoring cells")
  runSCENIC_3_scoreCells(scenicOptions, as.matrix(expr.mat.log),
                         skipTsne=TRUE, skipHeatmap=TRUE)
  
  message("\nProcess Complete!")
}

if(!interactive()){
  main(getArgs())
}







