# Andrew Elkins
# 2019-07-16
# scenic_pipeline.R
# Running the scenic pipeline

# Make sure to load:
# module load R/3.5.0
# module python/3.7.2
# Rscript scenic_pipeline.R expression_matrix_rds cell_info_rds --cis_path --nCores 10 --title ScenicSCBrain 

# Requirements
suppressPackageStartupMessages(library(Seurat))
library(argparse)
suppressPackageStartupMessages(library(AUCell))
library(SCENIC)
library(reticulate)
# Required for running on hoffman
reticulate::use_python("/u/local/apps/python/3.7.2/bin/python3", required=TRUE)

options(stringsAsFactors=FALSE)


getArgs <- function(){
  parser <- ArgumentParser(description='Rscript to run SCENIC on human single cell data')
  parser$add_argument('expression_mat_rds', help=paste0('Path to expression matrix rds file (genes by cells)',
                                  ' Genes as rownames and Cells as colnames'),
                      default=NULL )
  parser$add_argument('cell_info_rds', 
                      help=paste0('Path to cell info rds file (clusters with cells as rownames)'),
                      default=NULL)
  parser$add_argument('--cis_path', 
                      help=paste0('Path to cis Target pathway directory(e.g. ~/home/scenic/cisTarget_databases)'),
                      default='cisTarget_databases')
  parser$add_argument('--nCores', 
                      help=paste0('Number of cores to use'), type="integer", default=1)
  parser$add_argument('--title', help=paste0('add project title name'),default="SCENIC Exp")
  
  parser$add_argument('--min_counts_gene', help=paste0('Filter by the total number of reads per gene. ',
                                                       'By default it keeps only the genes with at least',
                                                       ' 6 UMI counts across all samples '),
                      default=NULL)
  parser$add_argument('--min_samples', 
                      help=paste0('Filter by the number of cells in which the gene is detected'),
                      default=NULL)
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
  message("\nInitializing...")
  createDirs()
  message("  ...reading in files")
  expr.mat <- readRDS(file=args$expression_mat_rds)
  cell.info <- readRDS(file=args$cell_info_rds)
  # Catch and formating
  if(ncol(expr.mat) != nrow(cell.info)){
    stop('Number of cells in two files does not match!')
  }
  if(grepl("\\.",colnames(expr.mat)[1])){
    colnames(expr.mat) <- sapply(colnames(expr.mat), function(x) gsub("\\.", "-", x))
  }
  if(grepl("\\.",rownames(cell.info)[1] )){
    rownames(cell.info) <- sapply(rownames(cell.info), function(x) gsub("\\.", "-", x) )
  }
  if( !all.equal(colnames(expr.mat),rownames(cell.info)) ){
    stop('cell names from expression matrix and cell info do not match!')
  }

  colnames(cell.info)[1] <- "cluster"
  cell.info$nGene <- apply(expr.mat, 2, function(c) sum(c!=0))
  cell.info$nUMI <- apply(expr.mat, 2, sum)
  saveRDS(cell.info, file="int/cellInfo.Rds") # save cell info for later

  data(defaultDbNames)
  dbs <- defaultDbNames[["hgnc"]]
  s.title <- ifelse(!is.null(args$title), "Scenic Analysis", args$title)

  scenic.options <- initializeScenic(org="hgnc", dbDir=args$cis_path,
                                     datasetTitle=s.title, nCores=args$nCores)
  scenic.options@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
  saveRDS(scenic.options, file="int/scenicOptions.Rds")

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
  genes.kept <- geneFiltering(expr.mat, scenicOptions=scenic.options,
                             minCountsPerGene=3*0.1*ncol(expr.mat),
                             minSamples=ncol(expr.mat) * 0.1)
  saveRDS(genes.kept, file="int/1.1_genesKept.Rds")
  message("   Total kept genes: ", genes.kept)

  message("   Running correlation matrix")
  expr.mat.filtered <- expr.mat[genes.kept, ]
  corr.mat <- cor(t(expr.mat.filtered), method="spearman")
  saveRDS(corr.mat, file=getIntName(scenic.options, "corrMat"))

  message("   exporting for GRNBoost!")
  exportsForArboreto(expr.mat.filtered, scenicOptions=scenic.options)
  message("\nGRNBoost")
  get.pyfile <- list.files(pattern = "grn_boost.py$", recursive = TRUE)
  if(is.null(get.pyfile)){
    stop("Cant find grn_boost.py!")
  }
  py_run_file(get.pyfile)

  grnboost.output <- importArboreto("int/1.2_grnoutput.txt")
  colnames(grnboost.output) <- c("TF", "Target", "weight")
  saveRDS(grnboost.output, file="int/1.4_GENIE3_linkList.Rds")
  message("\nCreating coexpression newtork modules")
  expr.mat.log <- log2(expr.mat + 1)
  runSCENIC_1_coexNetwork2modules(scenic.options)
  message("\nCreating regulons")
  runSCENIC_2_createRegulons(scenic.options)
  message("\nScoring cells")
  runSCENIC_3_scoreCells(scenic.options, as.matrix(expr.mat.log), skipTsne=TRUE,
                         skipHeatmap=TRUE, skipBinaryThresholds=TRUE)
  
  # make default thresholds (binarize mat)
  cells.auc.thres <- NULL
  regulon.auc <- readRDS(file=getIntName(scenic.options, "aucell_regulonAUC"))
  cells.auc.thres <- AUCell_exploreThresholds(regulon.auc,
                                              smallestPopPercent=getSettings(scenic.options,"aucell/smallestPopPercent"),
                                              assignCells=TRUE, plotHist=FALSE,
                                              verbose=FALSE, nCores=args$nCores)
  
  saveRDS(cells.auc.thres, file=getIntName(scenic.options, "aucell_thresholds"))
  message("\nProcess Complete!")
  
}

if(!interactive()){
  args <- getArgs()
  main(args)
}






