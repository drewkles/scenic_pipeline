## Running the scenic pipeline. 

For introduction and instructions on downloading cisTarget databases:
https://rawcdn.githack.com/aertslab/SCENIC/701cc7cc4ac762b91479b3bd2eaf5ad5661dd8c2/inst/doc/SCENIC_Setup.html

To get parameters and help run

`$ Rscript scenic_pipeline -h`

From command line

`$ Rscript scenic_pipeline.R /data/expr_mat.rds data/cell_info.rds --cis_path /cisTarget_databases --nCores 10 --title ScenicTestRun`

Rscript to run SCENIC on human single cell data

positional arguments:
  expression_mat_rds    Path to expression matrix rds file (genes by cells)
                        Genes as rownames and Cells as colnames
  cell_info_rds         Path to cell info rds file (clusters with cells as
                        rownames)

optional arguments:
  -h, --help            show this help message and exit
  --cis_path CIS_PATH   Path to cis Target pathway directory(e.g.
                        ~/home/scenic/cisTarget_databases)
  --nCores NCORES       Number of cores to use
  --title TITLE         add project title name
  --min_counts_gene MIN_COUNTS_GENE
                        Filter by the total number of reads per gene. By
                        default it keeps only the genes with at least 6 UMI
                        counts across all samples
  --min_samples MIN_SAMPLES
                        Filter by the 