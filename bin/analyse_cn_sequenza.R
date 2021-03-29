#!/usr/bin/env Rscript
# Load libs:
if (!require(sequenza)) stop("Package 'sequenza' missing\n.")

args <- commandArgs(TRUE)
print(args)
input <- args[1]
output_prefix <- args[2]
gender <- args[3]

params_list <- list("input" = input, "output_prefix" = output_prefix )
# Function:
sequenzaAnalysis <- function(input,
                             output_prefix,
                             window=1e6,
                             overlap=1,
                             gamma=80,
                             kmin=10,
                             min_reads=40,
                             min_reads_normal=10,
                             min_reads_baf=1,
                             max_mut_types=1,
                             breaks=NULL,
                             assembly="hg19",
                             weighted_mean=TRUE,
                             normalization_method="mean",
                             is_female=TRUE,
                             segment_filter=3e6,
                             ratio_priority=FALSE,
                             method="baf",
                             low_cell=0.1,
                             up_cell=1,
                             low_ploidy=1,
                             up_ploidy=7,
                             CNt_max=20) {
  
  # Define chromosomes to analyse (note these will subset to those that
  # are available for sequenza:
  chr_list <- c(1:22, "X")
  if (gender != "XX")
    chr_list <- c(chr_list, "Y")
  chr_list <- c(chr_list, paste0("chr", chr_list))
  
  
  # Extract sequenza data for model:
  cat(sprintf("- Starting analysis for %s\n", input))
  cat("- Calculating gc-stats\n")
  gc_stats   <- gc.sample.stats(input)
  
  
  cat("- Loading data\n")
  modDat <- sequenza.extract(input,
                             window=window,
                             overlap=overlap,
                             gamma=gamma,
                             kmin=kmin,
                             min.reads=min_reads,
                             min.reads.normal=min_reads_normal,
                             min.reads.baf=min_reads_baf,
                             max.mut.types=max_mut_types,
                             chromosome.list=chr_list,
                             breaks=breaks,
                             assembly=assembly,
                             weighted.mean=weighted_mean,
                             normalization.method=normalization_method,
                             parallel=8,
							 gc.stats=gc_stats)
  
  
  # Fit the model:
  cat("- Fitting the model\n")
  cells  <- seq(low_cell, up_cell, 0.01)
  plo    <- seq(low_ploidy, up_ploidy, 0.05)
  fit    <- sequenza.fit(modDat,
                         female=is_female,
                         segment.filter=segment_filter,
                         cellularity=cells,
                         ploidy=plo,
                         ratio.priority=ratio_priority,
                         method=method)
  
  # Export the data:
  cat("- Exporting results\n")
  outName <- basename(output_prefix)
  outDir  <- dirname(output_prefix)
  sequenza.results(modDat,
                   fit,
                   outName,
                   outDir,
                   female=is_female,
                   CNt.max=CNt_max,
                   ratio.priority=ratio_priority)
}

do.call(sequenzaAnalysis, params_list)

