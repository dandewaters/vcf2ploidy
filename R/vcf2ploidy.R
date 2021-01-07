#' @title Estimate ploidy directly from VCF files
#'
#' @description Read in a VCf file, convert to HAD format, and estimate ploidy in one function
#'
#' @param filename A character string of the data file path.
#' @param skip_lines A numeric of the number of metadata lines to skip over in the VCF file. If left null, metadata lines are skipped over automatically by the count_metadata_lines function. The count_metadata_lines function requires reading in the entire file, so if you have a large file and know the number of metadata lines in that file, you can save some run time by entering the number of metadata lines in this argument.
#' @param remove_double_hets Logical for determining if double heterozygous loci should be treated as missing information. Should fix issues with gbs2ploidy falsely labeling triploids.
#' @param props a vector containing valid allelic proportions given the expected cyotypes present in the sample.
#' @param mcmc.nchain number of chains for MCMC.
#' @param mcmc.steps number of post burnin iterations for each chain.
#' @param mcmc.burnin number of iterations to discard from each chain as a burnin.
#' @param mcmc.thin thinning interval for MCMC.
#' @param train a boolean specifying whether or not a training set with known ploidy should be used.
#' @param pl a vector of known ploidies with one entry per individual (use ‘NA’ for individuals with unknown ploidy); only used if train == TRUE.
#' @param set indixes for the training set; only used if train == TRUE.
#' @param nclasses the number of cyotypes expected.
#' @param pcs a vector giving the PC to use for DA.
#'
#' @return vcf2ploidy returns a list with three components:
#' @return pp A matrix with assignment probabilities for each individual (rows) to each group (columns); the first column gives the ids provided by the user. Only individuals that were not part of the training set are included.
#' @return pcwghts A matrix with the variable loadings (PC weights) from the ordination of residual heterozygosity and allelic proportions. Columns correspond with PCs in ascending order (i.e., the PC with the largest eigenvalue is first).
#' @return pcscrs A matrix of PC scores from the ordination of residual heterozygosity and allelic proportions. Columns correspond with PCs in ascending order (i.e., the PC with the largest eigenvalue is first).
#'
#' @examples
#' \dontrun{vcf2ploidy("./example.vcf")}
#' \dontrun{vcf2ploidy("./example.vcf", props=c(0.25, 0.5, 0.75))}
#'
#' @importFrom gbs2ploidy estprops
#' @importFrom gbs2ploidy estploidy
#'
#' @export
vcf2ploidy <- function(filename, skip_lines=NULL, remove_double_hets=FALSE,
                       props=c(0.25, 0.33, 0.5, 0.66, 0.75), mcmc.nchain=2,
                       mcmc.steps=10000, mcmc.burnin=1000, mcmc.thin=2,
                       train=FALSE, pl=NA, set=NA, nclasses=2, pcs=1:2){

  # Open file and convert to HAD format
  had_df <- vcf2had(filename, remove_double_hets=FALSE)

  # Separate the first and second alleles in the HAD data frame
  num_cols <- dim(had_df)[2]

  cov1_cols <- seq(1, num_cols, 2)
  cov2_cols <- seq(2, num_cols, 2)

  cov1 <- as.matrix(had_df[cov1_cols])
  cov2 <- as.matrix(had_df[cov2_cols])


  # Grab column names for "ids" argument in estploidy()
  ids_indeces <- seq(1, num_cols, 2)
  ids <- names(had_df)[ids_indeces]

  # Estimate allelic proportions
  prop <- estprops(cov1=cov1, cov2=cov2, props=props,
           mcmc.nchain=mcmc.nchain, mcmc.steps=mcmc.steps,
           mcmc.burnin=mcmc.burnin, mcmc.thin=mcmc.thin)

  # Calculate observed heterozygosity and depth of coverage from the allele count
  H <- apply(is.na(cov1)==FALSE, 2, mean)
  D <- apply(cov1+cov2, 2, mean, na.rm=TRUE)

  # Estimate ploidy
  estploidy(alphas=prop, het = H, depth = D, train=train, pl=pl, set=set,
            nclasses=nclasses, ids=ids, pcs=pcs)
}
