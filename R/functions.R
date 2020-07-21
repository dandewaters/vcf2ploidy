#' @title Read in VCF files
#'
#' @description Read in a Variant Calling Format file.
#'
#' @details Throws an error if file cannot be found.
#'
#' @param filename A character string of the data file path.
#' @param skip_lines Number of metadata lines to skip before reading in the VCF file. Metadata lines typically start with "##". Default is 10.
#'
#' @return This function returns a data frame table of the VCF data with the metadata and the first 9 columns removed.
#'
#' @examples
#' \dontrun{read_VCF("./inst/extdata/example.vcf")}
#'
#' @importFrom readr read_delim
#' @importFrom dplyr tbl_df
#' @importFrom dplyr rename
#' @importFrom dplyr select
#' @importFrom magrittr "%>%"
#'
#' @export
read_VCF <- function(filename, skip_lines=10){
  if(!file.exists(filename))
    stop("file '", filename, "' does not exist")

  # Read in VCF file
  VCF_file <- read_delim(file=filename, delim="\t", col_types = cols(.default = "c"), skip=skip_lines)

  # Remove first 9 columns
  VCF_file <-
    VCF_file %>%
    # Removing the pound symbol on the "#CHROM" column so the select() arguments don't get commented out
    rename("CHROM" = "#CHROM") %>%
    dplyr::select(-CHROM, -POS, -ID, -REF, -ALT, -QUAL, -FILTER, -INFO, -FORMAT)

  VCF_file
}



#' @title Extract Locus Information from VCF Data for HAD Format
#'
#' @description A function used to fill the HAD converted data frame with the number of reads on each allele for heterozygous loci. Not exported.
#'
#' @param locus A character string of the locus information from a column of the VCF data frame.
#' @param remove_double_hets Logical for determining if double heterozygous loci should be treated as missing information.
#'
#' @return A numeric vector of the number of reads for each allele at heterozygous loci.
#'
#' @examples
#' \dontrun{analyze_locus("1/0:8:0,6,0,2", remove_double_hets=TRUE)}
analyze_locus <- function(locus, remove_double_hets=FALSE){
  # If locus is heterozygous, extract the number of reads from
  if(startsWith(locus, "0/1") | startsWith(locus, "1/0")){
    # Isolate number of reads from locus data string
    locus <- strsplit(locus, ":")[[1]][3]
    locus <- strsplit(locus, ",")[[1]]
    # Cast to numeric vector and sort the numbers of reads
    locus <- as.numeric(locus)
    locus <- locus[ordered(-locus)]

    # Checks for double heterozygotes by looking for more than 2 nonzero reads
    if(remove_double_hets == TRUE & (locus[3] != 0 | locus[4] != 0)){return(c(NA, NA))}
    # Return NA if there are no reads for any allele or there are only reads for one allele
    # (iPyrad falsely labeling a homozygote or missing locus as a heterozygote)
    else if(locus[2] == 0){return(c(NA, NA))}
    # Return the 2 highest number of reads otherwise
    else{return(locus[1:2])}
  }
  # Return NA if locus is not heterozygous
  else{return(c(NA, NA))}
}


#' @title Extract Locus Information from VCF Data for Colony Format
#'
#' @description A function used to fill the Colony converted data frame. Not exported.
#'
#' @param locus A character string of the locus information from a column of the VCF data frame.
#'
#' @return A numeric vector of the number of reads for each allele at heterozygous loci.
#'
#' @examples
#' \dontrun{analyze_locus("1/0:8:0,6,0,2")}
analyze_locus_colony <- function(locus){
  # Return -9, -9 for missing locus
  if(startsWith(locus, "./.")){return(c(-9, -9))}

  else{
    locus <- strsplit(locus, ":")[[1]][1]
    locus <- strsplit(locus, "/")[[1]]
    locus <- as.numeric(locus)
    return(locus+1)
  }
}



#' @title Convert Raw VCF Data to HAD Format
#'
#' @description Convert tbl_df output from read_VCF() to a Heterozygous Allele Depth (HAD) format, to be read in by GBS2Ploidy.
#'
#' @param filename A character string of the data file path.
#' @param skip_lines Number of metadata lines to skip before reading in the VCF file. Metadata lines typically start with "##". Default is 10.
#' @param remove_double_hets Logical for determining if double heterozygous loci should be treated as missing information. Should fix issues with GBS2Ploidy falsely labeling triploids.
#'
#' @return A data frame in Heterozygous Allele Depth (HAD) format
#'
#' @examples
#' \dontrun{VCF2HAD(./inst/extdata/example.vcf)}
#' \dontrun{VCF2HAD(./inst/extdata/example.vcf, remove_dobule_hets=TRUE)}
#'
#' @export
VCF2HAD <- function(filename, skip_lines=10, remove_double_hets=FALSE){
  # Read in VCF file
  VCF_df <- read_VCF(filename, skip_lines)
  # Initialize an empty data frame
  HAD_df <- data.frame()
  # Get the number of rows and columns
  num_rows    <- dim(VCF_df)[1]
  num_columns <- dim(VCF_df)[2]

  for(i in 1:num_columns){
    # Grab column from VCF data frame
    col <- VCF_df[[i]]
    col_name <- colnames(VCF_df)[i]

    # Initialize an empty data frame
    new_col <- data.frame()

    # Loop through column and extract heterozygote information
    for(j in 1:num_rows){
      # Grab locus information from column
      locus <- col[j]
      # Analyze locus information and get numbers of reads for each allele to add to column
      locus <- analyze_locus(locus, remove_double_hets)
      # Add locus data to new column
      new_col <- rbind(new_col, locus)
    }

    # Get correct column names for new columns
    colnames(new_col) <- c(col_name, col_name)
    # Add new columns to final data frame
    if(dim(HAD_df)[1] == 0 & dim(HAD_df)[2] == 0){HAD_df = new_col}
    else{HAD_df <- cbind(HAD_df, new_col)}
  }
  return(HAD_df)
}



#' @title Convert Raw VCF Data to Colony Format
#'
#' @description Convert tbl_df output from read_VCF() to a format that can be read in by Colony.
#'
#' @param filename A character string of the data file path.
#' @param skip_lines Number of metadata lines to skip before reading in the VCF file. Metadata lines typically start with "##". Default is 10.
#' @param out_filename A character string of the resulting converted file path.
#'
#' @return NULL
#'
#' @examples
#' \dontrun{VCF2HAD(./inst/extdata/example.vcf)}
#'
#' @export
VCF2colony <- function(filename, skip_lines=10){
  # Read in VCF file
  VCF_df <- read_VCF(filename, skip_lines)
  # Initialize an empty data frame
  colony_df <- data.frame()

  # Get the number of rows and columns
  num_rows    <- dim(VCF_df)[1]
  num_columns <- dim(VCF_df)[2]

  for(i in 1:num_columns){
    # Grab column from VCF data frame
    col <- VCF_df[[i]]
    col_name <- colnames(VCF_df)[i]

    # Initialize an empty data frame
    new_col <- data.frame()

    # Loop through column and extract heterozygote information
    for(j in 1:num_rows){
      # Grab locus information from column
      locus <- col[j]
      # Analyze locus information and get numbers of reads for each allele to add to column
      locus <- analyze_locus_colony(locus)

      # Add locus data to new column
      new_col <- rbind(new_col, locus)
    }

    # Get correct column names for new columns
    colnames(new_col) <- c(col_name, col_name)
    # Add new columns to final data frame if data frame isn't empty and -
    if(dim(colony_df)[1] == 0 & dim(colony_df)[2] == 0){colony_df = new_col &
    # - if all loci aren't missing
      all(new_col[1] == 9) & all(new_col[2] == 9)}
    else{colony_df <- cbind(colony_df, new_col)}
  }
  # Write the converted data to a text file
  write.table(colony_df, file=out_filename, quote=FALSE, sep="\t", col.names=TRUE)
}



#' @title Estimate ploidy directly from VCF files
#'
#' @description Read in a VCf file, convert to HAD format, and estimate ploidy in one function
#'
#' @param filename A character string of the data file path.
#' @param skip_lines Number of metadata lines to skip before reading in the VCF file. Metadata lines typically start with "##". Default is 10.
#' @param remove_double_hets Logical for determining if double heterozygous loci should be treated as missing information. Should fix issues with GBS2Ploidy falsely labeling triploids.
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
#' @return VCF2Ploidy returns a list with three components:
#' @return pp A matrix with assignment probabilities for each individual (rows) to each group (columns); the first column gives the ids provided by the user. Only individuals that were not part of the training set are included.
#' @return pcwghts A matrix with the variable loadings (PC weights) from the ordination of residual heterozygosity and allelic proportions. Columns correspond with PCs in ascending order (i.e., the PC with the largest eigenvalue is first).
#' @return pcscrs A matrix of PC scores from the ordination of residual heterozygosity and allelic proportions. Columns correspond with PCs in ascending order (i.e., the PC with the largest eigenvalue is first).
#'
#' @examples
#' \dontrun{VCF2Ploidy("./example.vcf")}
#' \dontrun{VCF2Ploidy("./example.vcf", props=c(0.25, 0.5, 0.75))}
#'
#' @importFrom gbs2ploidy estprops
#' @importFrom gbs2ploidy estploidy
#'
#' @export
VCF2Ploidy <- function(filename, skip_lines=10, remove_double_hets=FALSE,
                       props=c(0.25, 0.33, 0.5, 0.66, 0.75), mcmc.nchain=2,
                       mcmc.steps=10000, mcmc.burnin=1000, mcmc.thin=2,
                       train=FALSE, pl=NA, set=NA, nclasses=2, pcs=1:2){

  # Open file and convert to HAD format
  had_df <- VCF2HAD(filename, skip_lines=10, remove_double_hets=FALSE)

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
