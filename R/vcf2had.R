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
    ## Isolate number of reads from locus data string
    # Split string by semicolon character
    locus <- strsplit(locus, ":")[[1]]
    # Get the index of the last element in the split locus vector,
    # which holds the number of nucleotide reads
    locus_length <- length(locus)
    locus <- locus[locus_length]
    # split up read numbers
    locus <- strsplit(locus, ",")[[1]]
    # Cast to numeric vector
    locus <- as.numeric(locus)
    # Keep only nonzero numbers of reads and sort remaining reads
    locus <- locus[locus != 0]
    locus <- locus[ordered(-locus)]


    # Return NA if there are no reads for any allele or there are only reads for one allele
    # (iPyrad falsely labeling a homozygote or missing locus as a heterozygote)
    if(length(locus) < 2){return(c(NA, NA))}

    # Return number of reads if there were exactly 2
    else if(length(locus) == 2){return(locus)}


    # Checks for double heterozygotes by looking for more than 2 nonzero reads
    else if((length(locus) > 2) & (remove_double_hets == TRUE))
    {return(c(NA, NA))}

    # Return the 2 highest number of reads otherwise
    else{return(locus[1:2])}
  }
  # Return NA if locus is not heterozygous
  else{return(c(NA, NA))}
}



#' @title Convert Raw VCF Data to HAD Format
#'
#' @description Convert tbl_df output from read_vcf() to a Heterozygous Allele Depth (HAD) format, to be read in by gbs2ploidy.
#'
#' @param filename A character string of the data file path.
#' @param skip_lines A numeric of the number of metadata lines to skip over in the VCF file. If left null, metadata lines are skipped over automatically by the count_metadata_lines function. The count_metadata_lines function requires reading in the entire file, so if you have a large file and know the number of metadata lines in that file, you can save some run time by entering the number of metadata lines in this argument.
#' @param remove_double_hets Logical for determining if double heterozygous loci should be treated as missing information. Should fix issues with gbs2ploidy falsely labeling triploids.
#'
#' @return A data frame in Heterozygous Allele Depth (HAD) format
#'
#' @examples
#' \dontrun{vcf2had("./inst/extdata/example.vcf")}
#' \dontrun{vcf2had("./inst/extdata/example.vcf", remove_dobule_hets=TRUE)}
#'
#' @export
vcf2had <- function(filename, skip_lines=NULL, remove_double_hets=FALSE){
  # Read in VCF file
  vcf_df <- read_vcf(filename, skip_lines)
  # Initialize an empty data frame
  had_df <- data.frame()
  # Get the number of rows and columns
  num_rows    <- dim(vcf_df)[1]
  num_columns <- dim(vcf_df)[2]

  for(i in 1:num_columns){
    # Grab column from VCF data frame
    col <- vcf_df[[i]]
    col_name <- colnames(vcf_df)[i]

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
    if(dim(had_df)[1] == 0 & dim(had_df)[2] == 0){had_df = new_col}
    else{had_df <- cbind(had_df, new_col)}
  }

  return(had_df)
}
