library(readr)
library(dplyr)

#' @title
#' @description
#' @details
#' @param
#' @return
#' @examples
#' \dontrun{}
#' @seealso \url{}
#' @importFrom
#' @export



#' @title Read in VCF files
#'
#' @description Read in a Variant Calling Format file.
#'
#' @details Throws an error if file cannot be found.
#'
#' @param filename A character string of the data file path.
#' @param skip_lines Number of lines to skip before reading in the VCF file. Default is 10.
#'
#' @return This function returns a data frame table of the VCF data with the metadata and the first 9 columns removed.
#'
#' @examples
#' \dontrun{read_VCF("./inst/extdata/example.vcf")}
#'
#' @seealso \url{}
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
    select(-CHROM, -POS, -ID, -REF, -ALT, -QUAL, -FILTER, -INFO, -FORMAT)

  VCF_file
}



#' @title Convert Raw VCF Data to HAD Format
#'
#' @description Convert tbl_df output from read_VCF() to a Heterozygous Allele Depth (HAD) format, to be read in by GBS2Ploidy.
#'
#' @details
#'
#' @param filename A character string of the data file path.
#' @param skip_lines Number of lines to skip before reading in the VCF. Default is 10.
#' @param remove_double_hets Logical
#'
#' @return a data frame in Heterozygous Allele Depth (HAD) format
#'
#' @examples
#' \dontrun{}
#'
#' @seealso \url{}
#'
#' @importFrom
#'
#' @export
VCF2HAD <- function(filename, skip_lines, remove_double_hets){

  # Read in VCF file
  VCF_df <- read_VCF(filename, skip_lines)
  # Initialize an empty data frame
  HAD_df <- data.frame()
  # Get the number of rows and columns
  num_rows    <- dim(VCF_df)[1]
  num_columns <- dim(VCF_df)[2]

  for(i in 1:num_columns){

    # Grab column from VCF data frame
    col = VCF_df[i]

    # Initialize an empty data frame
    new_col <- data.frame()

    # Loop through column and extract heterozygote information
    for(j in 1:num_rows){
      # Grab locus information from column
      locus <- col[j]

      # If locus is heterozygous, extract the number of reads from
      if(startsWith(locus, "0/1") | startsWith(locus, "1/0")){
        locus <- strsplit(locus, ":")[[1]][2]
        locus <- strsplit(locus, ",")
      }
      else{
        locus <- c(NA, NA)
      }
      # Add locus data to new column
      new_col <- rbind(new_col, locus)
    }

    # Add new columns to final data frame
    HAD_df <- cbind(HAD_df, new_col)
  }
}

