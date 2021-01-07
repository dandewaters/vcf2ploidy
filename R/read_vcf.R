#' @title Count number of metadata lines in VCF files
#'
#' @description Read in a Variant Calling Format file and count the number of metadata lines to skip over for read_vcf function. Not exported.
#'
#' @details Throws an error if file cannot be found. Metadata lines are identified by checking if the lines starts with "##". This function requires reading in the entire file, so if you have a large file and know the number of metadata lines in that file, you can save some run time in the read_vcf, vcf2had, vcf2ploidy, and vcf2colony functions by entering the number of metadata lines in their "skip_lines" argument.
#'
#' @param filename A character string of the data file path.
#'
#' @return Returns a numeric of the number of metadata lines in VCF file.
#'
#' @examples
#' \dontrun{count_metadata_lines("./inst/extdata/example.vcf")}
count_metadata_lines <- function(filename){
  if(!file.exists(filename))
    stop("file '", filename, "' does not exist")

  # Read in VCF File line by line
  vcf_file <- readLines(filename)

  # Initialize parameters to loop through lines of VCF file, keep counting until there are no more metadata lines
  line <- 1
  num_metadata_lines <- 0
  is_metadata_line <- TRUE

  # Loop through lines of VCF file,
  while(is_metadata_line == TRUE){
    # Count metadata lines
    if(startsWith(vcf_file[line], "##")){
      line <- line + 1
      num_metadata_lines <- num_metadata_lines + 1
    }
    # Break loop if line isn't a metadata line
    else{is_metadata_line <- FALSE}
  }

  return(num_metadata_lines)
}

#' @title Read in VCF files
#'
#' @description Read in a Variant Calling Format file.
#'
#' @details Throws an error if file cannot be found.
#'
#' @param filename A character string of the data file path.
#' @param skip_lines A numeric of the number of metadata lines to skip over in the VCF file. If left null, metadata lines are skipped over automatically by the count_metadata_lines function. The count_metadata_lines function requires reading in the entire file, so if you have a large file and know the number of metadata lines in that file, you can save some run time by entering the number of metadata lines in this argument.
#'
#' @return This function returns a data frame table of the VCF data with the metadata and the first 9 columns removed.
#'
#' @examples
#' \dontrun{read_vcf("./inst/extdata/example.vcf")}
#'
#' @importFrom readr cols
#' @importFrom readr read_delim
#' @importFrom dplyr tbl_df
#' @importFrom dplyr rename
#' @importFrom dplyr select
#' @importFrom magrittr "%>%"
#'
#' @export
read_vcf <- function(filename, skip_lines=NULL){
  if(!file.exists(filename))
    stop("file '", filename, "' does not exist")

  if(is.null(skip_lines) || is.na(skip_lines)){
    # Count the number of metadata lines to skip over in VCF file
    skip_lines <- count_metadata_lines(filename)
  }

  # Read in VCF file
  vcf_file <- read_delim(file=filename, delim="\t", col_types = cols(.default = "c"), skip=skip_lines)

  # Remove first 9 columns
  vcf_file <-
    vcf_file %>%
    # Removing the pound symbol on the "#CHROM" column so the select() arguments don't get commented out
    rename("CHROM" = "#CHROM") %>%
    dplyr::select(-CHROM, -POS, -ID, -REF, -ALT, -QUAL, -FILTER, -INFO, -FORMAT)

  return(vcf_file)
}
