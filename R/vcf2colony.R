#' @title Extract Locus Information from VCF Data for COLONY Format
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

  # Get heterozygosity numbers, add 1 to each
  else{
    locus <- strsplit(locus, ":")[[1]][1]
    locus <- strsplit(locus, "/")[[1]]
    locus <- as.numeric(locus)
    return(locus+1)
  }
}



#' @title Convert Raw VCF Data to COLONY Format
#'
#' @description Convert tbl_df output from read_vcf() to a format that can be read in by COLONY software.
#'
#' @param filename A character string of the data file path.
#' @param skip_lines A numeric of the number of metadata lines to skip over in the VCF file. If left null, metadata lines are skipped over automatically by the count_metadata_lines function. The count_metadata_lines function requires reading in the entire file, so if you have a large file and know the number of metadata lines in that file, you can save some run time by entering the number of metadata lines in this argument.
#' @param out_filename A character string of the resulting converted file path.
#'
#' @return NULL
#'
#' @examples
#' \dontrun{vcf2colony(filename = "./inst/extdata/example.vcf", out_filename = "./example.txt")}
#'
#' @export
vcf2colony <- function(filename, skip_lines=NULL, out_filename){
  # Read in VCF file
  vcf_df <- read_vcf(filename, skip_lines)
  # Initialize an empty data frame
  colony_df <- data.frame()

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
      locus <- analyze_locus_colony(locus)

      # Add locus data to new column
      new_col <- rbind(new_col, locus)
    }

    # Get correct column names for new columns
    colnames(new_col) <- c(col_name, col_name)


    # Add new columns if all loci aren't missing
    if(!all(new_col[1] == -9) & !all(new_col[2] == -9)){
      # Add new columns to final data frame if data frame isn't empty
      if(dim(colony_df)[1] == 0 & dim(colony_df)[2] == 0){colony_df = new_col}
      else{colony_df <- cbind(colony_df, new_col)}
    }
  }

  # Write the converted data to a text file
  write.table(colony_df, file=out_filename,
              quote=FALSE, sep="\t",
              row.names=FALSE, col.names=TRUE)
}
