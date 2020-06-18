#' @title Generate Random Locus Data
#' @description Generates heterozygous, homozygous, or missing locus data. Not exported.
#' @details Heterozygous or homozgous loci are given a random number of reads for each allele.
#' @param min_reads the minimum number of allele reads at each locus
#' @param max_reads the maximum number of allele reads at each locus
#' @return A character string of locus data
#' @examples
#' \dontrun{make_locus()}
#' \dontrun{make_locus(min_reads=16, max_reads=100)}
make_locus <- function(min_reads=4, max_reads=200){

  # Number of total reads
  total_num_reads <- sample(min_reads:max_reads, 1)

  # Determine if locus will be missing, heterozygous, or homozygous
  type_locus <- sample(1:3, 1)

  if(type_locus == 1){return("./.:0:0,0,0,0")}
  else if(type_locus == 2){
    # Choose beginning homozygote signal
    front <- sample(c("0/0", "1/1"), 1)

    return(paste0(front, ":", as.character(total_num_reads), ":",
                                         as.character(total_num_reads), ",0,0,0"))}
  else{
    # number of reads on each nucleotide
    num_reads <- c(0,0,0,0)

    # Add a number of reads to each nucleotide slot
    for(i in 1:4){
      add_to <- sample(1:4, 1)
      num_reads[add_to] <- num_reads[add_to] + (total_num_reads %/% 4)
    }
    # Fix total number of reads (integer division in previous line rounds down)
    sum_reads <- sum(num_reads)

    # Choose beginning heterozygote signal
    front <- sample(c("1/0", "0/1"), 1)
    # Combine number of reads
    num_reads <- paste(as.character(num_reads), collapse=",")
    return(paste0(front, ":", sum_reads, ":", num_reads))
  }
}



#' @title Generate specimen names
#' @description A very, very crude and confusing way to generate character strings of a number and a letter. Not exported.
#' @details The generated locus information will yield any meaningful information, this function is here only for testing the converter function.
#' @param n the number of specimen names to generate
#' @return A character vector of specimen names
#' @examples
#' \dontrun{make_specimen_names(8)}
make_specimen_names <- function(n){

  # initialize a blank vector of names and a prefix number
  spec_names <- vector()
  num=1

  for(i in 0:(n-1)){
    # Add 1 to the number prefix every time letters loop around
    if(i/26 == 1){num <- num+1}

    # Combine the number prefix to a letter
    new_name <- paste0(as.character(num), LETTERS[(i%%26)+1])

    # Add spec name to vector of names
    spec_names <- c(spec_names, new_name)
  }

  return(spec_names)
}


#' @title Generate a basic example VCF file
#' @details The generated file will hold no meaningful information, this function is here only for testing the converter function.
#' @param save_path the file path and file name to save the generated VCF file
#' @param min_reads the minimum of total number of reads per locus
#' @param max_reads the maximum of total number of reads per locus
#' @param num_loci the number of loci (rows) to generate in the file
#' @param num_specimens the number of specimen (columns) to generate in the file
#' @examples
#' \dontrun{make_example_VCF_file()}
#' \dontrun{make_example_VCF_file(save_path="./my_example.vcf")}
#' \dontrun{make_example_VCF_file(min_reads=24, max_reads=100)}
#' \dontrun{make_example_VCF_file(num_loci=100, num_specimens=4)}
#' @importFrom dplyr rename
#' @importFrom dplyr mutate
#' @importFrom magrittr "%>%"
#' @importFrom utils write.table
#' @export
make_example_VCF_file <- function(save_path="./generated_example.vcf",
                              min_reads=8, max_reads=200, num_loci=10, num_specimens=5){

  # Make initial column to make adding the other columns to data frame easier
  CHROM = paste0("locus_", 1:num_loci)

  df <- data.frame(CHROM)

  df <-
    df %>%
    dplyr::rename("#CHROM"="CHROM") %>%
    dplyr::mutate(POS = 1:num_loci) %>%
    dplyr::mutate(ID = ".") %>%
    dplyr::mutate(REF = sample(c("A","T","C","G"), 1)) %>%
    dplyr::mutate(ALT = sample(c("A","T","C","G"), 1)) %>%
    dplyr::mutate(QUAL = "13") %>%
    dplyr::mutate(FILTER = "PASS") %>%
    dplyr::mutate(INFO = paste0("NS=", sample(1:100, 1), ";DP=", sample(1000:9999, 1))) %>%
    dplyr::mutate(FORMAT = "GT:DP:CATG")

  ## Make a matrix of locus reads
  # Make an empty matrix
  loci <- matrix(NA, nrow=num_loci, ncol=num_specimens)

  # Fill empty matrix with locus reads
  for(row in 1:num_loci){
    for(specimen in 1:num_specimens){
      loci[row, specimen] <- make_locus(min_reads=min_reads, max_reads=max_reads)
    }
  }

  # Add specimen names
  colnames(loci) <- make_specimen_names(n=num_specimens)

  # Cast matrix to a data frame
  loci <- as.data.frame(loci)

  # Combine data frames of loci and loci information
  df <- cbind(df, loci)

  # Collapse rows into a single row of strings
  #paste_args <- c(df, sep="\t")
  #df <- data.frame(V1 = do.call(paste, paste_args))

  # Add metadata
  metadata <- data.frame(V1 = c('##fileformat=VCFv4.0',
                '##fileDate=2020/06/08',
                '##source=ipyrad_v.0.7.13',
                '##reference=pseudo-reference (most common base at site)',
                '##phasing=unphased',
                '##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">',
                '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">',
                '##FORMAT=<ID=CATG,Number=1,Type=String,Description="Base Counts (CATG)">'))

  # Save data frame to a file
  write.table(metadata, file=save_path, append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE)
  suppressWarnings(write.table(df, file=save_path, append=TRUE, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE))
}
