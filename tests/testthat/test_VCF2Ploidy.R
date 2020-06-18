context("Throwing Errors")

test_that("Throws errors", {
  throws_error(read_VCF("ThisFileNameDoesn'tExist.csv"))
})

context("Reading files are correct")

test_that("Reads files correctly", {
  filename <- system.file("extdata", "example1.vcf", package = "VCF2Ploidy", mustWork=TRUE)
  df <- read_VCF(filename)
  expect_that(df, is_a("data.frame"))
  expect_equal(length(df), 5)
})

context("Converting VCF Data is Correct")

test_that("Locus Analysis is correct", {
  # Missing locus
  locus1 <- "./.:0:0,0,0,0"
  # Homozygote
  locus2 <- "1/1:100:100,0,0,0"
  # Heterozygote
  locus3 <- "1/0:100:75,0,25,0"
  # Double heterozygote
  locus4 <- "1/0:150:75,0,25,50"

  locus_out1 <- analyze_locus(locus1)
  locus_out2 <- analyze_locus(locus2)
  locus_out3 <- analyze_locus(locus3)
  locus_out4 <- analyze_locus(locus4, remove_double_hets=TRUE)
  locus_out5 <- analyze_locus(locus4, remove_double_hets=FALSE)

  # Check that analyze_locus() returns a vector of length 2
  expect_that(length(locus_out1), equals(2))

  # Check that Missing loci and homozygotes returns NAs
  expect_that(locus_out1, is_equivalent_to(c(NA, NA)))
  expect_that(locus_out2, is_equivalent_to(c(NA, NA)))

  # Check that heterozygotes return numerics
  expect_that(locus_out3, is_a("numeric"))
  expect_that(locus_out3, is_equivalent_to(c(75,25)))

  # Check that double heterozygotes return NAs when remove_double_hets==TRUE
  expect_that(locus_out4, is_equivalent_to(c(NA, NA)))
  # Check that dobule hets return two largest reads when remove_double_hets==FALSE
  expect_that(locus_out5, is_equivalent_to(c(75, 50)))
})

test_that("VCF conversion is correct", {
  # Run conversion
  filename <- system.file("extdata", "example1.vcf", package = "VCF2Ploidy", mustWork=TRUE)
  df <- VCF2HAD(filename=filename, remove_double_hets=FALSE)

  # Check dimensions converted example data
  expect_that(dim(df)[1], equals(10))
  expect_that(dim(df)[2], equals(10))

  # Check that column names are correct
  correct_names <- c("1A", "1A", "1B", "1B", "1C", "1C", "1D", "1D", "1E", "1E")
  expect_that(colnames(df), is_identical_to(correct_names))

  expect_that(df[1, 1:2], is_equivalent_to(c(30, 10)))
})

context("Ploidy estimation is correct")

test_that("Ploidy estimation is correct", {
  filename <- system.file("extdata", "example1.vcf", package = "VCF2Ploidy", mustWork=TRUE)
  df <- VCF2Ploidy(filename=filename,
                   remove_double_hets=TRUE, props=c(0.25, 0.5, 0.75))

  expect_that(df, is_a("list"))
  expect_that(names(df), is_equivalent_to(c("pp", "pcwghts", "pcscrs")))
})
