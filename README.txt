This is the R package for the Hidden Ancestries project overseen by Dr Audrey Hendricks at CU Denver.
This package contains one function, ancestr, written by Ian. This function takes in a NxK+5 dataframe containing a chromosome column, a RSID column, reference and alternate allele columns, a column for each of the K reference ancestries, and a column for the observed homogeneous ancestry. It outputs a vector of hidden ancestry estimates for each of the K reference ancestries, minimization value, run time, and number of iterations.

#' @examples
#' data("ancestryData")
#' ref <-ancestryData[,c("ref_afr", "ref_eas", "ref_eur", "ref_nam", "ref_sas")]
#' obs <- cbind(ancestryData$gnomad_afr)
#' ancestr(ref,obs)

