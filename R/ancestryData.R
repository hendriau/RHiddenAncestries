#' Example allele frequency data
#'
#' reference data is 1000 Genomes and NAM. 
#' Observed data is from gnomAD:
#' african/african american, 
#' american/admixed latinx, 
#' and other  
#'
#' @docType data
#'
#' @usage data(ancestryData)
#' @format SNP, chromosome, refernce allele frequencies, observed allele frequencies
#'
#' @examples
#' data("ancestryData")
#' ref <-ancestryData[,c("ref_afr", "ref_eas", "ref_eur", "ref_nam", "ref_sas")]
#' obs <- cbind(ancestryData$gnomad_afr)
#' ancestr(ref,obs)
#'
"ancestryData"