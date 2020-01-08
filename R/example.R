data <- read.csv("/nfs/storage/math/gross-s2/projects/mixtures/RHiddenAncestries/data/packagedata.csv")
ref <-data[,c("ref_afr", "ref_eas", "ref_eur", "ref_nam", "ref_sas")]
obs <- cbind(data$gnomad_afr)
ancestr(ref,obs)