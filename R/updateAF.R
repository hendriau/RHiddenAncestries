#' update allele frequencies function
#' 
#' updates allele frequency data in heterogeneous genetic variant data
#'
#' @param D  is a dataframe: This dataframe must be in the same format as described
#' in the data format document.
#' @param ancestry dataframe: This dataframe is described below in the code.
#' It is also described in the README.
#' @param t is the string name of the observed ancestry
#' @return data.frame D with an additional column at the end containing updated allele frequencies.
#'
#' @author Gregory Matesi, \email{gregory.matesi@ucdenver.edu}
#'
#' @reference \url{http://www.ucdenver.edu/academics/colleges/PublicHealth/Academics/departments/Biostatistics/About/Faculty/Pages/HendricksAudrey.aspx}
#' @keywords genomics
#' 
#' @examples
#' When you make the second input dataframe, the name of the second ancestry 
#' (gnomad_afr in the example) must match the observed column in D that you want
#' to generate updated allele frequencies for.
#' 
#' # Load in the dataframe for the first argument in the function, D, as described below.
#' data(ancestryData)
#' 
#' # User must create a second input datafram containing the pi_hat and pi_star values
#' # for each ancestry. A is also described below.
#' A <-
#'    data.frame(eur_1000G  = c(pi_star = 0, pi_hat = .15),
#'       gnomad_afr = c(pi_star = 1, pi_hat = .85))
#'  
#' # Call the funtion using two inputs 
#' #   (1) ancestryData
#' #   (2) A
#' # and store the results in a new dataframe called E.
#' # The final dolumn of E will contain the updated allele frequencies.
#' E <- updateAF(ancestryData,A)
#' 
#' @export
######################
######################
# INPUT: 2 data.frames
#   D:
#       
#
#   CHR  RSID       bP       A1 A2  ref_eur     ...  ref_iam   obs_afr    obs_amr    obs_oth
#   1    rs2887286  1156131  C  T   0.173275495 ...  0.7093    0.4886100  0.52594300 0.22970500
#   1    rs41477744 2329564  A  G   0.001237745 ...  0.0000    0.0459137  0.00117925 0.00827206
#   1    rs9661525  2952840  G  T   0.168316089 ...  0.2442    0.1359770  0.28605200 0.15561700
#   1    rs2817174  3044181  C  T   0.428212624 ...  0.5000    0.8548790  0.48818000 0.47042500
#   1    rs12139206 3504073  T  C   0.204214851 ...  0.3372    0.7241780  0.29550800 0.25874800
#   1    rs7514979  3654595  T  C   0.004950604 ...  0.0000    0.3362490  0.01650940 0.02481620
#
#   ancestry:
#      data.frame(ref_eur  = c(pi_star = 0, pi_hat = .15),
#                 obs_afr  = c(pi_star = 1, pi_hat = .85))
#
#      data.frame(ref_eur  = c(pi_star = 0, pi_hat = NA),
#                 obs_afr  = c(pi_star = 1, pi_hat = NA))
#
#
#            ref_eur    ref_afr
#   pi_star  0.00       1.00
#   pi_hat   0.15       0.85
#
#            ref_eur    ref_afr
#   pi_star  0.00       1.00
#   pi_hat   NA         NA
######################
######################

updateAF <- function(D=NULL, ancestry=NULL, k=NULL, t="obs_afr"){

#   Check to see if pi_hats are NA. If the user entered NA for the pi_hats, use ancestr to obtain 
#      ancestry proportion estimates.
  if(all(is.na(ancestry[2,])) == TRUE){
     t1 <- which(colnames(D) == t) - (5+k)
  ancestr(D,k,t1)
  }
    
#   Normalize pi. We need the pi_hats and pi_stars to sum to 1
  ancestry[2,] <- ancestry[2,] / sum(ancestry[2,])
  ancestry[1,] <- ancestry[1,] / sum(ancestry[1,])
  
#   Grab the names of the reference ancestries and the single
#     observed ancestry.
#   num_ref_ancestries is the number of reference ancestries
#     being used for the method.
  ref_ancestries <- colnames(ancestry[1:dim(ancestry)[2]-1])
  obs_ancestry <- colnames(ancestry[dim(ancestry)[2]])
  num_ref_ancestries <- length(ref_ancestries)

#   We need to sum up the reference ancestries multiplied by pi_star. We will call this sum "starred"
#   We also need to sum up the reference ancestries multiplied by pi_hat. We will call this sum "hatted"
  
#   Initialize "hatted" and "starred"
  hatted  <- vector(mode = "double", length = dim(D)[1])
  starred <- vector(mode = "double", length = dim(D)[1])
  
#   Sum the K-1 reference ancestries multiplied by pi_hat.
#   Also the sum the k-1 reference ancestries multiplied by pi_star.
  for ( k in 1:length(ref_ancestries) ){
    hatted  <- hatted +  ancestry[2,ref_ancestries[k]] * D[ref_ancestries[k]]
    starred <- starred + ancestry[1,ref_ancestries[k]] * D[ref_ancestries[k]] #comment
  }
  
#   Add a new column to D:
#      Ancestry adjusted allele frequency time the ratio or pi_star/pi_hat
#      plus the sum "starred"
  D$updatedAF <- (
    
    ( ancestry[1,obs_ancestry] / ancestry[2,obs_ancestry] ) *
      
    ( D[obs_ancestry] - hatted )[[1]]
    
    + starred[[1]] )
    return(D)
}

