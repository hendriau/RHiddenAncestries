#' update allele frequencies function
#' 
#' updates allele frequency data in heterogeneous genetic variant data
#'
#' @param target
#' @param reference
#' @param pi_target
#' @param pi_reference
#' @param D  is a dataframe: This dataframe must be in the same format as described
#' in the data format document.
#' @param k number of reference ancestries in the model
#' Only to be used when pi_hats are unknown/NA.
#' Only used in the ancestr function.
#' @return data.frame D with an additional column at the end containing updated allele frequencies.
#'
#' @author Gregory Matesi, \email{gregory.matesi@ucdenver.edu}
#'
#' @reference \url{http://www.ucdenver.edu/academics/colleges/PublicHealth/Academics/departments/Biostatistics/About/Faculty/Pages/HendricksAudrey.aspx}
#' @keywords genomics
#' 
#' @examples
#' 
#' #########
#' Example 1
#' #########
#' 
#' The vignette is not updated yet
#' 
#' # Load in the dataframe for the first argument in the function, D, as described in the vignette.
#' data(ancestryData)
#'  
#' # Call the funtion using 5 inputs 
#'
#' #   (1) target       = "obs_afr
#' #   (2) reference    = c("ref_eur")
#' #   (3) pi_target    = c(1,0)
#' #   (4) pi_reference = c(.85, .15)
#' #   (5) data         = ancestryData 
#' #
#' # and store the results in a new dataframe called E.
#' # The final column of E will contain the updated allele frequencies.
#' E <- updateAF(ancestryData,A)
#' 
#' 
#' #########
#' Example 2
#' #########
#' 
#' The vignette is not updated yet
#' 
#' # Load in the dataframe for the first argument in the function, D, as described in the vignette.
#' data(ancestryData)
#'       
#' # Call the funtion using 6 inputs 
#'
#' #   (1) target       = "obs_afr
#' #   (2) reference    = c("ref_eur")
#' #   (3) pi_target    = c(1,0)
#' #   (4) pi_reference = c(NA)
#' #   (5) data         = ancestryData 
#' #   (6) k            = 5
#' 
#' # and store the results in a new dataframe called E.
#' # The final column of E will contain the updated allele frequencies.
#' E <- updateAF(ancestryData,A)
#' 
#' @export
######################
######################
# INPUTS:
#
#   target (string): The column name of the target ancestry
#
#   reference (string valued vector): 
#     column names of the reference
#     ancestries
#
#   pi_target (real valued vector):
#     target ancestry proportions starting with
#     the target ancestry and following in order of
#     reference ancestries
#
#   pi_reference (real valued vector):
#     reference ancestry proportions starting with
#     the target ancestry and following in order of
#     reference ancestries
#
#   data (dataframe):
#     
#
#     CHR  RSID       POS       A1 A2  ref_eur     ...  ref_iam   obs_afr    obs_amr    obs_oth
#     1    rs2887286  1156131  C  T   0.173275495 ...  0.7093    0.4886100  0.52594300 0.22970500
#     1    rs41477744 2329564  A  G   0.001237745 ...  0.0000    0.0459137  0.00117925 0.00827206
#     1    rs9661525  2952840  G  T   0.168316089 ...  0.2442    0.1359770  0.28605200 0.15561700
#     1    rs2817174  3044181  C  T   0.428212624 ...  0.5000    0.8548790  0.48818000 0.47042500
#     1    rs12139206 3504073  T  C   0.204214851 ...  0.3372    0.7241780  0.29550800 0.25874800
#     1    rs7514979  3654595  T  C   0.004950604 ...  0.0000    0.3362490  0.01650940 0.02481620
#

#
#   k (integer):
#     Number of reference ancestries. Only used if the user
#     does not enter a pi_reference vector
######################
######################

updateAF2 <- function(target       = "None", 
                      reference    = c("None"),
                      pi_target    = c(NA), 
                      pi_reference = c(NA), 
                      data=NULL, 
                      k=NULL){
#   Output the users ancestries and pi values in a dataframe as a check.
  input_display <- data.frame(matrix(ncol = length(pi_target), nrow = 2))
  
  # Set row and column names
  colnames(input_display) <- c(target,reference); rownames(input_display) <- c("target", "reference")
  
  # Set target and reference proportion values
  input_display[1,] <- pi_target; input_display[2,] <- pi_reference
  
  # Show the user thier entries
  print(input_display)


#   Check to see if pi_reference are NA. 
#   If yes, use the HA function.
  
    
#   Normalize pi. We need the pi_reference and pi_target to sum to 1
  pi_target <- pi_target / sum(pi_target)
  pi_reference <- pi_reference / sum(pi_reference)

 

#   We need to sum the reference ancestry allele frequencies multiplied by pi_target. 
#   We will call this sum "starred"
#   We also need to sum up the reference ancestry allele frequencies multiplied by pi_reference. 
#   We will call this sum "hatted"
  
#   Initialize "hatted" and "starred"
  hatted  <- vector(mode = "double", length = dim(D)[1])
  starred <- vector(mode = "double", length = dim(D)[1])
  
#   Sum the K-1 reference ancestries multiplied by pi_hat.
#   Also the sum the k-1 reference ancestries multiplied by pi_star.
  for ( k in 1:length(reference) ){
    hatted  <- hatted +  pi_target[k]    * D[reference[k]]
    starred <- starred + pi_reference[k] * D[reference[k]]
  }
  
#   Add a new column to D:
#      Ancestry adjusted allele frequency time the ratio or pi_star/pi_hat
#      plus the sum "starred"
  D$updatedAF <- (
    
    ( pi_target[1] / pi_reference[1] ) *
      
    ( D[target] - hatted )[[1]]
    
    + starred[[1]] )
    return(D)
}

