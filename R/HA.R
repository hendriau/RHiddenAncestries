#' Hidden Ancestries function
#'
#' Estimates ancestry proportions in heterogeneous allele frequency data
#'
#' @param refmatrix matrix of reference allele frequency data
#' @param obsvector observed heterogeneous allele frequency vector
#' @return Estimated ancestry proportions
#'
#' @author Gregory Matesi, \email{gregory.matesi@ucdenver.edu}
#' @reference \url{http://www.ucdenver.edu/academics/colleges/PublicHealth/Academics/departments/Biostatistics/About/Faculty/Pages/HendricksAudrey.aspx}
#' @keywords genomics
#' 
#' @examples
#' data <- data("packagedata")
#' ref <-data[,c("ref_afr", "ref_eas", "ref_eur", "ref_nam", "ref_sas")]
#' obs <- cbind(data$gnomad_afr)
#' ancestr(ref,obs)
#'
#' @export
#' @importFrom nloptr slsqp
#'
ancestr = function(refmatrix, obsvector){
  
  
  testmatrix = cbind(refmatrix, obsvector)
  
  
  starting = numeric(ncol(refmatrix))
  for (i in 1:(ncol(refmatrix))){
    starting[i] = 1/ncol(refmatrix)
  }
  
  
  fn.ancmix = function(x){
    minfunc = 0
    for (i in 1:ncol(refmatrix)){
      minfunc = minfunc + x[i]*testmatrix[i]
    }
    minfunc = minfunc - testmatrix[ncol(refmatrix) + 1]
    minfunc = sum((minfunc)**2)
    return(minfunc)
  }
  
  
  gr.ancmix <- function(x){
    
    gradvec = matrix(0,ncol(refmatrix),1)
    
    gradfunc = 0
    for (i in 1:ncol(refmatrix)){
      gradfunc = gradfunc + x[i]*testmatrix[i]
    }
    gradfunc = gradfunc - testmatrix[ncol(refmatrix) + 1]
    
    for (i in 1:ncol(refmatrix)){
      gradvec[i] = sum(2 * testmatrix[i] * gradfunc)
    }
    
    return(gradvec)
  }
  
  
  heq.ancmix = function(x){
    equality = 0
    for (i in 1:ncol(refmatrix)){
      equality = equality + x[i]
    }
    return(equality - 1)
  }
  
  
  hin.ancmix <- function(x){
    h = numeric(ncol(refmatrix))
    for (i in 1:ncol(refmatrix)){
      h[i] = x[i]
    }
    return(h)
  }
  
  
  start_time = Sys.time()
  S = slsqp(starting,
             fn = fn.ancmix,
             gr = gr.ancmix,
             hin = hin.ancmix,
             heq = heq.ancmix)
  end_time = Sys.time()
  ttime = end_time - start_time
  
  
  val = c( S$par,
           S$value,
           S$iter,
           ttime
           )
  
  
  return(val)
}