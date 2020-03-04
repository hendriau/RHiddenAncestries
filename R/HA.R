#' Hidden Ancestries function
#'
#' Estimates ancestry proportions in heterogeneous allele frequency data
#'
#' @param D dataframe: This dataframe should be in the format described in our data formatting document.
#' @param k integer: This parameter is the number of reference ancestries present in your dataframe D.
#' @param t integer: variable t is the number of the observed allele frequency
#' column the user wishes to run through SQP. If not specified, the function will just use the first
#' observed vector.
#' @param x_0 : x_0 is the starting guess for the SLSQP algorithm.
#' @return Estimated ancestry proportions
#'
#' @author Gregory Matesi, \email{gregory.matesi@ucdenver.edu}
#' @reference \url{http://www.ucdenver.edu/academics/colleges/PublicHealth/Academics/departments/Biostatistics/About/Faculty/Pages/HendricksAudrey.aspx}
#' @keywords genomics
#' 
#' @examples
#' data("ancestryData")
#' ancestr(ancestryData, 5, 1)
#'
#' @export
#' @importFrom nloptr slsqp
#'
ancestr = function(D=NULL, k=0, t=0, x_0 =NULL){  
  
  ##########################
  # Initial Checks
  ##########################
  # Check if D was specified
  if (length(D)==0){
    stop("ERROR: D is of length 0. 
         Make sure that your file format fits:
         chromosome, rsid, bP, a1, a2, reference ancesries, observed ancestries")
  }
  
  # Check if D has at least 5+k+t columns
  if(dim(D)[2] < 5+k+t){
    stop("Make sure there are k reference ancestries and at least t observed ancestries in your dataframe.
         Make sure that your file format fits:
         chromosome, rsid, bP, a1, a2, reference ancesries, observed ancestries")
  }
  
  # Check if k=0
  if (k==0){
    stop("ERROR: k=0. Please specify a number of reference ancestries for k")
  }
  
  # Check if user specified an observed ancestry
  if (t==0){
    stop("ERROR: t=0 Please specify an observed ancestry allele frequency column from your data.")
  }
  ################################################################################################
  ################################################################################################
  
  # The user enters a Nx... matrix containging 
  # chromosome, rsid, base pair, A1, A2, k reference ancestries, t observed ancestries.
  # This function processes this matrix and produces two outputs
  # A Nxk reference matrix
  # A N vector of the observed allele frequencis
  # Essentially, we are getting rid of the chromosome, rsid, position, and 
  # referance/alternate allele information as this is unnecessary for 
  # our calculation.
  data_processor <- function( D, k, t ){
    data_processor_list <- list(D[,c(6:(5+k))],
                                D[,c(5+k+t), drop = FALSE])
    return(data_processor_list)
  }
  
  # First we create a new matrix that concatenates the reference and observed columns
  # The first item in the list returned by the data_processor function is our
  # reference columns. The second item is the observed column.
  reference_observed  = cbind(data_processor(D,k,t)[[1]],
                              data_processor(D,k,t)[[2]])
  refmatrix <- data_processor(D,k,t)[[1]]
  
  # the Sequential Quadratic Programming algorithm our method uses (slsqp in the nloptr package) 
  #requires a starting guess at which to evaluate the objective function.
  # We need K starting guesses, one for each of the ancestry proportion values being estimated.
  # Here we set each starting guess to be 1/K.
  # It is necessary that these starting guesses are non negative and sum to one.
  
  # If user specified a starting guess, set to to the variable called "starting"
  if(length(x_0) != 0){
    
    #########################
    # Check if x_0 is numeric
    if (is.numeric(x_0)==FALSE){
      stop("ERROR: Please make sure x_0 is a positive numeric vector of length k that sums to one")
    }
    
    # Check if length(x_0)==k
    if (length(x_0)!= k){
      stop("ERROR: Please make sure x_0 is a positive numeric vector of length k that sums to one")
    }
    
    # Check that x_0 is positive
    if (all(x_0>0) == FALSE){
      stop("ERROR: Please make sure x_0 is a positive numeric vector of length k that sums to one")
    }
    
    # Check if sum(x_0) = 1
    if (sum(x_0)!=1){
      stop("ERROR: Please make sure x_0 is a positive numeric vector of length k that sums to one")
    }
    #########################
    
    ###############################
    # Set the starting guess to x_0
    ###############################
    starting = x_0
  }
  
  # If the user did not specify a starting guess then one is provided.
  # The starting guess will be 1/k for each of the k ancestry proportion values.
  else{
    starting = numeric(ncol(refmatrix))
    for (i in 1:(ncol(refmatrix))){
      starting[i] = 1/ncol(refmatrix)
    }
  }
  
  # Here we are defining the objective function. This function is evaluated at a 
  # k-dimensional point x . Each of our K reference allele frequencies are multiplied by 
  # our current best guess for the ancestry proportion.
  # We then subtract the allele frequency values from the observed homogeneous population.
  # And finally this sum is squared to achieve a least squares form.
  fn.ancmix = function(x){
    minfunc = 0
    for (i in 1:ncol(refmatrix)){
      minfunc = minfunc + x[i]*reference_observed[i]
    }
    minfunc = minfunc - reference_observed[ncol(refmatrix) + 1]
    minfunc = sum((minfunc)**2)
    return(minfunc)
  }
  
  # Here we are defining the gradient of the objective function.
  gr.ancmix <- function(x){
    gradvec = matrix(0,ncol(refmatrix),1)
    gradfunc = 0
    for (i in 1:ncol(refmatrix)){
      gradfunc = gradfunc + x[i]*reference_observed[i]
    }
    gradfunc = gradfunc - reference_observed[ncol(refmatrix) + 1]
    
    for (i in 1:ncol(refmatrix)){
      gradvec[i] = sum(2 * reference_observed[i] * gradfunc)
    }
    return(gradvec)
  }
  
  # H equality
  # This function returns the equality constraints for the nloptr slsqp algorithm
  # We sum up the K current proportion estimate values and subtract 1. If the estimated proportion values sum to 1 than this value should equal zero.
  heq.ancmix = function(x){
    equality = 0
    for (i in 1:ncol(refmatrix)){
      equality = equality + x[i]
    }
    return(equality - 1)
  }
  
  # H inequality
  # This function returns a K vector of 
  hin.ancmix <- function(x){
    h = numeric(ncol(refmatrix))
    for (i in 1:ncol(refmatrix)){
      h[i] = x[i]
    }
    return(h)
  }
  
  # We use the start_time function base function 
  #to record the run time for our convex optimization algorithm.
  # The output for the nloptr slsqp function is stored in the variable S.
  # We are inputing 5 values into this function
  # 1. fn.ancmix is the objective function.
  # 2. gr.ancmix is the gradient of the objective function.
  # 3. hin.ancmix 
  # 4. heq.ancmix defining the equality constraints
  start_time = Sys.time()
  S = slsqp(starting,
            fn = fn.ancmix,
            gr = gr.ancmix,
            hin = hin.ancmix,
            heq = heq.ancmix)
  end_time = Sys.time()
  ttime = end_time - start_time
  
  # par is the K optimal ancestry propotion estimates
  # value is the minimization function evaluated at par
  # iter is the number of iterations that the algorithm took to reach the optimal solution of par
  # finally ttime is the run time for the algorithm
  
  d <- data.frame() 
  for (i in 1:k){
    d[1,i] <- colnames(D)[5+i] 
    d[2,i] <- S$par[i]
  }
  val = c( d,
           S$value,
           S$iter,
           ttime
  )
  
  
  return(print(val))
}