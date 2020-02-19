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
#' data("ancestryData")
#' ancestr(ancestryData, 5, 1)
#'
#' @export
#' @importFrom nloptr slsqp
#'
#ancestr = function(refmatrix, obsvector){
ancestr = function(D=NULL, k=NULL, t=1){  
  
    

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
    reference_observed  = cbind(data_processor(ancestryData,k,obs)[[1]],
                                data_processor(ancestryData,k,obs)[[2]])
    refmatrix <- data_processor(ancestryData,k,obs)[[1]]
    
    # the Sequential Quadratic Programming algorithm our method uses (slsqp in the nloptr package) 
    #requires a starting guess at which to evaluate the objective function.
    # We need K starting guesses, one for each of the ancestry proportion values being estimated.
    # Here we set each starting guess to be 1/K.
    # It is necessary that these starting guesses non negative and sum to one.
    starting = numeric(ncol(refmatrix))
    for (i in 1:(ncol(refmatrix))){
    starting[i] = 1/ncol(refmatrix)
    }

    # Here we are defining the objective function. This function is evaluated at a point
    # x in R^K. Each of our K reference allele frequencies are multiplied by 
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
    val = c( S$par,
           S$value,
           S$iter,
           ttime
           )


    return(val)
}