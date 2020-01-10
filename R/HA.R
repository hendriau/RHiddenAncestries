## SLSQP Solving Function


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