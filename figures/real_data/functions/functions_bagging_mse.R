enlist <- function (...) 
{
  result <- list(...)
  if ((nargs() == 1) & is.character(n <- result[[1]])) {
    result <- as.list(seq(n))
    names(result) <- n
    for (i in n) result[[i]] <- get(i)
  }
  else {
    n <- sys.call()
    n <- as.character(n)[-1]
    if (!is.null(n2 <- names(result))) {
      which <- n2 != ""
      n[which] <- n2[which]
    }
    names(result) <- n
  }
  result
}


se <- function(x, margin =2){
  
  if (margin ==2) {
    out <- apply(x, margin, sd, na.rm = T)/
      sqrt(colSums(!is.na(x)))
  }else if(margin ==1){
    out <- apply(x, margin, sd, na.rm = T)/
      sqrt(rowSums(!is.na(x)))
  }
  
  
  return(out)
}

#function to generate lambda params
generate_params <- function(n, points = 25) {
  raw_vals <- unique(round(exp(seq(log(1), log(n/2), length.out = points))))
  params <- unique(c(0, raw_vals[raw_vals <= n/2]))
  return(params)
}