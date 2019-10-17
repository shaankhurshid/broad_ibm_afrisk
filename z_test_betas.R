z <- function(beta1,beta2,se1,se2){
  z <- abs(beta1 - beta2)/sqrt(se1^2+se2^2)
  p <- 2*pnorm(-abs(z))
  out <- c(z,p)
  return(out)
}

# Script to assess for differences in betas according to presence/absence of precipitant
hf <- z(4.21,2.75,)
