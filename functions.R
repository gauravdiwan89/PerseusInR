##these are custom functions that MUST be loaded at the beginning of the perseus pipeline

impute_norm <- function(x, down_shift = 1.8, width = 0.3, n = 1) {
  x <- x[x>0]
  if(length(x) > 1) {
    m <- mean(x, na.rm = T) - (down_shift*sd(x, na.rm = T))
    s <- width*sd(x, na.rm = T)
    sample(rnorm(n = 100, mean = m, sd = s), n)
  # } else if(length(x) == 1) {
  #   m <- mean(x, na.rm = T)
  #   s <- 0
  #   sample(rnorm(n = 100, mean = m, sd = s), n)
  } else 0
}

erfc <- function(x) 2 * pnorm(x * sqrt(2), lower = FALSE)
significanceA <- function(r, log = F, log.base=2){
  
  ###################################################
  #  transform to log scale
  ###################################################
  if(log) r <- log(r, log.base)
  
  quant.idx <- which(!is.na(r))
  ##################################################
  # robust and asymmetrical estimate of the sd
  #
  ##################################################
  estimate.sd <- quantile(r[quant.idx], c(0.1587, 0.5, 0.8413))
  
  # the notation used in the MQ supplement
  r.minus1 <- estimate.sd[1]
  r.zero <- estimate.sd[2]
  r.plus1 <- estimate.sd[3]
  
  ##################################################
  # z transformation
  ##################################################
  r.z <- ifelse(r >= r.zero, (r-r.zero)/(r.plus1-r.zero),  (r.zero-r)/(r.zero-r.minus1) )
  
  #sigA <- 0.5*erfc( r.z/sqrt(2) )
  sigA <- erfc( r.z/sqrt(2) )
  
  return(sigA)
}