rkwbs <- function(n, alpha, beta, a, b){
  V = rkumar(n, a, b)
  InversePhi = qnorm(V,0,1)
  return(beta*((alpha*InversePhi)/2 + (1+(alpha^(2))*(InversePhi^2))^(0.5))^(2) )
}