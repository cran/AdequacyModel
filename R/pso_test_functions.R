#' Test function for the naive_psoptim function.

booth_function <- function(par,x){
  x1 = par[1]
  x2 = par[2]
  (x1+2*x2-7)^2 + (2*x1+x2-5)^2
}

beale_function <- function(par,x){
  x1 = par[1]
  x2 = par[2]
  (1.5-x1+x1*x2)^2 + (2.25-x1+x1*x2^2)^2 + (2.625-x1+x1*x2^3)^2
}

bukin_function <- function(par,x){
  x1 = par[1]
  x2 = par[2]
  100*sqrt(abs(x2-0.01*x1^2))+0.01*abs(x1+10)
}

easom_function <- function(par,x){
  x1 = par[1]
  x2 = par[2]
  -cos(x1)*cos(x2)*exp(-((x1-pi)^2 + (x2-pi)^2))
}

eggholder_function <- function(par,x){
  x1 = par[1]
  x2 = par[2]
  -(x2+47)*sin(sqrt(abs(x2+x1/2+47)))-x1*sin(sqrt(abs(x1-(x2+47))))
} 

f1 <- function(par,x){
  theta=par[1]
  theta*cos(theta)
}

f2 <- function(par,x){
  theta=par[1]
  theta*cos(theta)
  6+(theta^2)*sin(14*theta)
}

schaffer_function <- function(par,x){
  x1 = par[1]
  x2 = par[2]
  0.5+(sin(x1^2-x2^2)^2-0.5)/(1+0.001*(x1^2+x2^2))^2
}

mccormick_function<-function(par,x){
  x1 = par[1]
  x2 = par[2]
  sin(x1+x2)+(x1-x2)^2 -1.5*x1+2.5*x2+1
}

holder_function <- function(par,x){
  x1 = par[1]
  x2 = par[2]
  -abs(sin(x1)*cos(x2)*exp(abs(1-sqrt(x1^2+x2^2)/pi)))
}

goldstein_function <- function(par,x){
  x1 = par[1]
  x2 = par[2]
  (1+((x1+x2+1)^2)*(19-14*x1+3*x1^2-14*x2+6*x1*x2+3*x2^2))*(30+((2*x1-3*x2)^2)*(18-32*x1+12*x1^2+48*x2-36*x1*x2+27*x2^2))
}

cross_function <- function(par,x){
  x1 = par[1]
  x2 = par[2]
  -0.0001*(abs(sin(x1)*sin(x2)*exp(abs(100-sqrt(x1^2+x2^2)/pi)+1)))^0.1
}