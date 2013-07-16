\name{goodness.fit}
\alias{goodness.fit}
\encoding{latin1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Adequacy of models
}
\description{
This function provides some useful statistics to assess the quality of fit of probabilistic models, including the statistics Cramer-von Mises and Anderson-Darling. These statistics are often used to compare models not fitted. You can also calculate other media goodness of fit such as AIC, CAIC, BIC and Kolmogorov-Smirnov test.
}
\usage{
goodness.fit(fdp, fda, starts, datas, method="L-BFGS-B", domain=c(0,Inf),
             emv=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fdp}{
        Probability density function. 
}
  \item{fda}{
        Cumulative distribution function.
}
  \item{starts}{
        Initial parameters to maximize the likelihood function.
}
  \item{datas}{
        Data vector.
}
  \item{method}{
        Method used for minimization of the function \code{-log(likelihood)}. The methods supported are: \code{L-BFGS-B} (default), \code{BFGS}, \code{Nelder-Mead}, \code{SANN}, \code{CG}. Can also be transmitted only the first letter of the methodology, ie, \code{L}, \code{B}, \code{N}, \code{S} or \code{C} respectively.
}
  \item{domain}{
        Domain of probability density function. By default the domain of probability density function is the open interval 0 to infinity.This option must be an vector with two values.  
}
  \item{emv}{
        Vector with the estimation maximum likelihood. This option should be used if you already have knowledge of the maximum likelihood estimates. The default is \code{NULL}, ie, the function will try to obtain the estimates of maximum likelihoods.   
}

}
\details{
The function \code{goodness.fit} returns statistics (Anderson-Darling), W (Cramer-von Misses). Are also calculated other measures of goodness of fit. These functions are: AIC (Akaike Information Criterion), CAIC (Consistent Akaikes Information Criterion), BIC (Bayesian Information Criterion) and KS (Kolmogorov-Smirnov).

The Kolmogorov-Smirnov test may return \code{NA} with a certain frequency. The return \code{NA} informs that the statistical \code{KS} is not reliable for the data set used. More details about this issue can be obtained from \code{\link{ks.test}}.

By default, the function calculates the maximum likelihood estimates. The errors of the estimates are also calculated. In cases that the function can not obtain the maximum likelihood estimates, the change of the values initial, in some cases, resolve the problem. You can also enter with the maximum likelihood estimation if there is already prior knowledge. 
}

\value{
\item{W}{Statistic Cramer-von Misses.}
\item{A}{Statistic Anderson Darling.}
\item{KS}{Kolmogorov Smirnov test.}
\item{EMV}{Maximum likelihood estimates.}
\item{AIC}{Akaike Information Criterion}
\item{CAIC}{Consistent Akaikes Information Criterion}
\item{BIC}{Bayesian Information Criterion}
\item{Erro}{Error in the maximum likelihood estimates.}
\item{Value}{Minimum value of the function -log(likelihood).}
\item{Convergence}{0 indicates successful completion and 1 indicates that the iteration limit maxit had been reached. More details at \code{\link{optim}}.}
}

\seealso{For details about the optimization methodologies may view the functions \code{\link{optim}}, \code{\link{ks.test}}, \code{\link{nlminb}}.}

\references{
Chen, G., Balakrishnan, N. (1995). A general purpose approximate goodness-of-fit test. Journal of Quality Technology, 27, 154-161.

Nocedal, J. and Wright, S. J. (1999) Numerical Optimization. Springer.
}

\note{
It is not necessary to define the likelihood function or log-likelihood. You only need to define the probability density function and distribution function. 
}

\examples{

# Example 1:

data(carbone)

# Exponentiated Weibull - Probability density function.
fdp_expweibull <- function(par,x){
  beta = par[1]
  c = par[2]
  a = par[3]
  a * beta * c * exp(-(beta*x)^c) * (beta*x)^(c-1) * (1 - exp(-(beta*x)^c))^(a-1)
}

# Exponentiated Weibull - Cumulative distribution function.
fda_expweibull <- function(par,x){
  beta = par[1]
  c = par[2]
  a = par[3]
  (1 - exp(-(beta*x)^c))^a
}

goodness.fit(fdp=fdp_expweibull, fda=fda_expweibull, 
             starts = c(1,1,1), datas = carbone,
             method="L-BFGS-B", domain=c(0,Inf),emv=NULL)
             
# Example 2:

data = c(0.08, 2.09, 3.48, 4.87, 6.94, 8.66, 13.11, 23.63,
              0.20, 2.23, 3.52, 4.98, 6.97, 9.02, 13.29, 0.40,
              2.26, 3.57, 5.06, 7.09, 9.22, 13.80, 25.74, 0.50,
              2.46, 3.64, 5.09, 7.26, 9.47, 14.24, 25.82, 0.51,
              2.54, 3.70, 5.17, 7.28, 9.74, 14.76, 26.31, 0.81,
              2.62, 3.82, 5.32, 7.32, 10.06, 14.77, 32.15, 2.64,
              3.88, 5.32, 7.39, 10.34, 14.83, 34.26, 0.90, 2.69,
              4.18, 5.34, 7.59, 10.66, 15.96, 36.66, 1.05, 2.69,
              4.23, 5.41, 7.62, 10.75, 16.62, 43.01, 1.19, 2.75,
              4.26, 5.41, 7.63, 17.12, 46.12, 1.26, 2.83, 4.33, 
              5.49, 7.66, 11.25, 17.14, 79.05, 1.35, 2.87, 5.62,
              7.87, 11.64, 17.36, 1.40, 3.02, 4.34, 5.71, 7.93,
              11.79, 18.10, 1.46, 4.40, 5.85, 8.26, 11.98, 19.13,
              1.76, 3.25, 4.50, 6.25, 8.37, 12.02, 2.02, 3.31, 4.51,
              6.54, 8.53, 12.03, 20.28, 2.02, 3.36, 6.76, 12.07,
              21.73, 2.07, 3.36, 6.93, 8.65, 12.63, 22.69)

# Kumaraswamy Weibull Poisson - Probability density function.
fdp_kwweibullpoisson <- function(par,x){
  a = par[1]
  b = par[2]
  c = par[3]
  lambda = par[4]
  beta = par[5]
  (a*b*c*lambda*(beta^c)*(x^(c-1))*((1-exp(-(x*beta)^c))^(a-1)) *
   ((1-(1-exp(-(beta*x)^c))^a)^(b-1)) * 
   exp(-lambda*(1-(1-(1-exp(-(beta*x)^c))^a)^b)
   -  (beta*x)^c))/(1-exp(-lambda))
}

# Kumaraswamy Weibull Poisson - Cumulative distribution function.
fda_kwweibullpoisson <- function(par,x){
  a = par[1]
  b = par[2]
  c = par[3]
  lambda = par[4]
  beta = par[5]
  (1 - exp(lambda*(-(1-(1-(1-exp(-(x*beta)^c))^a)^b))))/(1-exp(-lambda)) 
}

goodness.fit(fdp=fdp_kwweibullpoisson, fda=fda_kwweibullpoisson, 
             starts = c(0.120,1.010,1.000,1.000,0.100), datas = data,
             method="L-BFGS-B", domain=c(0,Inf),emv=NULL)

data = rweibull(250,2,3)
goodness.fit(fdp=fdp_kwweibullpoisson, fda=fda_kwweibullpoisson, 
             starts = c(0.120,1.010,1.000,1.000,0.100), datas = data,
             method="L-BFGS-B", domain=c(0,Inf),emv=NULL)

# Example 3:

# Kumaraswamy Beta - Probability density function.
fdp_kwbeta <- function(par,x){
  beta = par[1]
  a = par[2]
  alpha = par[3]
  b = par[4]
  (a*b*x^(alpha-1)*(1-x)^(beta-1)*(pbeta(x,alpha,beta))^(a-1)*
  (1-pbeta(x,alpha,beta)^a)^(b-1))/beta(alpha,beta) 
}
# Kumaraswamy Beta - Cumulative distribution function.
fda_kwbeta <- function(par,x){
  beta = par[1]
  a = par[2]
  alpha = par[3]
  b = par[4]
  1 - (1 - pbeta(x,alpha,beta)^a)^b
}

data = rbeta(1000,2,2.2)

goodness.fit(fdp = fdp_kwbeta, fda=fda_kwbeta, starts=c(0.5,1.02,2,1),
datas=data, method="L-BFGS-B", domain=c(0,1), emv=NULL)

}

             
\author{
Pedro Rafael Diniz Marinho \email{pedro.rafael.marinho@gmail.com}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Adequacy model}

