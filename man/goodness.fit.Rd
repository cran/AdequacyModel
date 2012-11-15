\name{goodness.fit}
\alias{goodness.fit}
\encoding{latin1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Adequacy of models
}
\description{
This package provides some useful functions for calculating measures of    goodness of fit, including statistics Cramer-von Mises and Anderson-Darling. These statistics are often used to compare models not fitted. You can also calculate other medias goodness of fit as the AIC, CAIC and BIC.
}
\usage{
goodness.fit(fdp, fda, starts, datas, method = "SANN")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fdp}{
        Specification of the probability density function. 
}
  \item{fda}{
        Specification of cumulative distribution function.
}
  \item{starts}{
        Initial parameters to maximize the likelihood function.
}
  \item{datas}{
        Data vector considered.
}
  \item{method}{
        This option refers to the choice of method to be used to maximize the log-likelihood function. The default method is the \code{SANN} (Simulated Annealing). Outas options accepted are: \code{BFGS} and \code{CG}. It is also accepted as the option argument nlminb (Unconstrained and constrained optimization using PORT routines). More details on these methods, see: \code{\link{optim}}, \code{\link{nlminb}}.
}
}
\details{
The function returns the maximum likelihood estimates as well as some statistical goodness of fit. The goodness of fit statistics are returned: Akaike Information Criterion (AIC), Bayesian Information Criterion (BIC) and Consistent Akaikes Information Criterion (CAIC). Such measures are commonly used to compare two models fitted. For models not fitted comparing statistics is returned Cramer-von Mises (W) and Anderson-Darling (A). These statistics are described in detail in Chen and Balakrishnan (1995).
}

\seealso{\code{\link{nlminb}}, \code{\link{optim}}}

\references{
Chen, G., Balakrishnan, N. (1995). A general purpose approximate goodness-of-fit test. Journal of Quality Technology 27, 154-161.
}

\note{
It is not necessary to define the likelihood function or log-likelihood. You only need to define the probability density function and distribution function. The \code{goodness.fit} needs argument \code{starts} for maximizing the log-likelihood function. Depending on the choice of the initial kick algorithm may not converge. It is advisable to inspect the graphical esolha kicks early.
}

\examples{

data(carbone)
carbone = c(3.70, 2.74, 2.73, 2.50, 3.60, 3.11, 3.27, 2.87, 1.47,
            3.11, 4.42, 2.4, 3.19, 3.22, 1.69, 3.28, 3.09, 1.87,
            3.15, 4.90, 3.75,2.43, 2.95, 2.97, 3.39, 2.96, 2.53,
           2.67, 2.93, 3.22, 3.39, 2.81, 4.20, 3.33, 2.55, 3.31,
           3.31, 2.85, 2.56, 3.56, 3.15, 2.35, 2.55, 2.59, 2.38,
           2.81, 2.77, 2.17, 2.83, 1.92, 1.41, 3.68, 2.97, 1.36,
           0.98, 2.76, 4.91, 3.68, 1.84, 1.59, 3.19, 1.57, 0.81,
           5.56, 1.73, 1.59, 2.00, 1.22, 1.12, 1.71, 2.17, 1.17, 
           5.08, 2.48, 1.18, 3.51, 2.17, 1.69, 1.25, 4.38, 1.84,
           0.39, 3.68, 2.48, 0.85, 1.61, 2.79, 4.70, 2.03, 1.80,
           1.57, 1.08, 2.03, 1.61, 2.12, 1.89, 2.88, 2.82, 2.05,
           3.65)
           
we <- function(par,x){
  c = par[1]
  theta = par[2]
  alpha = par[3]
  ((c*theta)/alpha)*((x/alpha)^(c-1))*exp(-(x/alpha)^c)*((1 - exp(-(x/alpha)^c))^(theta-1))
}

fda_we <- function(par,x){
  c = par[1]
  theta = par[2]
  alpha = par[3]
  (1-exp(-(x/alpha)^c))^theta
}

goodness.fit(fdp=we,fda = fda_we, start=c(1,2,1,1), method="Nelder-Mead",
             datas=carbone)
}

\author{
Pedro Rafael Diniz Marinho \email{pedro.rafael.marinho@gmail.com}

Cicero Rafael Barros Dias \email{cicerorafael@gmail.com}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Adequacy model}

