\name{rew}
\alias{rew}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pseudo-Random Numbers - Exponentiated Weibull
}
\description{
Generates pseudorandom numbers from Exponentiated Weibull Distribution.
}
\usage{
rew(n, beta, c, a)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
  Amount of generated numbers.
}
  \item{beta}{
  Scale parameter of the Weibull distribution exponentiated
}
  \item{c}{
  Shape parameter Exponentiated Weibull distribution.
}
  \item{a}{
  Shape parameter exponentiated Weibull distribution.
}
}
\details{
The exponentiated Weibull family of probability distributions was introduced by Mudholkar and Srivastava (1993) as an extension of the Weibull family obtained by adding a second shape parameter. The function is defined in positive real, ie \eqn{x>0}. 
    \deqn{G(x, \beta, c, a) = (1 - e^{-(\beta*x)^c})^{a},}
wherein \eqn{\beta>0}, \eqn{c>0} and \eqn{a>0}.
}

%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
\references{
Mudholkar, G.S.; Hutson, A.D. (1996). "The exponentiated Weibull family: some properties and a flood data application". Commun Stat Theory Meth 25: 3059-3083.

Mudholkar, G.S.; Srivastava, D.K. (1993). "Exponentiated Weibull family for analyzing bathtub failure-ratedata". IEEE Transactions on Reliability 42 (2): 299-302.

Mudholkar, G.S.; Srivastava, D.K.; Freimer, M. (1995). "The exponentiated Weibull family; a reanalysis of the bus motor failure data".
}
\author{
Pedro Rafael Diniz Marinho \email{pedro.rafael.marinho@gmail.com}
}
\note{
Making \eqn{a = 1}, we have the Weibull distribution. Making \eqn{a = 1}, we have the Weibull distribution. For \eqn{k = 1} we have the exponentiated exponential distribution.
}

\seealso{
\code{\link[stats]{rweibull}}
}
\examples{
rew(20,7,2,1.7)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{rweibull}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
