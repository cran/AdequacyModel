\name{kwp} 
\alias{kwp}
\alias{rkwp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pseudo-Random Numbers - Kumaraswamy Pareto
}
\description{
Generates pseudorandom numbers from Kumaraswamy Pareto Distribution.
}
\usage{
rkwp(n, k, beta, a, b)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
  Amount of generated numbers;
}
  \item{k}{
  Shape parameter of the Kumaraswamy Pareto distribution;
}
  \item{beta}{
  Scale parameter Kumaraswamy Pareto distribution;
}
  \item{a}{
  Shape parameter Kumaraswamy Pareto distribution;
}
  \item{b}{
  Shape parameter Kumaraswamy Pareto distribution.
}
}
\details{
The Kw-P (Kumaraswamy Pareto) distribution is not in fact very tractable. However, its heavy tail can adjust skewed data that cannot be properly fitted by existing distributions. Furthermore, the cumulative and hazard rate functions are simple. The Pareto distribution function Kumaraswamy function is given by

\deqn{
F(x; \beta, k, a, b) = 1 - \left\{1 - \left[1 - \left(\frac{\beta}{x}\right)^k\right]^a\right\}^b,
}
where \eqn{\beta>0} (\eqn{x>=\beta}) is a scale parameter and the other positive parameters \eqn{k}, \eqn{a} and \eqn{b} are shape parameters.

}

%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}

\references{
Bourguinon, M.; Silva, R.B.; Zea, L.M. and Cordeiro, G.M. (2013). The Kumaraswamy Pareto distribution. Journal of Statistical Theory and Applications 12, 129-144.
}
\author{
Pedro Rafael Diniz Marinho \email{pedro.rafael.marinho@gmail.com}

Marcelo Bourguignon \email{m.p.bourguignon@gmail.com}
}

\seealso{
\code{\link[stats]{rweibull}}
}
\examples{
beta = 1.5
k = 1.5
a = 1.5
b = 3.5

x<-seq(beta, 10, 0.01)

pdfKWP <- function(k,bet,a,b,x){
a*b*k*bet^(k)*(1-(beta/x)^k)^(a-1)*(1-(1-(beta/x)^k)^(a))^(b-1)/x^(k+1)
}

hist(rkwp(100,k,beta,a,b), freq = FALSE,main="",xlab="x",
     ylab="Density",xlim=c(1.5,7),ylim=c(0,1.5),lwd=1)
lines(x,pdfKWP(k,beta,a,b,x),lty=1,col="black",
      lwd=1,type="l",xlim=c(1.5,7),ylim=c(0,1.5))
box()  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{rweibull}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
