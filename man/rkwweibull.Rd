\name{rkwweibull}
\alias{rkwweibull}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pseudo-Random Numbers - Kumaraswamy Weibull
}
\description{
Generates pseudorandom numbers from Kumaraswamy Weibull Distribution.
}
\usage{
rkwweibull(n, a, b, c, beta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
  Amount of generated numbers;
}
  \item{a}{
  Shape parameter of the Kumaraswamy Weibull Distribution;
}
  \item{b}{
  Shape parameter of the Kumaraswamy Weibull Distribution;
}
  \item{c}{
  Shape parameter of the Kumaraswamy Weibull Distribution;
}
  \item{beta}{
  Scale parameter of the Kumaraswamy Weibull Distribution.
}

}
\details{
The Kumaraswamy Weibull Distribution function is given by:

\deqn{
F(x,a,b,c) = 1 - \{1 - [1-e^{-(\beta*x)^c}]^a\}^b
}
The distribution was proposed by Cordeiro, Ortega and Nadarajah in 2010 in the Journal of the Franklin Institute. The distribution is defined on the positive real (\eqn{x>0}). All parameters are also defined in the positive reals.
}
%\value{
%  ~Describe the value returned
%  If it is a LIST, use
%  \item{comp1 }{Description of 'comp1'}
%  \item{comp2 }{Description of 'comp2'}
% ...
%}

\references{
Cordeiro, G.M.; Ortega, E.M.M., Nadarajah, S. (2010). "The Kumaraswamy Weibull distribution with application to failure data". Journal of the Franklin Institute, 347, 1399-1429.
}
\author{
Pedro Rafael Diniz Marinho \email{pedro.rafael.marinho@gmail.com}
}

%\note{
%sdfg
%}

%\seealso{
%sadfgh
%}

\examples{
rkwweibull(20,2,2.3,1.5,2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
