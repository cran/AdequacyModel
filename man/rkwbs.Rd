\name{kwbs} 
\alias{rkwbs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pseudo-Random Numbers - Kumaraswamy Birnbaum-Saunders Distribution
}
\description{
Generates pseudorandom numbers from Kumaraswamy Birnbaum-Saunders Distribution.
}
\usage{
rkwbs(n, alpha, beta, a, b)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
  Amount of generated numbers;
}
  \item{alpha}{
  Shape parameter of the Kumaraswamy Birnbaum-Saunders distribution;
}
  \item{beta}{
  Scale parameter Kumaraswamy Birnbaum-Saunders distribution;
}
  \item{a}{
  Shape parameter Kumaraswamy Birnbaum-Saunders distribution;
}
  \item{b}{
  Shape parameter Kumaraswamy Birnbaum-Saunders distribution.
}
}
\details{
The additional parameters introduced by the Kumaraswamy generalization are sought as a
manner to furnish a more flexible distribution. The Kw-BS distribution is expected to have
immediate application in reliability and survival studies. The BS distribution arises as the
basic exemplar for \eqn{a = b = 1}. The exponentiated Birnbaum-Saunders (EBS) distribution
corresponds to \eqn{b = 1}. The Kw-BS is a particular case of the Mc-BS distribution when
\eqn{a = c}; see Cordeiro et al. (2011).

\deqn{
F(x; \alpha, \beta, a, b) = 1 - \{1 - \Phi(\nu)^a\}^b
}
where \eqn{\beta > 0} is a scale parameter and the other positive parameters \eqn{\alpha}, \eqn{a}, and \eqn{b} are shape parameters in that \eqn{\nu = \alpha^{-1}\rho(x/\beta)}, \eqn{\rho(z) = z^{1/2} - z^{-1/2}} and \eqn{\Phi(\cdotp)} denotes the standard normal distribution function.
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
}

\references{
Saulo, H.; Leao, J. and Bourguignon, M. The Kumaraswamy Birnbaum-Saunders Distribution. Journal of Statistical Theory and Practice 6, 745-759.

Cordeiro, G. M., and M. Castro. 2011. A new family of generalized distributions. J. Stat. Comput.
Simulation, 81(7), 883-898.

Cordeiro, G.M. Lemonte, A.J. and Ortega, E.M.M. (2013). An extended fatigue life distribution. Statistics 47, 626-653.
}
\author{
Pedro Rafael Diniz Marinho \email{pedro.rafael.marinho@gmail.com}

Marcelo Bourguignon \email{m.p.bourguignon@gmail.com}
}

%\seealso{
%\code{\link[stats]{rweibull}}, \code{\link[AdequacyModel]{rweibull}}, \code{\link[AdequacyModel]{rwp}}
%}
\examples{
a = 2.0
b = 2.5
alpha = 1.0
beta = 1.0
rkwbs(5, alpha, beta, a, b)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{rweibull}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
