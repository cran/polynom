\name{polynomial}
\alias{polynomial}
\alias{as.polynomial}
\alias{is.polynomial}
\alias{as.character.polynomial}
\alias{print.polynomial}
\alias{coef.polynomial}
\title{Polynomials}
\usage{
polynomial(coef = c(0, 1))
as.polynomial(p)
is.polynomial(p)
}
\arguments{
  \item{coef}{numeric vector, giving the polynomial coefficients in
    \emph{in}creasing order.}
  \item{p}{an arbitrary \R object.}
}
\description{
  Construct, coerce to and test for polynomial objects.
}
\details{
  \code{polynomial} constructs a polynomial from its coefficients,
  i.e., \code{p[1:k]} specifies the polynomial
  \deqn{p_1 + p_2 x + p_3 x^2 + \dots + p_k x^{k-1}}%
  {p[1] + p[2]* x + p[3]* x^2 + ... + p[k]* x^(k-1)}.
  Internally, polynomials are simply numeric coefficient vectors of
  class \code{"polynomial"}.  Several useful methods are available for
  this class, such as coercion to character and function, extraction of
  the coeeficients, printing and plotting.

  \code{as.polynomial} tries to coerce its arguments to a polynomial.

  \code{is.polynomial} tests whether its argument is a polynomial (in
  the sense that it has class \code{"polynomial"}.
}
\examples{
polynomial(1:4)
as.polynomial(c(1,0,3,0))
polynomial(c(2,rep(0,10),1))
}
\keyword{symbolmath}
