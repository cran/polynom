\name{deriv.polynomial}
\alias{deriv.polynomial}
\title{Differentiate a Polynomial}
\usage{
deriv.polynomial(p)
}
\arguments{
    \item{p}{Object of class polynomial}
}
\description{
    Calculates the derivative of a univariate polynomial.
}
\details{
    This is a method for the generic function \code{\link{deriv}}.
}
\value{
    Derivative of the polynomial.
}
\seealso{
    \code{\link{integral.polynomial}},
    \code{\link{deriv}}.
}
\examples{
pr <- poly.calc(1:5)
pr
## -120 + 274*x - 225*x^2 + 85*x^3 - 15*x^4 + x^5
deriv(pr)
##  274 - 450*x + 255*x^2 - 60*x^3 + 5*x^4
}
\keyword{symbolmath}
