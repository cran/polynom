\name{summary.polynomial}
\alias{summary.polynomial}
\alias{print.summary.polynomial}
\title{Summary of a Polynomial}
\usage{
summary.polynomial(p)
}
\arguments{
    \item{p}{A polynomial object.}
}
\value{
    A list of class \code{summary.polynomial} (which has its own
    \code{print} method) containing information on zeros, stationary and
    inflexion points.
}
\description{
    Summarize a polynomial by describing its ``key'' points.
}
\details{
    This is a method for the generic function \code{\link{summary}}.
}
}
\examples{
p <- polynomial(6:1)
p
## 6 + 5*x + 4*x^2 + 3*x^3 + 2*x^4 + x^5 
pz <- summary(p)
pz
## [1] -1.49180+0.0000i -0.80579-1.2229i -0.80579+1.2229i
## [4]  0.55169-1.2533i  0.55169+1.2533i
## To retrieve the original polynomial from the zeros:
poly.calc(pz)
## Warning: imaginary parts discarded in coercion
## 6 + 5*x + 4*x^2 + 3*x^3 + 2*x^4 + x^5
}
\keyword{symbolmath}
