\name{lines.polynomial}
\alias{lines.polynomial}
\title{Lines Method for Polynomials}
\usage{
lines.polynomial(p, length=100, \dots)
}
\arguments{
    \item{p}{An object of class ``polynomial''.}
    \item{length}{Size of vector at which evaluations are to be made.}
    \item{\dots}{Additional arguments as for the \code{\link{lines}}
	generic.}
}
\description{
    Add a polynomial to an existing plot usually as a line plot.
}
\details{
    This is a method for the generic function \code{\link{lines}}.
    
    Lines representing the given polynomial are added to an existing
    plot.  Values outside the current plot region are not shown.
}
\seealso{
    \code{\link{lines}},
    \code{\link{points}},
    \code{\link{points.polynomial}},
    \code{\link{plot}},
    \code{\link{plot.polynomial}}.
}
\examples{
plot  (poly.calc(-1:5))
lines (poly.calc( 2:4), lty = 2)
points(poly.calc(-2:6), pch = 4)
}
\keyword{symbolmath}
