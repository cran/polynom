\name{plot.polynomial}
\alias{plot.polynomial}
\title{Plot Method for Polynomials}
\description{
  Plots polynomials, optionally allowing the ``interesting'' region to
  be automatically determined.
}
\usage{
\method{plot}{polynomial}(p, xlim = 0:1, ylim = range(Px), type = "l",
                len = 100, \dots)
}
\arguments{
  \item{p}{an object of class \code{"polynomial"}.}
  \item{xlim}{the range to be encompased by the x axis.}
  \item{ylim}{the range to be encompased by the y axis.}    
  \item{type}{as for \code{\link{plot}()}.}
  \item{len}{number of x points drawn.}
  \item{\dots}{additional arguments as for \code{\link{plot}}.}
}
\details{
  This is a method for the generic function \code{\link{plot}}.
    
  A plot of the polynomial is produced on the currently active device.
  Unless otherwise specified, the domain is chosen to enclose the real
  parts of all zeros, stationary points and zero itself.
}
\seealso{
  \code{\link{plot}},
  \code{\link{lines}},
  \code{\link{points}},
  \code{\link{lines.polynomial}},
  \code{\link{points.polynomial}}.
}
\examples{
plot(p <- poly.calc(-1:5))
}
\keyword{symbolmath}
