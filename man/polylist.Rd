\name{polylist}
\alias{polylist}
\alias{deriv.polylist}
\alias{integral.polylist}
\alias{plot.polylist}
\alias{print.polylist}
\title{Lists of Polynomials}
\description{Create and manipulate lists of polynomials.}
\usage{polylist(\dots)}
\arguments{
  \item{\dots}{a list of \R objects.}
}
\details{
  \code{polylist} takes a list of arguments, tries to convert each into
  a polynomial (see \code{\link{polynomial}}), and sets the class of the
  list to \code{"polylist"}.

  This class has several useful methods, such as taking derivatives
  (\code{link{deriv}}) and antiderivatives (\code{\link{integral}}),
  printing and plotting.
}
\examples{
## Calculate orthogonal polynomials
pl <- poly.orth(rep(1:4, 1:4), 3)
pl
plot(pl)
deriv(pl)
integral(pl)
}
\keyword{symbolmath}
