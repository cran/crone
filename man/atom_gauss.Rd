% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/direct_space.R
\name{atom_gauss}
\alias{atom_gauss}
\title{Gaussian atom}
\usage{
atom_gauss(x, a, x0 = 0, Z = 1, B = 0, k = ksigma)
}
\arguments{
\item{x}{Point in the 1D cell at which this function is calculated.}

\item{a}{A real number. The width of the unit cell in which the gaussian
atom is placed.}

\item{x0}{A real number. The point corresponding to the atom's peak.}

\item{Z}{An integer number. Z is the atomic number of the atom (Z(H)=1,
Z(He)=2,Z(Li)=3,Z(B)=4, etc).}

\item{B}{A real number. This is the B factor characterizing the atom's
thermal agitation. It is given as B=8*pi^2*U, where U is the variance
of the position of the atoms' nucleus around the equilibrium position.}

\item{k}{A real number. It controls the standard deviation of the 
gaussian function describing the atom and, thus, the shape of the
associated peak. The standard deviation sigma is given by:
         \code{sigma = k * sqrt(Z)}}
}
\value{
A vector of length equal to the length of vector x, with
values equal to the evaluated gaussian atom.
}
\description{
Gaussian atom
}
\examples{
# Carbon gaussian atom in the middle of a cell
a <- 10
x0 <- 5
Z <- 6
x <- seq(0,a,length=1000)
rho <- atom_gauss(x,a,x0,Z)
plot(x,rho,type="l",xlab="x",ylab=expression(rho))
}
