\name{siarmcmc}
\alias{siarmcmc}
\title{ MCMC for stable isotope data }
\description{
  Runs an MCMC on stable isotope data from animals to determine their dietary habits.
}
\usage{
siarmcmc(data, sources, corrections = 0, iterations = 2e+05, burnin = 50000, howmany = 10000, thinby = 15, alpha = rep(1, nrow(sources)), beta = rep(1, nrow(sources)))
}
\arguments{
  \item{data}{ A matrix with each food source as a seperate row and each isotope as a seperate column.}
  \item{sources}{ A matrix containing the mean and standard deviations of the fractionated correction 
  values for each of the isotopes. Also allows corrections = 0 for pre-corrected data. }
  \item{corrections}{ A matrix containing the mean and standard deviations of the fractional correction
  values for each of the isotopes. Also allows corrections = 0 for pre-corrected data.}
  \item{iterations}{ The number of iterations to run. }
  \item{burnin}{ The size of the burnin }
  \item{howmany}{ How often to report the number of iterations. }
  \item{thinby}{ The amount of thinning of the iterations. }
  \item{alpha}{ The beta distribution prior alpha parameters, the default is the uniform for each proportion. }
  \item{beta}{ As above. }
}
\details{
  The model assumes that each target value comes from a Gaussian distribution with an unknown mean and standard 
  deviation. The structure of the mean is a weighted combination of the food sources' isotopic 
  values. The proportional weights are of key interest and are given independent Beta prior distributions. The 
  standard deviation is divided up between the uncertainty around the fractionation corrections (if 
  corrections are given) and the natural variability between target individuals within a defined group 
  (or between all individuals if no grouping structure is specified). The default iterations numbers work well 
  for the demo data sets, but advanced users will want to adjust them to suit their analysis. 

}
\value{
  A parameter matrix consisting of (iterations-burnin)/thinby rows with numgroups*(numsources+numiso) 
  columns, where numsources is the number of food sources, numiso is the number of isotopes, and numgroups 
  is the number of groups. The parameter matrix is structured so that, for each group, the first columns 
  are those of the proportions of each food source eaten, the next columns are the standard deviations for 
  each isotope. This format repeats across rows to each group. The parameters may then subsequently be 
  used for plotting, convergence checks, summaries, etc, etc.
}
\references{ }
\author{ Andrew Parnell }
\note{ }
\seealso{ \code{\link{siarmenu}} }
\examples{
# Should take around 10 seconds to run
#out <- siarmcmc(geese1demo,plantsdemo,correctionsdemo)
}
\keyword{ programming }
\keyword{ models }
\keyword{ multivariate }
