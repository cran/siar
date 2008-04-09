\name{siarplotdata}
\alias{siarplotdata}
\title{ Produces plots of target data and sources }
\description{
    Produces colourful scatter plots of siar target data and sources.
}
\usage{
siarplotdata(siardata, siarversion)
}
\details{
  Can be called at any time after running loadsiardata() or when running siarmenu()
}
\arguments{
  \item{siardata}{ A list containing some or all of the following parts: targets, sources, corrections,
  PATH, TITLE, numgroups, numdata, numsources, numiso, SHOULDRUN, GRAPHSONLY, EXIT, and output. For more details
  of these inputs see the \code{\link{siarmenu}} function. }
  \item{siarversion}{ The siar version number as a string. }
}
\author{ Andrew Parnell }
\keyword{ programming}
\keyword{ models}
