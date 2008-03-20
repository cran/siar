\name{siarproportionbysourceplot}
\alias{siarproportionbysourceplot}
\title{ siar proportion plots by source }
\description{
  Useful for siar data with multiple groups where the variability of each source over time is of interest.
}
\usage{
siarproportionbysourceplot(siardata, siarversion)
}
\arguments{
  \item{siardata}{ A list containing some or all of the following parts: targets, sources, corrections,
  PATH, TITLE, numgroups, numdata, numsources, numiso, SHOULDRUN, GRAPHSONLY, EXIT, and output. For more details
  of these inputs see the \code{\link{siarmenu}} function. }
  \item{siarversion}{ The siar version number as a string. }
}
\author{ Andrew Parnell }
\keyword{ programming}
