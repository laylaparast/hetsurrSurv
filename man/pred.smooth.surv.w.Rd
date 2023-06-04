\name{pred.smooth.surv.w}
\alias{pred.smooth.surv.w}
\title{
Estimates P(T>t|W=w)
}
\description{
Estimates P(T>t|W=w)
}
\usage{
pred.smooth.surv.w(x.ref, delta.ref, w.ref, w.apply, myt, extrapolate = T, h.use = NULL, warn.support = FALSE)
}
\arguments{
  \item{x.ref}{
x reference}
  \item{delta.ref}{
delta reference}
  \item{w.ref}{
w reference}
  \item{w.apply}{
w to apply estimate to}
  \item{myt}{
time of interest}
  \item{extrapolate}{
TRUE or FALSE}
  \item{h.use}{
bandwidth
}
  \item{warn.support}{
TRUE or FALSE}
}
\value{
Conditional survival estimates}
\keyword{internal}