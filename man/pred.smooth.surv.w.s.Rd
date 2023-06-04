\name{pred.smooth.surv.w.s}
\alias{pred.smooth.surv.w.s}
\title{
Estimates P(T>t|T>t_0, S=s, W=w)}
\description{
Estimates P(T>t|T>t_0, S=s, W=w)}
\usage{
pred.smooth.surv.w.s(x.ref, delta.ref, w.ref, s.ref, w.apply, s.apply, h.s, h.w, myt, extrapolate = T, kerni.ss.s, tmpind, tj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x.ref}{
x reference}
  \item{delta.ref}{
delta reference}
  \item{w.ref}{
w reference}
  \item{s.ref}{
s reference}
  \item{w.apply}{
w to apply estimate to}
  \item{s.apply}{
s to apply estimate to}
  \item{h.s}{
bandwidth}
  \item{h.w}{
bandwidth}
  \item{myt}{
time of interest}
  \item{extrapolate}{
TRUE or FALSE}
  \item{kerni.ss.s}{
kernel matrix}
  \item{tmpind}{
temporary indicator}
  \item{tj}{
temporary tj}
}
\value{
Conditional survival estimate
}
\keyword{internal}