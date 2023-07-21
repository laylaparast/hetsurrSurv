\name{R.surv.s.w}
\alias{R.surv.s.w}
\title{
Estimates the propoportion of treatment effect explained as a function of a baseline covariate, w
}
\description{
Estimates the propoportion of treatment effect explained as a function of a baseline covariate, w
}
\usage{
R.surv.s.w(xone, xzero, deltaone, deltazero, sone, szero, wone, wzero, w.grd, myt, landmark, extrapolate = T, h.0 = NULL, h.1 = NULL, h.w = NULL, h.s = NULL, h.w.1 = NULL)
}
\arguments{
  \item{xone}{
x1}
  \item{xzero}{
x0}
  \item{deltaone}{
delta1}
  \item{deltazero}{
delta0}
  \item{sone}{
s1}
  \item{szero}{
s0}
  \item{wone}{
w1}
  \item{wzero}{
w0}
  \item{w.grd}{
grid for w}
  \item{myt}{
t of interest}
  \item{landmark}{
landmark time}
  \item{extrapolate}{
TRUE or FALSE}
  \item{h.0}{
bandwidth}
  \item{h.1}{
bandwidth}
  \item{h.w}{
bandwidth}
  \item{h.s}{
bandwidth}
  \item{h.w.1}{
bandwidth}
}

\value{
A list with the following components:
  \item{R.s.w}{The propoportion of treatment effect explained as a function of the baseline covariate, w}
  \item{delta.w}{The treatment effect as a function of the baseline covariate, w}
  \item{delta.s.w}{The residual treatment effect as a function of the baseline covariate, w}

}
\keyword{internal}
