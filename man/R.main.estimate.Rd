\name{R.main.estimate}
\alias{R.main.estimate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimates the proportion of treatment effect explained by the surrogate marker as a function of a baseline covariate}
\description{
Assesses heterogeneity in the utility of a surrogate marker with respect to a baseline covariate using censored (survival data) by estimates the proportion of treatment effect explained by the surrogate marker as a function of a baseline covariate, w }
\usage{
R.main.estimate(xone, xzero, deltaone, deltazero, sone, szero, wone, wzero, w.grd, myt, 
landmark, type = "cont", var = FALSE, test = FALSE, extrapolate = T, h.0 = NULL, 
h.1 = NULL, h.w = NULL, h.s = NULL, h.w.1 = NULL)
}
\arguments{
  \item{xone}{
x1, observed event time in the treated group}
  \item{xzero}{
x0, observed event time in the control group}
  \item{deltaone}{
delta1, event indicator in the treated group}
  \item{deltazero}{
delta0, event indicator in the control group}
  \item{sone}{
s1, surrogate marker in the treated group}  
\item{szero}{
s0, surrogate marker in the control group}  
  \item{wone}{
w1, baseline covariate in the treated group}
  \item{wzero}{
w0, baseline covariate in the control group}
  \item{w.grd}{
grid for w where estimation will be provided}
  \item{myt}{
t of interest}
  \item{landmark}{
t0, landmark time}
 \item{type}{
options are "cont" or "discrete"; type of baseline covariate, default is "cont"}
  \item{var}{
TRUE or FALSE, if variance/standard error estimates are wanted}
  \item{test}{
TRUE or FALSE, if test for heterogeneity is wanted wanted
}
  \item{extrapolate}{
TRUE or FALSE
}
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
A list is returned:
  \item{w.values}{grid for w where estimation is provided}
  \item{R.s.w}{The propoportion of treatment effect explained as a function of the baseline covariate, w}
  \item{delta.w}{The treatment effect as a function of the baseline covariate, w}
  \item{delta.s.w}{The residual treatment effect as a function of the baseline covariate, w}
  \item{sd.R}{Standard error estimate of R.s.w}
  \item{sd.delta}{Standard error estimate of delta.w}
  \item{sd.delta.s}{Standard error estimate of delta.s.w}
  \item{pval.omnibus}{p-value from the omnibus test for heterogeneity}
  \item{pval.con}{p-value from the conservative omnibus test for heterogeneity}
}
\references{
Parast, L., Tian L, Cai, T. (2023) "Assessing Heterogeneity in Surrogacy Using Censored Data." Under Review.
}
\author{
Layla Parast
}
\examples{
	data(example.data)
	names(example.data)
	R.main.estimate(xone=example.data$x1, xzero=example.data$x0, deltaone=example.data$d1, 
	deltazero=example.data$d0, sone=log(example.data$s1), szero=log(example.data$s0), 
	wone=log(example.data$w1), wzero=log(example.data$w0), 
	w.grd=log(seq(0.1,0.9, length=25)), myt=1, landmark=0.5)
	R.main.estimate(xone=example.data$x1, xzero=example.data$x0, deltaone=example.data$d1, 
	deltazero=example.data$d0, sone=log(example.data$s1), szero=log(example.data$s0), 
	wone=example.data$w1_cat, wzero=example.data$w0_cat, 
	myt=1, landmark=0.5,type = "discrete", w.grd = c(1,2,3,4))
	#computationally intensive
	\donttest{
	R.main.estimate(xone=example.data$x1, xzero=example.data$x0, deltaone=example.data$d1,
	deltazero=example.data$d0, sone=log(example.data$s1), szero=log(example.data$s0), 
	wone=log(example.data$w1), wzero=log(example.data$w0), 
	w.grd=log(seq(0.1,0.9, length=25)), myt=1, landmark=0.5, test=TRUE)
	}
}
