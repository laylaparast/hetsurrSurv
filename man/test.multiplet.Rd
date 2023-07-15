\name{test.multiplet}
\alias{test.multiplet}
\title{
Tests for heterogeneity across multiple timepoints}
\description{
Tests for heterogeneity across multiple timepoints}
\usage{
test.multiplet(t.mult, xone, xzero, deltaone, deltazero, sone, szero, wone, 
wzero, w.grd, landmark, extrapolate = T, h.0 = NULL, h.1 = NULL, h.w = NULL, 
h.s = NULL, h.w.1 = NULL)
}
\arguments{
  \item{t.mult}{
Vector of time points}
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
  \item{landmark}{
t0, landmark time}
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
A list is returned:
\item{pval.multi}{p-value for omnibus test}
\item{pval.con.multi}{p-value for conservative omnibus test}
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
	#computationally intensive
	\donttest{
	test.multiplet(t.mult = c(1,1.25,1.5), xone=example.data$x1, xzero=example.data$x0, 
	deltaone=example.data$d1, deltazero=example.data$d0, sone=log(example.data$s1), 
	szero=log(example.data$s0), wone=log(example.data$w1), wzero=log(example.data$w0),
	 w.grd=log(seq(0.1,0.9, length=25)), landmark=0.5)
	}
}

