\name{delta.s.surv.estimate.new}
\alias{delta.s.surv.estimate.new}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Calculates robust residual treatment effect accounting for surrogate marker information measured at a specified time and primary outcome information up to that specified time
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function calculates the robust estimate of the residual treatment effect accounting for surrogate marker information measured at \eqn{t_0} and primary outcome information up to \eqn{t_0} i.e. the hypothetical treatment effect if both the surrogate marker distribution at \eqn{t_0} and survival up to \eqn{t_0} in the treatment group look like the surrogate marker distribution and survival up to \eqn{t_0} in the control group. Ideally this function is only used as a helper function and is not directly called. 

}
\usage{
delta.s.surv.estimate.new(xone, xzero, deltaone, deltazero, sone, szero, t, weight.perturb = NULL, landmark, extrapolate = FALSE, transform = FALSE, approx = T, warn.te = FALSE, warn.support = FALSE)
}
\arguments{
  \item{xone}{
%%     ~~Describe \code{xone} here~~
numeric vector, the observed event times in the treatment group, X = min(T,C) where T is the time of the primary outcome and C is the censoring time.
}
  \item{xzero}{
%%     ~~Describe \code{xzero} here~~
numeric vector, the observed event times in the control group, X = min(T,C) where T is the time of the primary outcome and C is the censoring time.
}
  \item{deltaone}{
%%     ~~Describe \code{deltaone} here~~
numeric vector, the event indicators for the treatment group, D = I(T<C) where T is the time of the primary outcome and C is the censoring time.

}
  \item{deltazero}{
%%     ~~Describe \code{deltazero} here~~
numeric vector, the event indicators for the control group, D = I(T<C) where T is the time of the primary outcome and C is the censoring time.
}
  \item{sone}{
%%     ~~Describe \code{sone} here~~
numeric vector; surrogate marker measurement at \eqn{t_0} for treated observations, assumed to be continuous. If \eqn{X_{1i}<t_0}, then the surrogate marker measurement should be NA. 

}
  \item{szero}{
%%     ~~Describe \code{szero} here~~
numeric vector; surrogate marker measurement at \eqn{t_0} for control observations, assumed to be continuous. If \eqn{X_{1i}<t_0}, then the surrogate marker measurement should be NA.
}
  \item{t}{
%%     ~~Describe \code{t} here~~
the time of interest.
}
  \item{weight.perturb}{
%%     ~~Describe \code{weight.perturb} here~~
weights used for perturbation resampling.
}
  \item{landmark}{
%%     ~~Describe \code{landmark} here~~
the landmark time \eqn{t_0} or time of surrogate marker measurement.
}
  \item{extrapolate}{
%%     ~~Describe \code{extrapolate} here~~
TRUE or FALSE; indicates whether the user wants to use extrapolation.
}
  \item{transform}{
%%     ~~Describe \code{transform} here~~
TRUE or FALSE; indicates whether the user wants to use a transformation for the surrogate marker. 
}
  \item{approx}{
TRUE or FALSE indicating whether an approximation should be used when calculating the probability of censoring; most relevant in settings where the survival time of interest for the primary outcome is greater than the last observed event but before the last censored case, default is TRUE.
}
  \item{warn.te}{
%%     ~~Describe \code{type} here~~
value passed from R.s.estimate function to control warnings; user does not need to specify.
}
  \item{warn.support}{
%%     ~~Describe \code{type} here~~
value passed from R.s.estimate function to control warnings; user does not need to specify.
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
	Details are included in the documentation for R.s.surv.estimate.

}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\eqn{\hat{\Delta}_S(t,t_0)}, the robust residual treatment effect estimate accounting for surrogate marker information measured at \eqn{t_0} and primary outcome information up to \eqn{t_0}.
}
\references{
%% ~put references to the literature/web site here ~
Parast, L., Cai, T., & Tian, L. (2017). Evaluating surrogate marker information using censored data. Statistics in Medicine, 36(11), 1767-1782.
}

\keyword{internal}