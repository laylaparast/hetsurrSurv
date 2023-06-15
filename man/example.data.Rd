\name{example.data}
\alias{example.data}
\docType{data}
\title{
Example data
}
\description{
Example data
}
\usage{data("example.data")}
\format{
 A list with 8 elements representing 2000 observations from a treatment group and 1500 observations from a control group:
  \describe{
    \item{\code{s1}}{the surrogate marker in the treatment group}
    \item{\code{s0}}{the surrogate marker in the control group}
    \item{\code{w1}}{the baseline covariate of interest in the treatment group}
    \item{\code{w0}}{the baseline covariate of interest in the control group}
    \item{\code{d1}}{the event indicator in the treatment group}
    \item{\code{d0}}{the event indicator in the control group}
    \item{\code{x1}}{the observed event time in the treatment group}
    \item{\code{x0}}{the observed evenet time in the control group}
}
}
\examples{
data(example.data)
names(example.data)
}

