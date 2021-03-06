% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SSpsVoigt1.R
\name{SSpsVoigt1}
\alias{SSpsVoigt1}
\title{Self-Starting nls Pseudo-Voigt model}
\usage{
SSpsVoigt1(input,y0, A, nu, xc, w)
}
\arguments{
\item{input}{a numeric vector of values at which to evaluate the model.}

\item{y0}{offset}

\item{A}{area}

\item{nu}{is the mixing parameter, 0 is Gauss and 1 is Lorentz, if nu>1 this is called a super Lorentzian.}

\item{xc}{center}

\item{w}{width}
}
\value{
a numeric vector of the same length as input. It is the value of the expression y0+A*((nu*2*w)/(pi*4*(x-xc)^2+w^2)+(1-nu)*(sqrt(4*log(2))*exp((-4*log(2)*(x-xc)^2)/w^2))/sqrt(pi*w)).
}
\description{
Pseudo-Voigt peak function.
}
\note{
This is primarily intended for use in formulae given to the \code{\link{nls}} function or similar functions for example \code{\link{mle2}}.
}
\examples{
data(smoking)
model <- nlsLM(USA ~ SSpsVoigt1(year,y0,A,nu,xc,w),data=smoking)
}
\references{
qtiplot: \url{https://www.qtiplot.com/doc/manual-en/x9885.html}.
}
\seealso{
\code{\link{selfStart}}, \code{\link{getInitial}}
}
\author{
Krzysztof Trajkowski
}
\keyword{function}
