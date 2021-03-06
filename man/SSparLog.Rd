% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FitCurve.R
\name{SSparLog}
\alias{SSparLog}
\title{Self-Starting nls parabola logarithmic model}
\usage{
SSparLog(input,a,b,c)
}
\arguments{
\item{input}{a numeric vector of values at which to evaluate the model.}

\item{a}{estimated model parameter}

\item{b}{estimated model parameter}

\item{c}{estimated model parameter}
}
\value{
a numeric vector of the same length as input. It is the value of the expression a*x^(b+c*log(x)).
}
\description{
Computes the parabola logarithmic function
}
\note{
This is primarily intended for use in formulae given to the \code{\link{nls}} function or similar functions for example \code{\link{mle2}}.
}
\examples{
x <- 1:20
y <- 5.75*x^(0.75+0.95*log(x))
model <- nls(y ~ SSparLog(x,a,b,c),algorithm = "port")
}
\references{
Goryl A., Jedrzejczyk Z., Kukula K. (Eds.), Osiewalski J., Walkosz A. (2004) "Wprowadzenie do ekonometrii w przykladach i zadaniach", Wydawnictwo Naukowe PWN, Warszawa.

Borkowski B., Dudek H., Szczesny W. (2004) "Ekonometria. Wybrane zagadnienia" , Wydawnictwo Naukowe PWN, Warszawa.
}
\seealso{
\code{\link{selfStart}}, \code{\link{getInitial}}
}
\author{
Krzysztof Trajkowski
}
\keyword{function}
