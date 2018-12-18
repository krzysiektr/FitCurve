#' @title Self-Starting nls Pseudo-Voigt model
#'
#' @description
#' Pseudo-Voigt peak function.
#' @param input a numeric vector of values at which to evaluate the model.
#' @param y0 offset
#' @param A area
#' @param nu is the mixing parameter, 0 is Gauss and 1 is Lorentz, if nu>1 this is called a super Lorentzian.
#' @param xc center
#' @param w width
#' @usage SSpsVoigt1(input,y0, A, nu, xc, w)
#' @return a numeric vector of the same length as input. It is the value of the expression y0+A*((nu*2*w)/(pi*4*(x-xc)^2+w^2)+(1-nu)*(sqrt(4*log(2))*exp((-4*log(2)*(x-xc)^2)/w^2))/sqrt(pi*w)).
#' @note This is primarily intended for use in formulae given to the \code{\link{nls}} function or similar functions for example \code{\link{mle2}}.
#' @keywords function
#' @author
#' Krzysztof Trajkowski
#' @examples
#' data(smoking)
#' model <- nlsLM(USA ~ SSpsVoigt1(year,y0,A,nu,xc,w),data=smoking)
#' @rdname SSpsVoigt1
#'
#' @references qtiplot: \url{http://www.qtiplot.com/doc/manual-en/x7261.html#sec-fit-psdvoigt1}.
#'
#' @export
#'
#' @seealso \code{\link{selfStart}}, \code{\link{getInitial}}

SSpsVoigt1 <- selfStart(
function (input, y0, A, nu, xc, w)
{
    .expr1 <- nu * 2
    .expr2 <- .expr1 * w
    .expr3 <- pi * 4
    .expr4 <- input - xc
    .expr5 <- .expr4^2
    .expr7 <- w^2
    .expr8 <- .expr3 * .expr5 + .expr7
    .expr10 <- 1 - nu
    .expr11 <- log(2)
    .expr13 <- sqrt(4 * .expr11)
    .expr15 <- -4 * .expr11
    .expr16 <- .expr15 * .expr5
    .expr18 <- exp(.expr16/.expr7)
    .expr19 <- .expr13 * .expr18
    .expr20 <- .expr10 * .expr19
    .expr21 <- pi * w
    .expr22 <- sqrt(.expr21)
    .expr24 <- .expr2/.expr8 + .expr20/.expr22
    .expr27 <- 2 * w
    .expr32 <- 2 * .expr4
    .expr35 <- .expr8^2
    .value <- y0 + A * .expr24
    .actualArgs <- as.list(match.call()[c("y0", "A", "nu", "xc", "w")])
    if(all(unlist(lapply(.actualArgs, is.name))))
    {
    .grad <- array(0, c(length(.value), 5L), list(NULL, c("y0",
        "A", "nu", "xc", "w")))
    .grad[, "y0"] <- 1
    .grad[, "A"] <- .expr24
    .grad[, "nu"] <- A * (.expr27/.expr8 - .expr19/.expr22)
    .grad[, "xc"] <- A * (.expr2 * (.expr3 * .expr32)/.expr35 -
        .expr10 * (.expr13 * (.expr18 * (.expr15 * .expr32/.expr7)))/.expr22)
    .grad[, "w"] <- A * (.expr1/.expr8 - .expr2 * .expr27/.expr35 -
        (.expr10 * (.expr13 * (.expr18 * (.expr16 * .expr27/.expr7^2)))/.expr22 +
            .expr20 * (0.5 * (pi * .expr21^-0.5))/.expr22^2))
    attr(.value, "gradient") <- .grad
    }
    .value
},
function(mCall, data, LHS)
          {
xy <- data.frame(sortedXyData(mCall[["input"]], LHS, data))
if(nrow(xy) < 2) {
stop("Too few distinct input values to fit a model")
}
pars <- as.vector(coef(minpack.lm::nlsLM( y ~ y0+A*((nu*2*w)/(pi*4*(x-xc)^2+w^2)+(1-nu)*(sqrt(4*log(2))*exp((-4*log(2)*(x-xc)^2)/w^2))/sqrt(pi*w)), start=list(y0=min(xy$y),A=(pi*((max(xy$x)-min(xy$x))/3)*max(xy$y)-pi*((max(xy$x)-min(xy$x))/3)*min(xy$y))/2,nu=0.5,xc=xy$x[which.max(xy$y)],w=(max(xy$x)-min(xy$x))/3), lower=c(-Inf,-Inf,-Inf,-Inf,0),upper=c(Inf,Inf,Inf,Inf,Inf), data = xy)))
    value <- c( pars[1], pars[2], pars[3], pars[4], pars[5] )
    names(value)= mCall[c("y0", "A", "nu", "xc", "w")]
    value
  },parameters=c("y0","A","nu","xc","w"))
