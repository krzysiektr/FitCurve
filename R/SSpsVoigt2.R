#' @title Self-Starting nls Pseudo-Voigt model
#'
#' @description
#' lorentz peak function.
#' @param input a numeric vector of values at which to evaluate the model.
#' @param y0 offset
#' @param A area
#' @param nu profile shape factor
#' @param xc center
#' @param wl width lorentz
#' @param wg width gauss
#' @usage SSpsVoigt2(input,y0, A, nu, xc, wl, wg)
#' @return a numeric vector of the same length as input. It is the value of the expression y0+A*((nu*2*wl)/(pi*4*(x-xc)^2+wl^2)+(1-nu)*(sqrt(4*log(2))*exp((-4*log(2)*(x-xc)^2)/wg^2))/sqrt(pi*wg)).
#' @note This is primarily intended for use in formulae given to the \code{\link{nls}} function or similar functions for example \code{\link{mle2}}.
#' @keywords function
#' @author
#' Krzysztof Trajkowski
#' @examples
#' data(smoking)
#' model <- nlsLM(USA~SSpsVoigt2(year,y0, A, nu, xc, wl, wg),data=smoking)
#' @rdname SSpsVoigt2
#'
#' @references qtiplot: \url{http://www.qtiplot.com/doc/manual-en/x7261.html#sec-fit-psdvoigt2}.
#'
#' @export
#'
#' @seealso \code{\link{selfStart}}, \code{\link{getInitial}}

SSpsVoigt2 <- selfStart(
function (input, y0, A, nu, xc, wl, wg)
{
    .expr1 <- nu * 2
    .expr2 <- .expr1 * wl
    .expr3 <- pi * 4
    .expr4 <- input - xc
    .expr5 <- .expr4^2
    .expr8 <- .expr3 * .expr5 + wl^2
    .expr10 <- 1 - nu
    .expr11 <- log(2)
    .expr13 <- sqrt(4 * .expr11)
    .expr15 <- -4 * .expr11
    .expr16 <- .expr15 * .expr5
    .expr17 <- wg^2
    .expr19 <- exp(.expr16/.expr17)
    .expr20 <- .expr13 * .expr19
    .expr21 <- .expr10 * .expr20
    .expr22 <- pi * wg
    .expr23 <- sqrt(.expr22)
    .expr25 <- .expr2/.expr8 + .expr21/.expr23
    .expr28 <- 2 * wl
    .expr33 <- 2 * .expr4
    .expr36 <- .expr8^2
    .value <- y0 + A * .expr25
    .actualArgs <- as.list(match.call()[c("y0", "A", "nu", "xc", "w")])
    if(all(unlist(lapply(.actualArgs, is.name))))
    {
    .grad <- array(0, c(length(.value), 6L), list(NULL, c("y0",
        "A", "nu", "xc", "wl", "wg")))
    .grad[, "y0"] <- 1
    .grad[, "A"] <- .expr25
    .grad[, "nu"] <- A * (.expr28/.expr8 - .expr20/.expr23)
    .grad[, "xc"] <- A * (.expr2 * (.expr3 * .expr33)/.expr36 -
        .expr10 * (.expr13 * (.expr19 * (.expr15 * .expr33/.expr17)))/.expr23)
    .grad[, "wl"] <- A * (.expr1/.expr8 - .expr2 * .expr28/.expr36)
    .grad[, "wg"] <- -(A * (.expr10 * (.expr13 * (.expr19 * (.expr16 *
        (2 * wg)/.expr17^2)))/.expr23 + .expr21 * (0.5 * (pi *
        .expr22^-0.5))/.expr23^2))
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
pars <- as.vector(coef(minpack.lm::nlsLM( y ~ y0+A*((nu*2*wl)/(pi*4*(x-xc)^2+wl^2)+(1-nu)*(sqrt(4*log(2))*exp((-4*log(2)*(x-xc)^2)/wg^2))/sqrt(pi*wg)), start=list(y0=min(xy$y),A=(pi*((max(xy$x)-min(xy$x))/3)*max(xy$y)-pi*((max(xy$x)-min(xy$x))/3)*min(xy$y))/2,nu=0.5,xc=xy$x[which.max(xy$y)],wl=(max(xy$x)-min(xy$x))/3,wg=0.5), lower=c(-Inf,-Inf,-Inf,-Inf,0,0),upper=c(Inf,Inf,Inf,Inf,Inf,Inf), data = xy)))
    value <- c( pars[1], pars[2], pars[3], pars[4], pars[5], pars[6] )
    names(value)= mCall[c("y0", "A", "nu", "xc", "wl", "wg")]
    value
  },parameters=c("y0","A","nu","xc","wl","wg"))
