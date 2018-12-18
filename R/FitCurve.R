#' @title Self-Starting nls Tornquist I model
#'
#' @description
#' Computes the Tornquist I function
#' @param input a numeric vector of values at which to evaluate the model.
#' @param a estimated model parameter
#' @param b estimated model parameter
#' @usage SStorn1(input,a,b)
#' @return a numeric vector of the same length as input. It is the value of the expression a*x/(x+b).
#' @note This is primarily intended for use in formulae given to the \code{\link{nls}} function or similar functions for example \code{\link{mle2}}.
#' @keywords function
#' @import bbmle
#' @import minpack.lm
#' @author
#' Krzysztof Trajkowski
#' @examples
#' x <- 1:15
#' y <- 254*x/(x+698)
#' model <- nls(y~SStorn1(x,a,b),algorithm = "port")
#' @rdname SStorn1
#'
#' @references Goryl A., Jedrzejczyk Z., Kukula K. (Eds.), Osiewalski J., Walkosz A. (2004) "Wprowadzenie do ekonometrii w przykladach i zadaniach", Wydawnictwo Naukowe PWN, Warszawa.
#' @references Borkowski B., Dudek H., Szczesny W. (2004) "Ekonometria. Wybrane zagadnienia" , Wydawnictwo Naukowe PWN, Warszawa.
#'
#' @export
#'
#' @seealso \code{\link{SSmicmen}}, \code{\link{selfStart}}, \code{\link{getInitial}}

SStorn1 <- selfStart(
  function (input, a, b)
  {
    .expr1 <- a * input
    .expr2 <- input + b
    .value <- .expr1/.expr2
    .actualArgs <- as.list(match.call()[c("a", "b")])
    if(all(unlist(lapply(.actualArgs, is.name))))
    {
      .grad <- array(0, c(length(.value), 2L), list(NULL, c("a",
                                                            "b")))
      .grad[, "a"] <- input/.expr2
      .grad[, "b"] <- -(.expr1/.expr2^2)
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
    pars <- as.vector(coef(lm(I(1/y) ~ I(1/x), data = xy)))
    value <- c(1/pars[1],pars[2]/pars[1])
    names(value)= mCall[c("a", "b")]
    value
  },parameters=c("a","b"))

#' @title Self-Starting nls Tornquist II model
#'
#' @description
#' Computes the Tornquist II function
#' @param input a numeric vector of values at which to evaluate the model.
#' @param a estimated model parameter
#' @param b estimated model parameter
#' @param c estimated model parameter
#' @usage SStorn2(input,a,b,c)
#' @return a numeric vector of the same length as input. It is the value of the expression a*(x-c)/(x+b).
#' @note This is primarily intended for use in formulae given to the \code{\link{nls}} function or similar functions for example \code{\link{mle2}}.
#' @keywords function
#' @author
#' Krzysztof Trajkowski
#' @examples
#' x <- 1:15
#' y <- 325*(x-145)/(x+875)
#' model <- nls(y~SStorn2(x,a,b,c),algorithm = "port")
#' @rdname SStorn2
#'
#' @references Goryl A., Jedrzejczyk Z., Kukula K. (Eds.), Osiewalski J., Walkosz A. (2004) "Wprowadzenie do ekonometrii w przykladach i zadaniach", Wydawnictwo Naukowe PWN, Warszawa.
#' @references Borkowski B., Dudek H., Szczesny W. (2004) "Ekonometria. Wybrane zagadnienia" , Wydawnictwo Naukowe PWN, Warszawa.
#'
#' @export
#'
#' @seealso \code{\link{selfStart}}, \code{\link{getInitial}}

SStorn2 <- selfStart(
  function (input, a, b, c)
  {
    .expr1 <- input - c
    .expr2 <- a * .expr1
    .expr3 <- input + b
    .value <- .expr2/.expr3
    .actualArgs <- as.list(match.call()[c("a", "b", "c")])
    if(all(unlist(lapply(.actualArgs, is.name))))
    {
      .grad <- array(0, c(length(.value), 3L), list(NULL, c("a",
                                                            "b", "c")))
      .grad[, "a"] <- .expr1/.expr3
      .grad[, "b"] <- -(.expr2/.expr3^2)
      .grad[, "c"] <- -(a/.expr3)
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
    pars <- as.vector(coef(lm( y ~ I(y/x)+I(1/x), data = xy)))
    value <- c(pars[1], -pars[2], -pars[3]/pars[1])
    names(value)= mCall[c("a", "b", "c")]
    value
  },parameters=c("a","b","c"))

#' @title Self-Starting nls Tornquist III model
#'
#' @description
#' Computes the Tornquist III function
#' @param input a numeric vector of values at which to evaluate the model.
#' @param a estimated model parameter
#' @param b estimated model parameter
#' @param c estimated model parameter
#' @usage SStorn3(input,a,b,c)
#' @return a numeric vector of the same length as input. It is the value of the expression a*x*(x-c)/(x+b).
#' @note This is primarily intended for use in formulae given to the \code{\link{nls}} function or similar functions for example \code{\link{mle2}}.
#' @keywords function
#' @author
#' Krzysztof Trajkowski
#' @examples
#' x <- 1:15
#' y <- 120*x*(x-254)/(x+654)
#' model <- nls(y~SStorn3(x,a,b,c),algorithm = "port")
#' @rdname SStorn3
#'
#' @references Goryl A., Jedrzejczyk Z., Kukula K. (Eds.), Osiewalski J., Walkosz A. (2004) "Wprowadzenie do ekonometrii w przykladach i zadaniach", Wydawnictwo Naukowe PWN, Warszawa.
#' @references Borkowski B., Dudek H., Szczesny W. (2004) "Ekonometria. Wybrane zagadnienia" , Wydawnictwo Naukowe PWN, Warszawa.
#'
#' @export
#'
#' @seealso \code{\link{selfStart}}, \code{\link{getInitial}}

SStorn3 <- selfStart(
  function (input, a, b, c)
  {
    .expr1 <- a * input
    .expr2 <- input - c
    .expr3 <- .expr1 * .expr2
    .expr4 <- input + b
    .value <- .expr3/.expr4
    .actualArgs <- as.list(match.call()[c("a", "b", "c")])
    if(all(unlist(lapply(.actualArgs, is.name))))
    {
      .grad <- array(0, c(length(.value), 3L), list(NULL, c("a",
                                                            "b", "c")))
      .grad[, "a"] <- input * .expr2/.expr4
      .grad[, "b"] <- -(.expr3/.expr4^2)
      .grad[, "c"] <- -(.expr1/.expr4)
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
    pars <- as.vector(coef(lm( y ~ x+I(y/x), data = xy)))
    value <- c(pars[2], -pars[3], pars[1]/(-pars[2]))
    names(value)= mCall[c("a", "b", "c")]
    value
  },parameters=c("a","b","c"))

#' @title Self-Starting nls Working model
#'
#' @description
#' Computes the Working function
#' @param input a numeric vector of values at which to evaluate the model.
#' @param a estimated model parameter
#' @param b estimated model parameter
#' @usage SSworking(input,a,b)
#' @return a numeric vector of the same length as input. It is the value of the expression exp(a+b*(1/x)).
#' @note This is primarily intended for use in formulae given to the \code{\link{nls}} function or similar functions for example \code{\link{mle2}} and \code{\link{nlsLM}}.
#' @keywords function
#' @author
#' Krzysztof Trajkowski
#' @examples
#' x <- 1:15
#' y <- exp(144+3*(1/x))
#' model <- nlsLM(y~SSworking(x,a,b))
#' @rdname SSworking
#'
#' @references Goryl A., Jedrzejczyk Z., Kukula K. (Eds.), Osiewalski J., Walkosz A. (2004) "Wprowadzenie do ekonometrii w przykladach i zadaniach", Wydawnictwo Naukowe PWN, Warszawa.
#' @references Borkowski B., Dudek H., Szczesny W. (2004) "Ekonometria. Wybrane zagadnienia" , Wydawnictwo Naukowe PWN, Warszawa.
#'
#' @export
#'
#' @seealso \code{\link{selfStart}}, \code{\link{getInitial}}

SSworking <- selfStart(
  function (input, a, b)
  {
    .expr1 <- 1/input
    .expr4 <- exp(a + b * .expr1)
    .value <- .expr4
    .actualArgs <- as.list(match.call()[c("a", "b", "c")])
    if(all(unlist(lapply(.actualArgs, is.name))))
    {
      .grad <- array(0, c(length(.value), 2L), list(NULL, c("a",
                                                            "b")))
      .grad[, "a"] <- .expr4
      .grad[, "b"] <- .expr4 * .expr1
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
    pars <- as.vector(coef(lm( log(y) ~ I(1/x), data = xy)))
    value <- c(pars[1], pars[2])
    names(value)= mCall[c("a", "b")]
    value
  },parameters=c("a","b"))

#' @title Self-Starting nls logistic model
#'
#' @description
#' Computes the logistic function
#' @param input a numeric vector of values at which to evaluate the model.
#' @param a estimated model parameter
#' @param b estimated model parameter
#' @param c estimated model parameter
#' @usage SSlogis1(input,a,b,c)
#' @return a numeric vector of the same length as input. It is the value of the expression a/(1+exp(b+c*x)).
#' @note This is primarily intended for use in formulae given to the \code{\link{nls}} function or similar functions for example \code{\link{mle2}}.
#' @keywords function
#' @author
#' Krzysztof Trajkowski
#' @examples
#' set.seed(46)
#' x <- 1:15
#' y <- 1/(1+exp(7-1.5*x))+rnorm(15,0,0.1)
#' model <- nls(y~SSlogis1(x,a,b,c))
#' @rdname SSlogis1
#'
#' @references Goryl A., Jedrzejczyk Z., Kukula K. (Eds.), Osiewalski J., Walkosz A. (2004) "Wprowadzenie do ekonometrii w przykladach i zadaniach", Wydawnictwo Naukowe PWN, Warszawa.
#' @references Borkowski B., Dudek H., Szczesny W. (2004) "Ekonometria. Wybrane zagadnienia" , Wydawnictwo Naukowe PWN, Warszawa.
#'
#' @export
#'
#' @seealso \code{\link{SSlogis}}, \code{\link{selfStart}}, \code{\link{getInitial}}

SSlogis1 <- selfStart(
  function (input, a, b, c)
  {
    .expr3 <- exp(b + c * input)
    .expr4 <- 1 + .expr3
    .expr8 <- .expr4^2
    .value <- a/.expr4
    .actualArgs <- as.list(match.call()[c("a", "b", "c")])
    if(all(unlist(lapply(.actualArgs, is.name))))
    {
      .grad <- array(0, c(length(.value), 3L), list(NULL, c("a",
                                                            "b", "c")))
      .grad[, "a"] <- 1/.expr4
      .grad[, "b"] <- -(a * .expr3/.expr8)
      .grad[, "c"] <- -(a * (.expr3 * input)/.expr8)
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
    pars <- as.vector(coef(nls( y ~ SSlogis(x,Asym,xmid,scal), data = xy)))
    value <- c(pars[1], pars[2]/pars[3], -1/pars[3])
    names(value)= mCall[c("a", "b", "c")]
    value
  },parameters=c("a","b", "c"))

#' @title Self-Starting nls logistic model
#'
#' @description
#' Computes the logistic function
#' @param input a numeric vector of values at which to evaluate the model.
#' @param alpha estimated model parameter
#' @param beta estimated model parameter
#' @param gamma estimated model parameter
#' @usage SSlogis2(input,alpha,beta,gamma)
#' @return a numeric vector of the same length as input. It is the value of the expression alpha/(1+beta*exp(-gamma*x)).
#' @note This is primarily intended for use in formulae given to the \code{\link{nls}} function or similar functions for example \code{\link{mle2}}.
#' @keywords function
#' @author
#' Krzysztof Trajkowski
#' @examples
#' set.seed(46)
#' x <- 1:15
#' y <- 1/(1+exp(7-1.5*x))+rnorm(15,0,0.1)
#' model <- nls(y~SSlogis2(x,alpha,beta,gamma))
#' @rdname SSlogis2
#'
#' @references Goryl A., Jedrzejczyk Z., Kukula K. (Eds.), Osiewalski J., Walkosz A. (2004) "Wprowadzenie do ekonometrii w przykladach i zadaniach", Wydawnictwo Naukowe PWN, Warszawa.
#' @references Borkowski B., Dudek H., Szczesny W. (2004) "Ekonometria. Wybrane zagadnienia" , Wydawnictwo Naukowe PWN, Warszawa.
#'
#' @export
#'
#' @seealso \code{\link{SSlogis}}, \code{\link{selfStart}}, \code{\link{getInitial}}

SSlogis2 <- selfStart(
  function (input, alpha, beta, gamma)
  {
    .expr3 <- exp(-gamma * input)
    .expr5 <- 1 + beta * .expr3
    .expr9 <- .expr5^2
    .value <- alpha/.expr5
    .actualArgs <- as.list(match.call()[c("alpha", "beta", "gamma")])
    if(all(unlist(lapply(.actualArgs, is.name))))
    {
      .grad <- array(0, c(length(.value), 3L), list(NULL, c("alpha",
                                                            "beta", "gamma")))
      .grad[, "alpha"] <- 1/.expr5
      .grad[, "beta"] <- -(alpha * .expr3/.expr9)
      .grad[, "gamma"] <- alpha * (beta * (.expr3 * input))/.expr9
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
    pars <- as.vector(coef(nls( y ~ SSlogis(x,Asym,xmid,scal), data = xy)))
    value <- c(pars[1], exp(pars[2]/pars[3]), 1/pars[3])
    names(value)= mCall[c("alpha", "beta", "gamma")]
    value
  },parameters=c("alpha","beta", "gamma"))

#' @title Self-Starting nls power model
#'
#' @description
#' Computes the logistic function
#' @param input a numeric vector of values at which to evaluate the model.
#' @param a estimated model parameter
#' @param b estimated model parameter
#' @usage SSpower(input,a,b)
#' @return a numeric vector of the same length as input. It is the value of the expression a*x^b.
#' @note This is primarily intended for use in formulae given to the \code{\link{nls}} function or similar functions for example \code{\link{mle2}}.
#' @keywords function
#' @author
#' Krzysztof Trajkowski
#' @examples
#' set.seed(36)
#' x <- 1:100
#' y <- 1.25*x^0.45+rnorm(100,0,0.5)
#' model <- nls(y~SSpower(x,a,b))
#' @rdname SSpower
#'
#' @references Goryl A., Jedrzejczyk Z., Kukula K. (Eds.), Osiewalski J., Walkosz A. (2004) "Wprowadzenie do ekonometrii w przykladach i zadaniach", Wydawnictwo Naukowe PWN, Warszawa.
#' @references Borkowski B., Dudek H., Szczesny W. (2004) "Ekonometria. Wybrane zagadnienia" , Wydawnictwo Naukowe PWN, Warszawa.
#'
#' @export
#'
#' @seealso \code{\link[vegan]{SSarrhenius}}, \code{\link{selfStart}}, \code{\link{getInitial}}

SSpower <- selfStart(
  function (input, a, b)
  {
    .expr1 <- input^b
    .value <- a * .expr1
    .actualArgs <- as.list(match.call()[c("a", "b")])
    if(all(unlist(lapply(.actualArgs, is.name))))
    {
      .grad <- array(0, c(length(.value), 2L), list(NULL, c("a",
                                                            "b")))
      .grad[, "a"] <- .expr1
      .grad[, "b"] <- a * (.expr1 * log(input))
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
    pars <- as.vector(coef(lm( log(y) ~ log(x), data = xy)))
    value <- c( exp(pars[1]), pars[2] )
    names(value)= mCall[c("a", "b")]
    value
  },parameters=c("a","b"))

#' @title Self-Starting nls exponential model
#'
#' @description
#' Computes the exponential function
#' @param input a numeric vector of values at which to evaluate the model.
#' @param a estimated model parameter
#' @param b estimated model parameter
#' @usage SSexpo(input,a,b)
#' @return a numeric vector of the same length as input. It is the value of the expression a*b^x.
#' @note This is primarily intended for use in formulae given to the \code{\link{nls}} function or similar functions for example \code{\link{mle2}}.
#' @keywords function
#' @author
#' Krzysztof Trajkowski
#' @examples
#' x <- 1:20
#' y <- 0.25*1.54^x+runif(20,0,100)
#' model <- nls(y ~ SSexpo(x,a,b))
#' plot(x,y)
#' p <- coef(model)
#' curve(SSexpo(x,p[1],p[2]),add=TRUE)
#' @rdname SSexpo
#'
#' @references Goryl A., Jedrzejczyk Z., Kukula K. (Eds.), Osiewalski J., Walkosz A. (2004) "Wprowadzenie do ekonometrii w przykladach i zadaniach", Wydawnictwo Naukowe PWN, Warszawa.
#' @references Borkowski B., Dudek H., Szczesny W. (2004) "Ekonometria. Wybrane zagadnienia" , Wydawnictwo Naukowe PWN, Warszawa.
#'
#' @export
#'
#' @seealso \code{\link{selfStart}}, \code{\link{getInitial}}

SSexpo <- selfStart(
  function (input, a, b)
  {
    .expr1 <- b^input
    .value <- a * .expr1
    .actualArgs <- as.list(match.call()[c("a", "b")])
    if(all(unlist(lapply(.actualArgs, is.name))))
    {
      .grad <- array(0, c(length(.value), 2L), list(NULL, c("a",
                                                            "b")))
      .grad[, "a"] <- .expr1
      .grad[, "b"] <- a * (b^(input - 1) * input)
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
    pars <- as.vector(coef(lm( log(y) ~ x, data = xy)))
    value <- c( exp(pars[1]), exp(pars[2]) )
    names(value)= mCall[c("a", "b")]
    value
  },parameters=c("a","b"))

#' @title Self-Starting nls exponential inverse model
#'
#' @description
#' Computes the exponential inverse function
#' @param input a numeric vector of values at which to evaluate the model.
#' @param a estimated model parameter
#' @param b estimated model parameter
#' @usage SSexpoInv(input,a,b)
#' @return a numeric vector of the same length as input. It is the value of the expression a*exp(b*(1/x)).
#' @note This is primarily intended for use in formulae given to the \code{\link{nls}} function or similar functions for example \code{\link{mle2}}.
#' @keywords function
#' @author
#' Krzysztof Trajkowski
#' @examples
#' x <- 1:20
#' y <- 0.25*exp(2.87*(1/x))
#' model <- nls(y ~ SSexpoInv(x,a,b),algorithm = "port")
#' @rdname SSexpoInv
#'
#' @references Goryl A., Jedrzejczyk Z., Kukula K. (Eds.), Osiewalski J., Walkosz A. (2004) "Wprowadzenie do ekonometrii w przykladach i zadaniach", Wydawnictwo Naukowe PWN, Warszawa.
#' @references Borkowski B., Dudek H., Szczesny W. (2004) "Ekonometria. Wybrane zagadnienia" , Wydawnictwo Naukowe PWN, Warszawa.
#'
#' @export
#'
#' @seealso \code{\link{selfStart}}, \code{\link{getInitial}}

SSexpoInv <- selfStart(
  function (input, a, b)
  {
    .expr1 <- 1/input
    .expr3 <- exp(b * .expr1)
    .value <- a * .expr3
    .actualArgs <- as.list(match.call()[c("a", "b")])
    if(all(unlist(lapply(.actualArgs, is.name))))
    {
      .grad <- array(0, c(length(.value), 2L), list(NULL, c("a",
                                                            "b")))
      .grad[, "a"] <- .expr3
      .grad[, "b"] <- a * (.expr3 * .expr1)
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
    pars <- as.vector(coef(lm( log(y) ~ I(1/x), data = xy)))
    value <- c( exp(pars[1]), pars[2] )
    names(value)= mCall[c("a", "b")]
    value
  },parameters=c("a","b"))

#' @title Self-Starting nls power-exponential model
#'
#' @description
#' Computes the power-exponential function
#' @param input a numeric vector of values at which to evaluate the model.
#' @param a estimated model parameter
#' @param b estimated model parameter
#' @param c estimated model parameter
#' @usage SSpowExpo(input,a,b,c)
#' @return a numeric vector of the same length as input. It is the value of the expression a*x^b*exp(c*x).
#' @note This is primarily intended for use in formulae given to the \code{\link{nls}} function or similar functions for example \code{\link{mle2}}.
#' @keywords function
#' @author
#' Krzysztof Trajkowski
#' @examples
#' x <- 1:20
#' y <- 7.25*x^2.85*exp(0.5*x)
#' model <- nls(y ~ SSpowExpo(x,a,b,c),algorithm = "port")
#' @rdname SSpowExpo
#'
#' @references Goryl A., Jedrzejczyk Z., Kukula K. (Eds.), Osiewalski J., Walkosz A. (2004) "Wprowadzenie do ekonometrii w przykladach i zadaniach", Wydawnictwo Naukowe PWN, Warszawa.
#' @references Borkowski B., Dudek H., Szczesny W. (2004) "Ekonometria. Wybrane zagadnienia" , Wydawnictwo Naukowe PWN, Warszawa.
#'
#' @export
#'
#' @seealso \code{\link{selfStart}}, \code{\link{getInitial}}

SSpowExpo <- selfStart(
  function (input, a, b, c)
  {
    .expr1 <- input^b
    .expr2 <- a * .expr1
    .expr4 <- exp(c * input)
    .value <- .expr2 * .expr4
    .actualArgs <- as.list(match.call()[c("a", "b", "c")])
    if(all(unlist(lapply(.actualArgs, is.name))))
    {
      .grad <- array(0, c(length(.value), 3L), list(NULL, c("a",
                                                            "b", "c")))
      .grad[, "a"] <- .expr1 * .expr4
      .grad[, "b"] <- a * (.expr1 * log(input)) * .expr4
      .grad[, "c"] <- .expr2 * (.expr4 * input)
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
    pars <- as.vector(coef(lm( log(y) ~ log(x)+x, data = xy)))
    value <- c( exp(pars[1]), pars[2], pars[3] )
    names(value)= mCall[c("a", "b", "c")]
    value
  },parameters=c("a","b", "c"))

#' @title Self-Starting nls power-exponential-inverse model
#'
#' @description
#' Computes the power-exponential-inverse function
#' @param input a numeric vector of values at which to evaluate the model.
#' @param a estimated model parameter
#' @param b estimated model parameter
#' @param c estimated model parameter
#' @usage SSpowExpoInv(input,a,b,c)
#' @return a numeric vector of the same length as input. It is the value of the expression a*x^b*exp(c*(1/x)).
#' @note This is primarily intended for use in formulae given to the \code{\link{nls}} function or similar functions for example \code{\link{mle2}}.
#' @keywords function
#' @author
#' Krzysztof Trajkowski
#' @examples
#' x <- 1:20
#' y <- 7.25*x^2.85*exp(0.5*(1/x))
#' model <- nls(y ~ SSpowExpoInv(x,a,b,c),algorithm = "port")
#' @rdname SSpowExpoInv
#'
#' @references Goryl A., Jedrzejczyk Z., Kukula K. (Eds.), Osiewalski J., Walkosz A. (2004) "Wprowadzenie do ekonometrii w przykladach i zadaniach", Wydawnictwo Naukowe PWN, Warszawa.
#' @references Borkowski B., Dudek H., Szczesny W. (2004) "Ekonometria. Wybrane zagadnienia" , Wydawnictwo Naukowe PWN, Warszawa.
#'
#' @export
#'
#' @seealso \code{\link{selfStart}}, \code{\link{getInitial}}

SSpowExpoInv <- selfStart(
  function (input, a, b, c)
  {
    .expr1 <- input^b
    .expr2 <- a * .expr1
    .expr3 <- 1/input
    .expr5 <- exp(c * .expr3)
    .value <- .expr2 * .expr5
    .actualArgs <- as.list(match.call()[c("a", "b", "c")])
    if(all(unlist(lapply(.actualArgs, is.name))))
    {
      .grad <- array(0, c(length(.value), 3L), list(NULL, c("a",
                                                            "b", "c")))
      .grad[, "a"] <- .expr1 * .expr5
      .grad[, "b"] <- a * (.expr1 * log(input)) * .expr5
      .grad[, "c"] <- .expr2 * (.expr5 * .expr3)
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
    pars <- as.vector(coef(lm( log(y) ~ log(x)+I(1/x), data = xy)))
    value <- c( exp(pars[1]), pars[2], pars[3] )
    names(value)= mCall[c("a", "b", "c")]
    value
  },parameters=c("a","b", "c"))

#' @title Self-Starting nls parabola logarithmic model
#'
#' @description
#' Computes the parabola logarithmic function
#' @param input a numeric vector of values at which to evaluate the model.
#' @param a estimated model parameter
#' @param b estimated model parameter
#' @param c estimated model parameter
#' @usage SSparLog(input,a,b,c)
#' @return a numeric vector of the same length as input. It is the value of the expression a*x^(b+c*log(x)).
#' @note This is primarily intended for use in formulae given to the \code{\link{nls}} function or similar functions for example \code{\link{mle2}}.
#' @keywords function
#' @author
#' Krzysztof Trajkowski
#' @examples
#' x <- 1:20
#' y <- 5.75*x^(0.75+0.95*log(x))
#' model <- nls(y ~ SSparLog(x,a,b,c),algorithm = "port")
#' @rdname SSparLog
#'
#' @references Goryl A., Jedrzejczyk Z., Kukula K. (Eds.), Osiewalski J., Walkosz A. (2004) "Wprowadzenie do ekonometrii w przykladach i zadaniach", Wydawnictwo Naukowe PWN, Warszawa.
#' @references Borkowski B., Dudek H., Szczesny W. (2004) "Ekonometria. Wybrane zagadnienia" , Wydawnictwo Naukowe PWN, Warszawa.
#'
#' @export
#'
#' @seealso \code{\link{selfStart}}, \code{\link{getInitial}}

SSparLog <- selfStart(
  function (input, a, b, c)
  {
    .expr1 <- log(input)
    .expr4 <- input^(b + c * .expr1)
    .value <- a * .expr4
    .actualArgs <- as.list(match.call()[c("a", "b", "c")])
    if(all(unlist(lapply(.actualArgs, is.name))))
    {
      .grad <- array(0, c(length(.value), 3L), list(NULL, c("a",
                                                            "b", "c")))
      .grad[, "a"] <- .expr4
      .grad[, "b"] <- a * (.expr4 * .expr1)
      .grad[, "c"] <- a * (.expr4 * (.expr1 * .expr1))
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
    pars <- as.vector(coef(lm( log(y) ~ log(x)+I(log(x)^2), data = xy)))
    value <- c( exp(pars[1]), pars[2], pars[3] )
    names(value)= mCall[c("a", "b", "c")]
    value
  },parameters=c("a","b", "c"))

#' @title Self-Starting nls hyperbolic model
#'
#' @description
#' Computes the hyperbolic function
#' @param input a numeric vector of values at which to evaluate the model.
#' @param a estimated model parameter
#' @param b estimated model parameter
#' @usage SShyper(input,a,b)
#' @return a numeric vector of the same length as input. It is the value of the expression (a*x^2)/(x+b).
#' @note This is primarily intended for use in formulae given to the \code{\link{nls}} function or similar functions for example \code{\link{mle2}}.
#' @keywords function
#' @author
#' Krzysztof Trajkowski
#' @examples
#' x <- 1:20
#' y <- (255*x^2)/(x+587)
#' model <- nls(y ~ SShyper(x,a,b),algorithm = "port")
#' @rdname SShyper
#'
#' @references Goryl A., Jedrzejczyk Z., Kukula K. (Eds.), Osiewalski J., Walkosz A. (2004) "Wprowadzenie do ekonometrii w przykladach i zadaniach", Wydawnictwo Naukowe PWN, Warszawa.
#' @references Borkowski B., Dudek H., Szczesny W. (2004) "Ekonometria. Wybrane zagadnienia" , Wydawnictwo Naukowe PWN, Warszawa.
#'
#' @export
#'
#' @seealso \code{\link{selfStart}}, \code{\link{getInitial}}

SShyper <- selfStart(
  function (input, a, b)
  {
    .expr1 <- input^2
    .expr2 <- a * .expr1
    .expr3 <- input + b
    .value <- .expr2/.expr3
    .actualArgs <- as.list(match.call()[c("a", "b")])
    if(all(unlist(lapply(.actualArgs, is.name))))
    {
      .grad <- array(0, c(length(.value), 2L), list(NULL, c("a",
                                                            "b")))
      .grad[, "a"] <- .expr1/.expr3
      .grad[, "b"] <- -(.expr2/.expr3^2)
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
    pars <- as.vector(coef(lm( I(1/y) ~ I(1/x)+I(1/x^2), data = xy)))
    value <- c( 1/pars[2], pars[3]/pars[2] )
    names(value)= mCall[c("a", "b")]
    value
  },parameters=c("a","b"))

#' @title Self-Starting nls general exponential model
#'
#' @description
#' Computes the general exponential function
#' @param input a numeric vector of values at which to evaluate the model.
#' @param a estimated model parameter
#' @param b estimated model parameter
#' @usage SSexpoGen(input,a,b)
#' @return a numeric vector of the same length as input. It is the value of the expression a*exp(b*x).
#' @note This is primarily intended for use in formulae given to the \code{\link{nls}} function or similar functions for example \code{\link{mle2}}.
#' @keywords function
#' @author
#' Krzysztof Trajkowski
#' @examples
#' x <- 1:20
#' y <- 2.5*exp(0.55*x)
#' model <- nls(y ~ SSexpoGen(x,a,b),algorithm = "port")
#' @rdname SSexpoGen
#'
#' @references Goryl A., Jedrzejczyk Z., Kukula K. (Eds.), Osiewalski J., Walkosz A. (2004) "Wprowadzenie do ekonometrii w przykladach i zadaniach", Wydawnictwo Naukowe PWN, Warszawa.
#' @references Borkowski B., Dudek H., Szczesny W. (2004) "Ekonometria. Wybrane zagadnienia" , Wydawnictwo Naukowe PWN, Warszawa.
#'
#' @export
#'
#' @seealso \code{\link{selfStart}}, \code{\link{getInitial}}

SSexpoGen <- selfStart(
  function (input, a, b)
  {
    .expr2 <- exp(b * input)
    .value <- a * .expr2
    .actualArgs <- as.list(match.call()[c("a", "b")])
    if(all(unlist(lapply(.actualArgs, is.name))))
    {
      .grad <- array(0, c(length(.value), 2L), list(NULL, c("a",
                                                            "b")))
      .grad[, "a"] <- .expr2
      .grad[, "b"] <- a * (.expr2 * input)
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
    pars <- as.vector(coef(lm( log(y) ~ x, data = xy)))
    value <- c( exp(pars[1]), pars[2] )
    names(value)= mCall[c("a", "b")]
    value
  },parameters=c("a","b"))

#' @title Self-Starting nls amplitude of gaussian model
#'
#' @description
#' Amplitude version of gaussian peak function.
#' @param input a numeric vector of values at which to evaluate the model.
#' @param y0 offset
#' @param A amplitude
#' @param xc center
#' @param w width
#' @usage SSgaussAmp(input,y0,A,xc,w)
#' @return a numeric vector of the same length as input. It is the value of the expression y0+A*exp(-(x-xc)^2/2*w^2).
#' @note This is primarily intended for use in formulae given to the \code{\link{nls}} function or similar functions for example \code{\link{mle2}}.
#' @keywords function
#' @author
#' Krzysztof Trajkowski
#' @examples
#' data(smoking)
#' model <- nlsLM(Germany ~ SSgaussAmp(year,y0,A,xc,w),data=smoking)
#' @rdname SSgaussAmp
#'
#' @references Origin: \url{http://www.originlab.com/doc/Origin-Help/GaussAmp-FitFunc}.
#'
#' @export
#'
#' @seealso \code{\link{selfStart}}, \code{\link{getInitial}}

SSgaussAmp <- selfStart(
  function (input, y0, A, xc, w)
  {
    .expr1 <- input - xc
    .expr4 <- -.expr1^2/2
    .expr5 <- w^2
    .expr7 <- exp(.expr4 * .expr5)
    .value <- y0 + A * .expr7
    .actualArgs <- as.list(match.call()[c("y0", "A", "xc", "w")])
    if(all(unlist(lapply(.actualArgs, is.name))))
    {
      .grad <- array(0, c(length(.value), 4L), list(NULL, c("y0",
                                                            "A", "xc", "w")))
      .grad[, "y0"] <- 1
      .grad[, "A"] <- .expr7
      .grad[, "xc"] <- A * (.expr7 * (2 * .expr1/2 * .expr5))
      .grad[, "w"] <- A * (.expr7 * (.expr4 * (2 * w)))
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
    pars <- as.vector(coef(minpack.lm::nlsLM( y ~ y0+A*exp(-(x-xc)^2/2*w^2),start=list(y0=min(xy$y),A=max(xy$y)-min(xy$y),xc=xy$x[which.max(xy$y)],w=1),lower=c(-Inf,-Inf,-Inf,0),upper=c(Inf,Inf,Inf,Inf), data = xy)))
    value <- c( pars[1], pars[2], pars[3], pars[4] )
    names(value)= mCall[c("y0", "A", "xc", "w")]
    value
  },parameters=c("y0","A","xc","w"))

#' @title Self-Starting nls area of gaussian model
#'
#' @description
#' Area version of gaussian function.
#' @param input a numeric vector of values at which to evaluate the model.
#' @param y0 offset
#' @param A area
#' @param xc center
#' @param w width
#' @usage SSgaussAre(input,y0,A,xc,w)
#' @return a numeric vector of the same length as input. It is the value of the expression y0+(A/(w*sqrt(pi/2)))*exp(-2*(x-xc)^2/(w^2)).
#' @note This is primarily intended for use in formulae given to the \code{\link{nls}} function or similar functions for example \code{\link{mle2}}.
#' @keywords function
#' @author
#' Krzysztof Trajkowski
#' @examples
#' data(smoking)
#' model <- nlsLM(Germany ~ SSgaussAre(year,y0,A,xc,w),data=smoking)
#' @rdname SSgaussAre
#'
#' @references Origin: \url{http://www.originlab.com/doc/Origin-Help/Gauss-FitFunc}.
#'
#' @export
#'
#' @seealso \code{\link{selfStart}}, \code{\link{getInitial}}

SSgaussAre <- selfStart(
  function (input, y0, A, xc, w)
  {
    .expr2 <- sqrt(pi/2)
    .expr3 <- w * .expr2
    .expr4 <- A/.expr3
    .expr6 <- input - xc
    .expr8 <- -2 * .expr6^2
    .expr9 <- w^2
    .expr11 <- exp(.expr8/.expr9)
    .value <- y0 + .expr4 * .expr11
    .actualArgs <- as.list(match.call()[c("y0", "A", "xc", "w")])
    if(all(unlist(lapply(.actualArgs, is.name))))
    {
      .grad <- array(0, c(length(.value), 4L), list(NULL, c("y0",
                                                            "A", "xc", "w")))
      .grad[, "y0"] <- 1
      .grad[, "A"] <- 1/.expr3 * .expr11
      .grad[, "xc"] <- .expr4 * (.expr11 * (2 * (2 * .expr6)/.expr9))
      .grad[, "w"] <- -(.expr4 * (.expr11 * (.expr8 * (2 * w)/.expr9^2)) +
                          A * .expr2/.expr3^2 * .expr11)
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
    pars <- as.vector(coef(minpack.lm::nlsLM( y ~ y0+(A/(w*sqrt(pi/2)))*exp(-2*(x-xc)^2/(w^2)),
                                              start=list(y0=min(xy$y),A=(sqrt(pi)*((max(xy$x)-min(xy$x))/3)*max(xy$y)-sqrt(pi)*((max(xy$x)-min(xy$x))/3)*min(xy$y))/sqrt(2),xc=xy$x[which.max(xy$y)],w=(max(xy$x)-min(xy$x))/3),lower=c(-Inf,0,-Inf,0),upper=c(Inf,Inf,Inf,Inf), data = xy)))
    value <- c( pars[1], pars[2], pars[3], pars[4] )
    names(value)= mCall[c("y0", "A", "xc", "w")]
    value
  },parameters=c("y0","A","xc","w"))

#' @title Self-Starting nls lorentz model
#'
#' @description
#' lorentz peak function.
#' @param input a numeric vector of values at which to evaluate the model.
#' @param y0 offset
#' @param A area
#' @param xc center
#' @param w width
#' @usage SSlorentz(input,y0,A,xc,w)
#' @return a numeric vector of the same length as input. It is the value of the expression y0+A*(2*w)/(pi*(4*(x-xc)^2+w^2)).
#' @note This is primarily intended for use in formulae given to the \code{\link{nls}} function or similar functions for example \code{\link{mle2}}.
#' @keywords function
#' @author
#' Krzysztof Trajkowski
#' @examples
#' data(smoking)
#' smoking$t <- 1:nrow(smoking)
#' model <- nlsLM(Germany ~ SSlorentz(t,y0,A,xc,w),data=smoking)
#' @rdname SSlorentz
#'
#' @references Origin: \url{http://www.originlab.com/doc/Origin-Help/Lorentz-FitFunc}.
#'
#' @export
#'
#' @seealso \code{\link{selfStart}}, \code{\link{getInitial}}

SSlorentz <- selfStart(
  function (input, y0, A, xc, w)
  {
    .expr1 <- 2 * w
    .expr2 <- A * .expr1
    .expr3 <- input - xc
    .expr8 <- pi * (4 * .expr3^2 + w^2)
    .expr16 <- .expr8^2
    .value <- y0 + .expr2/.expr8
    .actualArgs <- as.list(match.call()[c("y0", "A", "xc", "w")])
    if(all(unlist(lapply(.actualArgs, is.name))))
    {
      .grad <- array(0, c(length(.value), 4L), list(NULL, c("y0",
                                                            "A", "xc", "w")))
      .grad[, "y0"] <- 1
      .grad[, "A"] <- .expr1/.expr8
      .grad[, "xc"] <- .expr2 * (pi * (4 * (2 * .expr3)))/.expr16
      .grad[, "w"] <- A * 2/.expr8 - .expr2 * (pi * .expr1)/.expr16
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
    pars <- as.vector(coef(minpack.lm::nlsLM( y ~ y0+A*(2*w)/(pi*(4*(x-xc)^2+w^2)),
                                              start=list(y0=min(xy$y),A=(pi*((max(xy$x)-min(xy$x))/3)*max(xy$y)-pi*((max(xy$x)-min(xy$x))/3)*min(xy$y))/2,xc=xy$x[which.max(xy$y)],w=(max(xy$x)-min(xy$x))/3), lower=c(-Inf,0,-Inf,0),upper=c(Inf,Inf,Inf,Inf), data = xy)))
    value <- c( pars[1], pars[2], pars[3], pars[4] )
    names(value)= mCall[c("y0", "A", "xc", "w")]
    value
  },parameters=c("y0","A","xc","w"))



#' @title Smoking data
#'
#' @description Sales of cigarettes per adult per day in early industrialized countries, 1920-2010.
#' @name smoking
#' @docType data
#' @usage smoking
#' @format This data frame contains the observations of 91.
#' @source Max Roser (2016) – 'Smoking'. Published online at OurWorldInData.org. Retrieved from: \url{https://ourworldindata.org/smoking} [Online Resource].
#' @source The data for the single countries from the 'International Smoking Statistics' \url{http://www.pnlee.co.uk/ISS.htm}.
#' @keywords datasets
#' @examples
#' head(smoking)
NULL

#' @title Population in USA 1790-2010
#'
#' @description United States population data 1790-2010
#' @name population
#' @docType data
#' @usage population
#' @format This data frame contains the observations of 23.
#' @source US Bureau of the Census. Web: \url{http://www.census.gov}.
#' @keywords datasets
#' @examples
#' head(population)
NULL
