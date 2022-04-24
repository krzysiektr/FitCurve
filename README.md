# FitCurve

## Install

From GitHub:

* dependencies

```
install.packages("bbmle")
install.packages("minpack.lm")
```

* FitCurve package:

```
devtools::install_github("krzysiektr/FitCurve", build_vignettes = TRUE)
```

## Example

Tornquist II model:
$$f(x)=\frac{a(x-c)}{x+b}$$

```{r}
library("FitCurve")

set.seed(2305)
x <- 1:15
y <- log(x)+rnorm(1)

model <- nls(y~SStorn2(x,a,b,c))
summary(model)

plot(x,y,pch=20)
init <- getInitial( y~SStorn2( x, a,b,c), data.frame(x,y) )
p <- coef(model)
curve(SStorn2(x,p[1],p[2],p[3]),add=TRUE,col="orange")
```
