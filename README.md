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
devtools::install_github("krzysiektr/FitCurve")
```

## Example

Tornquist III model:
$$f(x)=\frac{ax}{x+b}$$

```{r}
library("FitCurve")

x <- 1:15
y <- 120*x*(x-254)/(x+654)

model <- nls(y~SStorn3(x,a,b,c),algorithm = "port")
summary(model)
```
