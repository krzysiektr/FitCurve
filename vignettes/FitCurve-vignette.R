## ----setup, include=FALSE------------------------------------------------
library(knitr)
opts_chunk$set(fig.align="center",
               fig.width=6, fig.height=4,
               dev.args=list(pointsize=8),
               par=TRUE)

knit_hooks$set(par = function(before, options, envir)
  { if(before && options$fig.show != "none") 
       par(mar=c(4,4,1,1)+0.1, mgp=c(3,0.6,0),las=1)
})

## ------------------------------------------------------------------------
library("FitCurve")

## ------------------------------------------------------------------------
expo <- nlsLM(USA~SSexpo(year,a,b),
              data=subset(smoking,year<="1949"))
summary(expo)

## ------------------------------------------------------------------------
with(plot(y=USA,x=year,type="l",lwd=2,
          xlim=c(1900,2030),ylim=c(0,12)),
     data=subset(smoking,year<="1949"))
title("Sales of cigarettes per adult per day in USA 1927-1949")
p <- coef(expo)
curve(SSexpo(x,p[1],p[2]),add=TRUE,col="orange3",lwd=2)

## ------------------------------------------------------------------------
logi <- nls(USA~SSlogis(year,Asym,xmid,scal),
            algorithm = "port",
            data=subset(smoking,year<="1966"))
summary(logi)

## ------------------------------------------------------------------------
with(plot(y=USA,x=year,type="l",lwd=2,
          xlim=c(1900,2030),ylim=c(0,12)),
     data=subset(smoking,year<="1966"))
title("Sales of cigarettes per adult per day in USA 1927-1966")
p <- coef(logi)
curve(SSlogis(x,p[1],p[2],p[3]),add=TRUE,col="orange3",lwd=2)

## ------------------------------------------------------------------------
amp <- nls(USA ~ SSgaussAmp(year,y0,A,xc,w),
           data=smoking)
summary(amp)

## ----fig.height=4,fig.width=6--------------------------------------------
with(plot(year,USA,type="l",lwd=2,
          xlim=c(1900,2030),ylim=c(0,12)),
     data=smoking)
title("Sales of cigarettes per adult per day in USA 1927-2010")
p <- coef(amp)
curve(SSgaussAmp(x,p[1],p[2],p[3],p[4]),add=TRUE,col="orange3",lwd=2)

## ------------------------------------------------------------------------
p <- getInitial(population~SSfpl(year,A,B,xmid,scal),data=population)

## ------------------------------------------------------------------------
fpl_nor <- mle2(population~dnorm(mean= SSfpl(year,A,B,xmid,scal),sd= sd),
                start=list(A= p[1], B= p[2], xmid= p[3], scal= p[4], sd= 1),
                method= "L-BFGS-B", data= population,
                lower=c(A=-Inf,B=-Inf,xmid=-Inf,scal=-Inf,sd=0))
summary(fpl_nor)

## ----fig.height=4,fig.width=6--------------------------------------------
with(plot(year,population,type="l",lwd=3,
          xlim=c(1800,2100),ylim=c(0,45e+07)),
     data=population)
title("Population in USA 1790-2010")
p <- coef(fpl_nor)
curve(SSfpl(x,p[1],p[2],p[3],p[4]),add=TRUE,col="orange3",lwd=1.5)

## ------------------------------------------------------------------------
p <- getInitial(population~SSlogis(year,Asym,xmid,scal),data=population)

## ------------------------------------------------------------------------
logi_nor <- mle2(population~dnorm(mean= SSlogis(year,Asym,xmid,scal),sd= sd),
                 start=list(Asym= p[1], xmid= p[2], scal= p[3], sd= 1),
                 method= "L-BFGS-B", data= population,
                 lower=c(Asym=-Inf,xmid=-Inf,scal=-Inf,sd=0))
summary(logi_nor)

## ----fig.height=4,fig.width=6--------------------------------------------
with(plot(year,population,type="l",lwd=3,
          xlim=c(1800,2100),ylim=c(0,45e+07)),
     data=population)
title("Population in USA 1790-2010")
p <- coef(logi_nor)
curve(SSlogis(x,p[1],p[2],p[3]),add=TRUE,col="orange3",lwd=1.5)

