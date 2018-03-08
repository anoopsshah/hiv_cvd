## this R code assumes that the working directory is set to same directory as that where the R code is saved####


# load packages####
library(tableone)
library(caret)
library(reshape2)
library(tidyr)
library(dplyr)
library(binom)
library(survival)
library(gridExtra)
library(grid)
library(gdata)
library(survMisc)
library(scales)
library(ggplot2)
library(tidyr)
library(metafor)
library(xlsx)
library(xlsxjars)
library(rJava)
library(readxl)
library(stringr)
library(reshape)
library(data.table)
library(plyr)
library(RColorBrewer)
library(rworldmap)
library(classInt)
library(gdata)

# load beta pharm function ####

beta.parms.from.quantiles <- function(q, p=c(0.025,0.975),
                                      precision=0.001, derivative.epsilon=1e-3, start.with.normal.approx=T, start=c(1, 1), plot=F)
{
  # Version 1.2.2 (December 2012)
  #
  # Function developed by 
  # Lawrence Joseph and Patrick Belisle
  # Division of Clinical Epidemiology
  # Montreal General Hospital
  # Montreal, Qc, Can
  #
  # patrick.belisle@clinepi.mcgill.ca
  # http://www.medicine.mcgill.ca/epidemiology/Joseph/PBelisle/BetaParmsFromQuantiles.html
  #
  # Please refer to our webpage for details on each argument.
  
  f <- function(x, theta){dbeta(x, shape1=theta[1], shape2=theta[2])}
  F.inv <- function(x, theta){qbeta(x, shape1=theta[1], shape2=theta[2])}
  f.cum <- function(x, theta){pbeta(x, shape1=theta[1], shape2=theta[2])}
  f.mode <- function(theta){a <- theta[1]; b <- theta[2]; mode <- ifelse(a>1, (a-1)/(a+b-2), NA); mode}
  theta.from.moments <- function(m, v){a <- m*m*(1-m)/v-m; b <- a*(1/m-1); c(a, b)}
  plot.xlim <- c(0, 1)
  
  dens.label <- 'dbeta'
  parms.names <- c('a', 'b')
  
  if (length(p) != 2) stop("Vector of probabilities p must be of length 2.")
  if (length(q) != 2) stop("Vector of quantiles q must be of length 2.")
  p <- sort(p); q <- sort(q)
  
  #_____________________________________________________________________________________________________
  
  print.area.text <- function(p, p.check, q, f, f.cum, F.inv, theta, mode, cex, plot.xlim, M=30, M0=50)
  {
    par.usr <- par('usr')
    par.din <- par('din')
    
    p.string <- as.character(round(c(0,1) + c(1,-1)*p.check, digits=4))
    str.width <- strwidth(p.string, cex=cex)
    str.height <- strheight("0", cex=cex)
    
    J <- matrix(1, nrow=M0, ncol=1)
    
    x.units.1in <- diff(par.usr[c(1,2)])/par.din[1]
    y.units.1in <- diff(par.usr[c(3,4)])/par.din[2]
    aspect.ratio <- y.units.1in/x.units.1in
    
    # --- left area  -----------------------------------------------------------
    
    scatter.xlim <- c(max(plot.xlim[1], par.usr[1]), q[1])
    scatter.ylim <- c(0, par.usr[4])
    x <- seq(from=scatter.xlim[1], to=scatter.xlim[2], length=M)
    y <- seq(from=scatter.ylim[1], to=scatter.ylim[2], length=M)
    x.grid.index <- rep(seq(M), M)
    y.grid.index <- rep(seq(M), rep(M, M))
    
    grid.df <- f(x, theta)
    
    # Estimate mass center
    tmp.p <- seq(from=0, to=p[1], length=M0)
    tmp.x <- F.inv(tmp.p, theta)
    h <- f(tmp.x, theta)
    mass.center <- c(mean(tmp.x), sum(h[-1]*diff(tmp.x))/diff(range(tmp.x)))
    
    # Identify points under the curve
    # (to eliminate them from the list of candidates)
    gridpoint.under.the.curve <- y[y.grid.index] <= grid.df[x.grid.index]
    w <- which(gridpoint.under.the.curve)
    x <- x[x.grid.index]; y <- y[y.grid.index]
    if (length(w)){x <- x[-w]; y <- y[-w]}
    
    # Eliminate points to the right of the mode, if any
    w <- which(x>mode)
    if (length(w)){x <- x[-w]; y <- y[-w]}
    
    # Eliminate points for which the text would fall out of the plot area
    w <- which((par.usr[1]+str.width[1]) <= x & (y + str.height) <= par.usr[4])
    x <- x[w]; y <- y[w]
    
    # For each height, eliminate the closest point to the curve
    # (we want to stay away from the curve to preserve readability)
    w <- which(!duplicated(y, fromLast=T))
    if (length(w)){x <- x[-w]; y <- y[-w]}
    
    # For each point, compute distance from mass center and pick the closest point
    d <- ((x-mass.center[1])^2) + ((y-mass.center[2])/aspect.ratio)^2
    w <- which.min(d)
    x <- x[w]; y <- y[w]
    
    if (length(x))
    {
      text(x, y, labels=p.string[1], adj=c(1,0), col='gray', cex=cex)
    }
    else
    {
      text(plot.xlim[1], mean(par.usr[c(3,4)]), labels=p.string[1], col='gray', cex=cex, srt=90, adj=c(1,0))
    }
    
    # --- right area  ----------------------------------------------------------
    
    scatter.xlim <- c(q[2], plot.xlim[2])
    scatter.ylim <- c(0, par.usr[4])
    x <- seq(from=scatter.xlim[1], to=scatter.xlim[2], length=M)
    y <- seq(from=scatter.ylim[1], to=scatter.ylim[2], length=M)
    x.grid.index <- rep(seq(M), M)
    y.grid.index <- rep(seq(M), rep(M, M))
    grid.df <- f(x, theta)
    
    # Estimate mass center
    tmp.p <- seq(from=p[2], to=f.cum(plot.xlim[2], theta), length=M0)
    tmp.x <- F.inv(tmp.p, theta)
    h <- f(tmp.x, theta)
    mass.center <- c(mean(tmp.x), sum(h[-length(h)]*diff(tmp.x))/diff(range(tmp.x)))
    
    # Identify points under the curve
    # (to eliminate them from the list of candidates)
    gridpoint.under.the.curve <- y[y.grid.index] <= grid.df[x.grid.index]
    w <- which(gridpoint.under.the.curve)
    x <- x[x.grid.index]; y <- y[y.grid.index]
    if (length(w)){x <- x[-w]; y <- y[-w]}
    
    # Eliminate points to the left of the mode, if any
    w <- which(x<mode)
    if (length(w)){x <- x[-w]; y <- y[-w]}
    
    # Eliminate points for which the text would fall out of the plot area
    w <- which((par.usr[2]-str.width[2]) >= x & (y + str.height) <= par.usr[4])
    x <- x[w]; y <- y[w]
    
    # For each height, eliminate the closest point to the curve
    # (we want to stay away from the curve to preserve readability)
    w <- which(!duplicated(y))
    if (length(w)){x <- x[-w]; y <- y[-w]}
    
    # For each point, compute distance from mass center and pick the closest point
    d <- ((x-mass.center[1])^2) + ((y-mass.center[2])/aspect.ratio)^2
    w <- which.min(d)
    x <- x[w]; y <- y[w]
    
    if (length(x))
    {
      text(x, y, labels=p.string[2], adj=c(0,0), col='gray', cex=cex)
    }
    else
    {
      text(plot.xlim[2], mean(par.usr[c(3,4)]), labels=p.string[2], col='gray', cex=cex, srt=-90, adj=c(1,0))
    }
  }
  
  # ......................................................................................................................................
  
  Newton.Raphson <- function(derivative.epsilon, precision, f.cum, p, q, theta.from.moments, start.with.normal.approx, start)
  {
    Hessian <- matrix(NA, 2, 2)
    
    if (start.with.normal.approx)
    {
      # Probably not a very good universal choice, but proved good in most cases in practice
      m <-  diff(q)/diff(p)*(0.5-p[1]) + q[1]
      v <- (diff(q)/diff(qnorm(p)))^2
      theta <- theta.from.moments(m, v)
    }
    else theta <- start
    
    
    change <- precision + 1
    niter <- 0
    # Newton-Raphson multivariate algorithm
    while (max(abs(change)) > precision)
    {
      Hessian[,1] <- (f.cum(q, theta) - f.cum(q, theta - c(derivative.epsilon, 0))) / derivative.epsilon
      Hessian[,2] <- (f.cum(q, theta) - f.cum(q, theta - c(0, derivative.epsilon))) / derivative.epsilon
      
      f <- f.cum(q, theta) - p
      change <- solve(Hessian) %*% f
      last.theta <- theta
      theta <- last.theta - change
      
      # If we step out of limits, reduce change
      
      if (any(theta<0))
      {
        k <- min(last.theta/change)
        theta <- last.theta - k/2*change
      }
      
      niter <- niter + 1
    }
    
    list(theta=as.vector(theta), niter=niter, last.change=as.vector(change))
  }
  
  # ...............................................................................................................
  
  plot.density <- function(p, q, f, f.cum, F.inv, mode, theta, plot.xlim, dens.label, parms.names, cex)
  {
    if (length(plot.xlim) == 0)
    {
      plot.xlim <- F.inv(c(0, 1), theta)
      
      if (is.infinite(plot.xlim[1]))
      {
        tmp <- min(c(0.001, p[1]/10))
        plot.xlim[1] <- F.inv(tmp, theta)
      }  
      
      if (is.infinite(plot.xlim[2]))
      {
        tmp <- max(c(0.999, 1 - (1-p[2])/10))
        plot.xlim[2] <- F.inv(tmp, theta)
      }
    }
    plot.xlim <- sort(plot.xlim)
    
    
    x <- seq(from=min(plot.xlim), to=max(plot.xlim), length=1000)
    h <- f(x, theta)
    x0 <- x; f0 <- h
    ylab <- paste(c(dens.label, '(x, ', parms.names[1], ' = ', round(theta[1], digits=5), ', ', parms.names[2], ' = ', round(theta[2], digits=5), ')'), collapse='')
    plot(x, h, type='l', ylab=ylab)
    
    # fill in area on the left side of the distribution
    x <- seq(from=plot.xlim[1], to=q[1], length=1000)
    y <- f(x, theta)
    x <- c(x, q[1], plot.xlim[1]); y <- c(y, 0, 0)
    polygon(x, y, col='lightgrey', border='lightgray')
    # fill in area on the right side of the distribution
    x <- seq(from=max(plot.xlim), to=q[2], length=1000)
    y <- f(x, theta)
    x <- c(x, q[2], plot.xlim[2]); y <- c(y, 0, 0)
    polygon(x, y, col='lightgrey', border='lightgray')
    # draw distrn again
    points(x0, f0, type='l')
    h <- f(q, theta)
    points(rep(q[1], 2), c(0, h[1]), type='l', col='orange')
    points(rep(q[2], 2), c(0, h[2]), type='l', col='orange')
    # place text on both ends areas
    print.area.text(p, p.check, q, f, f.cum, F.inv, theta, mode, cex, plot.xlim)  
    
    xaxp <- par("xaxp")
    x.ticks <- seq(from=xaxp[1], to=xaxp[2], length=xaxp[3]+1)
    q2print <- as.double(setdiff(as.character(q), as.character(x.ticks)))
    
    mtext(q2print, side=1, col='orange', at=q2print, cex=0.6, line=2.1)
    points(q, rep(par('usr')[3]+0.15*par('cxy')[2], 2), pch=17, col='orange')
  }
  
  #________________________________________________________________________________________________________________
  
  
  parms <- Newton.Raphson(derivative.epsilon, precision, f.cum, p, q, theta.from.moments, start.with.normal.approx, start=start)
  p.check <- f.cum(q, parms$theta)
  
  if (plot) plot.density(p, q, f, f.cum, F.inv, f.mode(parms$theta), parms$theta, plot.xlim, dens.label, parms.names, 0.8)
  
  list(a=parms$theta[1], b=parms$theta[2], last.change=parms$last.change, niter=parms$niter, q=q, p=p, p.check=p.check)
}


### Part 1 -  estimate pooled rate and risk ratio####


# load files

dat.rate <- read.csv("hiv.rate.csv",header = T,na.strings = "")
dat.ratio <- read.csv("hiv.ratio.csv",header = T,na.strings = "")

##=============================================
## calculate meta estimates of rate
##=============================================####

names(dat.rate) <- tolower(names(dat.rate))
dat.rate$outcome.2 <- as.character(dat.rate$outcome.2)

dat.rate$patient.years <- dat.rate$patient.years/10000 

dat.rate$author.year <- paste0(substr(dat.rate$authors..primary,1,regexpr(",",dat.rate$authors..primary)-1)," ","et al, ",dat.rate$pub.year)
dat.rate <- dat.rate[!is.na(dat.rate$master.id),]
dat.rate[dat.rate$master.id==52,"author.year"]<- "INITIO cohort, 2006"
dat.rate[dat.rate$master.id==5115,"author.year"]<- "ATCC cohort, 2006"
dat.rate[dat.rate$master.id==5117,"author.year"]<- "DAD cohort, 2006"


# no of estimates for rate
length(dat.rate$outcome.2)

# no of studies for rate
length(unique(dat.rate$master.id))

# no of patients and total follow up years
x <- dat.rate[!duplicated(dat.rate$master.id),]
sum(x$n)
sum(x$patient.years)*10000


# meta estimate for incident CVD disease (rate)
res.cvd.data <- dat.rate[ dat.rate$outcome.2=="Incidence - CVD",]
res.cvd <- rma.glmm(measure = "IRLN", xi = events, ti = patient.years, data = res.cvd.data,slab = author.year)
forest.cvd <- forest(res.cvd, at=log(c(1,5, 10, 25,50,100,250,500)),atransf =exp,order=order(res.cvd.data$pub.year))
res.cvd <- predict(res.cvd,transf = exp,digits=2)

# meta estimate for incident stroke (rate)
res.stroke.data <- dat.rate[ dat.rate$outcome.2=="Incidence - Stroke",]
res.stroke <- rma.glmm(measure = "IRLN", xi = events, ti = patient.years, data = res.stroke.data,slab = author.year)
forest.stroke <- forest(res.stroke, at=log(c(1,5, 10, 25,50,100,250,500)),atransf =exp,order=order(res.stroke.data$pub.year))
res.stroke <- predict(res.stroke,transf = exp,digits=2)

# meta estimate for incident mi (rate)
res.mi.data <- dat.rate[ dat.rate$outcome.2=="Incidence - MI",]
res.mi <- rma.glmm(measure = "IRLN", xi = events, ti = patient.years, data = res.mi.data,slab = author.year)
forest.mi <- forest(res.mi, at=log(c(1,5, 10, 25,50,100,250,500)),atransf =exp,order=order(res.mi.data$pub.year))
res.mi <- predict(res.mi,transf = exp,digits=2)

# meta estimate for mortality from cvd  (rate)
res.mort.cvd.data <- dat.rate[ dat.rate$outcome.2=="Mortality - CVD",]
res.mort.cvd <- rma.glmm(measure = "IRLN", xi = events, ti = patient.years, data = res.mort.cvd.data,slab = res.mort.cvd.data$author.year)
forest.res.mort.cvd <- forest(res.mort.cvd, at=log(c(1,5, 10, 25,50,100,250,500)),atransf =exp,order=order(res.mort.cvd.data$pub.year))
res.mort.cvd <-predict(res.mort.cvd,transf = exp,digits=2)

# meta estimate for mortality from stroke  (rate)
res.mort.stroke.data <- dat.rate[ dat.rate$outcome.2=="Mortality - Stroke",]
res.mort.stroke <- rma.glmm(measure = "IRLN", xi = events, ti = patient.years, data = res.mort.stroke.data,slab=res.mort.stroke.data$author.year)
forest.res.mort.stroke <- forest(res.mort.stroke, at=log(c(1,5, 10, 25,50,100,250,500)),atransf =exp,order=order(res.mort.stroke.data$pub.year))
res.mort.stroke <-predict(res.mort.stroke,transf = exp,digits=2)

# meta estimate for mortality from mi  (rate)
res.mort.mi.data <- dat.rate[ dat.rate$outcome.2=="Mortality - MI",]
res.mort.mi <- rma.glmm(measure = "IRLN", xi = events, ti = patient.years, data = res.mort.mi.data,slab=res.mort.mi.data$author.year)
forest.res.mort.mi <- forest(res.mort.mi, at=log(c(1,5, 10, 25,50,100,250,500)),atransf =exp,order=order(res.mort.mi.data$pub.year))
res.mort.mi <-predict(res.mort.mi,transf = exp,digits=2)

# summary rates
list <- c(res.cvd,res.mi,res.stroke,res.mort.cvd,res.mort.mi,res.mort.stroke)

summary.rates <- as.matrix(do.call("rbind",list))
summary.rates <- matrix(summary.rates,ncol=6)
summary.rates <- as.data.frame(t(summary.rates))
names(summary.rates) <- c("pred","V2","ll","ul")
summary.rates <-summary.rates[,c("pred","ll","ul")] 
summary.rates$variable <- c("res.cvd","res.mi","res.stroke","res.mort.cvd","res.mort.mi","res.mort.stroke")
summary.rates[c("pred","ll","ul")] <- lapply(summary.rates[c("pred","ll","ul")], function(x){as.numeric(as.character((x)))})

summary.rates$factor <- letters[6:1]

theme_set(theme_bw())
theme_update(
  axis.line = element_line(colour = "black"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  plot.margin = unit(c(0,0,0,0), "lines")
)

# plots summary rates and table

p1 <- ggplot(summary.rates,aes(pred,factor)) + 
  geom_errorbarh(aes(xmax = ul, xmin = ll), height = 0.15) +
  geom_point(size=5, shape=21,colour="black",fill="red") +
  geom_hline(yintercept = 1.5, linetype = "longdash") +
  geom_hline(yintercept=2.5,linetype = "longdash")+
  geom_hline(yintercept=3.5,linetype = "solid")+
  geom_hline(yintercept=4.5,linetype = "longdash")+
  geom_hline(yintercept=5.5,linetype = "longdash")+
  scale_x_continuous(breaks = seq(0,80,10), labels = seq(0,80,10), trans="log2") +
  labs(x="Crude rates per 10,000 persons (95% Confidence Intervals)", y="")

p1

no.estimates <- as.vector(table(dat.rate$outcome.2))
events <- tapply(dat.rate$events,dat.rate$outcome.2,sum)
person.years <- (tapply(dat.rate$patient.years,dat.rate$outcome.2,sum)*10000)/1000
person.years <- round(person.years,1)
prop <- paste0(round(summary.rates$pred,1)," ","(",round(summary.rates$ll,1),"-",round(summary.rates$ul,1),")")

a <- rep(c(letters[6:1]),6)
b <- rep(c(1:6),each=6)
c <- c("Incidence","","","Mortality","","","All CVD","Myocardial infarction","Stroke","All CVD","Myocardial infarction","Stroke",no.estimates,events,person.years,prop)
d <- as.data.frame(cbind(a,b,c))
row.names(d) <- 1:42

data_table.rate <- ggplot(d, aes(x = d$b, y = d$a, label = format(d$c, nsmall = 1))) +
  geom_text(size = 4, hjust=0, vjust=0) + theme_bw()  + 
  geom_hline(yintercept=1.5,linetype="longdash")+
  geom_hline(yintercept=2.5,linetype="longdash")+
  geom_hline(yintercept=3.5,linetype="solid")+
  geom_hline(yintercept=4.5,linetype="longdash")+
  geom_hline(yintercept=5.5,linetype="longdash")+
  theme(panel.grid.major = element_blank(), 
        legend.position = "none",
        panel.border = element_blank(), 
        axis.text.x = element_text(colour="white"),#element_blank(),
        axis.text.y = element_blank(), 
        axis.ticks = element_line(colour="white"),#element_blank(),
        plot.margin = unit(c(0,0,0,0), "lines")) +
  labs(x="",y="") +
  coord_cartesian(xlim=c(1,6))

data_table.rate

p1 <- grid.arrange(data_table.rate, p1, ncol=2)


##============================================
## calculate meta-estimates for risk ratio
##============================================

names(dat.ratio) <- tolower(names(dat.ratio))
dat.ratio <- subset(dat.ratio,!is.na(dat.ratio$outcome2))
dat.ratio <- subset(dat.ratio,!is.na(dat.ratio$ratio.ll))

dat.ratio$lgrr <- logb(dat.ratio$ratio.ce)
dat.ratio$lgll <- logb(dat.ratio$ratio.ll)
dat.ratio$lgul <- logb(dat.ratio$ratio.ul)
dat.ratio$se <- (dat.ratio$lgul-dat.ratio$lgll)/(2*1.96)

rr <- rma(yi = dat.ratio$lgrr,sei = dat.ratio$se, method = "ML") 
rr.exp <- predict(rr, transf = exp, digits=2)

# sampling from probability distributions to derive RR and uncertainty estimate with 95% CI
num_iter <- 10000
rr.sim <- exp(rnorm(num_iter, as.vector(rr$b), rr$se))
(quantile(rr.sim, probs = c(0.025,0.5, 0.975)))
hist((rr.sim))

# plot for ratio

dat.ratio$author <- paste0(substr(dat.ratio$authors..primary,1,regexpr(",",dat.ratio$authors..primary)-1)," ","et al, ",dat.ratio$pub.year)

dat.ratio.plot <- dat.ratio[c("author","pub.year","ratio.ce","ratio.ll","ratio.ul","se")]
dat.ratio.plot <- dat.ratio.plot[order(dat.ratio.plot$pub.year),]
dat.ratio.plot[,c("ratio.ce","ratio.ll","ratio.ul")] <- lapply(dat.ratio.plot[,c("ratio.ce","ratio.ll","ratio.ul")], as.numeric)
dat.ratio.plot[,c("ratio.ce","ratio.ll","ratio.ul")] <- lapply(dat.ratio.plot[,c("ratio.ce","ratio.ll","ratio.ul")], function(x){format(round(as.numeric(x),digits=2),nsmall=2)})
dat.ratio.plot$risk <-paste0(dat.ratio.plot$ratio.ce," ","(",dat.ratio.plot$ratio.ll,"-",dat.ratio.plot$ratio.ul,")") 

dat.ratio.plot <- rbind(dat.ratio.plot,
                        c(rep(NA,7)),
                        c("Pooled estimate",NA,round(rr.exp[[1]],2),round(rr.exp[[3]],2),round(rr.exp[[4]],2),0.05,NA))

dat.ratio.plot$risk <-paste0(dat.ratio.plot$ratio.ce," ","(",dat.ratio.plot$ratio.ll,"-",dat.ratio.plot$ratio.ul,")") 

dat.ratio.plot$factor <- c(rep("a",19),"b")

dat.ratio.plot$factor2 <- letters[1:20]

dat.ratio.plot[,c("ratio.ce","ratio.ll","ratio.ul")] <- lapply(dat.ratio.plot[,c("ratio.ce","ratio.ll","ratio.ul")], function(x){round(as.numeric(x),digits=1)})


## forest plot for ratio
p2 <- ggplot(dat.ratio.plot,aes(ratio.ce,rev(factor2),)) + 
  geom_errorbarh(aes(xmax = ratio.ul, xmin = ratio.ll), height = 0.15) +
  geom_point(aes(size=1/as.numeric(se)^2), shape=21,colour="black",fill="red") +
  scale_x_continuous(trans="log2",breaks = c(0.5,1,2,seq(4,36,4)), labels = c(0.5,1,2,seq(4,36,4))) +
  labs(x="Risk ratio (95% Confidence intervals)", y="")+
  geom_hline(yintercept=0.5,linetype="solid")+
  geom_vline(xintercept = 1,linetype="dashed")+
  theme(legend.position="none")

p2

a <- rep(c(letters[20:1]),2)
b <- rep(c(1:2),each=20)
c <- c(dat.ratio.plot$author[1:18],"",dat.ratio.plot$author[20],dat.ratio.plot$risk[1:18],"",dat.ratio.plot$risk[20])
d <- as.data.frame(cbind(a,b,c))


data_table.ratio <- ggplot(d, aes(x = d$b, y = d$a, label = format(d$c, nsmall = 1))) +
  geom_text(size = 3, hjust=0, vjust=0) + theme_bw()+
  theme(panel.grid.major = element_blank(), 
        legend.position = "none",
        panel.border = element_blank(), 
        axis.text.x = element_text(colour="white"),#element_blank(),
        axis.text.y = element_blank(), 
        axis.ticks = element_line(colour="white"),#element_blank(),
        plot.margin = unit(c(0,0,0,0), "lines")) +
  labs(x="",y="") +
  coord_cartesian(xlim=c(1,2))

data_table.ratio

p2 <- grid.arrange(data_table.ratio, p2, ncol=2)

# subgroup analysis of ratio


# by outcome
levels(dat.ratio$outcome2)
rr.stroke <- rma(yi = lgrr,sei = se, data = dat.ratio[dat.ratio$outcome2=="Stroke",],method = "ML") 
predict(rr.stroke, transf = exp, digits=2)

rr.cvd <- rma(yi = lgrr,sei = se, data = dat.ratio[dat.ratio$outcome2=="CVD",],method = "ML") 
predict(rr.cvd, transf = exp, digits=2)

rr.mi <- rma(yi = lgrr,sei = se, data = dat.ratio[dat.ratio$outcome2=="MI",],method = "ML") 
predict(rr.mi, transf = exp, digits=2)

# by adjustment - low
rr.adjust.low <- rma(yi = lgrr,sei = se, data = dat.ratio[dat.ratio$risk.of.bias=="low",],method = "ML")
predict(rr.adjust.low,transf = exp,digits=2)

# by adjustment - moderate/high
rr.adjust.mod.high <- rma(yi = lgrr,sei = se, data = dat.ratio[!dat.ratio$risk.of.bias=="low",],method = "ML")
predict(rr.adjust.mod.high,transf = exp,digits=2)

# by median publication year
rr.recent <- rma(yi = lgrr,sei = se, data = dat.ratio[dat.ratio$pub.year>=median(dat.ratio$pub.year),],method = "ML")
predict(rr.recent,transf = exp,digits=2)

rr.old <- rma(yi = lgrr,sei = se, data = dat.ratio[dat.ratio$pub.year<median(dat.ratio$pub.year),],method = "ML")
predict(rr.old,transf = exp,digits=2)


### metaregression by age

rr.age <- rma(yi = lgrr,sei = se, 
              data = dat.ratio, 
              mods = dat.ratio[, "age.at.baseline"], 
              method = "ML")


## calculate predicted risk ratios for 0 to 60 degrees absolute latitude
preds <- predict(rr.age, newmods=c(30:70), transf=exp)

## calculate point sizes by rescaling the standard errors
wi    <- 1/sqrt(dat.ratio$se^2)
size  <- 0.5 + 3.0 * (wi - min(wi))/(max(wi) - min(wi))

## plot the risk ratios against absolute latitude
plot(dat.ratio$age.at.baseline, exp(dat.ratio$lgrr), pch=19, cex=size, 
     xlab="Age at baseline", ylab="Risk Ratio",
     las=1, bty="l", log="y")

## add predicted values (and corresponding CI bounds)
lines(30:70, preds$pred)
lines(30:70, preds$ci.lb, lty="dashed")
lines(30:70, preds$ci.ub, lty="dashed")

### dotted line at RR=1 (no difference between groups)
abline(h=1, lty="dashed")

### metaregression by time
dat.ratio <- subset(dat.ratio,dat.ratio$pub.year>2005)
rr.time <- rma(yi = lgrr,sei = se, 
              data = dat.ratio, 
              mods = dat.ratio[, "pub.year"], 
              method = "ML")


## calculate predicted risk ratios for 0 to 60 degrees absolute latitude
preds <- predict(rr.time, newmods=c(1990:2015), transf=exp)

## calculate point sizes by rescaling the standard errors
wi    <- 1/dat.ratio$se^2
size  <- 0.5 + 3.0 * (wi - min(wi))/(max(wi) - min(wi))

## plot the risk ratios against absolute latitude
plot(dat.ratio$pub.year, exp(dat.ratio$lgrr), pch=19, cex=size, 
     xlab="year", ylab="Risk Ratio",
     las=1, bty="l", log="y")

## add predicted values (and corresponding CI bounds)
lines(1990:2015, preds$pred)
lines(1990:2015, preds$ci.lb, lty="dashed")
lines(1990:2015, preds$ci.ub, lty="dashed")

### dotted line at RR=1 (no difference between groups)
abline(h=1, lty="dashed")


# publication bias
regtest(rr)
predict(trimfill(rr))
funnel(rr)
funnel(trimfill(rr))



###=========================================================
### Part 2: estimating burden using daly and prevalence data
###=========================================================

##==============================
## Part 2a - sort IHME CVD data 
##==============================

# load country data into R with ISO country code

country.un <- read.csv("unaids.regions.csv", header = T)
country.un$countryid <- 1:182


# load daly data - country specific but in one file

dat <- read.csv("IHME-GBD_2015_DATA.csv",header = T)
dat$location_name <- as.character(dat$location_name)

#match country names to UN and relabel unmatched names 

combined <- merge(x = country.un,y = dat,by.x = "country",by.y ="location_name" ,all.y=T)
str(combined)
combined[,c("country","unaids.region","unaids.region.id")] <- lapply(combined[,c("country","unaids.region","unaids.region.id")], as.character)
names(table(combined$unaids.region.id))

combined$unaids.region.id <- ifelse(combined$country=="Global","global",combined$unaids.region.id)
combined$unaids.region <- ifelse(combined$country=="Global","global",combined$unaids.region)

combined <- subset(combined, !is.na(combined$unaids.region.id))
names(table(combined$country))

combined.countrylevel <- combined


# file for all patients

names(combined) <- c("country", "unaids.region", "unaids.region.id", "countryid", 
                       "measure_id", "measure_name", "location_id", "sex_id", "sex_name", 
                       "age_id", "age_name", "cause_id", "cause_name", "metric_id", 
                       "metric_name", "year_id", "nm_mean", "nm_upper", "nm_lower")

# file for all patients
x <- function(x){rbind(x, c("total",names(table(y$unaids.region)),rep(NA,3),names(table(y$year_id)),colSums(y[7:9])))}
combined.all <- combined[(!combined$age_name=="Age-standardized"&combined$sex_name=="Both"),c("country", "unaids.region","cause_name","sex_name","age_name","year_id","nm_mean", "nm_lower", "nm_upper")]
combined.all <- split(combined.all,list(combined.all$unaids.region,combined.all$year_id))
combined.all <- lapply(combined.all,function(x){rbind(x, c("total",names(table(x$unaids.region)),rep(NA,3),names(table(x$year_id)),colSums(x[7:9])))})

combined.all <- rbind.fill((combined.all))

total.both <- combined.all[combined.all$country=="total",]

# file for male patients
x <- function(x){rbind(x, c("total",names(table(y$unaids.region)),rep(NA,3),names(table(y$year_id)),colSums(y[7:9])))}
combined.male <- combined[(!combined$age_name=="Age-standardized"&combined$sex_name=="Male"),c("country", "unaids.region","cause_name","sex_name","age_name","year_id","nm_mean", "nm_lower", "nm_upper")]
combined.male <- split(combined.male,list(combined.male$unaids.region,combined.male$year_id))
combined.male <- lapply(combined.male,function(x){rbind(x, c("total.male",names(table(x$unaids.region)),rep(NA,3),names(table(x$year_id)),colSums(x[7:9])))})
combined.male <- rbind.fill((combined.male))
total.males <- combined.male[combined.male$country=="total.male",]


# file for female patients
x <- function(x){rbind(x, c("total",names(table(y$unaids.region)),rep(NA,3),names(table(y$year_id)),colSums(y[7:9])))}
combined.female <- combined[(!combined$age_name=="Age-standardized"&combined$sex_name=="Female"),c("country", "unaids.region","cause_name","sex_name","age_name","year_id","nm_mean", "nm_lower", "nm_upper")]
combined.female <- split(combined.female,list(combined.female$unaids.region,combined.female$year_id))
combined.female <- lapply(combined.female,function(x){rbind(x, c("total.female",names(table(x$unaids.region)),rep(NA,3),names(table(x$year_id)),colSums(x[7:9])))})
combined.female <- rbind.fill((combined.female))
total.females <- combined.female[combined.female$country=="total.female",]

# file for all patients, stratified by gender

total.all <- as.data.frame(rbind(total.both,total.males,total.females))
total.all <- total.all[,c("country", "unaids.region", "year_id", "nm_mean", "nm_lower", "nm_upper")]

total.all$compound.key <- paste0(total.all$country,total.all$unaids.region)
total.all$nm_mean <- as.numeric(total.all$nm_mean)
total.all$nm_lower <- as.numeric(total.all$nm_lower)
total.all$nm_upper <- as.numeric(total.all$nm_upper)


# daly.ce - estimating central daly estimate startified by sex, year and region

daly.ce <-total.all[,c(3,4,7)] 

daly.ce <- split(daly.ce, daly.ce$compound.key)

fun1 <- function(x){
  z <- approx(x[["year_id"]],x[["nm_mean"]],xout=1990:2015)
  z <- unlist(z)
  as.data.frame(matrix(z,ncol=2))
}

daly.ce <- lapply(daly.ce,fun1)
daly.ce <- do.call(rbind.data.frame, daly.ce)
daly.ce$V3 <- row.names(daly.ce)
daly.ce$V3 <- gsub('[[:digit:]]+', '', daly.ce$V3)
daly.ce$V3 <- str_sub(daly.ce$V3,-nchar(daly.ce$V3),-2)
names(daly.ce) <- c("year","daly.ce","subgroup")

# daly.ll - estimating daly.ll estimate startified by sex, year and region

daly.ll <-total.all[,c(3,5,7)] 

daly.ll <- split(daly.ll, daly.ll$compound.key)

fun1 <- function(x){
  z <- approx(x[["year_id"]],x[["nm_lower"]],xout=1990:2015)
  z <- unlist(z)
  as.data.frame(matrix(z,ncol=2))
}

daly.ll <- lapply(daly.ll,fun1)
daly.ll <- do.call(rbind.data.frame, daly.ll)
daly.ll$V3 <- row.names(daly.ll)
daly.ll$V3 <- gsub('[[:digit:]]+', '', daly.ll$V3)
daly.ll$V3 <- str_sub(daly.ll$V3,-nchar(daly.ll$V3),-2)
names(daly.ll) <- c("year","daly.ll","subgroup")

# daly.ul - estimating daly.ul estimate startified by sex, year and region

daly.ul <-total.all[,c(3,6,7)] 

daly.ul <- split(daly.ul, daly.ul$compound.key)

fun1 <- function(x){
  z <- approx(x[["year_id"]],x[["nm_upper"]],xout=1990:2015)
  z <- unlist(z)
  as.data.frame(matrix(z,ncol=2))
}

daly.ul <- lapply(daly.ul,fun1)
daly.ul <- do.call(rbind.data.frame, daly.ul)
daly.ul$V3 <- row.names(daly.ul)
daly.ul$V3 <- gsub('[[:digit:]]+', '', daly.ul$V3)
daly.ul$V3 <- str_sub(daly.ul$V3,-nchar(daly.ul$V3),-2)
names(daly.ul) <- c("year","daly.ul","subgroup")

# combined daly

daly <- as.data.frame(cbind(daly.ce,daly.ll,daly.ul))

daly <- daly[,c("year", "daly.ce", "daly.ll", "daly.ul", "subgroup")]


rownames(daly) <- 1:dim(daly)[1]

daly$area <- gsub("[.]",'',daly$subgroup)
daly$area <- gsub("totalmale|totalfemale|total",'',daly$area)

pattern1=c("total")
pattern2=c("totalmale")
pattern3=c("totalfemale")

daly$sex1 <- gsub("[.]",'',daly$subgroup)
daly$sex <- str_extract(daly$sex1,paste(rev(mget(ls(pattern = "pattern\\d+"))), collapse="|"))


daly$sex <- factor(daly$sex,
                  levels=c("total","totalmale","totalfemale"),
                  labels=c("adults","male","female"))

daly$area.id <- factor(daly$area,
                            levels=c("global", "UNAIDS Region - Asia and the Pacific", "UNAIDS Region - East and Southern Africa", 
                              "UNAIDS Region - Eastern Europe and Central Asia", "UNAIDS Region - Latin America and the Caribbean", 
                              "UNAIDS Region - Middle East and North Africa", "UNAIDS Region - West and Central Africa", 
                              "UNAIDS Region - Western & Central Europe and North America"),
                            labels=c("All countries","UNAIDS Region - Asia and the Pacific", "UNAIDS Region - East and Southern Africa", 
                                     "UNAIDS Region - Eastern Europe and Central Asia", "UNAIDS Region - Latin America and the Caribbean", 
                                     "UNAIDS Region - Middle East and North Africa", "UNAIDS Region - West and Central Africa", 
                                     "UNAIDS Region - Western & Central Europe and North America"))

##=================================
## Part 2b - load prevalence files
##=================================


hiv.prev <- read.csv("hiv.prev.unaids.210716.csv",header = T,stringsAsFactors = F)

names(hiv.prev) <- tolower(names(hiv.prev))


### scrub prevalence data

## remove estimates relating to young people

hiv.prev <- hiv.prev[grepl("Young",hiv.prev$subgroup)==F,]

## data with all countries 

hiv.prev.country <- hiv.prev
hiv.prev.country.all <- hiv.prev

## remove estimates from countries only leaving global and regional data

pattern <- "03M49WLD|UNAAP|UNAEECA|UNAESA|UNALAC|UNAMENA|UNAWCA|UNAWCENA"

hiv.prev <- hiv.prev[grepl(pattern,hiv.prev$area.id)==T,]

hiv.prev$subgroup2 <- factor(hiv.prev$subgroup,
                             levels=c("Adults (15-49) estimate modelled", "Adults (15-49) lower estimate modelled", 
                                      "Adults (15-49) upper estimate modelled", "Females Adults (15-49) estimate", 
                                      "Females Adults (15-49) lower estimate", "Females Adults (15-49) upper estimate", 
                                      "Males Adults (15-49) estimate", "Males Adults (15-49) lower estimate", 
                                      "Males Adults (15-49) upper estimate"),
                             labels = c("ce","ll","ul","ce","ll","ul","ce","ll","ul"))

hiv.prev$subgroup <- factor(hiv.prev$subgroup,
                            levels=c("Adults (15-49) estimate modelled", "Adults (15-49) lower estimate modelled", 
                                     "Adults (15-49) upper estimate modelled", "Females Adults (15-49) estimate", 
                                     "Females Adults (15-49) lower estimate", "Females Adults (15-49) upper estimate", 
                                     "Males Adults (15-49) estimate", "Males Adults (15-49) lower estimate", 
                                     "Males Adults (15-49) upper estimate"),
                            labels = c("adults","adults","adults","female","female","female","male","male","male"))

hiv.prev <- hiv.prev[c("subgroup", "area", "area.id", "time.period", "data.value", "subgroup2")]

hiv.prev_wide <- spread(hiv.prev, subgroup2, data.value)


##=========================================
## Part 2c - merge prevelance data to daly
##=========================================

final <- merge(x=hiv.prev_wide,y=daly,by.x = c("subgroup","area","time.period"), by.y=c("sex","area.id","year"),all = T)

## derive alpha and beta using beta pharm function for prevalence data 

# convert percentage into decimals

final[,c("ce","ll","ul")] <- lapply(final[,c("ce","ll","ul")],function(x){x/100})

alpha <- function(x){
  alpha <- as.vector(beta.parms.from.quantiles(c(as.numeric(x[6]),as.numeric(x[7])))[1])
}

beta <- function(x){
  beta <- as.vector(beta.parms.from.quantiles(c(as.numeric(x[6]),as.numeric(x[7])))[2])
}

final$alpha <- as.vector(unlist(apply(final, 1,alpha)))
final$beta <- as.vector(unlist(apply(final, 1,beta)))

# check beta results
final$uniq_id <- seq(1000, 1000+nrow(final)-1, 1)

check_res <- by(final, final$uniq_id, function(x) qbeta(shape1 = x$alpha, shape2= x$beta, 
                                           p = c(0.025, 0.5, 0.975)))
check_res <- do.call(rbind, check_res)
colnames(check_res) <- c("est", "lci", "uci")
check_res <- cbind(final, check_res)
plot(check_res$ll, check_res$lci)
plot(check_res$ul, check_res$uci)
plot(check_res$ce, check_res$est)
rm(check_res)
final$uniq_id <- NULL

a <- as.list(NULL)
num_iter <- 10000

## derive central estimate and 95% uncertainty term for CVD daly using log distribution

final[,c("daly.ce","daly.ll","daly.ul")] <- lapply(final[,c("daly.ce","daly.ll","daly.ul")],function(x){log(x,base = exp(1))})


final$nm_se <- (as.numeric(final$daly.ul)-as.numeric(final$daly.ll))/3.92
final$nm_mean <- as.numeric(final$daly.ce)

# check if log distribution and sampling method appropriate

final2 <- final

for(i in 1:nrow(final2)){
  y <- rnorm(10000,final2$nm_mean[i],final2$nm_se[i])
  
  final2$daly.ce.test[i] <- quantile(y,0.5)
  final2$daly.ll.test[i] <- quantile(y,0.025)
  final2$daly.ul.test[i] <- quantile(y,0.975)
}



hist(exp(final2$daly.ce.test)-exp(final2$daly.ce),breaks = 100)
hist(exp(final2$daly.ll.test)-exp(final2$daly.ll),breaks = 100)
hist(exp(final2$daly.ul.test)-exp(final2$daly.ul),breaks = 100)

###=======================================================
### Part 3 = calculate paf and corresponding CVD burden due to HIV 
###=======================================================
newvars <- data.frame(id = 1:nrow(final), 
                      burdence = NA,
                      burdenll = NA,
                      burdenul = NA,
                      pafce= NA,
                      pafll = NA,
                      paful = NA)
newvars$id <- NULL


final <- cbind(final, newvars)
for(i in 1:nrow(final)){
  x <- rbeta(num_iter, shape1 = final$alpha[i], shape2 = final$beta[i])
  y <- rnorm(num_iter,final$nm_mean[i],final$nm_se[i])
  y <- exp(y)
  z <- rr.sim
  p <- (x*(rr.sim-1))/(1+x*(rr.sim-1))
  burden <- p*y
  
  
  
  final$burdence[i] <- quantile(burden,0.5)
  final$burdenll[i] <- quantile(burden,0.025)
  final$burdenul[i] <- quantile(burden,0.975)
  final$pafce[i] <- quantile(p,0.5)
  final$pafll[i] <- quantile(p,0.025)
  final$paful[i] <- quantile(p,0.975)
}

# in thousands

final$burdence <-final$burdence/1000
final$burdenll <-final$burdenll/1000
final$burdenul <-final$burdenul/1000



# plot of all daly
dat <-  final[final$subgroup=="adults"&final$area=="All countries",]
p1 <- ggplot(data <- dat) 
p1 <- p1+
  geom_line(aes(y=dat$burdence, x=dat$time.period), colour = "blue")+
  geom_line(aes(y=dat$burdenll, x=dat$time.period), colour = "red",linetype="dashed")+
  geom_line(aes(y=dat$burdenul, x=dat$time.period), colour = "red",linetype="dashed")+
  scale_x_continuous(name="\nTime (year)")+
  scale_y_continuous(name="DALYs in thousands\n", limits=c(0,4000),breaks=seq(0,4000,500),labels=c("0","500","1,000","1,500","2,000","2,500","3,000","3,500","4,000"))+
  theme(legend.position="none")+
  theme_light()
  
p1

# plot of all daly
dat <-  final[!final$subgroup=="adults" & final$area.id=="03M49WLD",]

p2 <- ggplot(data <- dat, aes(y=dat$burdence, x=dat$time.period, ymin=dat$burdenll, ymax=dat$burdenul))+
      geom_ribbon(aes(fill=dat$subgroup), alpha = 0.2)+
      geom_line(aes(color=dat$subgroup))+
      scale_x_continuous(name="\nTime (year)")+
      scale_y_continuous(name="DALYs in thousands\n", limits=c(0,2500),breaks=seq(0,2500,500),labels=c("0","500","1,000","1,500","2,000","2,500"))+
      theme_light()+theme(legend.position='none')

p2

# plot by region
final$area.id <- factor(final$area.id,levels=c("03M49WLD",
                                               "UNAMENA",
                                               "UNAWCENA",
                                               "UNALAC",
                                               "UNAEECA",
                                               "UNAWCA",
                                               "UNAAP",
                                               "UNAESA"))


p3 <- ggplot(arrange(final[!final$area=="All countries"&final$subgroup=="adults",],area.id), aes(x=as.numeric(as.character(time.period)), y=burdence,order=area.id,fill=area.id))
p3 <- p3 + geom_area(alpha=0.8)+scale_fill_brewer(palette="YlOrRd", direction = 1,labels=c("Middle East and North Africa","Western & Central Europe and North America","Latin America and the Caribbean","Eastern Europe and Central Asia","West and Central Africa","Asia and the Pacific","East and Southern Africa")) +
  scale_x_continuous(name="\nTime (year)")+
  scale_y_continuous(name="DALYs in thousands\n", limits=c(0,2000),breaks=seq(0,2000,500),labels=c("0","500","1,000","1,500","2,000"))+
  theme_light()+
  theme(legend.position = c(0.2, 0.83),legend.text = element_text(size=10),legend.title = element_blank())
 

p3

# plot by 22 high HIV burden states


priority.countries <- "\\bAngola\\b|\\bBotswana|Burundi|Cameroon|Chad|Cote d'Ivoire|Democratic Republic of the Congo|Ethiopia|Ghana|India|Kenya|Lesotho|Malawi|Mozambique|Namibia|Nigeria|South Africa|Swaziland|Tanzania|Uganda|Zambia|Zimbabwe"
priority.countries1 <- "\\bAGO\\b|\\bBWA\\b|\\bBDI\\b|\\bCMR\\b|\\bTCD\\b|\\bCIV\\b|\\bCOD\\b|\\bETH\\b|\\bGHA\\b|\\bIndia\\b|\\bKEN\\b|\\bLSO\\b|\\bMWI\\b|\\bMOZ\\b|\\bNAM\\b|\\bNGA\\b|\\bZAF\\b|\\bSWZ\\b|\\bTZA\\b|\\bUGA\\b|\\bZMB\\b|\\bZWE\\b"

# file for all patients
x <- function(x){rbind(x, c("total",names(table(y$country)),rep(NA,3),names(table(y$year_id)),colSums(y[7:9])))}
combined.pc <- combined[(!combined$age_name=="Age-standardized"&combined$sex_name=="Both"),c("country", "unaids.region","cause_name","sex_name","age_name","year_id","nm_mean", "nm_lower", "nm_upper")]
combined.pc <- split(combined.pc,list(combined.pc$country,combined.pc$year_id))
combined.pc <- lapply(combined.pc,function(x){rbind(x, c("total",names(table(x$country)),rep(NA,3),names(table(x$year_id)),colSums(x[7:9])))})

combined.pc <- rbind.fill((combined.pc))

combined.pc <- combined.pc[combined.pc$country=="total",]

combined.pc$pc <- grepl(priority.countries,combined.pc$unaids.region)
combined.pc <- subset(combined.pc,combined.pc$pc==T)[,-c(3:5,10)]

# estimating the daly by country for 1990 to 2013

daly.ce.pc <- combined.pc[,c(2:4)]  
daly.ce.pc <- split(daly.ce.pc, daly.ce.pc$unaids.region)

fun1 <- function(x){
  z <- approx(x[["year_id"]],x[["nm_mean"]],xout=1990:2015)
  z <- unlist(z)
  as.data.frame(matrix(z,ncol=2))
}

daly.ce.pc <- lapply(daly.ce.pc,fun1)
daly.ce.pc <- do.call(rbind.data.frame, daly.ce.pc)
daly.ce.pc$V3 <- row.names(daly.ce.pc)
daly.ce.pc$V3 <- gsub('[[:digit:]]+', '', daly.ce.pc$V3)
daly.ce.pc$V3 <- str_sub(daly.ce.pc$V3,-nchar(daly.ce.pc$V3),-2)
names(daly.ce.pc) <- c("year","daly.ce","country")

# estimating the daly ll by country for 1990 to 2015

daly.ll.pc <- combined.pc[,c(2,3,5)]  
daly.ll.pc <- split(daly.ll.pc, daly.ll.pc$unaids.region)

fun2 <- function(x){
  z <- approx(x[["year_id"]],x[["nm_lower"]],xout=1990:2015)
  z <- unlist(z)
  as.data.frame(matrix(z,ncol=2))
}

daly.ll.pc <- lapply(daly.ll.pc,fun2)
daly.ll.pc <- do.call(rbind.data.frame, daly.ll.pc)
daly.ll.pc$V3 <- row.names(daly.ll.pc)
daly.ll.pc$V3 <- gsub('[[:digit:]]+', '', daly.ll.pc$V3)
daly.ll.pc$V3 <- str_sub(daly.ll.pc$V3,-nchar(daly.ll.pc$V3),-2)
names(daly.ll.pc) <- c("year","daly.ll","country")


# estimating the daly ll by country for 1990 to 2015

daly.ul.pc <- combined.pc[,c(2,3,6)]  
daly.ul.pc <- split(daly.ul.pc, daly.ul.pc$unaids.region)

fun3 <- function(x){
  z <- approx(x[["year_id"]],x[["nm_upper"]],xout=1990:2015)
  z <- unlist(z)
  as.data.frame(matrix(z,ncol=2))
}

daly.ul.pc <- lapply(daly.ul.pc,fun3)
daly.ul.pc <- do.call(rbind.data.frame, daly.ul.pc)
daly.ul.pc$V3 <- row.names(daly.ul.pc)
daly.ul.pc$V3 <- gsub('[[:digit:]]+', '', daly.ul.pc$V3)
daly.ul.pc$V3 <- str_sub(daly.ul.pc$V3,-nchar(daly.ul.pc$V3),-2)
names(daly.ul.pc) <- c("year","daly.ul","country")

# combined daly

daly.pc <- as.data.frame(cbind(daly.ce.pc,daly.ll.pc,daly.ul.pc))
daly.pc <- daly.pc[,c("year", "country","daly.ce", "daly.ll", "daly.ul")]
rownames(daly.pc) <- 1:dim(daly.pc)[1]
names(daly.pc)


daly.pc$country <- ifelse(daly.pc$country=="Cote d'Ivoire","Cote dvoire",daly.pc$country)
daly.pc$country <- ifelse(daly.pc$country=="United Republic of Tanzania","Tanzania",daly.pc$country)

# prevalence data for 21 priority countries - prevalence rates for india not available

hiv.prev.country <- hiv.prev.country[grepl(priority.countries1,hiv.prev.country$area.id)==T,]

hiv.prev.country$area.id <- noquote(hiv.prev.country$area.id)

hiv.prev.country$subgroup2 <- factor(hiv.prev.country$subgroup,
                             levels=c("Adults (15-49) estimate modelled", "Adults (15-49) lower estimate modelled", 
                                      "Adults (15-49) upper estimate modelled", "Females Adults (15-49) estimate", 
                                      "Females Adults (15-49) lower estimate", "Females Adults (15-49) upper estimate", 
                                      "Males Adults (15-49) estimate", "Males Adults (15-49) lower estimate", 
                                      "Males Adults (15-49) upper estimate"),
                             labels = c("ce","ll","ul","ce","ll","ul","ce","ll","ul"))

hiv.prev.country$subgroup <- factor(hiv.prev.country$subgroup,
                            levels=c("Adults (15-49) estimate modelled", "Adults (15-49) lower estimate modelled", 
                                     "Adults (15-49) upper estimate modelled", "Females Adults (15-49) estimate", 
                                     "Females Adults (15-49) lower estimate", "Females Adults (15-49) upper estimate", 
                                     "Males Adults (15-49) estimate", "Males Adults (15-49) lower estimate", 
                                     "Males Adults (15-49) upper estimate"),
                            labels = c("adults","adults","adults","female","female","female","male","male","male"))

hiv.prev.country <- hiv.prev.country[c("subgroup", "area", "area.id", "time.period", "data.value", "subgroup2")]

hiv.prev.country_wide <- spread(hiv.prev.country, subgroup2, data.value)
hiv.prev.country_wide$area <- ifelse(hiv.prev.country_wide$area.id=="CIV","Cote dvoire",hiv.prev.country_wide$area)
hiv.prev.country_wide$area <- ifelse(hiv.prev.country_wide$area.id=="TZA","Tanzania",hiv.prev.country_wide$area)


hiv.prev.country_wide <- subset(hiv.prev.country_wide,hiv.prev.country_wide$subgroup=="adults")

# combine data for daly and hiv prev for 21 hgh burden countries

daly.pc <- merge(daly.pc,hiv.prev.country_wide,by.x = c("country", "year"), by.y=c("area","time.period"),all.y = T)

# convert percentage into decimals

daly.pc[,c("ce","ll","ul")] <- lapply(daly.pc[,c("ce","ll","ul")],function(x){x/100})

alpha <- function(x){
  alpha <- as.vector(beta.parms.from.quantiles(c(as.numeric(x[9]),as.numeric(x[10])))[1])
}

beta <- function(x){
  beta <- as.vector(beta.parms.from.quantiles(c(as.numeric(x[9]),as.numeric(x[10])))[2])
}

daly.pc$alpha <- as.vector(unlist(apply(daly.pc, 1,alpha)))
daly.pc$beta <- as.vector(unlist(apply(daly.pc, 1,beta)))

a <- as.list(NULL)
num_iter <- 10000

daly.pc[,c("daly.ce","daly.ul","daly.ll")] <- lapply(daly.pc[,c("daly.ce","daly.ul","daly.ll")],function(x){log(x,base=exp(1))})

daly.pc$nm_se <- (as.numeric(daly.pc$daly.ul)-as.numeric(daly.pc$daly.ll))/3.92
daly.pc$nm_mean <- as.numeric(daly.pc$daly.ce)


for(i in 1:nrow(daly.pc)){
  x <- rbeta(num_iter, shape1 = daly.pc$alpha[i], shape2 = daly.pc$beta[i])
  y <- rnorm(num_iter,daly.pc$nm_mean[i],daly.pc$nm_se[i])
  y <- exp(y)
  z <- rr.sim
  p <- (x*(rr.sim-1))/(1+x*(rr.sim-1))
  burden <- p*y
  
  
  
  daly.pc$burdence[i] <- quantile(burden,0.5,na.rm = T)
  daly.pc$burdenll[i] <- quantile(burden,0.025,na.rm = T)
  daly.pc$burdenul[i] <- quantile(burden,0.975,na.rm = T)
  daly.pc$pafce[i] <- quantile(p,0.5,na.rm = T)
  daly.pc$pafll[i] <- quantile(p,0.025,na.rm = T)
  daly.pc$paful[i] <- quantile(p,0.975,na.rm = T)
}

ce <- (by(daly.pc,daly.pc$year,function(x){sum(x$burdence)})[1:26])
ll <- (by(daly.pc,daly.pc$year,function(x){sum(x$burdenll)})[1:26])
ul <- (by(daly.pc,daly.pc$year,function(x){sum(x$burdenul)})[1:26])

daly.pc.plot <- as.data.frame(cbind(1990:2015,ce,ll,ul))


pc.plot <- ggplot(data = daly.pc.plot)+
  geom_line(aes(y=daly.pc.plot$ce/1000, x=daly.pc.plot$V1), colour = "blue")+
  geom_line(aes(y=daly.pc.plot$ll/1000, x=daly.pc.plot$V1), colour = "red",linetype="dashed")+
  geom_line(aes(y=daly.pc.plot$ul/1000, x=daly.pc.plot$V1), colour = "red",linetype="dashed")+
  scale_x_continuous(name="\nTime (year)")+
  scale_y_continuous(name="DALYs in thousands\n", limits=c(0,1500),breaks=seq(0,1000,500))+
  theme(legend.position="none")+
  theme_light()

pc.plot

# in thousands

daly.pc$burdence <-daly.pc$burdence/1000
daly.pc$burdenll <-daly.pc$burdenll/1000
daly.pc$burdenul <-daly.pc$burdenul/1000

# ============================================================
# paf by country from UNAIDS data - has data on all countries
# ============================================================


all1 <- read.csv("unaids.country.csv")
all <- all1


## note to me - need to download rate data
# extract standardised rates for all countries to calculate country specific rates

all.daly <- combined.countrylevel
names(all.daly) <- c("country", "unaids.region", "unaids.region.id", "countryid", 
                     "measure_id", "measure_name", "location_id", "sex_id", "sex_name", 
                     "age_id", "age_name", "cause_id", "cause_name", "metric_id", 
                     "metric_name", "year_id", "nm_mean", "nm_upper", "nm_lower")

all.daly <- all.daly[(all.daly$age_name=="Age-standardized"&all.daly$sex_name=="Both"&all.daly$year_id==2015),c("country","unaids.region","cause_name","sex_name","age_name","year_id","nm_mean", "nm_lower", "nm_upper")]
all.daly <- split(all.daly,list(all.daly$country,all.daly$year_id))
all.daly <- lapply(all.daly,function(x){rbind(x, c(names(table(x$country))
                                                   ,names(table(x$unaids.region))
                                                   ,NA
                                                   ,NA
                                                   ,NA
                                                   ,NA
                                                   ,colSums(x[7:9])))})

all.daly <- rbind.fill((all.daly))

all.daly[,c("cause_name","sex_name","age_name","year_id")] <- lapply(all.daly[,c("cause_name","sex_name","age_name","year_id")],as.character)
all.daly$cause_name[is.na(all.daly$cause_name)] <- "all cvd"
all.daly$sex_name[is.na(all.daly$sex_name)] <- "both"
all.daly$age_name[is.na(all.daly$age_name)] <- "age standardised"
all.daly$year_id[is.na(all.daly$year_id)] <- 2015
all.daly[,c("nm_mean","nm_upper","nm_lower")] <- lapply(all.daly[,c("nm_mean","nm_upper","nm_lower")],as.numeric)

all.daly <- all.daly[all.daly$cause_name=="all cvd",]

all.daly <- merge(x=all.daly,y=country.un,by.x="country",by.y="country",all.y=T)

all.daly$rt_mean <- all.daly$nm_mean*100000 
all.daly$rt_lower <- all.daly$nm_lower*100000 
all.daly$rt_upper <- all.daly$nm_upper*100000 

all.daly <- all.daly[,c("country","unaids.region.id", "rt_mean", "rt_lower", "rt_upper")]

# remove countries where no daly available and highlight countries where either prevalence or daly not available

all.daly[is.na(all.daly$rt_mean),"country"]


all.daly <- subset(all.daly, !is.na(all.daly$rt_mean))

x <- merge(x=all1,y=all.daly,by.x="country.code",by.y="unaids.region.id",all=T)
x <- x[,c("country.x","country.y")]

names(x) <- c("paf.avail","daly.avail")
x <- x[is.na(x$paf.avail)|is.na(x$daly.avail),]

# combine with data for paf

all.paf <- all
all.paf <- all.paf[c("country","country.code","HIV.prevalence..15..","g..HIV.population..15...Total.Region..Male.Female","g..HIV.population..15...Total.Region..Male.Female...Lower","g..HIV.population..15...Total.Region..Male.Female...Upper","g..Population.aged.15.64.Total.Region..Male.Female", "g..Population.aged.65..Total.Region..Male.Female")]

all.combined <- merge(x=all.daly,y=all.paf,by.x="unaids.region.id",by.y="country.code",all.x=T) 

all.combined <- all.combined[!is.na(all.combined$country.y),]

names(all.combined) <- c("unaids.region.id", "country.x", "rt_mean", "rt_lower", "rt_upper", 
                         "country.y", "hiv.prev.15+.", "hiv.pop.all.adults.ce","hiv.pop.all.adults.ll","hiv.pop.all.adults.ul", 
                         "pop.15-64", "pop.64+")
all.combined$pop.all <- all.combined$`pop.15-64`+all.combined$`pop.64+`


# prevalece, ce, ll and ul from the population

all.combined$prev.ce <- as.numeric(all.combined$hiv.pop.all.adults.ce/all.combined$pop.all)
all.combined$prev.ll <- as.numeric(all.combined$hiv.pop.all.adults.ll/all.combined$pop.all)
all.combined$prev.ul <- as.numeric(all.combined$hiv.pop.all.adults.ul/all.combined$pop.all)

# exclude montenegro - data anomaly

all.combined <- subset(all.combined, !all.combined$unaids.region.id=="MNE")
all.combined <- subset(all.combined, !all.combined$unaids.region.id=="OMN")

alpha <- function(x){
  alpha <- as.vector(beta.parms.from.quantiles(c(as.numeric(x[15]),as.numeric(x[16])))[1])
}

beta <- function(x){
  beta <- as.vector(beta.parms.from.quantiles(c(as.numeric(x[15]),as.numeric(x[16])))[2])
}


all.combined$alpha <- as.vector(unlist(apply(all.combined, 1,alpha)))

all.combined$beta <- as.vector(unlist(apply(all.combined, 1,beta)))

all.combined[,c("rt_mean","rt_upper","rt_lower")] <-lapply(all.combined[,c("rt_mean","rt_upper","rt_lower")],function(x){log(x,base=exp(1))}) 

all.combined$rt_se <- (as.numeric(all.combined$rt_upper)-as.numeric(all.combined$rt_lower))/3.92

for(i in 1:nrow(all.combined)){
  x <- rbeta(num_iter, shape1 = all.combined$alpha[i], shape2 = all.combined$beta[i])
  y <- rnorm(num_iter,as.numeric(all.combined$rt_mean)[i],all.combined$rt_se[i])
  y <- exp(y)
  p <- (x*(rr.sim-1))/(1+x*(rr.sim-1))
  burden <- p*y
  all.combined$burdence[i] <- quantile(burden,0.5,na.rm = T)
  all.combined$burdenll[i] <- quantile(burden,0.025,na.rm = T)
  all.combined$burdenul[i] <- quantile(burden,0.975,na.rm = T)
  all.combined$paf.ce[i] <- quantile(p,0.5,na.rm = T)*100
  all.combined$paf.ll[i] <- quantile(p,0.025,na.rm = T)*100
  all.combined$paf.ul[i] <- quantile(p,0.975,na.rm = T)*100
}

all.combined$paf.test <- (all.combined$prev.ce*(2.07-1))/(1+all.combined$prev.ce*(2.07-1))
plot(all.combined$paf.test,all.combined$paf.ce/100)

all.combined[,c("rt_mean","rt_upper","rt_lower")] <-lapply(all.combined[,c("rt_mean","rt_upper","rt_lower")],function(x){exp(x)}) 

all.combined$rate <- paste0(round(all.combined$rt_mean,0)," ","(",round(all.combined$rt_lower,0),"-",round(all.combined$rt_upper,0),")")
all.combined$rate <- paste0(round(all.combined$rt_mean,0)," ","(",round(all.combined$rt_lower,0),"-",round(all.combined$rt_upper,0),")")

all.combined$burdence <- format(round(all.combined$burdence,2),nsmall=2)
all.combined$burdenll <- format(round(all.combined$burdenll,2),nsmall=2)
all.combined$burdenul <- format(round(all.combined$burdenul,2),nsmall=2)

all.combined$burden <- paste0(format(all.combined$burdence,nsmall=2)," ","(",format(all.combined$burdenll,nsmall=2),"-",format(all.combined$burdenul,nsmall=2),")")

all.combined$burdence <- as.numeric(all.combined$burdence)

#plot for paf
paf <- joinCountryData2Map( all.combined
                            ,joinCode = "ISO3"
                            ,nameJoinColumn = "unaids.region.id")

classInt <- classIntervals( paf[["paf.ce"]]
                            ,n=7, style = "quantile")
catMethod = classInt[["brks"]]

colourPalette <- brewer.pal(7,'YlOrRd')

#plot map
mapDevice() #create world map shaped window
mapParams <- mapCountryData(paf
                            ,nameColumnToPlot="paf.ce"
                            ,addLegend=T
                            ,catMethod = catMethod
                            ,colourPalette=colourPalette)

classInt

#plot for daly
daly <- joinCountryData2Map( all.combined
                             ,joinCode = "ISO3"
                             ,nameJoinColumn = "unaids.region.id")

classInt <- classIntervals( daly[["burdence"]]
                            ,n=7, style = "quantile")
catMethod = classInt[["brks"]]

colourPalette <- brewer.pal(7,'YlOrRd')
#plot map
mapDevice() #create world map shaped window
mapParams <- mapCountryData(daly
                            ,nameColumnToPlot="burdence"
                            ,addLegend=F
                            ,catMethod = catMethod
                            ,colourPalette=colourPalette
                            ,missingCountryCol = 'white')

classInt

# comparing HIV to other risk factors
setwd("/Users/anoopshah/Documents/Anoop/Anoop Medical/HIV and CV Disease sys review/manuscript/hiv manuscript/final/bmj/r code for publication/compare to other RF")
f <- list.files()

dat <- lapply(f,function(i){
  x <- read.xls(i)
})


dat <- as.data.frame(rbindlist(dat))
dat <- dat[dat$Cause.of.death.or.injury=="Cardiovascular diseases"&!dat$Risk.factor=="",]
dat$Location <- as.character(dat$Location)
dat$Location <- ifelse(dat$Location=="Cote d'Ivoire","Cote dvoire",dat$Location)
dat <- dat[,c("Location","Risk.factor","Value","Lower.bound","Upper.bound") ]

dat[,c("Value","Lower.bound","Upper.bound")] <- lapply(dat[,c("Value","Lower.bound","Upper.bound")],function(x){round(x*100,2)})
dat$paf <- paste0(format(dat$Value,nsmall=2)," ","(",format(dat$Lower.bound,nsmall=2),"-",format(dat$Upper.bound,nsmall=2),")")

dat <- dat[,c("Location","Risk.factor","paf")]

dat_wf <- spread(dat,"Risk.factor","paf")  

dat2 <- daly.pc[daly.pc$year==2015,]
dat2 <- dat2[,c("country","pafce","pafll","paful")]
dat2[,c("Value","Lower.bound","Upper.bound")] <- lapply(dat2[,c("pafce","pafll","paful")],function(x){round(x*100,2)})
dat2$paf <- paste0(format(dat2$Value,nsmall=2)," ","(",format(dat2$Lower.bound,nsmall=2),"-",format(dat2$Upper.bound,nsmall=2),")")

dat2 <- dat2[,c("country","paf")]
names(dat2) <- c("Location", "Value", "Lower.bound", "Upper.bound","Risk.factor")
names(dat2) <- c("country","HIV")

dat3 <- merge(x=dat_wf,y=dat2,by.x="Location",by.y="country")

write.csv(dat3, "/Users/anoopshah/Documents/Anoop/Anoop Medical/HIV and CV Disease sys review/manuscript/hiv manuscript/final/circulation/response/supplment.table5.csv")


