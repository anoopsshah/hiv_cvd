#### This R code assumes that the working directory is set to same directory as that where the R code is saved ####
#### Version used is R version 3.4.4 (2018-03-15) ####
#### The analysis approach described here uses sampling methods from various distributions. As such the final results estimates may vary slightly from those published.####

###===============================================####
### Part 1 -  estimate pooled rate and risk ratio ####
###===============================================####

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
res.si.data <- dat.rate[ dat.rate$outcome.2=="Incidence - MI",]
res.si <- rma.glmm(measure = "IRLN", xi = events, ti = patient.years, data = res.si.data,slab = author.year)
forest.mi <- forest(res.si, at=log(c(1,5, 10, 25,50,100,250,500)),atransf =exp,order=order(res.si.data$pub.year))
res.si <- predict(res.si,transf = exp,digits=2)

# meta estimate for mortality from cvd  (rate)
res.sort.cvd.data <- dat.rate[ dat.rate$outcome.2=="Mortality - CVD",]
res.sort.cvd <- rma.glmm(measure = "IRLN", xi = events, ti = patient.years, data = res.sort.cvd.data,slab = res.sort.cvd.data$author.year)
forest.res.sort.cvd <- forest(res.sort.cvd, at=log(c(1,5, 10, 25,50,100,250,500)),atransf =exp,order=order(res.sort.cvd.data$pub.year))
res.sort.cvd <-predict(res.sort.cvd,transf = exp,digits=2)

# meta estimate for mortality from stroke  (rate)
res.sort.stroke.data <- dat.rate[ dat.rate$outcome.2=="Mortality - Stroke",]
res.sort.stroke <- rma.glmm(measure = "IRLN", xi = events, ti = patient.years, data = res.sort.stroke.data,slab=res.sort.stroke.data$author.year)
forest.res.sort.stroke <- forest(res.sort.stroke, at=log(c(1,5, 10, 25,50,100,250,500)),atransf =exp,order=order(res.sort.stroke.data$pub.year))
res.sort.stroke <-predict(res.sort.stroke,transf = exp,digits=2)

# meta estimate for mortality from mi  (rate)
res.sort.mi.data <- dat.rate[ dat.rate$outcome.2=="Mortality - MI",]
res.sort.mi <- rma.glmm(measure = "IRLN", xi = events, ti = patient.years, data = res.sort.mi.data,slab=res.sort.mi.data$author.year)
forest.res.sort.mi <- forest(res.sort.mi, at=log(c(1,5, 10, 25,50,100,250,500)),atransf =exp,order=order(res.sort.mi.data$pub.year))
res.sort.mi <-predict(res.sort.mi,transf = exp,digits=2)

# summary rates
list <- c(res.cvd,res.si,res.stroke,res.sort.cvd,res.sort.mi,res.sort.stroke)

summary.rates <- as.matrix(do.call("rbind",list))
summary.rates <- matrix(summary.rates,ncol=6)
summary.rates <- as.data.frame(t(summary.rates))
names(summary.rates) <- c("pred","V2","ll","ul")
summary.rates <-summary.rates[,c("pred","ll","ul")] 
summary.rates$variable <- c("res.cvd","res.si","res.stroke","res.sort.cvd","res.sort.mi","res.sort.stroke")
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
dat.ratio$author <- paste0(substr(dat.ratio$authors..primary,1,regexpr(",",dat.ratio$authors..primary)-1)," ","et al")


# remove tripathi et al as paper provided two estimates with use of the same control group. have chosen the ART+ group to reflect current HIV population

dat.ratio <- dat.ratio[!dat.ratio$order.number==19,]
  
dat.ratio$lgrr <- logb(dat.ratio$ratio.ce)
dat.ratio$lgll <- logb(dat.ratio$ratio.ll)
dat.ratio$lgul <- logb(dat.ratio$ratio.ul)
dat.ratio$se <- (dat.ratio$lgul-dat.ratio$lgll)/(2*1.96)

rr <- rma(yi = dat.ratio$lgrr,sei = dat.ratio$se, method = "ML", slab=dat.ratio$author) 
rr.exp <- predict(rr, transf = exp, digits=2)

dat.ratio$weights <- weights(rr)

dat.ratio$rr <- paste(format(dat.ratio$ratio.ce,digits=3),"","[",format(dat.ratio$ratio.ll,digits=3),",","",format(dat.ratio$ratio.ul,digits=3),"]" )


# sampling from probability distributions to derive RR and uncertainty estimate with 95% CI
num_iter <- 10000
rr.sim <- exp(rnorm(num_iter, as.vector(rr$b), rr$se))
(quantile(rr.sim, probs = c(0.025,0.5, 0.975)))
hist((rr.sim))

# plot for ratio

#dat.ratio <- dat.ratio[order(dat.ratio$pub.year),]
#dat.ratio <- dat.ratio[order(dat.ratio$outcome2),]

#rownames(dat.ratio) <- seq(length=nrow(dat.ratio))
#dat.ratio$author <- format(dat.ratio$author,justify="left")

### decrease margins so the full space is used
par(mar=c(4,4,1,2), font=1)

forest(rr, xlim=c(-3, 4), at=log(c(0.75,1,2, 4, 8, 12,16,20)), atransf=exp,
      ilab=cbind(dat.ratio$pub.year,format(dat.ratio$weights,digits=2),dat.ratio$rr),
      ilab.xpos=c(-2,-1.5,-1.0),cex=0.75, ylim=c(-1, 34),
      order=rev(order(dat.ratio$outcome2,dat.ratio$pub.year)),
      rows=c(3:9,15:19,26:30),
      xlab="Risk Ratio", mlab="",
      col = "red", border="red")

### add text with Q-value, dfs, p-value, and I^2 statistic
text(-3, -1, pos=4, cex=0.70, bquote(paste("RE Model for All Studies (Q = ",
                                            .(formatC(rr$QE, digits=2, format="f")), ", df = ", .(rr$k - rr$p),
                                            "; ", I^2, " = ",.(formatC(rr$I2, digits=1, format="f")), "%"," [95% CI ",.(formatC(confint(rr)[[1]][[7]], digits=1, format="f")),
                                           " - ",
                                           .(formatC(confint(rr)[[1]][[11]], digits=1, format="f")),
                                           "])")))

### set font expansion factor (as in forest() above) and use bold italic
### font and save original settings in object 'op'
op <- par(cex=0.70, font=4)

### add text for the subgroups
text(-3, c(31,20,10), pos=4, c("Cardiovascular events",
                               "Myocardial infarction",
                               "Stroke"))

### switch to bold font
par(font=2)

### add column headings to the plot
text(-3,34, "Author(s)",  pos=4)
text(-2,34,"Year")
text(-1.5,34,"Weights")
text(-1,34, "Risk Ratio [95% CI]")

### set par back to the original settings
par(op)


### fit random-effects model in the three subgroups
res.c <- rma(yi = lgrr,sei = se, data = dat.ratio[dat.ratio$outcome2=="CVD",],method = "ML") 
res.m <- rma(yi = lgrr,sei = se, data = dat.ratio[dat.ratio$outcome2=="MI",],method = "ML") 
res.s <- rma(yi = lgrr,sei = se, data = dat.ratio[dat.ratio$outcome2=="Stroke",],method = "ML") 


### add summary polygons for the three subgroups
addpoly(res.c, row=23.5, cex=0.75, atransf=exp, mlab="",col = "blue", border="blue")
addpoly(res.m, row= 12.5, cex=0.75, atransf=exp, mlab="",col = "blue", border="blue")
addpoly(res.s, row= 1.5, cex=0.75, atransf=exp, mlab="",col = "blue", border="blue")

### add text with Q-value, dfs, p-value, and I^2 statistic for subgroups
text(-3, 24.5, pos=4, cex=0.70, bquote(paste("RE Model for All Studies (Q = ",
                                             .(formatC(res.c$QE, digits=2, format="f")), ", df = ", .(res.c$k - res.c$p),
                                             "; ", I^2, " = ",.(formatC(res.c$I2, digits=1, format="f")), "%"," [95% CI ",.(formatC(confint(res.c)[[1]][[7]], digits=1, format="f")),
                                             " - ",
                                             .(formatC(confint(res.c)[[1]][[11]], digits=1, format="f")),
                                             "])")))
res.c <- predict(res.c,transf = exp)
res.c <- paste(format(round(res.c[[1]],2),nsmall=2),"[",format(round(res.c[[3]],2),nsmall=2),",",format(round(res.c[[4]],2),nsmall=2),"]")

text(-1, 24.5, cex=0.70, res.c)
     
text(-3, 13.5, pos=4, cex=0.70, bquote(paste("RE Model for All Studies (Q = ",
                                             .(formatC(res.m$QE, digits=2, format="f")), ", df = ", .(res.m$k - res.m$p),
                                             "; ", I^2, " = ",.(formatC(res.m$I2, digits=1, format="f")), "%"," [95% CI ",.(formatC(confint(res.m)[[1]][[7]], digits=1, format="f")),
                                             " - ",
                                             .(formatC(confint(res.m)[[1]][[11]], digits=1, format="f")),
                                             "])")))

res.m <- predict(res.m,transf = exp)
res.m <- paste(format(round(res.m[[1]],2),nsmall=2),"[",format(round(res.m[[3]],2),nsmall=2),",",format(round(res.m[[4]],2),nsmall=2),"]")

text(-1, 13.5, cex=0.70, res.m)

text(-3, 1.5, pos=4, cex=0.70, bquote(paste("RE Model for All Studies (Q = ",
                                            .(formatC(res.s$QE, digits=2, format="f")), ", df = ", .(res.s$k - res.s$p),
                                            "; ", I^2, " = ",.(formatC(res.s$I2, digits=1, format="f")), "%"," [95% CI ",.(formatC(confint(res.s)[[1]][[7]], digits=1, format="f")),
                                            " - ",
                                            .(formatC(confint(res.s)[[1]][[11]], digits=1, format="f")),
                                            "])")))

res.s <- predict(res.s,transf = exp)
res.s <- paste(format(round(res.s[[1]],2),nsmall=2),"[",format(round(res.s[[3]],2),nsmall=2),",",format(round(res.s[[4]],2),nsmall=2),"]")

text(-1, 1.5, cex=0.70, res.s)


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

# by follow-uo
rr.long <- rma(yi = lgrr,sei = se, data = dat.ratio[dat.ratio$follow.up>=median(dat.ratio$follow.up,na.rm=T),],method = "ML")
predict(rr.long,transf = exp,digits=2)

rr.short <- rma(yi = lgrr,sei = se, data = dat.ratio[dat.ratio$follow.up<median(dat.ratio$follow.up,na.rm=T),],method = "ML")
predict(rr.short,transf = exp,digits=2)

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

funnel(rr, atransf = exp,xlab = "Risk Ratio",at=log(c(0.25,0.5, 1,2,2,4,8,16)),ylim = c(0,0.6))
  
funnel(trimfill(rr), atransf = exp,xlab = "Risk Ratio",at=log(c(0.25,0.5, 1,2,2,4,8,16)),ylim = c(0,0.6))

