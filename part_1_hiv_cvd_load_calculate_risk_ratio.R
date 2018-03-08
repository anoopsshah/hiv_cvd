#### This R code assumes that the working directory is set to same directory as that where the R code is saved ####

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

