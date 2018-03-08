#### This R code assumes that the working directory is set to same directory as that where the R code is saved ####

#### The analysis approach described here uses sampling methods from various distributions. As such the final results estimates may vary slightly from those published.####

#### At present we do not have permission to share country level HIV prevalence data. As such the code below is used to reproduce our HIV attributable burden estimates at a regional and global level. ####

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

