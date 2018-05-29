#### This R code assumes that the working directory is set to same directory as that where the R code is saved ####
#### Version used is R version 3.4.4 (2018-03-15) ####
#### The analysis approach described here uses sampling methods from various distributions. As such the final results estimates may vary slightly from those published.####


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

