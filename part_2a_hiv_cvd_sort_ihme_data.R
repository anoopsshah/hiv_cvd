#### This R code assumes that the working directory is set to same directory as that where the R code is saved ####
#### Version used is R version 3.4.4 (2018-03-15) ####
#### The analysis approach described here uses sampling methods from various distributions. As such the final results estimates may vary slightly from those published.####
#### This part will require part 1 to be run prior 

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



