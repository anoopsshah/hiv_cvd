#### This R code assumes that the working directory is set to same directory as that where the R code is saved ####
#### Version used is R version 3.4.4 (2018-03-15) ####
#### The analysis approach described here uses sampling methods from various distributions. As such the final results estimates may vary slightly from those published.####

###=========================================================
### Part 2: estimating burden using daly and prevalence data
###=========================================================

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

hiv.prev$subgroup <- factor(hiv.prev$subgroup,
                            levels=c("Adults (15-49) estimate modelled", "Adults (15-49) lower estimate modelled", 
                                     "Adults (15-49) upper estimate modelled", "Females Adults (15-49) estimate", 
                                     "Females Adults (15-49) lower estimate", "Females Adults (15-49) upper estimate", 
                                     "Males Adults (15-49) estimate", "Males Adults (15-49) lower estimate", 
                                     "Males Adults (15-49) upper estimate"),
                            labels = c("adults.ce","adults.ll","adults.ul","female.ce","female.ll","female.ul","male.ce","male.ll","male.ul"))

hiv.prev <- hiv.prev %>% separate(subgroup,c("subgroup","subgroup2"))

hiv.prev <- hiv.prev[c("subgroup", "area", "area.id", "time.period", "data.value", "subgroup2")]

hiv.prev_wide <- spread(hiv.prev, subgroup2, data.value)
