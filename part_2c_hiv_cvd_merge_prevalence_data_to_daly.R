#### This R code assumes that the working directory is set to same directory as that where the R code is saved ####
#### Version used is R version 3.4.4 (2018-03-15) ####
#### The analysis approach described here uses sampling methods from various distributions. As such the final results estimates may vary slightly from those published.####
#### This part will require part 1 and 2a and 2b to be run prior 

###=========================================================
### Part 2: estimating burden using daly and prevalence data
###=========================================================


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