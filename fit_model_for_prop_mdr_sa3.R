### Code to fit smooth spline to WHO data for chosen countries

###********** Libraries and ggplot theme ************************************************************************************************#######
library('rstan')
library('dplyr')
library('ggmcmc')
library("boot")
library("reshape")
library("ggforce")
library("brms")
library("data.table")
theme_set(theme_bw(base_size = 24))

###********** Home ************************************************************************************************#######
home <- "~/Documents/LTBI_MDR/"
setwd(home)

###********** Load code and data ************************************************************************************************#######
# Country list and WHO data
## NOT SUPPLIED WITH REPOSITORY
who0 <- as.data.frame(read.csv("data/new_who_edited_sub.csv")[,-1]) 
#who0$year <- who0$year_new

# Add in sigma (standard deviation) to data
who0$mdr_new <- who0$new_mdr_prop
who0$sigma <- (who0$mhi - who0$mdr_new)/1.96  

## JUST doing for 3 countries now: India / China / USA
uu <- c("USA", "CHN","IND")
luu <- 3 

# Store
pp_store <- c() # only three countries, not so big

# Spline function
sample_spline <- function(lambdav){
  # output
  ans <- list()
  # 1,000 samples from the data: with normal distributed by SD given by mhi - mlo
  N <- 1e3
  # Fit to 1,000 sets of data
  for(k in 1:N){
    ss <- smooth.spline(x=DR$year,
                        y= rnorm(n=nrow(DR),
                                 mean=DR$mdr_new,
                                 sd = (DR$mhi-DR$mlo)/3.92),
                        lambda = lambdav)
    so <- predict(ss,x=1970:2018) # predict for 1970:2018
    ans[[k]] <- as.data.table(so)
  }
  ans <- rbindlist(ans)
  ans[,iter:=rep(1:N,each=length(1970:2018))]
  colnames(ans) <- c("year","value","iter")
  return(ans)
}

# Cycle through all countries
for(ii in 1:luu) {
  # country for this fit 
  cnn <- uu[ii]
  
  # data for this country with extra 0 in 1970
  who_cnn <- subset(who0, iso3 == cnn)[,c("year_new","mdr_new","mlo","mhi")]
  
  DR <- rbind(data.frame(year_new=1970,mhi=0,mlo=0,mdr_new=0),
              who_cnn[,c("year_new","mhi","mlo","mdr_new")])
  
  ## Lamda value (internal smooth parameter) varied for each country - fitted by eye 
  if(cnn == "USA"){lambdav = 5e-5}
  if(cnn == "CHN"){lambdav = 5e-5}
  if(cnn == "IND"){lambdav = 1e-3}
  
  ans <- sample_spline(lambdav) # run fits to 1,000 data points
  
  # plot
  ggplot(data=DR,aes(year_new,mdr_new)) + geom_point(col="red") +
    geom_errorbar(data=subset(DR,year_new>1970),aes(ymin=mlo,ymax=mhi),col="red") +
    geom_line(data=ans,aes(year,value,group=iter),alpha=.1) + ggtitle(cnn) 
  ggsave(paste0("output_spline/",cnn,"_splinefit.pdf"))
  

  # Remove those negative
  ans2 <- ans %>% group_by(iter) %>% filter(!any(value < 0))
  ans2$new_iter <- NA # renumber iterations
  ans2$new_iter[order(ans2$iter)] <- rep(1:(length(unique(ans2$iter))),each = 49)
  
  # just take first 200 fits
  ppcm <- subset(ans2,new_iter %in% seq(1,200,1)) 
  
  ggplot(data=DR,aes(year_new,mdr_new)) + geom_point(col="red") +
    geom_errorbar(data=subset(DR,year_new>1970),aes(ymin=mlo,ymax=mhi),col="red") +
    geom_line(data=ppcm,aes(year,value,group=iter),alpha=.1) + 
    scale_x_continuous("Year") + scale_y_continuous("Prop. new with MDR")
  ggsave(paste0("output_spline/",cnn,"chosen_fits.pdf"))
  
  # save 
  ppcm$iso3 <- cnn
  
  pp_store <- rbind(pp_store, ppcm)
}

write.csv(pp_store,"output_spline/all_samples_p.csv")

# relabel iterations (removed those that were anywhere negative)
pp_store$iter <- pp_store$new_iter

##### make correct for adding to other data
# and merge
load("../MDR-LTBI-paper/data_final/rundata_ari_1000.Rdata") # From running the code in Houben & Dodd in 2016.
rundatar$ari <- exp(rundatar$lari)

# Add in 0 before 1970
ns <- max(as.numeric(unique(pp_store$iter)))
sample_v <- seq(1,ns,1)
years <- unique(rundatar$year) 
years_v <- years[which(years < 1970)]
ny <- length(years_v)
nc <- length(unique(pp_store$iso3))
pre_1970_mdr <- as.data.frame(cbind(rep(sample_v,ny*nc),rep(years_v, each = ns*nc)))
colnames(pre_1970_mdr) <- c("replicate", "year")   
pre_1970_mdr$iso3 <-rep(uu, each = ns)
pre_1970_mdr$value <- 0

pp_store$replicate <- pp_store$iter
all_samples <- rbind(pp_store[,c("replicate","year","iso3","value")], pre_1970_mdr)

#ggplot(all_samples, aes(x=year,  y = value, group = replicate)) + geom_line() + facet_wrap(~iso3)

all_samples_1970_2014 <-all_samples
write.csv(all_samples_1970_2014, "output_spline/all_samples_p_1970_2014.csv") 

## MERGE DS AND MDR
all0n <- merge(all_samples, rundatar, by = c("year","iso3","replicate"))
# MDR-ARI
all0n$mdr_ari <- all0n$value * all0n$ari

# TOTAL AND MDR
all0n$ds_ari <- all0n$ari - all0n$mdr_ari # ds ARI is the remainder of the ari

save(all0n,file="output_spline/all0n_p_ds_mdr.Rdata") 


### PLOT
# READ IN
#load("output_spline/all0n_p_ds_mdr.Rdata")
theme_set(theme_bw(base_size = 24))
pp <- "smooth.spline"
for(cci in 1:length(uu)){
  print(uu[cci])
  ### WHO data
  d <-subset(who0, iso3 == as.character(uu[cci]) )
  
  ### ARI for both DS and mDR in all0n
  rdata <- all0n[which(all0n$iso3 == as.character(uu[cci])),]
  
  a1 <- ggplot(d, aes(x=year_new, y = new_mdr_prop),col="red",pch = 10) + # points won't plot over lines unless do points first?!
    geom_point() +
    geom_line(data = rdata, aes(x=year, y = value, group = factor(replicate)),alpha = 0.2) +
    scale_y_continuous("Prop. new with MDR") + scale_x_continuous("Year",lim=c(1970,2016)) +
    geom_point(data = d, aes(x=year_new, y = new_mdr_prop),col="red",pch = 10, size = 3) + geom_errorbar(data = d, aes(ymin = mlo, ymax = mhi), col = "red")
  ggsave(paste0("output_spline/",uu[cci],"_mdr_trends_with_data_",pp,".pdf"),width=11, height=11)
  
  a2 <- ggplot(rdata, aes(x=year, y = mdr_ari, group = factor(replicate))) + geom_line(alpha = 0.2) +
    scale_y_continuous("MDR ARI") + scale_x_continuous("Year",lim=c(1970,2015))
  ggsave(paste0("output_spline/",uu[cci],"_mdr_ari_",pp,".pdf"),width=11, height=11)
  
}  


### New all0
load("../MDR-LTBI-paper/data_final/all0_p_ds_mdr.Rdata") # From running the code in Houben & Dodd in 2016.
w<-which(all0$iso3 %in% c("CHN","IND","USA"))
all0_new <- all0[-w,]
colnames(all0n)[which(colnames(all0n) == "value")] <- "prediction"
all0_new <- rbind(all0_new, all0n)

## Plots to explore / compare old and new behaviour

ggplot(subset(all0, iso3 %in% c("CHN","IND","USA")), aes(x=year, y = mdr_ari, group = replicate)) + geom_line() + facet_wrap(~iso3, scales = 'free') + 
  geom_line(data = subset(all0_new, iso3 %in% c("CHN","IND","USA")), aes(x=year, y = mdr_ari, group = replicate), col = "red", alpha = 0.2)  + ggtitle("MDR-ARI")
ggsave("output_spline/new_mdr_ARI_vs_old.pdf")       

ggplot(subset(all0, iso3 %in% c("CHN","IND","USA")), aes(x=year, y = ds_ari, group = replicate)) + geom_line() + facet_wrap(~iso3, scales = "free") + 
  geom_line(data = subset(all0_new, iso3 %in% c("CHN","IND","USA")), aes(x=year, y = ds_ari, group = replicate), col = "red", alpha = 0.2) + ggtitle("DS-ARI")
ggsave("output_spline/new_ds_ARI_vs_old.pdf")

ggplot(subset(all0, iso3 %in% c("CHN","IND","USA")), aes(x=year, y = ari, group = replicate)) + geom_line() + facet_wrap(~iso3, scales = "free") + 
  geom_line(data = subset(all0_new, iso3 %in% c("CHN","IND","USA")), aes(x=year, y = ari, group = replicate), col = "red", alpha = 0.2) + ggtitle("ARI")
ggsave("output_spline/new_ARI_vs_old.pdf")

## Restrict time
ggplot(subset(subset(all0, year > 1980), iso3 %in% c("CHN","IND","USA")), aes(x=year, y = mdr_ari, group = replicate)) + geom_line() + facet_wrap(~iso3, scales = 'free') + 
  geom_line(data = subset(subset(all0_new, year > 1980), iso3 %in% c("CHN","IND","USA")), aes(x=year, y = mdr_ari, group = replicate), col = "red", alpha = 0.2)  + ggtitle("MDR-ARI")
ggsave("output_spline/new_mdr_ARI_vs_old_recent.pdf")       

ggplot(subset(subset(all0, year > 1980), iso3 %in% c("CHN","IND","USA")), aes(x=year, y = ds_ari, group = replicate)) + geom_line() + facet_wrap(~iso3, scales = "free") + 
  geom_line(data = subset(subset(all0_new, year > 1980), iso3 %in% c("CHN","IND","USA")), aes(x=year, y = ds_ari, group = replicate), col = "red", alpha = 0.2) + ggtitle("DS-ARI")
ggsave("output_spline/new_ds_ARI_vs_old_recent.pdf")

ggplot(subset(subset(all0, year > 1980), iso3 %in% c("CHN","IND","USA")), aes(x=year, y = ari, group = replicate)) + geom_line() + facet_wrap(~iso3, scales = "free") + 
  geom_line(data = subset(subset(all0_new, year > 1980), iso3 %in% c("CHN","IND","USA")), aes(x=year, y = ari, group = replicate), col = "red", alpha = 0.2) + ggtitle("ARI")
ggsave("output_spline/new_ARI_vs_old_recent.pdf")

# control... checking same. Yep! 
ggplot(subset(all0, iso3 %in% c("BRA","HKG","BWA")), aes(x=year, y = mdr_ari, group = replicate)) + geom_line() + facet_wrap(~iso3) +
  geom_line(data = subset(all0_new, iso3 %in% c("BRA","HKG","BWA")), aes(x=year, y = mdr_ari, group = replicate), col = "red")

## Save updated all0_new
save(all0_new,file="output_spline/all0new_p_ds_mdr.Rdata") 

