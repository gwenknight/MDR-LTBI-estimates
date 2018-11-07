### Code to fit y = bt - ct^2 to WHO data for all countries

###********** Libraries and ggplot theme ************************************************************************************************#######
library('rstan')
library('dplyr')
library('ggmcmc')
library("boot")
library("reshape")
library("ggforce")
theme_set(theme_bw(base_size = 24))

###********** Home ************************************************************************************************#######
home <- "~/Documents/LTBI_MDR/"
setwd(home)

###********** Load code and data ************************************************************************************************#######
# Country list and WHO data
## NOT SUPPLIED WITH REPOSITORY
#who0 <- as.data.frame(read.csv("data_final/new_who_edited_sub.csv")[,-1]) 
#who0$year <- who0$year_new

# Add in sigma (standard deviation) to data
who0$mdr_new <- who0$new_mdr_prop
who0$sigma <- (who0$mhi - who0$mdr_new)/1.96  

uu <- read.csv("data_final/138_final_list_included_countries.csv",stringsAsFactors = FALSE)[,-1]
luu <- length(unique(who0$iso3)) # 138

# Store
pred_samples_p <- list()

# Cycle through all countries
for(ii in 1:luu) {
  # country for this fit 
  country <- uu[ii]
  
  # data for this country
  who_l <- who0 %>% filter(iso3==country) 
  who_l <- who_l[,c("year", "mdr_new", "mlo", "mhi", "sigma")]
  
  years_to_predict <- seq(1970, 2018)
  print(c(ii,as.character(country),nrow(who_l)))
  
  # if only one data point - have to make some inputs to stan "as.array"s
  if(nrow(who_l) == 1){ 
    # model fit in stan
    m.p <- stan(file="model_trend.stan",
                data = list(N=nrow(who_l), 
                            N2=length(years_to_predict), 
                            q=as.array(who_l$mdr_new),
                            years_obs=as.array(who_l$year), 
                            sigma=as.array(who_l$sigma), 
                            years=years_to_predict),
                control = list(adapt_delta = 0.99),
                chains=2, iter=20000, warmup=10000, thin=10)
    
  } else {
    m.p <- stan(file="quadratic_priors_years.stan",
                data = list(N=nrow(who_l), 
                            N2=length(years_to_predict), 
                            q=who_l$mdr_new, 
                            years_obs=who_l$year, 
                            sigma=who_l$sigma, 
                            years=years_to_predict),
                control = list(adapt_delta = 0.99),
                chains=2, iter=20000, warmup=10000, thin=10)
  }
  
  ## With informative Priors
  posterior.p<-As.mcmc.list(m.p,pars=c("b", "t_m","rho"))
  ggmcmc(ggs(posterior.p), file=paste0("output_final/",country, "_inforprior_mcmc.pdf"))

  # 200 samples for predicted levels 
  samples.p.p <- rstan::extract(m.p, pars="p_pred", permuted = FALSE)[801:1000,1,] # pick first chain
  
  colnames(samples.p.p) <-NULL
  mcmc.samples.p.p <- data.frame(samples.p.p) %>%
    tbl_df %>%
    dplyr::mutate(sample=1:n()) %>%
    gather(col, prediction, starts_with("X")) %>%
    dplyr::mutate(year=as.integer(sub("^X", "", col)) + 1969) 

  pred_samples_p[[country]] <- mcmc.samples.p.p %>%
    mutate(country = country)
  
}

all_samples_p <- pred_samples_p[] %>% bind_rows

write.csv(pred_samples_p,"output_final/pred_samples_p.csv")
write.csv(all_samples_p,"output_final/all_samples_p.csv")

##### make correct for adding to other data
# and merge
load("data_final/rundata_ari_1000.Rdata") # From running the code in Houben & Dodd in 2016.
rundatar$ari <- exp(rundatar$lari)

# Rename
all_samples <- all_samples_p %>%
  mutate(replicate = sample, iso3 = country) %>%
  dplyr::select(-sample) %>% dplyr::select(-country)
# Set negative to 0
all_samples[which(all_samples$prediction < 0), "prediction"] <- 0

# Set greater than 1 to 1! Data on proportion (0-1) so cannot go above 1. 
all_samples[which(all_samples$prediction > 1), "prediction"] <- 1

# Add in 0 before 1970
ns <- max(all_samples_p$sample)
sample_v <- seq(1,ns,1)
years <- unique(rundatar$year) 
years_v <- years[which(years < 1970)]
ny <- length(years_v)
nc <- llu
pre_1970_mdr <- as.data.frame(cbind(rep(sample_v,ny*nc),rep(years_v, each = ns*nc)))
colnames(pre_1970_mdr) <- c("replicate", "year")   
pre_1970_mdr$iso3 <-rep(uu, each = ns)
pre_1970_mdr$prediction <- 0

all_samples <- rbind(all_samples[,c("replicate","year","iso3","prediction")], pre_1970_mdr)

all_samples_1970_2014 <-all_samples
write.csv(all_samples_1970_2014, "output_final/all_samples_p_1970_2014.csv") 

## MERGE DS AND MDR
all0 <- merge(all_samples, rundatar, by = c("year","iso3","replicate"))
# MDR-ARI
all0$mdr_ari <- all0$prediction * all0$ari

# TOTAL AND MDR
all0$ds_ari <- all0$ari - all0$mdr_ari # ds ARI is the remainder of the ari

save(all0,file="output_final/all0_p_ds_mdr.Rdata") 
