#### Run cohort_ltbi_mdr for sensitivity analysis 3: just the countries with peak and crash
# Take in ARI trends, run cohort model
# Calculate proportion of MDR-LTBI acquired in different time points
# Code for some additional plots given

###********** Libraries and ggplot theme ************************************************************************************************#######
library(ggplot2)
theme_set(theme_bw(base_size = 24))
library(plyr)
library(dplyr)
library(cowplot)
library(data.table)
library(reshape2)

library(maps)
library(rworldmap)
library(RColorBrewer)


###********** Home ************************************************************************************************#######
home <- "~/Documents/LTBI_MDR/"
setwd(home)

###********** Load code and data ************************************************************************************************#######
source("../MDR-LTBI-paper/code_final/cohort_ltbi_mdr.R") # loads function for underyling cohort model

## Population size 2014
load('../MDR-LTBI-paper/data_final/POP2014.Rdata')  

## WHO data
w_data <- read.csv("data/new_who_edited_sub.csv")[,-1]

## Which countries? 
cni <- c("USA", "CHN","IND")

# READ IN
load("output_spline/all0new_p_ds_mdr.Rdata")
all0 <- all0_new

###********** Run for different countries ************************************************************************************************************************#######
# How many countries?
length(cni) # 3 now
llu <- length(cni)

cni_rem <- c() # blank to store what else to remove

# Number of runs
nari = 200

# Store all? 
store_all <- as.data.frame(matrix(0,length(cni)*4*81*100,10))
runn <- 1
level2014 <- c(); #breakdown proportions infected by age
s_level <- c(); #sum proportions infected 


# Run for all countries
for(cci in 1:llu){
  sa <- c() # store for this country
  print(c(cci,cni[cci]))
  
  ### WHO data
  d <-subset(w_data, iso3 == as.character(cni[cci]) )
  
  ### ARI for both DS and mDR in all0
  rdata <- all0[which(all0$iso3 == as.character(cni[cci])),]
  
  a1 <- ggplot(d, aes(x=year_new, y = new_mdr_prop),col="red",pch = 10) + # points won't plot over lines unless do points first?!
    geom_point() + 
    geom_line(data = rdata, aes(x=year, y = prediction, group = factor(replicate)),alpha = 0.2) + 
    scale_y_continuous("Prop. new with MDR") + scale_x_continuous("Year") +
    geom_point(data = d, aes(x=year_new, y = new_mdr_prop),col="red",pch = 10)
  save_plot(paste0("output_spline/",cni[cci],"_mdr_trends_with_data.pdf"), a1, base_aspect_ratio = 2 )
  
  a2 <- ggplot(rdata, aes(x=year, y = mdr_ari, group = factor(replicate))) + geom_line() + 
    scale_y_continuous("MDR ARI") + scale_x_continuous("Year") 
  save_plot(paste0("output_spline/",cni[cci],"_mdr_ari.pdf"), a2, base_aspect_ratio = 2 )
  
  for(i in 1:nari){
    print(c(i,"ari rep"))
    ari <- rdata[which(rdata$replicate == i),c("ds_ari","mdr_ari")]
    colnames(ari) <- c("ds","mdr")
    pop <- as.data.frame(POP2014[which(as.character(POP2014$iso3) == as.character(cni[cci])),"value"])
    
    cc <- cohort_ltbi(ari, pop)
    
    combs <- cc$combs
    
    # by age
    combs$mdr_rep <- i
    combs$age_group <- seq(1:17)
    combs$popf <- cci
    level2014 <- rbind(level2014,cbind(cc$c_2014,combs[1,c("mdr_rep","popf")],row.names = NULL))
    
    # total percentage infected sums
    ltbi_dr <- sum(combs$perc_dr) # percentage infected
    ltbi_ds <- sum(combs$perc_ds)
    pltbi_dr <- sum(combs$size_dr) # number of people infected
    pltbi_ds <- sum(combs$size_ds)
    
    ltbi_dr_kids <- sum(combs$perc_dr[1:3]) # percentage infected
    ltbi_ds_kids <- sum(combs$perc_ds[1:3])
    pltbi_dr_kids <- sum(combs$size_dr[1:3]) # number of people infected
    pltbi_ds_kids <- sum(combs$size_ds[1:3])
    
    # Bind together.
    s_level <- rbind(s_level,c(i,ltbi_dr,ltbi_ds,cci,pltbi_dr, pltbi_ds,ltbi_dr_kids,ltbi_ds_kids, pltbi_dr_kids, pltbi_ds_kids, sum(pop), sum(pop[1:3,1])))
    
    ssc <- cc$store_c
    lowi <- ((runn-1)*(dim(ssc)[1])+1)
    uppi <- ((runn)*(dim(ssc)[1]))
    store_all[lowi:uppi,1] <- i;
    store_all[lowi:uppi,2] <- cni[cci];
    store_all[lowi:uppi,3:10] <- ssc
    
    runn <- runn + 1
    sa <- rbind(sa,cbind(i,cni[cci], ssc)) # just for this country
  }
  
  # store all for this country
  sa <- as.data.frame(sa)
  colnames(sa) <- c(c("mdr_rep","cn"),colnames(cc$store_c)) 
  write.csv(sa, paste0("output_spline/",cni[cci],"_sa_",nari,".csv"))
  
  # Just recent infection 
  w<-which(sa$year > 2012)
  write.csv(sa[w,], paste0("output_spline/",cni[cci],"_rec_infec_",nari,".csv"))
  
}
# Totals
s_level <- as.data.frame(s_level)
colnames(s_level) <- c("rep","ltbir","ltbis","popf","pltbir", "pltbis","ltbir_kids","ltbis_kids",
                       "pltbir_kids", "pltbis_kids","pop","pop_kids")

dim(level2014) #nari * 107
level2014$cn <- cni[as.numeric(level2014$popf)]
# Store
write.csv(level2014, paste0("output_spline/level2014_",nari,".csv"))

s_level0 <- s_level
s_level$pop_name <- cni[as.numeric(s_level0$popf)]
s_level <- as.data.table(s_level)
# Store
write.csv(s_level, paste0("output_spline/s_level_",nari,".csv"))

### OR READ IN
#level2014 <- read.csv("output_spline/level2014_10_lin.csv",stringsAsFactors = FALSE)[,-1]
#s_level <- read.csv("output_spline/s_level_10_lin.csv", stringsAsFactors = FALSE)[,-1]

###********** PLOTS *****************************************************************#####
a2r<-ggplot(s_level, aes(x=pop_name, y = ltbir, col=factor(rep) )) + geom_point() + 
  guides(colour=FALSE) + 
  scale_x_discrete("Country") + scale_y_continuous("LTBI-MDR\n(% population infected)") +
  theme(axis.text.y = element_text(size = 6)) + 
  coord_flip() + aes(x=reorder(pop_name,ltbir),y=ltbir) 
save_plot(paste0("output_spline/ltbi_all_countries_r",nari,".pdf"), a2r, base_aspect_ratio = 1.5 )

a2s<-ggplot(s_level, aes(x=pop_name, y = ltbis, col=factor(rep) )) + geom_point() + 
  guides(colour=FALSE) + 
  scale_x_discrete("Country") + scale_y_continuous("LTBI-DS\n(% population infected)") +
  theme(axis.text.y = element_text(size = 6)) + 
  coord_flip() + aes(x=reorder(pop_name,ltbis),y=ltbis) 
save_plot(paste0("output_spline/ltbi_all_countries_s",nari,".pdf"), a2s, base_aspect_ratio = 1.5 )

s_level_mean <- s_level %>%
  group_by(pop_name) %>% 
  summarise_at(c("ltbir","ltbis","pltbir","pltbis"),funs(mean)) 

# functions for uncertainty interval
ub <- function(x)quantile(x,probs = .975)
lb <- function(x)quantile(x,probs = .025)

s_level_lo <- s_level %>%
  group_by(pop_name) %>% 
  summarise_at(c("ltbir","ltbis","pltbir","pltbis"),funs(lb))
colnames(s_level_lo) <- c("pop_name","ltbir.lo","ltbis.lo","pltbir.lo","pltbis.lo")

s_level_hi <- s_level %>%
  group_by(pop_name) %>% 
  summarise_at(c("ltbir","ltbis","pltbir","pltbis"),funs(ub))
colnames(s_level_hi) <- c("pop_name","ltbir.hi","ltbis.hi","pltbir.hi","pltbis.hi")

s_level_mean <- merge(s_level_mean, s_level_lo, by = "pop_name")
s_level_mean <- merge(s_level_mean, s_level_hi, by = "pop_name")

# plots
a2r<-ggplot(s_level_mean, aes(x=pop_name, y = ltbir )) + geom_point() + 
  scale_x_discrete("") + scale_y_continuous("") + #LTBI-MDR\n(% population infected)") +
  theme(axis.text.y = element_text(size = 6)) + 
  coord_flip() + aes(x=reorder(pop_name,ltbir),y=ltbir) +
  geom_errorbar(aes(min = ltbir.lo, max = ltbir.hi))
save_plot(paste0("output_spline/ltbi_all_countries_r_mean",nari,".pdf"), a2r, base_aspect_ratio = 0.8)

a2s<-ggplot(s_level_mean, aes(x=pop_name, y = ltbis )) + geom_point() +  
  scale_x_discrete("Country") + scale_y_continuous("LTBI-DS\n(% population infected)") +
  theme(axis.text.y = element_text(size = 6)) + 
  coord_flip() + aes(x=reorder(pop_name,ltbis),y=ltbis)  +
  geom_errorbar(aes(min = ltbis.lo, max = ltbis.hi))
save_plot(paste0("output_spline/ltbi_all_countries_s_mean",nari,".pdf"), a2s, base_aspect_ratio = 1.5 )

#### *** MAP **** ###

# define color buckets
cols = c("#F1EEF6", "#D4B9DA", "#C994C7", "#DF65B0", "#DD1C77", "#980043","grey")
s_level_mean$ltbis_c <- as.numeric(cut(s_level_mean$ltbis, c(0, 10, 20, 30, 40, 50,100)))
s_level_mean$ltbis_r <- as.numeric(cut(s_level_mean$ltbis, c(0, 1, 2, 3, 4, 5, 6)))

mapped_data <- joinCountryData2Map(s_level_mean, joinCode = "ISO3", nameJoinColumn = "pop_name")
par(mai=c(0,0,0.2,0),xaxs="i",yaxs="i")
cols <- colorRampPalette(brewer.pal(11,"Reds"), bias = 2)(13)
mapParams <- mapCountryData(mapped_data, nameColumnToPlot = "ltbis", catMethod = seq(0,50,10),
                            colourPalette = cols,
                            addLegend = FALSE)

pdf(paste0("output_spline/map_ltbis_",nari,".pdf"))
mapParams <- mapCountryData(mapped_data, nameColumnToPlot = "ltbis", catMethod = seq(0,50,10),
                            colourPalette = cols,
                            addLegend = FALSE)
do.call( addMapLegend, c( mapParams
                          , legendLabels="all"
                          , legendWidth=0.5
                          
))
dev.off()

pdf(paste0("output_spline/map_ltbir_",nari,".pdf"))
mapParams <- mapCountryData(mapped_data, nameColumnToPlot = "ltbir", catMethod = seq(0,3,0.25),
                            colourPalette = cols,
                            addLegend = FALSE)
do.call( addMapLegend, c( mapParams
                          , legendLabels="all"
                          , legendWidth=0.5
                          
))
dev.off()


################**######################################################################################################################################################################################################
################**######################################################################################################################################################################################################
################**######################################################################################################################################################################################################
###******* AGE *****################################################################################################################
################**######################################################################################################################################################################################################
################**######################################################################################################################################################################################################
################**######################################################################################################################################################################################################

ss_mean <- c()

# Run for all countries
for(cci in 1:llu){
  print(c(cci,cni[cci]))
  
  # Read in all data for this country
  sa <- read.csv(paste0("output_spline/",cni[cci],"_sa_",nari,".csv"))[,-1]
  
  # 2014 level
  level2014 <- read.csv(paste0("output_spline/level2014_",nari,".csv"))[,-1]
  
  sa_rec <- sa[which(sa$year > 1965),] # recent
  
  ####*********** When contributes most to latent burden? *********######################################################################################################################################################
  age_groups <- cbind(seq(1,85,5),seq(5,85,5))
  age_groups[17,2] <- 100 # last one 81 - 100 years old
  
  # This countries data for each age and year
  s <- sa
  # level at 2014 for this country
  l14 <- level2014[which(as.character(level2014$cn) == as.character(cni[cci])),]
  
  # where store?
  s_new <- c()
  
  for(k in 1:nari){# for each mdr_rep
    print(c("mdr_rep",k))
    # subset the yearly data
    s_k <- subset(s, mdr_rep == k)
    
    for(j in 1:100){# all ages
      # print(c("age",j))
      
      # subset the proportion infected by mdr_rep and age group
      l14_k_j <- subset(l14, mdr_rep == k)[j,] # proportion with dr in 2014 for mdr_rep = k and age = j
      
      s_temp <- c()
      # for each year get the data for those that are that age in 2014
      for(yr in 2014:1934){
        
        agen = j - (2014-yr) # age in that year
        
        if(agen > 0){ # need age > 0!
          s_temp <- rbind(s_temp,subset(s_k, year == yr)[agen,]) # remeber this is age in 2014
        }
      }
      
      ## Cumulative change in proportion with DR or DS
      ## OK for this to go negative... as suggests proportion is decreasing.
      s_temp$cumr_py <- s_temp$new_dr - s_temp$rei_rs + s_temp$rei_sr
      s_temp$cums_py <- s_temp$new_ds - s_temp$rei_sr + s_temp$rei_rs
      
      w<-which(s_temp$year == 1934) 
      if(length(w) > 0 ){s_temp[w,"cums_py"] <- 0} # remove initial infections prior to 1934}
      
      
      ### Gives the right proportion as in l14_k_j
      # s_temp <- as.data.frame(s_temp)
      # colwise(sum)(s_temp)[,"new_dr"] - colwise(sum)(s_temp)[,"rei_rs"] + colwise(sum)(s_temp)[,"rei_sr"]
      # tail(s_temp,1)[,"pr_ds"] + colwise(sum)(s_temp)[,"new_ds"] - colwise(sum)(s_temp)[,"rei_sr"] + colwise(sum)(s_temp)[,"rei_rs"]
      # tail(s_temp,1)[,"pr_ds"] + sum(s_temp$new_ds) - sum(s_temp$rei_sr) + sum(s_temp$rei_rs)
      # sum(s_temp$cumr_py)
      # sum(s_temp$cums_py)
      # l14_k_j
      
      ## proportion of amount in 2014 that is from this cumulative change
      s_propr <- matrix(0,81,1);   s_props <- matrix(0,81,1)
      if(l14_k_j$pr_dr > 0){
        w<-which(s_temp$cumr_py > 0) # only take those that add to the LTBI burden
        s_propr[w,] <- s_temp[w,"cumr_py"] / sum(s_temp[w,"cumr_py"])} # Divide by added LTBI not final total. #s_temp$cumr_py /l14_k_j$pr_dr}
      
      w<-which(s_temp$cums_py > 0) # only take those that add to the LTBI burden
      s_props[w,] <- s_temp[w,"cums_py"] / sum(s_temp[w,"cums_py"]) # Divide by added LTBI not final total. #s_temp$cumr_py /l14_k_j$pr_dr}
      
      s_npropr <- matrix(0,81,1);
      s_nprops <- matrix(0,81,1);
      for(kk in 1:length(s_props)){
        s_npropr[kk] <-  s_propr[kk]
        s_nprops[kk] <-  s_props[kk]
      }
      
      ## should be 1
      #sum(s_props)
      #sum(s_propr)
      
      s_new <- rbind(s_new,(cbind(
        s_npropr, s_nprops, # proportion of amount in 2014 that is from this cumulative change
        seq(2014,1934,-1),j,k # years, age, rep
      ))) 
    }
    
  }
  
  ## All data
  s_new <- as.data.frame(s_new)
  colnames(s_new)<-c("pr_r","pr_s","year","age","mdr_rep")
  
  ## s_new: each row has the age in 2014 the year from which some contribution may come 
  # and the size of the contribution
  write.csv(s_new,paste0("output_spline/s_all_",cni[cci],".csv"))
  
  ### s_new => has the proportion from each year at each age for each mdr_rep
  w<-intersect(which(s_new$age == 32),which(s_new$mdr_rep == 1))
  sum(s_new[w,"pr_r"]) # = 1
  sum(s_new[w,"pr_s"]) # = 1
  
  ## Want: of all LTBI, when was it gained? So need not just by age... but by total population.
  s_new$pr_ltbir <- 0
  s_new$pr_ltbis <- 0
  s_new$yearcat<-cut(s_new$year, seq(1929,2018,5))

  ## For each country get population distribution
  pop <-  POP2014[which(as.character(POP2014$iso3) == as.character(cni[cci])),"value"]
  
  ## For the population in 2014. Calculate the proportion of the total population in each yearly age group.
  pr_2014_age = pop / sum(pop) / 5 ## Divided by 5 to make per subunit (think works as equivalent to averaging proportions and multiplying by total)
  pr_2014_age[17] = pop[17] / sum(pop) / 20 ## Apart from last which is 20 yrs long
  
  ss_here <- s_new
  
  for(i in 1:100){
    w <- which(ss_here$age == i)
    m <- intersect(which(age_groups[,1] <= i), which(age_groups[,2] >= i))
    ss_here[w,"pr_ltbir"] = ss_here[w,"pr_r"] * as.numeric(pr_2014_age[m])
    ss_here[w,"pr_ltbis"] = ss_here[w,"pr_s"] * as.numeric(pr_2014_age[m])
  }
  
  ## Grouped by mdr_rep
  # w<-which(ss_here$mdr_rep == 5)
  # sum(ss_here[w,"pr_ltbis"]) # = 1
  # sum(ss_here[w,"pr_ltbir"]) # = 1
  
  #setwd(output)
  ## This says: by age, when were they infected. The proportion of their % infected that can be
  # # allocated to past times.
  # ggplot(ss_here[which(ss_here$mdr_rep < 5),], aes(age, pr_s, fill = factor(year))) +
  #   geom_bar(position = "fill", stat = "identity") +
  #   scale_y_continuous() + facet_wrap(~mdr_rep) + ggtitle(paste0("DS-TB, ",cni[cci])) +
  #   scale_fill_hue("clarity")
  # ggsave(paste0("output_spline/eg_DS_age_",cni[cci],"_ltbis_when.pdf"), height = 10, width = 10)
  # 
  # ggplot(ss_here[which(ss_here$mdr_rep < 5),], aes(age, pr_r, fill = factor(year))) +
  #   geom_bar(position = "fill", stat = "identity") +
  #   scale_y_continuous() + facet_wrap(~mdr_rep) + ggtitle(paste0("MDR-TB, ",cni[cci])) +
  #   scale_fill_hue("clarity")
  # ggsave(paste0("output_spline/eg_DR_age_",cni[cci],"_ltbir_when.pdf"), height = 10, width = 10)
  # 
  # 
  # ## This says: by mdr_rep, when is the time window that contributes most
  # ## Tried to highlight 1980 period... but not working
  # #w <- which(ss_here$yearcat == "(1989,1994]")
  # #ss_here$extra_label_fill <- 0
  # #ss_here[w,"extra_label_fill"] <- 1
  # #scale_colour_manual( values = c( "1"="black","0" = "white"), guide = FALSE )
  # 
  # ggplot(ss_here, aes(mdr_rep, pr_ltbis, fill = factor(yearcat))) +
  #   geom_bar(position = "fill", stat = "identity") +
  #   scale_y_continuous("Percentage of LTBI DS\nfrom this 5 year time interval") + 
  #   ggtitle("DS-TB") +
  #   scale_fill_hue("Year")
  # ggsave(paste0("output_spline/DS_",cni[cci],"_ltbis_when.pdf"), height = 10, width = 20)
  # 
  # ggplot(ss_here,aes(mdr_rep, pr_ltbir, fill = factor(yearcat))) +
  #   geom_bar(position = "fill", stat = "identity") +
  #   scale_y_continuous("Percentage of LTBI MDR\nfrom this 5 year time interval") + 
  #   ggtitle("MDR-TB") + 
  #   scale_fill_hue("Year")
  # ggsave(paste0("output_spline/DR_",cni[cci],"_ltbir_when.pdf"), height = 10, width = 20)
  # 
  # Average over yearcat
  meanv <- ss_here %>% group_by(yearcat,mdr_rep) %>%
    dplyr::summarise(sum_prltbir = sum(pr_ltbir), sum_prltbis = sum(pr_ltbis)) 
  
  write.csv(meanv, paste0("output_spline/",cni[cci],"_props_ltbi_when.csv"))
  
  # Average over all reps
  means <- ss_here %>% group_by(yearcat,mdr_rep) %>%
    dplyr::summarise(totals = sum(pr_ltbis)) %>%
    ungroup %>%
    group_by(yearcat) %>%
    dplyr::summarise(mean=mean(totals), min=quantile(totals, 0.025), max=quantile(totals, 0.975))
  
  meanr <- ss_here %>% group_by(yearcat,mdr_rep) %>%
    dplyr::summarise(totals = sum(pr_ltbir)) %>%
    ungroup %>%
    group_by(yearcat) %>%
    dplyr::summarise(mean=mean(totals), min=quantile(totals, 0.025), max=quantile(totals, 0.975))
  
  means$type <- 0
  meanr$type <- 1
  mean_both <- rbind(means,meanr)
  
  ggplot(mean_both, aes(x= yearcat, y= mean, fill = factor(type) )) + 
    geom_bar(stat = "identity", position = "dodge") + geom_errorbar(aes(ymin = min, ymax = max), position = "dodge") + 
    scale_fill_discrete("TB type",labels = c("DS","MDR")) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    scale_x_discrete("Year grouping") + scale_y_continuous("Contribution of this year\nto LTBI burden")
  ggsave(paste0("output_spline/DR_mean_",cni[cci],"_ltbir_when.pdf"), height = 10, width = 20)
  
  ss_mean <- rbind(ss_mean, cbind(mean_both,cni[cci],ii))
  
  write.csv(paste0("output_spline/",cni[cci],"_ss_mean_",nari,".csv"))[,-1]
  
}



