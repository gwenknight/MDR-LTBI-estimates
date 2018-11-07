### Combine results to give output

###********** Libraries ************************************************************************************************#######
library(plyr)
library(xtable)
library(magrittr)
library(dplyr)
library(ggplot2)
library("ggforce")
library(reshape)

library(maps)
library(rworldmap)
library(RColorBrewer)

###********** Home ************************************************************************************************#######
home <- "~/Documents/LTBI_MDR/"
setwd(home)

###********** Load code and data ************************************************************************************************#######
load('data_final/whokey.Rdata') # WHOkey has global region and iso3 codes

# number of replicates 
nari <- 200

# functions for 95% uncertainty
ub <- function(x)quantile(x,probs = .975, na.rm = TRUE)
lb <- function(x)quantile(x,probs = .025, na.rm = TRUE)

# level at 2014
s_level <- read.csv(paste0("output_final/s_level_",nari,".csv"))[,-1]
s_level$iso3 <- s_level$pop_name

### ADDITIONAL CALCULATIONS
# Kid levels
s_level$perc_ds_kids <- 100*s_level$pltbis_kids / s_level$pltbis
s_level$perc_dr_kids <- 100*s_level$pltbir_kids / s_level$pltbir

s_level$perc_ds_kids_all <- 100*s_level$pltbis_kids / (s_level$pltbis_kids + s_level$pltbir_kids)
s_level$perc_dr_kids_all <- 100*s_level$pltbir_kids / (s_level$pltbis_kids + s_level$pltbir_kids)

# Percentage of LTBI that is MDR
s_level$perc_ltbi_mdr <- 100*s_level$pltbir/(s_level$pltbir + s_level$pltbis)

#### By country
med.ltbir <- aggregate(s_level[,c("ltbir","ltbis","pltbir","pltbis","ltbir_kids","ltbis_kids",
                                  "pltbir_kids","pltbis_kids",
                                  "perc_ds_kids","perc_dr_kids","perc_ltbi_mdr","perc_ds_kids_all","perc_dr_kids_all")], 
                       list(s_level$iso3), median, na.rm = TRUE)

ub.ltbir <- aggregate(s_level[,c("ltbir","ltbis","pltbir","pltbis","ltbir_kids","ltbis_kids","pltbir_kids","pltbis_kids",
                                 "perc_ds_kids","perc_dr_kids","perc_ltbi_mdr","perc_ds_kids_all","perc_dr_kids_all")], 
                      list(s_level$iso3), ub)

lb.ltbir <- aggregate(s_level[,c("ltbir","ltbis","pltbir","pltbis","ltbir_kids","ltbis_kids","pltbir_kids","pltbis_kids",
                                 "perc_ds_kids","perc_dr_kids","perc_ltbi_mdr","perc_ds_kids_all","perc_dr_kids_all")], 
                      list(s_level$iso3), lb)

med.ltbir$iso3 <- med.ltbir$Group.1
ub.ltbir$iso3 <- ub.ltbir$Group.1
lb.ltbir$iso3 <- lb.ltbir$Group.1

#### By global region
# Add in global regions
ltbi <- merge(s_level,WHOkey[,c('iso3','g_whoregion')],by='iso3',all.x=TRUE)

# Sum over countries
ltbi_global <- aggregate(ltbi[,c("pltbir","pltbis","pltbir_kids","pltbis_kids","pop","pop_kids")], 
                         list(ltbi$g_whoregion,ltbi$rep), sum)

ltbi_global$perc_ds_kids <- 100*ltbi_global$pltbis_kids / ltbi_global$pltbis
ltbi_global$perc_dr_kids <- 100*ltbi_global$pltbir_kids / ltbi_global$pltbir

ltbi_global$perc_ds_kids_all <- 100*ltbi_global$pltbis_kids / (ltbi_global$pltbis_kids + ltbi_global$pltbir_kids)
ltbi_global$perc_dr_kids_all <- 100*ltbi_global$pltbir_kids / (ltbi_global$pltbis_kids + ltbi_global$pltbir_kids)

ltbi_global$perc_ltbi_mdr <- 100*ltbi_global$pltbir / (ltbi_global$pltbir + ltbi_global$pltbis)

### Risk ratio
ltbi_global$rr <- (ltbi_global$pltbir_kids / (ltbi_global$pltbir_kids + ltbi_global$pltbis_kids) ) / (ltbi_global$pltbir / (ltbi_global$pltbir + ltbi_global$pltbis))

rr.med.g <- aggregate(ltbi_global[,"rr"],list(ltbi_global$Group.1),median)
rr.ub.g <- aggregate(ltbi_global[,"rr"],list(ltbi_global$Group.1),ub)
rr.lb.g <- aggregate(ltbi_global[,"rr"],list(ltbi_global$Group.1),lb)

# Total global level
ltbi_total <- aggregate(ltbi[,c("pltbir","pltbis","pltbir_kids","pltbis_kids",
                                "pop","pop_kids")], 
                        list(ltbi$rep), sum) # still got replicates in there
ltbi_total$perc_ds_kids <- 100*ltbi_total$pltbis_kids / ltbi_total$pltbis
ltbi_total$perc_dr_kids <- 100*ltbi_total$pltbir_kids / ltbi_total$pltbir
ltbi_total$perc_ds_kids_all <- 100*ltbi_total$pltbis_kids / (ltbi_total$pltbis_kids + ltbi_total$pltbir_kids)
ltbi_total$perc_dr_kids_all <- 100*ltbi_total$pltbir_kids / (ltbi_total$pltbis_kids + ltbi_total$pltbir_kids)

ltbi_total$perc_ltbi_mdr <- 100*ltbi_total$pltbir / (ltbi_total$pltbir + ltbi_total$pltbis)

med.ltbir.g <- aggregate(ltbi_global[,c("pltbir","pltbis","pltbir_kids","pltbis_kids",
                                        "pop","pop_kids","perc_ds_kids","perc_dr_kids","perc_ltbi_mdr","perc_ds_kids_all","perc_dr_kids_all")],
                         list(ltbi_global$Group.1),median)

ub.ltbir.g <- aggregate(ltbi_global[,c("pltbir","pltbis","pltbir_kids","pltbis_kids",
                                       "pop","pop_kids","perc_ds_kids","perc_dr_kids","perc_ltbi_mdr","perc_ds_kids_all","perc_dr_kids_all")],
                        list(ltbi_global$Group.1),ub)

lb.ltbir.g <- aggregate(ltbi_global[,c("pltbir","pltbis","pltbir_kids","pltbis_kids",
                                       "pop","pop_kids","perc_ds_kids","perc_dr_kids","perc_ltbi_mdr","perc_ds_kids_all","perc_dr_kids_all")],
                        list(ltbi_global$Group.1),lb)

med.total <- colwise(median)(ltbi_total[,c("pltbir","pltbis","pltbir_kids","pltbis_kids","pop","pop_kids","perc_ds_kids","perc_dr_kids","perc_ltbi_mdr","perc_ds_kids_all","perc_dr_kids_all")])
lb.total <- colwise(lb)(ltbi_total[,c("pltbir","pltbis","pltbir_kids","pltbis_kids","pop","pop_kids","perc_ds_kids","perc_dr_kids","perc_ltbi_mdr","perc_ds_kids_all","perc_dr_kids_all")])
ub.total <- colwise(ub)(ltbi_total[,c("pltbir","pltbis","pltbir_kids","pltbis_kids","pop","pop_kids","perc_ds_kids","perc_dr_kids","perc_ltbi_mdr","perc_ds_kids_all","perc_dr_kids_all")])

### Risk ratio: global
ltbi_total$rr <- (ltbi_total$pltbir_kids / (ltbi_total$pltbir_kids + ltbi_total$pltbis_kids) ) / (ltbi_total$pltbir / (ltbi_total$pltbir + ltbi_total$pltbis))

paste0(round(median(ltbi_total$rr),2), " [",round(lb(ltbi_total$rr),2),"-",round(ub(ltbi_total$rr),2),"]")



rr.med.g <- aggregate(ltbi_total[,"rr"],median)
rr.ub.g <- aggregate(ltbi_total[,"rr"],list(ltbi_total$Group.1),ub)
rr.lb.g <- aggregate(ltbi_total[,"rr"],list(ltbi_total$Group.1),lb)


## Global level
med.ltbir.g$ltbir <- 100*med.ltbir.g$pltbir / med.ltbir.g$pop
ub.ltbir.g$ltbir <- 100*ub.ltbir.g$pltbir / ub.ltbir.g$pop
lb.ltbir.g$ltbir <- 100*lb.ltbir.g$pltbir / lb.ltbir.g$pop

med.ltbir.g$ltbis <- 100*med.ltbir.g$pltbis / med.ltbir.g$pop
ub.ltbir.g$ltbis <- 100*ub.ltbir.g$pltbis / ub.ltbir.g$pop
lb.ltbir.g$ltbis <- 100*lb.ltbir.g$pltbis / lb.ltbir.g$pop

med.ltbir.g$ltbir_kids <- 100*med.ltbir.g$pltbir_kids / med.ltbir.g$pop_kids
ub.ltbir.g$ltbir_kids <- 100*ub.ltbir.g$pltbir_kids / ub.ltbir.g$pop_kids
lb.ltbir.g$ltbir_kids <- 100*lb.ltbir.g$pltbir_kids / lb.ltbir.g$pop_kids

med.ltbir.g$ltbis_kids <- 100*med.ltbir.g$pltbis_kids / med.ltbir.g$pop_kids
ub.ltbir.g$ltbis_kids <- 100*ub.ltbir.g$pltbis_kids / ub.ltbir.g$pop_kids
lb.ltbir.g$ltbis_kids <- 100*lb.ltbir.g$pltbis_kids / lb.ltbir.g$pop_kids

## Total levels
med.total$ltbir <- 100*med.total$pltbir / med.total$pop
ub.total$ltbir <- 100*ub.total$pltbir / ub.total$pop
lb.total$ltbir <- 100*lb.total$pltbir / lb.total$pop

med.total$ltbis <- 100*med.total$pltbis / med.total$pop
ub.total$ltbis <- 100*ub.total$pltbis / ub.total$pop
lb.total$ltbis <- 100*lb.total$pltbis / lb.total$pop

med.total$ltbir_kids <- 100*med.total$pltbir_kids / med.total$pop_kids
ub.total$ltbir_kids <- 100*ub.total$pltbir_kids / ub.total$pop_kids
lb.total$ltbir_kids <- 100*lb.total$pltbir_kids / lb.total$pop_kids

med.total$ltbis_kids <- 100*med.total$pltbis_kids /med.total$pop_kids
ub.total$ltbis_kids <- 100*ub.total$pltbis_kids / ub.total$pop_kids
lb.total$ltbis_kids <- 100*lb.total$pltbis_kids / lb.total$pop_kids

#### Output tables
med.ltbir <- merge(med.ltbir,WHOkey[,c('iso3','g_whoregion')],by='iso3',all.x=TRUE) # add global region to country output
## All countries 
table1_countries <- as.data.frame(cbind(as.character(med.ltbir$iso3),
                                        as.character(med.ltbir$g_whoregion),
                                        paste0(sprintf('%.1f',med.ltbir$ltbis), " [", 
                                               sprintf('%.1f',lb.ltbir$ltbis), "-", 
                                               sprintf('%.1f',ub.ltbir$ltbis),"]"),
                                        paste0(sprintf('%.2f',med.ltbir$ltbir), " [", 
                                               sprintf('%.2f',lb.ltbir$ltbir), "-", 
                                               sprintf('%.2f',ub.ltbir$ltbir),"]"),
                                        paste0(sprintf('%.1f',med.ltbir$perc_ltbi_mdr), " [", 
                                               sprintf('%.1f',lb.ltbir$perc_ltbi_mdr), "-", 
                                               sprintf('%.1f',ub.ltbir$perc_ltbi_mdr),"]"),
                                        paste0(sprintf('%.1f',med.ltbir$perc_dr_kids_all), " [", 
                                               sprintf('%.1f',lb.ltbir$perc_dr_kids_all), "-", 
                                               sprintf('%.1f',ub.ltbir$perc_dr_kids_all),"]")
))

colnames(table1_countries) <- c("iso3","WHO region", "Perc. with LTBIS", "Perc. with LTBIR", 
                                "Perc. of LTBI that is MDR","Perc. of LTBI in <15yos that is MDR")

## All countries: NUMBER
table1_countries_numb <- as.data.frame(cbind( as.character(med.ltbir$iso3),
                                              as.character(med.ltbir$g_whoregion),
                                              paste0(sprintf('%.1f',1000 * med.ltbir$pltbis), " [", 
                                                     sprintf('%.1f',1000 * lb.ltbir$pltbis), "-", 
                                                     sprintf('%.1f',1000 * ub.ltbir$pltbis),"]"),
                                              paste0(sprintf('%.1f',1000 * med.ltbir$pltbir), " [", 
                                                     sprintf('%.1f',1000 * lb.ltbir$pltbir), "-", 
                                                     sprintf('%.1f',1000 * ub.ltbir$pltbir),"]"),
                                              paste0(sprintf('%.1f',1000 * (med.ltbir$pltbir+med.ltbir$pltbis)), " [", 
                                                     sprintf('%.1f',1000 * (lb.ltbir$pltbir+lb.ltbir$pltbis)), "-", 
                                                     sprintf('%.1f',1000 * (ub.ltbir$pltbir+ub.ltbir$pltbis)),"]")
))

colnames(table1_countries_numb) <- c("iso3","WHO region", "Number with LTBIS","Number with LTBIR", "Total with LTBI")


## Global
table1_global <- as.data.frame(cbind( as.character(med.ltbir.g$Group.1),
                                      paste0(sprintf('%.1f',med.ltbir.g$ltbis), " [", 
                                             sprintf('%.1f',lb.ltbir.g$ltbis), "-", 
                                             sprintf('%.1f',ub.ltbir.g$ltbis),"]"),
                                      paste0(sprintf('%.2f',med.ltbir.g$ltbir), " [", 
                                             sprintf('%.2f',lb.ltbir.g$ltbir), "-", 
                                             sprintf('%.2f',ub.ltbir.g$ltbir),"]"),
                                      paste0(sprintf('%.1f',med.ltbir.g$perc_ltbi_mdr), " [", 
                                             sprintf('%.1f',lb.ltbir.g$perc_ltbi_mdr), "-", 
                                             sprintf('%.1f',ub.ltbir.g$perc_ltbi_mdr),"]"),
                                      paste0(sprintf('%.1f',med.ltbir.g$perc_dr_kids_all), " [", 
                                             sprintf('%.1f',lb.ltbir.g$perc_dr_kids_all), "-", 
                                             sprintf('%.1f',ub.ltbir.g$perc_dr_kids_all),"]")
))
## Total
total <- as.data.frame(cbind("GLOBAL",
                             paste0(sprintf('%.1f',med.total$ltbis), " [", 
                                    sprintf('%.1f',lb.total$ltbis), "-", 
                                    sprintf('%.1f',ub.total$ltbis),"]"),
                             paste0(sprintf('%.2f',med.total$ltbir), " [", 
                                    sprintf('%.2f',lb.total$ltbir), "-", 
                                    sprintf('%.2f',ub.total$ltbir),"]"),
                             paste0(sprintf('%.1f',med.total$perc_ltbi_mdr), " [", 
                                    sprintf('%.1f',lb.total$perc_ltbi_mdr), "-", 
                                    sprintf('%.1f',ub.total$perc_ltbi_mdr),"]"),
                             paste0(sprintf('%.1f',med.total$perc_dr_kids_all), " [", 
                                    sprintf('%.1f',lb.total$perc_dr_kids_all), "-", 
                                    sprintf('%.1f',ub.total$perc_dr_kids_all),"]")
))

table1_global <- rbind(table1_global, total)


colnames(table1_global) <- c("WHO Region","Perc. with LTBIS", "Perc. with LTBIR", "Perc. of LTBI that is MDR", "Perc. of LTBI in <15yos that is MDR")
table1_global <- table1_global[c(1,2,5,3,6,4,7),]

##Global: NUMBER
table1_global_numb <- as.data.frame(cbind( as.character(med.ltbir.g$Group.1),
                                           paste0(prettyNum(signif(med.ltbir.g$pltbis,3),big.mark=","), " [", 
                                                  prettyNum(signif(lb.ltbir.g$pltbis,3),big.mark=","), "-", 
                                                  prettyNum(signif(ub.ltbir.g$pltbis,3),big.mark=","),"]"),
                                           paste0(prettyNum(signif(med.ltbir.g$pltbir,3),big.mark=","), " [", 
                                                  prettyNum(signif(lb.ltbir.g$pltbir,3),big.mark=","), "-", 
                                                  prettyNum(signif(ub.ltbir.g$pltbir,3),big.mark=","),"]")
))
## Total
total_numb <- as.data.frame(cbind("GLOBAL",
                                  paste0(prettyNum(signif(med.total$pltbis,3),big.mark=","), " [", 
                                         prettyNum(signif(lb.total$pltbis,3),big.mark=","), "-", 
                                         prettyNum(signif(ub.total$pltbis,3),big.mark=","),"]"),
                                  paste0(prettyNum(signif(med.total$pltbir,3),big.mark=","), " [", 
                                         prettyNum(signif(lb.total$pltbir,3),big.mark=","), "-", 
                                         prettyNum(signif(ub.total$pltbir,3),big.mark=","),"]")
))

table1_global_numb <- rbind(table1_global_numb, total_numb)

# percentage in each
med.ltbir.g$pop_r_each <- med.ltbir.g$pltbir / med.total$pltbir
med.ltbir.g$pop_each <- med.ltbir.g$pop / med.total$pop


colnames(table1_global_numb) <- c("WHO Region","# with LTBIS", "# with LTBIR")
table1_global_numb <- table1_global_numb[c(1,2,5,3,6,4,7),]

## Store
write.csv(table1_countries, paste0("output_final/table1_countries.csv"))
write.csv(table1_countries_numb, paste0("output_final/table1_countries_number.csv"))
write.csv(table1_global, paste0("output_final/table1_global.csv"))
write.csv(table1_global_numb, paste0("output_final/table1_global_numb.csv"))

## Which country biggest burden? 
ww <- which.max(med.ltbir$perc_ltbi_mdr)
med.ltbir[ww,]

#### *** MAP **** ###
### Map of kids ltbi prevalence
# define color buckets
# if run above this will be the best one
med.ltbir$ltbir_kids_c <- as.numeric(cut(med.ltbir$ltbir_kids, c(0, 10, 20, 30, 40, 50, 100)))
med.ltbir$ltbis_kids_c <- as.numeric(cut(med.ltbir$ltbis_kids, c(0, 10, 20, 30, 40, 50, 60)))
cols <- colorRampPalette(brewer.pal(11,"Reds"), bias = 2)(13)

mapped_data <- joinCountryData2Map(med.ltbir, joinCode = "ISO3", nameJoinColumn = "Group.1")

pdf(paste0("output_final/map_ltbis_kids_",nari,"_infor_prior.pdf"))
mapParams <- mapCountryData(mapped_data, nameColumnToPlot = "ltbis_kids", catMethod = seq(0,10,0.25),
                            colourPalette = cols,
                            addLegend = FALSE,missingCountryCol = gray(.8))
do.call( addMapLegend, c( mapParams
                          , legendLabels="all"
                          , legendWidth=0.5))
dev.off()

pdf(paste0("output_final/map_ltbir_kids_",nari,"_infor_prior.pdf"))
mapParams <- mapCountryData(mapped_data, nameColumnToPlot = "ltbir_kids", catMethod = seq(0,10,0.25),
                            colourPalette = cols,
                            addLegend = FALSE,missingCountryCol = gray(.8))
do.call( addMapLegend, c( mapParams
                          , legendLabels="all"
                          , legendWidth=0.5))
dev.off()

pdf(paste0("output_final/map_ltbir_",nari,"_infor_prior.pdf"))
mapParams <- mapCountryData(mapped_data, nameColumnToPlot = "ltbir", catMethod = seq(0,4,0.25),
                            colourPalette = cols,
                            addLegend = FALSE,missingCountryCol = gray(.8))
do.call( addMapLegend, c( mapParams
                          , legendLabels="all"
                          , legendWidth=0.5))
dev.off()


#### *** LEVEL in 2035 / 2050 **** ###
load('data_final/POP2035.Rdata')
load('data_final/POP2050.Rdata')

cni <- read.csv("data_final/138_final_list_included_countries.csv", stringsAsFactors = FALSE)[,-1]
# store
datam <- as.data.frame(matrix(0,138*20000, 11))

### Gather data
for(i in 1:138){
  cnn <- cni[i]
  d <- read.csv(paste0("output_final/",cnn,"level2014_200.csv"))[,-1]
  d$age <- seq(1,100,1)
  d$iso3 <- cnn
  datam[(1 + (i-1)*20000):(i*20000),] <- d
}
colnames(datam) <- colnames(d)

# 2035 population
N2035 <- sum(POP2035[,"value"])
POP2035$acatn <- 0:16
POP2035$acat <- POP2035$acatn-4 # 20 yrs earlier these were acat ... 
POP2035 <- subset(POP2035, POP2035$acat>=0)
POP2035$pop35 = POP2035$value
POP2035 <- POP2035[,c("iso3","acat","pop35")]

# 2050 population 
N2050 <- sum(POP2050[,"value"])
POP2050$acatn <- 0:16
POP2050$acat <- POP2050$acatn-7 # 35 year earlier 
POP2050 <- subset(POP2050, POP2050$acat>=0)
POP2050$pop50 = POP2050$value
POP2050 <- POP2050[,c("iso3","acat","pop50")]

# Group by age categories
datam$acat <- c(rep(seq(0,15,1), each = 5), rep(16,20))
new_data <- datam[,c("pr_dr","pr_ds","mdr_rep","acat","popf","iso3")] %>% 
  group_by(mdr_rep,acat,popf,iso3) %>% 
  summarise_at(c("pr_dr","pr_ds"),funs(mean)) 
dim(new_data) # 138*200*17 = 469200

new_data2 <- merge(new_data,WHOkey[,c('iso3','g_whoregion')],by='iso3',all.x=TRUE)

## 2035
# need iso, replicate, act, prev, g_whoregion 
fruns <- merge(new_data2[,c("iso3","mdr_rep","acat","pr_dr","pr_ds","g_whoregion")],
               POP2035,by=c('iso3','acat'),all=TRUE)
fruns <- fruns[!is.na(fruns$mdr_rep),]     # ditch those countries not in the 138
fruns <- fruns[!is.na(fruns$pop35),]       # ditch the dead - only acat to 12 included
dim(fruns) ## 200*138*13 = 358800

fruns$pLTBIR <- fruns$pop35 * fruns$pr_dr 
fruns$pLTBIS <- fruns$pop35 * fruns$pr_ds


# each number in POP2035 is the actual number divided by 1,000
fruns.total.35 <- aggregate(fruns[,c("pLTBIR","pLTBIS")], list(fruns$mdr_rep), sum)
med.fruns.total.35 <- colwise(median)(fruns.total.35) 
ub.fruns.total.35 <- colwise(ub)(fruns.total.35) 
lb.fruns.total.35 <- colwise(lb)(fruns.total.35) 

# each number in POP2035 is the actual number divided by 1,000
print(c(paste0(signif(1e3*med.fruns.total.35$pLTBIR,2)," [",
               signif(1e3*lb.fruns.total.35$pLTBIR,2),", ",
               signif(1e3*ub.fruns.total.35$pLTBIR,2),"]")))
print(c(paste0(signif(1e3*med.fruns.total.35$pLTBIS,2)," [",
               signif(1e3*lb.fruns.total.35$pLTBIS,2),", ",
               signif(1e3*ub.fruns.total.35$pLTBIS,2),"]")))


med.rate.35 <- (0.15*0.01) * med.fruns.total.35/N2035 * 1e5 # 0.15% reactivation x total  / total population x 100,000
ub.rate.35 <- (0.15*0.01) * ub.fruns.total.35/N2035 * 1e5 
lb.rate.35 <- (0.15*0.01) * lb.fruns.total.35/N2035 * 1e5 

print(c(paste0(sprintf('%.2f',med.rate.35$pLTBIR)," [",
               sprintf('%.2f',lb.rate.35$pLTBIR),", ",
               sprintf('%.2f',ub.rate.35$pLTBIR),"]")))
print(c(paste0(sprintf('%.1f',med.rate.35$pLTBIS)," [",
               sprintf('%.1f',lb.rate.35$pLTBIS),", ",
               sprintf('%.1f',ub.rate.35$pLTBIS),"]")))

100*med.fruns.total.35/N2035

100*med.fruns.total.35$pLTBIR/(med.fruns.total.35$pLTBIR+med.fruns.total.35$pLTBIS)

## 2050
# need iso, replicate, act, prev, g_whoregion 
fruns <- merge(new_data2[,c("iso3","mdr_rep","acat","pr_dr","pr_ds","g_whoregion")],
               POP2050,by=c('iso3','acat'),all=TRUE)
fruns <- fruns[!is.na(fruns$mdr_rep),]     # ditch those countries not in the 138
fruns <- fruns[!is.na(fruns$pop50),]       # ditch the dead - only acat to 12 included
dim(fruns) ## 200*138*13 = 358800

fruns$pLTBIR <- fruns$pop50 * fruns$pr_dr 
fruns$pLTBIS <- fruns$pop50 * fruns$pr_ds

# each number in POP2050 is the actual number divided by 1,000
fruns.total.50 <- aggregate(fruns[,c("pLTBIR","pLTBIS")], 
                            list(fruns$mdr_rep), sum)
med.fruns.total.50 <- colwise(median)(fruns.total.50) 
ub.fruns.total.50 <- colwise(ub)(fruns.total.50) 
lb.fruns.total.50 <- colwise(lb)(fruns.total.50) 

# each number in POP2050 is the actual number divided by 1,000
print(c(paste0(signif(1e3*med.fruns.total.50$pLTBIR,2)," [",
               signif(1e3*lb.fruns.total.50$pLTBIR,2),", ",
               signif(1e3*ub.fruns.total.50$pLTBIR,2),"]")))
print(c(paste0(signif(1e3*med.fruns.total.50$pLTBIS,2)," [",
               signif(1e3*lb.fruns.total.50$pLTBIS,2),", ",
               signif(1e3*ub.fruns.total.50$pLTBIS,2),"]")))

100*med.fruns.total.50/N2050

100*med.fruns.total.50$pLTBIR/(med.fruns.total.50$pLTBIR+med.fruns.total.50$pLTBIS)

med.rate.50 <- (0.15*0.01) * med.fruns.total.50/N2050 * 1e5 # 0.15% reactivation x total  / total population x 100,000
ub.rate.50 <- (0.15*0.01) * ub.fruns.total.50/N2050 * 1e5 
lb.rate.50 <- (0.15*0.01) * lb.fruns.total.50/N2050 * 1e5 

print(c(paste0(sprintf('%.2f',med.rate.50$pLTBIR)," [",
               sprintf('%.2f',lb.rate.50$pLTBIR),", ",
               sprintf('%.2f',ub.rate.50$pLTBIR),"]")))
print(c(paste0(sprintf('%.1f',med.rate.50$pLTBIS)," [",
               sprintf('%.1f',lb.rate.50$pLTBIS),", ",
               sprintf('%.1f',ub.rate.50$pLTBIS),"]")))

med.num.50 <- (0.15*0.01) * 1e3 * med.fruns.total.50
ub.num.50 <- (0.15*0.01) * 1e3 * ub.fruns.total.50
lb.num.50 <- (0.15*0.01) * 1e3 * lb.fruns.total.50

print(c(paste0(sprintf('%.2f',med.num.50$pLTBIR)," [",
               sprintf('%.2f',lb.num.50$pLTBIR),", ",
               sprintf('%.2f',ub.num.50$pLTBIR),"]")))

print(c(paste0(sprintf('%.2f',med.num.50$pLTBIS)," [",
               sprintf('%.2f',lb.num.50$pLTBIS),", ",
               sprintf('%.2f',ub.num.50$pLTBIS),"]")))

#### SENSITIVITY ANALYSIS: lower reactivation
### 
react <- 0.15*0.6

med.rate.35 <- (react*0.01) * med.fruns.total.35/N2035 * 1e5 # 0.15% reactivation x total  / total population x 100,000
ub.rate.35 <- (react*0.01) * ub.fruns.total.35/N2035 * 1e5 
lb.rate.35 <- (react*0.01) * lb.fruns.total.35/N2035 * 1e5 

print(c(paste0(sprintf('%.2f',med.rate.35$pLTBIR)," [",
               sprintf('%.2f',lb.rate.35$pLTBIR),", ",
               sprintf('%.2f',ub.rate.35$pLTBIR),"]")))


med.rate.50 <- (react*0.01) * med.fruns.total.50/N2050 * 1e5 # 0.15% reactivation x total  / total population x 100,000
ub.rate.50 <- (react*0.01) * ub.fruns.total.50/N2050 * 1e5 
lb.rate.50 <- (react*0.01) * lb.fruns.total.50/N2050 * 1e5 

print(c(paste0(sprintf('%.2f',med.rate.50$pLTBIR)," [",
               sprintf('%.2f',lb.rate.50$pLTBIR),", ",
               sprintf('%.2f',ub.rate.50$pLTBIR),"]")))


####**** Group trend figures *****######

## WHO data
w_data <- read.csv("~/Dropbox/MDR/new_who_edited_sub.csv")[,-1]

nari = 200 # up to 200

# DS and MDR data
# Label for plots 
pp <- "infor_prior"

# READ IN
load("output_final/all0_p_ds_mdr.Rdata")
theme_set(theme_bw(base_size = 24))
for(cci in 1:length(cni)){
  print(cni[cci])
  ### WHO data
  d <-subset(w_data, iso3 == as.character(cni[cci]) )
  
  ### ARI for both DS and mDR in all0
  rdata <- all0[which(all0$iso3 == as.character(cni[cci])),]
  
  a1 <- ggplot(d, aes(x=year_new, y = new_mdr_prop),col="red",pch = 10) + # points won't plot over lines unless do points first?!
    geom_point() +
    geom_line(data = rdata, aes(x=year, y = prediction, group = factor(replicate)),alpha = 0.2) +
    scale_y_continuous("Prop. new with MDR") + scale_x_continuous("Year",lim=c(1970,2016)) +
    geom_point(data = d, aes(x=year_new, y = new_mdr_prop),col="red",pch = 10, size = 3) + geom_errorbar(data = d, aes(ymin = mlo, ymax = mhi), col = "red")
  ggsave(paste0("output_final/",cni[cci],"_mdr_trends_with_data_",pp,".pdf"),width=11, height=11)
  
  a2 <- ggplot(rdata, aes(x=year, y = mdr_ari, group = factor(replicate))) + geom_line(alpha = 0.2) +
    scale_y_continuous("MDR ARI") + scale_x_continuous("Year",lim=c(1970,2015))
  ggsave(paste0("output_final/",cni[cci],"_mdr_ari_",pp,".pdf"),width=11, height=11)
  
}    

## All trends
theme_set(theme_bw(base_size = 10))
pdf("output_final/trends_all.pdf")
for(i in 1:9){
  print(ggplot(w_data, aes(x=year_new, y = new_mdr_prop),col="red",pch = 10) + # points won't plot over lines unless do points first?!
          geom_point() +
          geom_line(data = all0, aes(x=year, y = prediction, group = factor(replicate)),alpha = 0.2) +
          scale_y_continuous("Prop. new with MDR") + scale_x_continuous("Year",lim=c(1970,2020)) +
          geom_point(data = w_data, aes(x=year_new, y = new_mdr_prop),col="red",pch = 10) +
          geom_errorbar(data = w_data, aes(ymin = mlo, ymax = mhi), col = "red") +
          facet_wrap_paginate(~iso3, scales = "free",ncol = 4, nrow = 4, page = i)
  )
}
dev.off()

# All ARI trends
theme_set(theme_bw(base_size = 10))
pdf("output_final/trends_ari_all.pdf")
for(i in 1:9){
  print(ggplot(all0, aes(x=year, y = mdr_ari, group = factor(replicate))) + geom_line(alpha = 0.2) +
          scale_y_continuous("ARI with MDR-Mtb") + scale_x_continuous("Year",lim=c(1970,2015)) +
          facet_wrap_paginate(~iso3, scales = "free",ncol = 4, nrow = 4, page = i)
  )
}
dev.off()



### By region
wwm <- merge(w_data, WHOkey, by = "iso3")
wwr <- merge(all0, WHOkey, by = "iso3")

ggplot(wwm[which(wwm$g_whoregion == "SEA"),], aes(x=year_new, y = new_mdr_prop),col="red",pch = 10) + # points won't plot over lines unless do points first?!
  geom_point() + facet_wrap(~country, scales = "free") + 
  theme(strip.text.x = element_text(size = 10)) + 
  geom_line(data = wwr[which(wwr$g_whoregion == "SEA"),], aes(x=year, y = prediction, group = factor(replicate)),alpha = 0.2) +
  scale_y_continuous("Prop. new with MDR") + scale_x_continuous("Year",lim=c(1970,2020), breaks = c(1970, 1990, 2010)) +
  geom_point(data = wwm[which(wwm$g_whoregion == "SEA"),], aes(x=year_new, y = new_mdr_prop),col="red",pch = 10, size = 3) + 
  geom_errorbar(data = wwm[which(wwm$g_whoregion == "SEA"),], aes(ymin = mlo, ymax = mhi), col = "red")
ggsave(paste0("output_final/SEA_mdr_trends_with_data_",pp,".pdf"),width=13, height=11)


# ribbon effect
sea.med  <- aggregate(wwr[which(wwr$g_whoregion == "SEA"),"prediction"], list(wwr[which(wwr$g_whoregion == "SEA"),"year"],wwr[which(wwr$g_whoregion == "SEA"),"country"]), median)
colnames(sea.med) <- c("year","country","med")
sea.ub  <- aggregate(wwr[which(wwr$g_whoregion == "SEA"),"prediction"], list(wwr[which(wwr$g_whoregion == "SEA"),"year"],wwr[which(wwr$g_whoregion == "SEA"),"country"]), ub)
colnames(sea.ub) <- c("year","country","ub")
sea.lb  <- aggregate(wwr[which(wwr$g_whoregion == "SEA"),"prediction"], list(wwr[which(wwr$g_whoregion == "SEA"),"year"],wwr[which(wwr$g_whoregion == "SEA"),"country"]), lb)
colnames(sea.lb) <- c("year","country","lb")

sea.rib <- merge(sea.ub, sea.lb, by = c("year","country"))
sea.rib <- merge(sea.rib, sea.med, by = c("year","country"))

ggplot(wwm[which(wwm$g_whoregion == "SEA"),], aes(x=year_new)) + # points won't plot over lines unless do points first?!
  geom_point( aes(x=year_new, y = new_mdr_prop),col="red",pch = 10, size = 3) + 
  facet_wrap(~country, scales = "free") + 
  theme(strip.text.x = element_text(size = 10)) + 
  geom_line(data = sea.rib, aes(x=year, y = med)) +
  scale_y_continuous("Prop. new with MDR") + scale_x_continuous("Year",lim=c(1970,2020), breaks = c(1970, 1990, 2010)) +
  geom_errorbar(data = wwm[which(wwm$g_whoregion == "SEA"),], aes(ymin = mlo, ymax = mhi), col = "red") + 
  geom_ribbon(data = sea.rib, aes(x = year, ymin = lb, ymax = ub), alpha = 0.3, fill = "blue")
ggsave(paste0("output_final/SEA_mdr_trends_with_data_ribbon_",pp,".pdf"),width=13, height=11)

ggplot(wwr[which(wwr$g_whoregion == "SEA"),], aes(x=year, y = mdr_ari, group = replicate)) + # points won't plot over lines unless do points first?!
  geom_line(alpha = 0.2) + facet_wrap(~country, scales = "free") + 
  theme(strip.text.x = element_text(size = 10)) + 
  scale_y_continuous("ARI with MDR-M.tb") + 
  scale_x_continuous("Year",lim=c(1970,2020), breaks = c(1970, 1990, 2010)) 
ggsave(paste0("output_final/SEA_ari_trends_",pp,".pdf"),width=13, height=11)


### By age
## All ARI trends
theme_set(theme_bw(base_size = 24))
## 2014
# need iso, replicate, act, prev, g_whoregion 
POP2014$acat <- 0:16
fruns <- merge(new_data2[,c("iso3","mdr_rep","acat","pr_dr","pr_ds","g_whoregion")],
               POP2014,by=c('iso3','acat'),all=TRUE)
fruns <- fruns[!is.na(fruns$mdr_rep),]     # ditch those countries not in the 138
dim(fruns) ## 200*138*17 = 469200

fruns$pLTBIR <- fruns$value * fruns$pr_dr 
fruns$pLTBIS <- fruns$value * fruns$pr_ds

fruns.total.14 <- aggregate(fruns[,c("pLTBIR","pLTBIS")], list(fruns$mdr_rep,fruns$g_whoregion,fruns$age), sum)
dim(fruns.total.14) ## 200 * 6 * 17 = 20400
colnames(fruns.total.14) <- c("mdr_rep","g_whoregion","age","pLTBIR","pLTBIS")
med.fruns.total.14 <- aggregate(fruns.total.14[,c("pLTBIR","pLTBIS")], list(fruns.total.14$g_whoregion,fruns.total.14$age), median)
colnames(med.fruns.total.14) <- c("g_who","age","r","s")
mmed <- melt(med.fruns.total.14,  id.vars = c("g_who","age"))
ub.fruns.total.14 <-  aggregate(fruns.total.14[,c("pLTBIR","pLTBIS")], list(fruns.total.14$g_whoregion,fruns.total.14$age), ub)
colnames(ub.fruns.total.14) <- c("g_who","age","r","s")
mub <- melt(ub.fruns.total.14,  id.vars = c("g_who","age"))
lb.fruns.total.14 <-  aggregate(fruns.total.14[,c("pLTBIR","pLTBIS")], list(fruns.total.14$g_whoregion,fruns.total.14$age),lb)
colnames(lb.fruns.total.14) <- c("g_who","age","r","s")
mlb <- melt(lb.fruns.total.14,  id.vars = c("g_who","age"))

plot.14 <- merge(mmed, mub, by = c('g_who', 'age','variable'))
plot.14 <- merge(plot.14, mlb,by = c('g_who', 'age','variable'))

total_by_region_age <- aggregate(fruns[,c("value")], list(fruns$g_whoregion,fruns$age), sum)
colnames(total_by_region_age) <- c("g_who","age","pop_size")# /200 as 200 mdr runs in fruns and summed over all of these
total_by_region_age$pop_size <- total_by_region_age$pop_size / 200

mplot14 <- merge(plot.14, total_by_region_age, by = c('g_who', 'age'))
mplot14$med.perc <- 100*mplot14$value.x / (mplot14$pop_size) 
mplot14$ub.perc <- 100*mplot14$value.y / (mplot14$pop_size)
mplot14$lb.perc <- 100*mplot14$value / (mplot14$pop_size)

ggplot(mplot14, aes(x=age, y = med.perc/100)) + geom_bar(aes(fill = variable),stat='identity', pos = "dodge") + 
  facet_wrap(~g_who) + scale_fill_discrete("LTBI",labels = c("MDR-","DS-")) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_errorbar(aes(ymin = lb.perc/100, ymax = ub.perc/100)) +
  scale_y_continuous("Percentage infected",labels = scales::percent_format(accuracy = 1)) +
  scale_x_discrete("Age")
ggsave(paste0("output_final/LTBI_by_ds_mdr_region_",pp,".pdf"),width=14, height=11)

mplot14r <- mplot14[which(mplot14$variable == "r"),]
ggplot(mplot14r, aes(x=age, y = med.perc/100)) + geom_bar(aes(fill = variable),stat='identity') + 
  geom_errorbar(aes(ymin = lb.perc/100, ymax = ub.perc/100)) + 
  facet_wrap(~g_who)  + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + guides(fill=FALSE) + 
  scale_y_continuous("Percentage infected",labels = scales::percent_format(accuracy = .1)) +
  scale_x_discrete("Age")
ggsave(paste0("output_final/LTBI_mdr_region_",pp,".pdf"),width=14, height=11)


### ARI by region
arim <- merge(all0, WHOkey, by = "iso3")
arimm <- aggregate(arim[,c("mdr_ari")], list(arim$g_whoregion,arim$year,arim$replicate), median)
colnames(arimm) <- c("g_reg","year","mdr_rep","mdr_ari")

# functions
ub <- function(x)quantile(x,probs = .975, na.rm = TRUE)
lb <- function(x)quantile(x,probs = .025, na.rm = TRUE)
arim2 <- aggregate(arim[,c("mdr_ari")], list(arim$g_whoregion,arim$year), median)
arim2_ub <- aggregate(arim[,c("mdr_ari")], list(arim$g_whoregion,arim$year), ub)
arim2_lb <- aggregate(arim[,c("mdr_ari")], list(arim$g_whoregion,arim$year), lb)
colnames(arim2) <- c("g_reg","year","mdr_ari")
colnames(arim2_lb) <- c("g_reg","year","mdr_ari_lb")
colnames(arim2_ub) <- c("g_reg","year","mdr_ari_ub")

arimplot <- merge(arim2, arim2_lb)
arimplot <- merge(arimplot, arim2_ub)

ggplot(arimplot, aes(x=year)) + geom_line(aes(y = mdr_ari)) + 
  geom_ribbon(aes(ymin=mdr_ari_lb, ymax=mdr_ari_ub), alpha=0.3, fill = "red") +
  facet_wrap(~g_reg, scales = "free") + 
  scale_y_continuous("ARI with MDR M.tb") + scale_x_continuous(lim = c(1960, 2020))
ggsave(paste0("output_final/ARI_region_ribbon_",pp,".pdf"),width=14, height=11)

# functions - 50% range
ub7 <- function(x)quantile(x,probs = .75, na.rm = TRUE)
lb7 <- function(x)quantile(x,probs = .25, na.rm = TRUE)
arim2_ub7 <- aggregate(arim[,c("mdr_ari")], list(arim$g_whoregion,arim$year), ub7)
arim2_lb7 <- aggregate(arim[,c("mdr_ari")], list(arim$g_whoregion,arim$year), lb7)
colnames(arim2_lb7) <- c("g_reg","year","mdr_ari_lb")
colnames(arim2_ub7) <- c("g_reg","year","mdr_ari_ub")

arimplot7 <- merge(arim2, arim2_lb7)
arimplot7 <- merge(arimplot7, arim2_ub7)

ggplot(arimplot7, aes(x=year)) + geom_line(aes(y = mdr_ari)) + 
  geom_ribbon(aes(ymin=mdr_ari_lb, ymax=mdr_ari_ub), alpha=0.3, fill = "red") +
  facet_wrap(~g_reg, scales = "free") + 
  scale_y_continuous("ARI with MDR M.tb") + scale_x_continuous(lim = c(1960, 2020))
ggsave(paste0("output_final/ARI_region_ribbon_50_",pp,".pdf"),width=14, height=11)

ggplot(arimm, aes(x=year)) + geom_line(aes(y = mdr_ari, group = mdr_rep), alpha = 0.2) +
  facet_wrap(~g_reg) + scale_x_continuous(lim = c(1960, 2020))
ggsave(paste0("output_final/ARI_region_",pp,".pdf"),width=14, height=11)
