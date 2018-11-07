#### Cohort model for MDR-LTBI estimates

cohort_ltbi <- function(ari,pop,prot = 0){
  ### inputs 
  # ari: 2 columns: ds and mdr ARI
  # pop: population in 2014
  # prot: if switching on fitness cost to resistance, prot = 1
  
  #### Protection from re-infecton
  ## Andrews: 0.79 .7-.86 
  pm <- 0.79                      # Mean        
  pv <- (0.86-0.7)^2/3.92^2
  apb <- pm*(1-pm)/pv-1
  # Shape parameters
  pa <- pm*apb                            #77.88
  pb <- (1-pm)*apb                        #20.70
  
  ### Random sample of level of protection - beta distribution
  # same for dr/ds
  alphs <- rbeta(100*81,shape1=pb,shape2=pa)
  dim(alphs) <-  c(81, 100)
  
  alphr <- alphs # standard assumption that same level of protection
  
  # if prot == 1 then Scenario analysis = 1
  if(prot == 1){
    print("Fitness affects protection ON")
    rf = 0.6
    alphr <- (1-rf*pm) / (1-pm) * alphs ### alphr = 1 - relfit * pm = y * (1 - pm) = y * alphs # reorganise to give factor multiplying alphs
    # mean(alphr) # ~ 0.53 = 1 - 0.6*0.78
  } 
  
  
  # Matrix framework
  c_last <- as.data.frame(matrix(0,100,8))
  colnames(c_last) <- c("pr_ds","pr_dr","new_ds","new_dr","rei_sr","rei_rs","rei_ss","rei_rr")
  c_now <- as.data.frame(matrix(0,100,8))
  colnames(c_now) <- c("pr_ds","pr_dr","new_ds","new_dr","rei_sr","rei_rs","rei_ss","rei_rr")
  store_c <- as.data.frame(matrix(0,100*81,8))
  colnames(store_c) <- colnames(c_last)
  
  # Number of years (1934:2014)
  nyrs <- 81
  
  ## Initial conditions
  # hazard = 
  H = seq(0,99,1) * as.numeric(ari[1,"ds"]) # constant ARI pre-1934
  c_last[,"pr_ds"] <- 1 - exp(-H)
  
  # Cycle thru from 1934 to 2014, updating the cohort
  for(i in 1:nyrs){
    c_now[1,] <- 0
    c_now[2:100,] <- c_last[1:99,] # New time, 100+yos not included
    
    ## Add together existing with new proportions ds / dr
    ## remove reinfecteds move to the other... 
    c_now$pr_ds <- c_now$pr_ds + c_now$new_ds + c_now$rei_rs - c_now$rei_sr
    c_now$pr_dr <- c_now$pr_dr + c_now$new_dr + c_now$rei_sr - c_now$rei_rs 
    # Set new to zero
    c_now[,c("new_ds","rei_rs","new_dr","rei_sr")] <- 0
    
    # New calculations
    # ARI for this year
    ari_s <- as.numeric(ari[i,"ds"]) 
    ari_s[ari_s > 1] <- 1 # cap at everyone!
    ari_r <- as.numeric(ari[i,"mdr"]) 
    ari_r[ari_r > 1] <- 1
    
    # Calculations
    c_now$new_ds <- ari_s * (1 - c_now$pr_ds - c_now$pr_dr) # currently none, new infection DS
    c_now$new_dr <- ari_r * (1 - c_now$pr_ds - c_now$pr_dr - c_now$new_ds) # currently none, new infection DR
    c_now[c_now$new_ds < 0, "new_ds"] <- 0 # no new if go negative
    c_now[c_now$new_dr < 0, "new_dr"] <- 0 # no new if go negative
    c_now$rei_sr <- c_now$pr_ds * ari_r * alphs[i,] 
    c_now$rei_rs <- c_now$pr_dr * ari_s * alphr[i,]
    c_now$rei_rr <- c_now$pr_dr * ari_r * alphr[i,] # reinfection but no change to overall proportions: for recent infection calcs
    c_now$rei_ss <- c_now$pr_ds * ari_s * alphs[i,] # reinfection but no change to overall proportions: for recent infection calcs
    if(length(which(c_now < 0))>0){print(c(i,"negative c_now", ari_s, ari_r,c_now[i,],"which",which(c_now < 0)))} # error catching
    
    ### Store and update
    c_last <- c_now  
    store_c[((i-1)*100 + 1):(i*100),] <- c_now
  }
  
  # Add in year
  store_c$year <- rep(1934:2014,each = 100)
  # Add in age
  store_c$age <- rep(1:100,each = 1)
  
  ## Capture last one = 2014
  c_2014 <- c_last
  c_2014$pr_ds <- c_2014$pr_ds + c_2014$new_ds + c_2014$rei_rs - c_2014$rei_sr
  c_2014$pr_dr <- c_2014$pr_dr + c_2014$new_dr + c_2014$rei_sr - c_2014$rei_rs 
  
  #Number in each age group 
  # group by 5 yr to 80+ 
  c_ds_age <- colwise(mean)(as.data.frame(matrix(c_2014$pr_ds, 5)))
  ds_age <- as.numeric(c(c_ds_age[1:16], mean(as.numeric(c_ds_age[17:20]))))
  c_dr_age <- colwise(mean)(as.data.frame(matrix(c_2014$pr_dr, 5)))
  dr_age <- as.numeric(c(c_dr_age[1:16], mean(as.numeric(c_dr_age[17:20]))))
  
  # size of the population infected in 2014
  size_ds <- pop * ds_age
  size_dr <- pop * dr_age
  
  # percentage of the population
  perc_ds <- 100*size_ds / sum(pop)
  perc_dr <- 100*size_dr / sum(pop)
  
  combs <- as.data.frame(cbind(size_ds,size_dr,perc_ds,perc_dr))
  colnames(combs)<- c("size_ds","size_dr","perc_ds","perc_dr")
  
  return(list(store_c = store_c, c_2014 = c_2014,combs = combs ))
  
}

