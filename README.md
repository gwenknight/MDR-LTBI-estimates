# MDR-LTBI-paper
Code for MDR-LTBI burden estimates (Knight et al.)

# Code description
cohort_ltbi_mdr.R:
Main cohort model

Run_cohort_ltbi_mdr_final.R:
Run the cohort model 

Fit_model_for_prop_mdr.R:
Runs Stan model

model_trend.stan:
Stan code for model fit

Output_production.R:
Code for gathering the data and plotting figures for the paper

# Data description
138_final_list_included_countries.csv: the final 138 included countries

all0_p_ds_mdr.Rdata: the 200 posterior estimates of MDR and DS ARI over time for each country

POP2014.Rdata / POP2035.Rdata / POP2050.Rdata: UN population estimates of population size

whokey.Rdata: WHO key for translating iso3 to country names and grouping by regions
