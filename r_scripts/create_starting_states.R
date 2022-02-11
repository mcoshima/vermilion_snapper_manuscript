### Create 100 unique starting scenarios for VS MSE
##

library(r4ss)
library(dplyr)
library(stringr)

setwd('..')
dir. <- paste0(getwd(), "/with_comp")
dir.mcmc <- file.path(dir., "mcmc")
dir.create(dir.mcmc)
set.seed <- read.csv(file.path(getwd(), "TidyData", "setseed.csv"))
source(file.path(getwd(), "r_scripts", "utils.R"))

## Read in starter file and change settings for running MCMC 
start <- SS_readstarter(file = paste0(dir., "/initial_pop/starter.ss"))
start$init_values_src <- 0
start$run_display_detail <- 0
start$MCMCburn <- 100
start$MCMCthin <- 1
SS_writestarter(start, dir = paste0(dir., "/initial_pop/"), overwrite = TRUE)

## Copy files over to mcmc directory
files. <- c("forecast.ss", "vs.ctl", "vs_envM2.dat", "starter.ss", "ss_osx")
files.path <- file.path(dir., "initial_pop", files.)
file.copy(from = files.path, to = dir.mcmc, overwrite = TRUE, recursive = TRUE)


#Run ss with mcmc and mcsave first then run ss -mceval to estimate posteriors
run_ss(dir.mcmc, os = "osx", xtras = "-mcmc 11000 -mcsave 10")
run_ss(dir.mcmc, os = "osx", xtras = "-mceval")
# system(paste("cd", dir.mcmc, "&& ./ss_osx -mcmc 11000 -mcsave 10 > /dev/null 2>&1", sep = " "))
# system(paste("cd", dir.mcmc, "&& ./ss_opt_osx -mceval > /dev/null 2>&1", sep = " "))

## for windows
# shell(paste("cd/d", dir.mcmc, "&& ss -mcmc 11000 -mcsave 10 >NUL 2>&1", sep = " "))
# shell(paste("cd/d", dir.mcmc, "&& ss -mceval >NUL 2>&1", sep = " "))

files. <- list("vs_envM2.dat", "vs.ctl", "forecast.ss", "starter.ss", "ss_osx", "ss.par")

## Create starting states folders
iterations <- seq(1, 10, by = 1)
dir.it <- paste0(dir., "/start_pop_", iterations)
sapply(dir.it, dir.create)

for(i in 1:length(dir.it)){
  files.path <- file.path(dir.mcmc, files.)
  file.copy(files.path, file.path(dir.it[i]), overwrite = TRUE)
}


rep.mcmc <- SS_output(dir., dir.mcmc = "mcmc")
mcmc <- rep.mcmc$mcmc

#### All of this needs to be repeated for 100 directories (using i to indicate which iteration(row of mcmc) to use)
####################################################################################
#Read in pars file from iteration directory
for(i in 1:length(iterations)){
pars <- SS_readpar_3.30(parfile = file.path(dir.it[i], "ss.par"), 
                        datsource = file.path(dir.it[i], "vs_envM2.dat"), 
                        ctlsource = file.path(dir.it[i], "vs.ctl"))

#Create IDs for pars by group to match up to mcmc labels
mg <- row.names(pars$MG_parms)
sr <- row.names(pars$SR_parms)
s <- row.names(pars$S_parms)
recdevs <- paste0("Main_RecrDev_", pars$recdev1[,1])
recdevs.fore <- c(paste0("Late_RecrDev_", 
                         pars$recdev_forecast[which(pars$recdev_forecast[,1] <= 2022),1]),
                  paste0("ForeRecr_", pars$recdev_forecast[which(pars$recdev_forecast[,1] > 2022),1]))
q <- row.names(pars$Q_parms)
f <- paste("F_fleet", 
           pars$F_rate$fleet, 
           "YR",
           pars$F_rate$year, 
           "s", 
           pars$F_rate$seas, 
           sep = "_")
mcmc.cols <- colnames(rep.mcmc$mcmc)

## Subset mcmc.cols
mcmc.sr <- mcmc[which(str_detect(mcmc.cols, "SR"))]
mcmc.mg <- mcmc[which(str_detect(mcmc.cols, 
                            "NatM
                            |L_at_
                            |VonBert_k
                            |CV_
                            |Wtlen_
                            |Mat
                            |Eggs/kg_
                            |CohortGrow
                            |FracFemale
                            |M2"))]
mcmc.s <- mcmc[which(str_detect(mcmc.cols, "Age|Size"))]
mcmc.f <- mcmc[which(str_detect(mcmc.cols, "F_fleet_"))]
mcmc.q <- mcmc[which(str_detect(mcmc.cols, "LnQ"))]
mcmc.recdevs <- mcmc[which(str_detect(mcmc.cols, paste0(recdevs, collapse = "|")))]
mcmc.recfore <- mcmc[which(str_detect(mcmc.cols, paste0(recdevs.fore, collapse = "|")))]

#check to make sure you have the right number of values to replace
if(nrow(pars$F_rate) == ncol(mcmc.f)){
  
  pars$F_rate$F <- t(mcmc.f[i,])
  
}else{
  message("Mismatched number of F_rate values to replace original values")
}

if(nrow(pars$recdev1) == ncol(mcmc.recdevs)){
  
  pars$recdev1[,2] <- t(mcmc.recdevs[i,])
  
}else{
  message("Mismatched number of rec.dev values to replace original values")
}

if(nrow(pars$recdev_forecast) == ncol(mcmc.recfore)){
  
  pars$recdev_forecast[,2] <- t(mcmc.recfore[i,])
}else{
  message("Mismatch of forecast recdevs to replace original values")
}

if(ncol(mcmc.q) == 1){
  
  pars$Q_parms[which(row.names(pars$Q_parms) == colnames(mcmc.q)), 2] <- mcmc.q[i,]
  
}else{
  message("More than 1 Q value is present")
}

if(ncol(mcmc.sr) == 2){
  
  pars$SR_parms[which(row.names(pars$SR_parms) == colnames(mcmc.sr)[1]),2] <- mcmc.sr[i,1]
  pars$SR_parms[which(row.names(pars$SR_parms) == colnames(mcmc.sr)[2]),2] <- mcmc.sr[i,2]
  
}else{
  message("More than 2 SR variables present")
}

S_parms <- pars$S_parms %>% 
  filter(str_detect(row.names(.), "BLK", negate = TRUE)) %>% 
  filter(str_detect(row.names(.), "VIDEO|SEAMAP|Age")) %>% 
  filter(str_detect(row.names(.), "SMP_BYC|COMP", negate = TRUE)) %>% 
  tibble::rownames_to_column(var = "ID")

if(nrow(S_parms) == ncol(mcmc.s)){
  
  S_parms$ESTIM <- t(mcmc.s[i,])
  S_parms <- pars$S_parms %>% 
    tibble::rownames_to_column(var = "ID") %>% 
    mutate(order = seq(1, nrow(.))) %>% 
    merge(S_parms, by = "ID", all.x = TRUE) %>% 
    mutate(INIT.x = ifelse(is.na(INIT.y), INIT.x, INIT.y),
           ESTIM.x = ifelse(is.na(ESTIM.y), ESTIM.x, ESTIM.y)) %>% 
    arrange(order) %>% 
    select(c(ID, INIT.x, ESTIM.x)) %>% 
    tibble::column_to_rownames(var = "ID") %>% 
    rename(INIT = INIT.x,
      ESTIM = ESTIM.x)
  
  pars$S_parms[,1:2] <- S_parms[,1:2]
}else{
  message("Mismatch number of selectivity variables to replace originial values")
}

if(ncol(mcmc.mg) == 2){
  
  pars$MG_parms[which(str_detect(row.names(pars$MG_parms), "M2")), 2] <- t(mcmc.mg[i,])
  
}else{
  message("More than 2 MG values are present")
}

#Save new pars to directory
SS_writepar_3.30(pars, outfile = file.path(dir.it[i], "ss.par"), overwrite = TRUE)

### rerun model to get starting states
starter.i <- SS_readstarter(file = file.path(dir.it[i], "starter.ss"))
starter.i$init_values_src <- 1                             
starter.i$last_estimation_phase <- 0
starter.i$seed <- set.seed[i,2]
SS_writestarter(starter.i, dir = file.path(dir.it[i]), overwrite = TRUE)

run_ss(dir.it[i], os = "osx", xtras = "-nohess")

}

### Compare starting states ######
dir.create(path = file.path(dir., "starting_states_comp_plots"))
big.reps <- SSgetoutput(dirvec = c(dir.,  dir.it), verbose = FALSE) ###TODO fix this so it gets all rep files from any start_pop_ folder
sum.reps <- SSsummarize(big.reps)
SSplotComparisons(sum.reps, print = TRUE, plotdir = file.path(dir., "starting_states_comp_plots"), legend = FALSE)
