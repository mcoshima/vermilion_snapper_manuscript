### Creating dummy data used for running bootstrap
## Run this after creating starting states - starting states are in the start_pop_i folder so run bootstrap in there for creating first 5 years of data then copy files to OM and EM folders and continue simulations from there

library(r4ss)
library(tidyverse)

## args is a vector that will be specified on work script but for now is a vector that contains: OM_scenario name, iteration (1-100)

args <- c("perfect", 2)
om.scenario <- args[1]
iteration <- args[2]

dir. <- paste0(getwd(), "/vermilion_snapper_manuscript/with_comp")

dir.OMscen <- paste0(getwd(), "/vermilion_snapper_manuscript/with_comp/start_pop_", iteration, "/", om.scenario)
dir.create(dir.OMscen)
dir.OM <- paste0(dir.OMscen, "/OM")
dir.create(dir.OM)
#dir.boot <- dir.OM

#to get functions like splt.recombine
source(file.path(getwd(), "vermilion_snapper_manuscript", "r_scripts", "utils.R"))
source(file.path(getwd(), "vermilion_snapper_manuscript", "r_scripts", "init_OM.R"))

## Copy files into OM directory
files. <- c("forecast.ss", "ss.exe", "ss.par", "vs.ctl", "starter.ss", "vs_envM2.dat")
file.copy(paste0(dir., "/start_pop_", iteration, "/", files.), dir.OM, overwrite = TRUE)


## Create sample structure data
colnames(boot_dat$catch) <- c("year", "seas", "fleet", "catch_se")
colnames(boot_dat$CPUE) <- c("year", "seas", "index", "se_log")
colnames(boot_dat$lencomp) <- c("Yr", "Seas", "FltSvy", "Gender", "Part", "Nsamp")
colnames(boot_dat$agecomp) <- c("Yr", "Seas", "FltSvy", "Gender", "Part", "Ageerr", "Lbin_lo", "Lbin_hi", "Nsamp")

proj.index.full <-
  read.csv(file.path(dir., "../TidyData/RS.biomass.report.csv"))

proj.index <-
  proj.index.full %>% 
  select(Yr, RS_relative) %>% 
  filter(Yr >= 2018) %>% 
  rename(Year = Yr)


#number of years to generate data for
n.years.fwd <- 5

dat. <- SS_readdat_3.30(file = file.path(dir.OM, "vs_envM2.dat"))

#year starting data generation from 
start.yr <- dat.$endyr + 1

## Updating DAT file
endyr <- dat.$endyr
#dat.$endyr <- dat.$endyr + n.years.fwd 

## Find patterns in historical data
cpue.patterns <- dat.$CPUE %>% 
  group_by(index) %>% 
  summarise(term.yr = max(year),
            mean_se = mean(se_log)) %>% 
  filter(term.yr == endyr)

catch.patterns <- dat.$catch %>% 
  group_by(fleet) %>% 
  summarise(term.yr = max(year),
            mean.se = mean(catch_se)) %>% 
  filter(term.yr == endyr)

discard.patterns <- dat.$discard_data %>% 
  group_by(Flt) %>% 
  summarise(term.yr = max(Yr),
            mean.se = mean(Std_in)) %>% 
  filter(term.yr == endyr)

lencomp.patterns <- dat.$lencomp %>% 
  group_by(FltSvy) %>% 
  summarise(term.yr = max(Yr),
            meanN = mean(Nsamp))
lbins <- dat.$lbin_vector

agecomp.patterns <- dat.$agecomp %>% 
  group_by(FltSvy) %>% 
  summarise(term.yr = max(Yr),
            meanN = mean(Nsamp))
abins <- dat.$agebin_vector

## Build dummy data dfs
### Catch df
catch <- data.frame(
  year = rep(seq(start.yr, start.yr + n.years.fwd - 1), nrow(catch.patterns)),
  seas = rep(1, n.years.fwd*nrow(catch.patterns)),
  fleet = rep(unique(catch.patterns$fleet), each = n.years.fwd),
  catch_se = rep(catch.patterns$mean.se, each = n.years.fwd)
)


### CPUE df
CPUE <- data.frame(
  year = rep(seq(start.yr, start.yr + n.years.fwd - 1), nrow(cpue.patterns)),
  seas = rep(7, n.years.fwd*nrow(cpue.patterns)),
  index = rep(cpue.patterns$index, each = n.years.fwd),
  se_log = rep(cpue.patterns$mean_se, each = n.years.fwd)
)

discard <- data.frame(
  Yr = rep(seq(start.yr, start.yr + n.years.fwd - 1), nrow(discard.patterns)),
  Seas = rep(7, n.years.fwd*nrow(discard.patterns)),
  Flt = rep(discard.patterns$Flt, each = n.years.fwd), 
  Discard = rep(1, n.years.fwd*nrow(discard.patterns)), 
  Std_in = rep(discard.patterns$mean.se, each = n.years.fwd)
)
### length comp
lencomp <- data.frame(
  Yr = rep(seq(start.yr, start.yr + n.years.fwd - 1), nrow(lencomp.patterns)), 
  Seas = rep(7, n.years.fwd*nrow(lencomp.patterns)), 
  FltSvy = rep(lencomp.patterns$FltSvy, each = n.years.fwd), 
  Gender = rep(0, n.years.fwd*nrow(lencomp.patterns)), 
  Part = rep(0, n.years.fwd*nrow(lencomp.patterns)), 
  Nsamp = rep(75, n.years.fwd*nrow(lencomp.patterns))
)

### age comp
agecomp <- data.frame(
  Yr = rep(seq(start.yr, start.yr + n.years.fwd - 1), nrow(agecomp.patterns)), 
  Seas = rep(7, n.years.fwd*nrow(agecomp.patterns)), 
  FltSvy = rep(agecomp.patterns$FltSvy, each = n.years.fwd), 
  Gender = rep(0, n.years.fwd*nrow(agecomp.patterns)), 
  Part = rep(0, n.years.fwd*nrow(agecomp.patterns)), 
  Ageerr = rep(1, n.years.fwd*nrow(agecomp.patterns)),
  Lbin_lo = rep(-1, n.years.fwd*nrow(agecomp.patterns)),
  Lbin_hi = rep(-1, n.years.fwd*nrow(agecomp.patterns)),
  Nsamp = rep(75, n.years.fwd*nrow(agecomp.patterns))
)


### Envrionmental data (competition index)
### Note right now it is just the perfectly known
envdat <- data.frame(
  Yr = seq(start.yr, start.yr + n.years.fwd - 1),
  Variable = rep(1, n.years.fwd),
  Value = proj.index[which(proj.index$Year >= start.yr & proj.index$Year <= start.yr + n.years.fwd - 1),2]
)


boot_dat <- list(
  "catch" = catch,
  "CPUE" = CPUE,
  "discard_data" = discard,
  "lencomp" = lencomp,
  "agecomp" = agecomp,
  "envdat" = envdat
)

### Create OM 
create_OM(OM_out_dir = dir.OM, overwrite = TRUE, nyrs = 5, nyrs_assess = 4, scen_name = "perfect", sample_struct = boot_dat, os = "win", seed = seed)

## Create bootstrapped data and update OM and EM
start <- r4ss::SS_readstarter(paste0(dir.OM, "/starter.ss"))
start$N_bootstraps <- 3
r4ss::SS_writestarter(start, dir = dir.OM, overwrite = TRUE)

run_ss(dir.OM, os = "win", "-maxfn 0 -phase 50 -nohess")

OM.dat <- r4ss::SS_readdat_3.30(file = paste0(dir.OM, "/data.ss_new"), section = 2)


EM.dat <- r4ss::SS_readdat_3.30(file = paste0(dir.OM, "/data.ss_new"), section = 3)

### Why is recreational catch so different in bootstrapped years?
EM.dat$catch %>% 
  filter(fleet ==3) %>% 
  filter(year > 2000 & year < 2018) %>% 
  ggplot(aes(x = year, y = catch)) +
  geom_point() +
  geom_line()

OM.dat$catch %>% 
  filter(year > 1980) %>% 
  ggplot(aes(x = year, y = catch)) +
  geom_point() +
  geom_line() +
  geom_point(data = EM.dat$catch[which(EM.dat$catch$year > 1980),], 
             aes(x = EM.dat$catch[which(EM.dat$catch$year > 1980), 1],
                 y = EM.dat$catch[which(EM.dat$catch$year > 1980), 4]), color = "blue") +
  facet_wrap(~fleet, scales = "free") 

report <- SS_output(dir = dir.OM)
report$timeseries
