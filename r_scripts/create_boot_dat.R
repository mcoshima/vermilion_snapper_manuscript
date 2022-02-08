### Creating dummy data used for running bootstrap
## Run this after creating starting states - starting states are in the start_pop_i folder so run bootstrap in there for creating first 5 years of data then copy files to OM and EM folders and continue simulations from there

library(r4ss)
library(tidyverse)

## args is a vector that will be specified on work script but for now is a vector that contains: OM_scenario name, iteration (1-100)

args <- c("perfect", 2)
om.scenario <- args[1]
iteration <- args[2]

dir. <- paste0(getwd(), "/vermilion_snapper_manuscript/with_comp")
dir.boot <- file.path(dir., paste0("start_pop_", iteration))
#to get functions like splt.recombine
source(file.path(getwd(), "vermilion_snapper_manuscript", "r_scripts", "utils.R"))


proj.index.full <-
  read.csv(file.path(dir., "../TidyData/RS.biomass.report.csv"))

proj.index <-
  proj.index.full %>% 
  select(Yr, RS_relative) %>% 
  filter(Yr >= 2018) %>% 
  rename(Year = Yr)
# proj.index$low <- NA
# proj.index$med <- NA
# proj.index$high <- NA
# 
# for(i in 1:nrow(proj.index)){
#   proj.index$low[i]<- rnorm(1, proj.index$RS_relative[i], .)
#   proj.index$med[i]<- rnorm(1, proj.index$RS_relative[i], .01)
#   proj.index$high[i]<- rnorm(1, proj.index$RS_relative[i], .1)
#   
# }


#number of years to generate data for
n.years.fwd <- 5

dat. <- SS_readdat_3.30(file = file.path(dir.boot, "vs_envM2.dat"))

#year starting data generation from 
#start.yr <- dat.$endyr + 1

## Updating DAT file
# old.endyr <- dat.$endyr
# dat.$endyr <- dat.$endyr + n.years.fwd 

## Find patterns in historical data
cpue.patterns <- dat.$CPUE %>% 
  group_by(index) %>% 
  summarise(term.yr = max(year),
    mean_se = mean(se_log)) %>% 
  filter(term.yr == old.endyr)

catch.patterns <- dat.$catch %>% 
  group_by(fleet) %>% 
  summarise(term.yr = max(year),
            mean.se = mean(catch_se)) %>% 
  filter(term.yr == old.endyr)

discard.patterns <- dat.$discard_data %>% 
  group_by(Flt) %>% 
  summarise(term.yr = max(Yr),
            mean.se = mean(Std_in)) %>% 
  filter(term.yr == old.endyr)

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
  catch = rep(0.001, n.years.fwd*nrow(catch.patterns)),
  catch_se = rep(catch.patterns$mean.se, each = n.years.fwd)
)

### CPUE df
CPUE <- data.frame(
  year = rep(seq(start.yr, start.yr + n.years.fwd - 1), nrow(cpue.patterns)),
  seas = rep(7, n.years.fwd*nrow(cpue.patterns)),
  index = rep(cpue.patterns$index, each = n.years.fwd),
  obs = rep(1, n.years.fwd*nrow(cpue.patterns)),
  se_log = rep(cpue.patterns$mean_se, each = n.years.fwd)
)

### discard df
discard_data <- data.frame(
  Yr = rep(seq(start.yr, start.yr + n.years.fwd - 1), nrow(discard.patterns)),
  Seas = rep(7, n.years.fwd*nrow(discard.patterns)),
  Flt = rep(discard.patterns$Flt, each = n.years.fwd),
  Discard = rep(0, n.years.fwd*nrow(discard.patterns)),
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

## make a matrix of 0s for dummy len/age comp then set colnames to paste0("l", lbins)
dum.len <- as.data.frame(matrix(data = 1, nrow = nrow(lencomp), ncol = length(lbins)))
colnames(dum.len) <- paste0("l", lbins)
lencomp <- cbind(lencomp, dum.len)

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

## make a matrix of 0s for dummy len/age comp then set colnames to paste0("l", lbins)
dum.age <- as.data.frame(matrix(data = 1, nrow = nrow(agecomp), ncol = length(abins)))
colnames(dum.age) <- paste0("a", abins)
agecomp <- cbind(agecomp, dum.age)

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
                "discard_data" = discard_data,
                "lencomp" = lencomp,
                "agecomp" = agecomp,
                "envdat" = envdat
                 )

### Combine new and old data
dat.$catch <- as.data.frame(splt.recombine(dat.$catch, 
                            boot_dat$catch, 
                            fleet, 
                            length(unique(boot_dat$catch$fleet))))

dat.$CPUE <- as.data.frame(splt.recombine(dat.$CPUE, 
                           boot_dat$CPUE, 
                           index, 
                           length(unique(boot_dat$CPUE$index))))

discard.splt <- boot_dat$discard_data %>% 
  group_by(Flt) %>% 
  group_split() 

dat.discard.splt <- dat.$discard_data %>% 
  group_by(Flt) %>% 
  group_split()
dis <- list()
dis[[1]] <- rbind(dat.discard.splt[[2]], discard.splt[[2]])
dis[[2]] <- rbind(dat.discard.splt[[3]], discard.splt[[3]])
dis[[3]] <- rbind(dat.discard.splt[[4]], discard.splt[[4]])
dis[[4]] <- rbind(dat.discard.splt[[5]], dat.discard.splt[[1]], discard.splt[[1]])
dat.$discard_data <- as.data.frame(do.call(rbind, dis))
dat.$discard_data <- dat.$discard_data %>% 
  mutate(Seas = ifelse(Yr == 1972 & Flt == 4 | Yr == dat.$endyr & Flt == -4, -7, 7))

dat.$lencomp <- as.data.frame(splt.recombine(dat.$lencomp, 
                          boot_dat$lencomp, 
                          FltSvy, 
                          2))

dat.$agecomp <- as.data.frame(splt.recombine(dat.$agecomp, 
                          boot_dat$agecomp, 
                          FltSvy, 
                          3))

dat.$envdat <- rbind(dat.$envdat, boot_dat$envdat)

SS_writedat_3.30(dat., outfile = file.path(dir.boot, "vs_envM2.dat"), overwrite = T)

## Update Starter for bootstrap files
starter <- SS_readstarter(file = file.path(dir.boot, "starter.ss"))
starter$N_bootstraps <- 3
starter$init_values_src <- 1
starter$last_estimation_phase <- 0
SS_writestarter(starter, dir = dir.boot, overwrite = TRUE)

shell(paste("cd/d", dir.boot, "&& ss -nohess >NUL 2>&1", sep = " "))
