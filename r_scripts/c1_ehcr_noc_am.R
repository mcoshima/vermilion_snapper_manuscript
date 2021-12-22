##This code is the framework for the MSE with VS competition index with the empirical rules setting catch limit
#Running the model but ignoring the competition in the assessment model

# library(devtools)
# devtools::install_github("mcoshima/moMSE")
# install.packages("C:/Users/w986430/Documents/GitHub/moMSE/moMSE_0.1.0.tar.gz",
#repos = NULL,
#type = "source")

library(moMSE)
library(r4ss)
library(dplyr)
library(ggplot2)
library(stringr)
library(psych) #for geometric mean
library(tidyr)
library(purrr)
library(truncnorm)

args = commandArgs(trailingOnly = TRUE)
setwd('..')
dir. <- paste0(getwd(), "/one_plus")
report. <- MO_SSoutput(dir = file.path(dir., "initial_pop", "initial_w_compF"))

switch_sr_phases <- function(directory, ln_ro = TRUE, sig_r = FALSE){
  
  ctl = readLines(directory)  
  sr_ln <- ctl %>% str_detect(., "SR_LN") %>% which(isTRUE(.))
  
  sr_sigr <- ctl %>% str_detect(., "SR_sigmaR") %>% which(isTRUE(.))
  
  if(isTRUE(ln_ro)){
    
    ro.vec <- strsplit(ctl[sr_ln], split = " ") 
    ro.vec[[1]][8] <- as.character(as.numeric(ro.vec[[1]][8]) * -1)
    ctl[sr_ln] <- paste(ro.vec[[1]], collapse = " ")
    
    sirg.vec <- strsplit(ctl[sr_sigr], split = " ") 
    neg <- sign(as.numeric(sirg.vec[[1]][8]))
    if(neg == 1){
      
      sirg.vec[[1]][8] <- as.character(as.numeric(sirg.vec[[1]][8]) * -1)
      ctl[sr_sigr] <- paste(sirg.vec[[1]], collapse = " ")
    }
    
  }
  
  if(isTRUE(sig_r)){
    
    sirg.vec <- strsplit(ctl[sr_sigr], split = " ") 
    sirg.vec[[1]][8] <- as.character(as.numeric(sirg.vec[[1]][8]) * -1)
    ctl[sr_sigr] <- paste(sirg.vec[[1]], collapse = " ")
    
    ro.vec <- strsplit(ctl[sr_ln], split = " ") 
    neg <- sign(as.numeric(ro.vec[[1]][8]))
    if(neg == 1){
      
      ro.vec[[1]][8] <- as.character(as.numeric(ro.vec[[1]][8]) * -1)
      ctl[sr_ln] <- paste(ro.vec[[1]], collapse = " ")
    }
  }
  
  writeLines(ctl, directory)
  
}


dat.init <- SS_readdat(file.path(dir., "initial_pop", "initial_w_compF", "/VS.dat"))
seed <- read.csv("setseed.csv")
seed <- as.vector(seed[,2])
start.year <- 2014
#Data and parameter setup
Nyrs <- 100 #Double the number of projection years
Year.vec <- seq(2, 100, by = 2)
year.seq <- seq(start.year, start.year + 50, by = .5) #for function
Nfleet <- 7 #ComEIFQ/ComWIFQ/Rec/ShrimpBC/Comp/Video/Groundfish
Nages <- 15 #0 to 14+
Nsize <- report.$nlbins
assessYrs <- seq(2019, 2069, by = 5)

#Recruitment
recruit <- report.$recruit
colnames(recruit) <-
  c(
    "Year",
    "spawn_bio",
    "exp_recr",
    "with_env",
    "adjusted",
    "pred_recr",
    "dev",
    "biasedaj",
    "era"
  )
head(recruit)
rhat <-
  geometric.mean(recruit$pred_recr[c(39:66)]) #recruitment estimates from 1993 to 2014, prior years were just based on equation
rsig <-
  report.$parameters[which(report.$parameters$Label == "SR_sigmaR"), 3]
rec.dev <- rsig * rnorm(1000, 0, 1) - (rsig ^ 2 / 2)
rec.dis <- c()
for (i in 1:100) {
  rec.dis[i] <- rhat * exp(sample(rec.dev, 1))
}

B0 <- report.$timeseries %>%
  slice(1L) %>%
  select(Bio_all) %>%
  pull()
##Population dynamics
##Weight-at-length for midpoint lengths
a <- 2.19E-5
b <- 2.916
wtatlen <- a * report.$lbins ^ b
wtatage <- tail(report.$wtatage[, -c(1:6)], 1)
linf <- report.$Growth_Parameters$Linf
k <- report.$Growth_Parameters$K
t0 <- report.$Growth_Parameters$A_a_L0

##Mortality
###Fishing mortality
catage <- report.$catage[,-c(4,5,9)]
E_4 <-
  report.$discard %>% filter(Yr > 2011 &
                               Yr <= start.year) %>% select(F_rate) %>% pull()
E_4_mu <- mean(log(E_4))
E_4_sd <- sd(log(E_4))

c. <- dat.init$catch %>% 
  filter(year > 2006) %>% 
  select(1:3) %>% 
  mutate(EP = CM_E/rowSums(.), 
         WP = CM_W/rowSums(.), 
         RP = REC/rowSums(.)) %>% 
  mutate(ce.sd = sd(CM_E), 
         cw.sd = sd(CM_W), 
         rec.sd = sd(REC)) %>% 
  colMeans(.) 

c_1 <- c.[c(1,7)]
c_2 <- c.[c(2,8)]
c_4 <- dat.init$discard_data %>% tail(n = 4) %>% select(Discard) %>% range()

c_3 <- catage %>%
  filter(Yr > 2010 & Fleet == 3) %>%
  select(-c(1:7)) %>%
  rename(
    "Zero" = "0",
    "One" = "1",
    "Two" = "2",
    "Three" = "3",
    "Four" = "4",
    "Five" = "5",
    "Six" = "6",
    "Seven" = "7",
    "Eight" = "8",
    "Nine" = "9",
    "Ten" = "10",
    "Eleven" = "11",
    "Tweleve" = "12",
    "Thirteen" = "13",
    "Fourteen" = "14",
  ) %>%
  transmute(Zero = Zero * wtatage[1, 1],
            One = One * wtatage[1, 2],
            Two = Two * wtatage[1, 3],
            Three = Three * wtatage[1, 4],
            Four = Four * wtatage[1, 5],
            Five = Five * wtatage[1, 6],
            Six = Six * wtatage[1, 7],
            Seven = Seven * wtatage[1, 8],
            Eight = Eight * wtatage[1, 9],
            Nine = Nine * wtatage[1, 10],
            Ten = Ten * wtatage[1, 11],
            Eleven = Eleven * wtatage[1, 12],
            Tweleve = Tweleve * wtatage[1, 13],
            Thirteen = Thirteen * wtatage[1, 14],
            Fourteen = Fourteen * wtatage[1, 15] * wtatage[1, 1]) %>%
  mutate(total = rowSums(.)) %>%
  select(total) %>%
  summarise(mean = mean(total), sd = sd(total)) 

c_3 <- c(c_3[1,1], c_3[1,2])
c_4 <- dat.init$discard_data %>% tail(n = 3) %>% select(Discard) %>% range()

#mean proportion in each age group
mpatage <- dat.init$agecomp %>% 
  filter(FltSvy == 3) %>% 
  select(-c(1:9)) %>% 
  pivot_longer(everything(), 
               names_to = c("p", "a"), 
               names_pattern = "(.)(.)") %>% 
  select(-p) %>% 
  group_by(a) %>% 
  summarise(mp = mean(value)) %>% 
  filter(a > 1 & a < 8) %>% 
  pull(mp)
###Natural mortality
M <- report.$M_at_age %>%
  filter(Year == start.year) %>%
  select(-c(Bio_Pattern, Gender, Year))
M <- apply(M, 2, rep, 101)
M[, 15] <- M[, 14]

#Selectivities for each fleet by age
#Used selectivity for post commercial fleets post IFQ
sel <- report.$ageselex %>%
  filter(Factor == "Asel") %>%
  filter(Yr == start.year) %>%
  select(-c(Factor, Fleet, Yr, Seas, Sex, Morph, Label))
asel <- sel[c(8, 9, 3, 4, 5), ]
asel[4, 2] <- .75
#FI surveys are length based selectivity
lsel <- report.$sizeselex %>%
  filter(Fleet > 10) %>%
  filter(Yr == start.year) %>%
  select(-c(Factor, Fleet, Yr, Sex, Label))

#Get age-length key and flip bc it is inverted
ALK <- as.data.frame(report.$ALK)
ALK <- ALK[, c(1:15)]
ALK <- apply(ALK, 2, rev)
trans.prob <- c(.011, 4, 2)
ageerror <- report.$age_error_sd[-1, 2]

q_order <- c(8,9,3,4,11,12)
q <- report.$index_variance_tuning_check %>%
  slice(match(q_order, Fleet))%>%
  select(Q) %>%
  pull() %>%
  as.numeric()

SE <- c(0.369, 0.138, 0.200, 0.2, 0.200, 0.200)
se_log <- matrix(data = NA, nrow = Nyrs, ncol = Nfleet)
CPUE.se <- report.$cpue %>%
  group_by(Fleet) %>%
  select(Fleet, SE) %>%
  filter(Fleet == 8 |
           Fleet == 9 |
           Fleet == 3 | Fleet == 4 | Fleet == 11 | Fleet == 12) %>%
  summarise_all(mean)
CPUE.se <- CPUE.se[c(3, 4, 1, 2, 5, 6), ]

fecund <- report.$ageselex %>%
  filter(Factor == "Fecund") %>%
  filter(Yr == start.year) %>%
  select(-c(1:7))

#Number-at-age matrix
N <- data.frame(
  Year = seq(start.year, start.year + 50, by = 0.5),
  Zero = rep(NA, 101),
  One = rep(NA, 101),
  Two = rep(NA, 101),
  Three = rep(NA, 101),
  Four = rep(NA, 101),
  Five = rep(NA, 101),
  Six = rep(NA, 101),
  Seven = rep(NA, 101),
  Eight = rep(NA, 101),
  Nine = rep(NA, 101),
  Ten = rep(NA, 101),
  Eleven = rep(NA, 101),
  Tweleve = rep(NA, 101),
  Thirteen = rep(NA, 101),
  Fourteen = rep(NA, 101)
)
N[c(1), 2:16] <- report.$natage %>% filter(Time == 2014.0) %>% select(-c(1:11))

f.by.fleet <- c()

harvest.rate <- c()

catch.se <- c(report.$catch_error[1:3], .1, report.$catch_error[5])

f.list <- list()

z.list <- list()
#projected competition index
proj.index.full <-
  read.csv(file.path(dir., "TidyData/RS.biomass.report.csv"))

proj.index <-
  proj.index.full %>% select(Yr, RS_relative) %>% filter(Yr >= start.year & Yr <= start.year+102) %>% rename(Year = Yr)
proj.index$low <- NA
proj.index$med <- NA
proj.index$high <- NA

for(i in 1:nrow(proj.index)){
  proj.index$low[i]<- rnorm(1, proj.index$RS_relative[i], .005)
  proj.index$med[i]<- rnorm(1, proj.index$RS_relative[i], .01)
  proj.index$high[i]<- rnorm(1, proj.index$RS_relative[i], .05)
  
}

rs.scen = args[3]
col.num <- which(colnames(proj.index) == rs.scen)
fixed.comp.f <-  proj.index[,col.num] * report.$exploitation %>% filter(Yr == start.year) %>% select(COMP) %>% pull() 
comp.f <- proj.index[,col.num] * report.$exploitation %>% filter(Yr == start.year) %>% select(COMP) %>% pull()


smp <- data.frame(Year = seq(start.year, start.year+102), Seas = rep(1, 103), Fleet = rep(4, 103), f = rep(0.07356127, 103))
cmp <- data.frame(Year = seq(start.year, start.year+102), Seas = rep(1, 103), Fleet = rep(5, 103), f = fixed.comp.f)
full.forecast.f <- smp %>% arrange(Year)
#full.forecast.f <- bind_rows(smp, cmp) %>% arrange(Year)

##Catch dataframes
catch <- matrix(NA, nrow = 4, ncol = 15)
###Maybe try reintroducing this but change the shrimp bycatch value from -2
catch.err <- report.$catch_error[1:4]
catch.by.year <- matrix(data = NA, nrow = Nyrs, ncol = Nages)
catch.by.fleet <- matrix(data = NA, nrow = 4, ncol = Nages)
catch.fleet.year <- matrix(data = NA, nrow = Nyrs, ncol = 4)
.datcatch <- matrix(NA, nrow = 100, ncol = 4)
comp.discard <- matrix(NA, nrow = 100, ncol = Nages)
init.catch <- 1555
ind.w.mat <- read.csv(file.path(dir., "TidyData/eHCR.csv"), row.names = 1)
init.slopes <- dat.init$CPUE %>%
  filter(year >2009) %>%
  filter(index == 3 |
           index == 4 |
           index == 8 |
           index == 9 |
           index > 10) %>%
  mutate(index = factor(index,
                        levels = c(8, 9, 3, 4, 5, 11, 12))) %>%
  group_by(index) %>%
  nest() %>%
  mutate(model = map(data, ~ lm(obs ~ year, data = .) %>%
                       coef %>%
                       as.list %>%
                       as_tibble())) %>%
  unnest(model) %>%
  select(index, year) %>% mutate(year = year) %>% 
  ungroup() %>%
  mutate(index = recode_factor(index,
                             `8` = 7,
                             `9` = 8,
                             `3` = 3,
                             `4` = 4,
                             `11` =10,
                             `12` = 11))


##Numbers at length matrix
cvdist <-
  rnorm(1000, 0, .2535) #create a sample of error for each length with mean = CV of growth .2535 and sd .1
length.list <- list(vector("list", 16))
lbins <- c(report.$lbins, 100)
natlen <- matrix(data = NA, nrow = ncol(ALK), ncol = nrow(ALK))
Nlen <- matrix(data = NA, nrow = Nyrs, ncol = nrow(ALK))

##Age and length comp data
agecomp.list <- list()
len.comps.list <- list()
survey.len <- matrix(NA, nrow = 2, ncol = 12)

##Biomass matrices and set up for IOA (index matrix and qs)
b.age <- matrix(data = NA, nrow = 3, ncol = Nages)
b.len <- matrix(data = NA, nrow = 2, ncol = 12)
Index.list <- list()
I <- matrix(data = NA, nrow = Nyrs, ncol = 6)

tot.catch.b <- c()
ssb_tot <- c()
b.list <- list()
hist.catch.list <- list()

OFL <- list()
f.reduced <- c()
dat.list <- list(
  Nages = Nages,
  age_selectivity = asel,
  length_selectivity = lsel,
  M = M,
  CPUE_se = CPUE.se,
  catch_se = catch.se,
  wtatage = wtatage,
  N_fishfleet = 3,
  N_survey = 2,
  ageerror = ageerror,
  q = q,
  year_seq = seq(start.year, start.year+50, by = 0.5),
  N_areas = 1,
  N_totalfleet = 4,
  catch_proportions = c(0.35, 0.17, 0.48),
  RS_projections = proj.index[,c(1:2)],
  full_forecast = full.forecast.f,
  files.to.copy = list("forecast.ss",
                       "starter.ss",
                       "VS.dat",
                       "VS.ctl",
                       "Report.sso",
                       "Forecast-report.sso",
                       "ss3.PAR",
                       "CompReport.sso",
                       "covar.sso",
                       "wtatage.ss_new")
  
)
dat.list$RS_projections$RS_relative <- fixed.comp.f

## Reference points
ref.points <- list()
jit <- final_try <- FALSE
rebuild.mat <- c()
ref.points <- list()
iteration = as.numeric(args[1])
T.target <- 0
recommend_catch = F
x <- 1
scenario = args[2]
set.seed(seed[iteration])
escen = args[4]

if(length(which(str_detect(list.files(dir.), "Assessments_2") == TRUE)) == 0){ dir.create(file.path(dir., "Assessments_2"))}
dir.it <- paste0(dir., "/Assessments_2/", scenario, "/", iteration)
dir.create(file.path(dir., "Assessments_2", scenario))
dir.create(file.path(dir., "Assessments_2", scenario, paste(iteration)))
files. <- list("VS.dat", "VS.ctl", "forecast.ss", "starter.ss")
file.copy(file.path(dir.,"/initial_pop/initial_no_comp", files.), file.path(dir.it), overwrite = T)
ss3.file <- list("SS3") 
file.copy(file.path(dir., ss3.file), file.path(dir.it))
#unlink(here("one_plus", "Assessments"), recursive = T)
writeLines(c(rs.scen, scenario), paste0("one_plus/Assessments_2/", scenario, "/", iteration, "/README.txt"))


for (year in Year.vec) {
  
  if (recommend_catch == F) {
    .datcatch[year,1] <- rtruncnorm(1, a = 0, b = Inf, c_1[1], c_1[2])
    .datcatch[year,2] <- rtruncnorm(1, a = 0, b = Inf, c_2[1], c_1[2])
    .datcatch[year,3] <- rtruncnorm(1, a = 0, b = Inf, c_3[1], c_3[2])
    .datcatch[year,4] <- runif(1, c_4[1], c_4[2])
    
    f.by.fleet[1] <-
      H_rate(.datcatch[year,1], 1, N, year, dat.list, F)
    f.by.fleet[2] <-
      H_rate(.datcatch[year,2], 2, N, year, dat.list, F)
    f.by.fleet[3] <-
      H_rate(.datcatch[year,3], 3, N, year, dat.list, F)
    f.by.fleet[4] <-
      H_rate(.datcatch[year,4], 4, N, year, dat.list, F, F)
    f.by.fleet[5] <- comp.f[year/2]
    
  }
  
  if (recommend_catch == T) {
    
    .datcatch[year,1] <- rtruncnorm(1, a = 0, b = Inf, mean = c_1[x], 341.489)
    .datcatch[year,2] <- rtruncnorm(1, a = 0, b = Inf, mean = c_2[x], 60.46)
    .datcatch[year,3] <- rtruncnorm(1, a = 0, b = Inf, mean = c_3[x], 126.1625)
    .datcatch[year,4] <- runif(1, c_4[1], c_4[2])
    
    f.by.fleet[1] <-
      H_rate(.datcatch[year,1], 1, N, year, dat.list, F)
    f.by.fleet[2] <-
      H_rate(.datcatch[year,2], 2, N, year, dat.list, F)
    f.by.fleet[3] <-
      H_rate(.datcatch[year,3], 3, N, year, dat.list, F)
    f.by.fleet[4] <-
      H_rate(.datcatch[year,4], 4, N, year, dat.list, F, F)
    f.by.fleet[5] <- comp.f[year/2]
  }
  
  ##SSB caluclated at beginning of each year
  ssb_tot[year - 1] <- sum(N[year - 1, -1] * fecund)
  
  f.list[[year]] <- f.by.fleet
  
  # z <- zatage(dat.list, year, f.by.fleet)
  # z.list[[year]] <- z
  
  #Catches
  catch.by.fleet[,] <- exploit_catch(f.by.fleet, N, year, dat.list)
  ncatch.by.fleet <- rbind(num_to_bio(catch.by.fleet[c(1,2),], dat.list, Natage = F, age0 = T), catch.by.fleet[c(3:4),])
  catch.by.year[year,] <- colSums(ncatch.by.fleet)  #catch from 2014.0 N
  catch.fleet.year[year,] <- rowSums(catch.by.fleet)
  
  hist.catch.list[[year]] <- catch.by.fleet
  #.datcatch[year,5] <- sum(catch.by.fleet[5,])
  
  ### simulate age comps for fleets 1 - 3
  if (sum(catch.fleet.year[year, c(1:3)]) > 0.001) {
    agecomp.list[[year]] <- simAgecomp(catch.by.fleet, year, dat.list)
  }
  
  ## N-at-age for the next year
  #Age-0 year.0
  N[year + 1, 2] <- sample(rec.dis, 1, replace = F)
  
  #Age 1 to 14
  for (age in 3:ncol(N)) {
    if (age < ncol(N)) {
      N[year + 1, age] <- N[year-1,age-1]*exp(-M[year-1,age-2])-catch.by.year[year,age-2]
    } else{
      N[year + 1, age] <- N[year-1,age-1]*exp(-M[year-1,age-2])-catch.by.year[year,age-2] + N[year-1,age]*exp(-M[year-1,age-1])-catch.by.year[year,age-1] #plus age group
    }
  }
  
  ## Numbers in mid-year
  for(age in 3:ncol(N)){
    N[year, age-1] <- N[year-1,age-1] - (N[year-1,age-1]-N[year+1,age])/2
    if(age == 16){
      N[year, age-1] <- N[year-1,age-1] - (N[year-1,age-1]*exp(-M[year-1,age-2])-catch.by.year[year,age-2])/2
      N[year, age] <- N[year-1,age] - (N[year-1,age]*exp(-M[year-1,age-1])-catch.by.year[year,age-1])/2
    }
  }
  
  ## Proportion of numbers in each age and length bins
  
  for (i in 1:ncol(ALK)) {
    natlen[i, ] <- N[year - 1, i + 1] * ALK[, i]
  }
  
  Nlen[year - 1, ] <- colSums(natlen)
  
  new.num <- Nlen[year - 1, c(1, 2, 3)] * trans.prob
  Nlen[year, ] <- c(new.num, Nlen[year - 1, -c(1:3)])
  
  
  #### IOA ####
  
  ## Calculate biomass avaiable for surveys
  
  b.age <- simBatage(N, dat.list, year)
  b.len <- simBatlen(Nlen, dat.list, year)
  byc.bio <- rlnorm(1, E_4_mu, E_4_sd)
  
  b <- c(apply(b.age, 1, sum), byc.bio, apply(b.len, 1, sum))
  b.list[[year]] <- sum(N[year, -1] * wtatage[1, ])
  I[year, ] <- simIndex(dat.list, b)
  
  x <- x+1

  #### Add data into .dat file and run assessment ####
  if (year %% 10 == 0) {
    dat. <- SS_readdat(paste0(dir.it, "/VS.dat"), version = "3.24")
    dat.update(year = year,
               dat.list = dat.list,
               dat. = dat.,
               agecomp.list = agecomp.list,
               I = I,
               .datcatch = .datcatch,
               comp.I = NULL,
               dir. = dir.it,
               write = T)
    dat. <- SS_readdat(paste0(dir.it, "/VS.dat"), version = "3.24")
    #check stock status compared to ref points.
    yr.end <- floor(dat.list$year_seq[year])
    yr.start <- floor(dat.list$year_seq[year-9])
    
    bratio <- rowSums(N[year,-1]*wtatage[1,])/B0
    if(bratio < 0.3){
      OFL[[year]] <- 0
    }else{
      ind.slopes <- dat.$CPUE %>%
        filter(year >= yr.start & year <= yr.end) %>%
        filter(index == 3 |
                 index == 4 |
                 index == 7 |
                 index == 8 |
                 index > 9) %>%
        mutate(index = factor(index,
                              levels = c(7, 8, 3, 4, 10, 11))) %>%
        group_by(index) %>%
        nest() %>%
        mutate(model = map(data, ~ lm(obs ~ year, data = .) %>%
                             coef %>%
                             as.list %>%
                             as_tibble())) %>%
        unnest(model) %>%
        select(index, year) %>%
        arrange(index)
      
      rows <- seq(year-10,year)
      
      past.catch <- .datcatch %>%
        as.data.frame() %>%
        slice(rows)  %>%
        na.omit() %>%
        rowSums() %>%
        mean()
      
      OFL[[year]] <- past.catch * ((ind.w.mat[escen,1] * 1+ind.slopes$year[1]) +
                                            (ind.w.mat[escen,2] * 1+ind.slopes$year[2]) +
                                            (ind.w.mat[escen,3] * 1+ind.slopes$year[3]) +
                                            (ind.w.mat[escen,4] * 1+ind.slopes$year[4]) +
                                            #(ind.w.mat[escen,5] * 1+ind.slopes$year[5]) +
                                            (ind.w.mat[escen,6] * 1+ind.slopes$year[5]) +
                                            (ind.w.mat[escen,7] * 1+ind.slopes$year[6]))
    }
    
    # if(OFL[[year]] > 1555){
    #   OFL[[year]] <- 1555
    # }
    #calculate slopes for each fleet based on previous 5 years
    
    c_1[x:(x+5)] <- OFL[[year]]  * dat.list$catch_proportions[1]
    c_2[x:(x+5)] <- OFL[[year]]  * dat.list$catch_proportions[2]
    c_3[x:(x+5)] <- OFL[[year]]  * dat.list$catch_proportions[3]
    recommend_catch <- T
    
    ss <- run_ss(dir.it, lin = TRUE)
    
    if (isFALSE(ss)) {
      rep.file <-
        MO_SSoutput(dir = dir.it,
                    verbose = F,
                    printstats = F)
    } else{
      
      message("The model did not converge.")
      jit <- try_jit(dir.it, lin = TRUE)
      
      if (isFALSE(jit)) {
        rep.file <-
          MO_SSoutput(dir = dir.it,
                      verbose = F,
                      printstats = F)
      } else{
        starter <- SS_readstarter(paste0(dir.it, "/starter.ss"))
        starter$init_values_src <- 0
        SS_writestarter(starter,
                        dir = dir.it,
                        file = "starter.ss",
                        overwrite = T)
        final_try <- run_ss(dir.it, lin = TRUE)
      }
      if (isTRUE(final_try)) {
        stop("Failed to converge")
      } else{
        rep.file <-
          MO_SSoutput(dir = dir.it,
                      verbose = F,
                      printstats = F)
      }
      
    }
    
    spr <- rep.file$derived_quants %>%
      filter(str_detect(Label, "Bratio")) %>%
      slice(tail(row_number(), 10)) %>%
      summarise(mean(Value)) %>%
      pull() %>% round(3)
    
    
    if (spr >= .299 & spr <=.31) {
      opt = "A"
    } else{
      opt = "B"
    }
    
    if (opt == "B") {
      find_spr(dir.it, dat.list$N_totalfleet, notifications = F, lin = TRUE)
      
    }
    
    rep.file <- MO_SSoutput(dir = dir.it,
                            verbose = F,
                            printstats = F)
    ref.points[[year]] <- getRP(rep.file, dat.list, year)
    
    #check stock status compared to ref points.
    MSST_rel <- ref.points[[year]]$status_cur
    
    copy_files(dat.list, dir.it, paste0(dir.it, "/", dat.list$year_seq[year]))
    
    save.list <- list(N = N,
                      I= I,
                      catch.data = catch.fleet.year,
                      F. = do.call(rbind, compact(f.list)),
                      OFL = do.call(rbind, compact(OFL)),
                      ref.points = do.call(rbind, compact(ref.points)))
    
    for(i in names(save.list)){
      write.csv(save.list[[i]], paste0(dir.it, "/", i, ".csv"))
    }
    start <- SS_readstarter(file = paste0(dir.it, "/starter.ss"))
    start$init_values_src <- 0
    SS_writestarter(start, dir = dir.it, overwrite = T)
    rm(rep.file, dat.)
    jit <- final_try <- FALSE
    
  }
}
