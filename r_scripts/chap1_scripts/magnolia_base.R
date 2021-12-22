### No competition OM
# install.packages(c("tidyverse", "r4ss", "psych", "truncnorm"))
# devtools::install_github("mcoshima/moMSE")
library(moMSE)
library(r4ss)
library(dplyr)
#library(here)
library(ggplot2)
library(stringr)
#library(RPushbullet)
library(psych) #for geometric mean
library(tidyr)
library(purrr)
library(truncnorm)

args = commandArgs(trailingOnly = TRUE)
setwd('..')
dir. <- paste0(getwd(), "/one_plus")
report. <- MO_SSoutput(dir = file.path(dir., "initial_pop", "initial_no_comp"))

dat.init <- SS_readdat(file.path(dir., "initial_pop", "initial_no_comp", "/VS.dat"))
seed <- read.csv("setseed.csv")
seed <- as.vector(seed[,2])

#Data and parameter setup
Nyrs <- 100 #Double the number of projection years
Year.vec <- seq(2, 100, by = 2)
year.seq <- seq(2014, 2014 + 50, by = .5) #for function
Nfleet <- 6 #ComEIFQ/ComWIFQ/Rec/ShrimpBC/Video/Groundfish
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

#Number-at-age matrix

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
###Catch and effort

E_4 <-
  report.$discard %>% filter(Yr > 2011 &
                               Yr < 2015) %>% select(F_rate) %>% pull()
E_4_mu <- mean(log(E_4))
E_4_sd <- sd(log(E_4))

c_1 <- dat.init$catch %>% filter(year > 2012) %>% select(c(1)) %>% range()
c_2 <- dat.init$catch %>% filter(year > 2006) %>% select(c(2)) %>% range() 
c_3 <- dat.init$catch %>% filter(year > 2010) %>% select(c(3)) %>% range() 
c_4 <- dat.init$discard_data %>% tail(n = 3) %>% select(Discard) %>% range()

###Natural mortality
M <- report.$M_at_age %>%
  filter(Year == 2013) %>%
  select(-c(Bio_Pattern, Gender, Year))
M <- apply(M, 2, rep, 101)
M[, 15] <- M[, 14]

#Selectivities for each fleet by age
#Used selectivity for post commercial fleets post IFQ
sel <- report.$ageselex %>%
  filter(Factor == "Asel") %>%
  filter(Yr == 2014) %>%
  select(-c(Factor, Fleet, Yr, Seas, Sex, Morph, Label))
asel <- sel[c(7, 8, 3, 4), ]
asel[4, 2] <- .75
#FI surveys are length based selectivity
lsel <- report.$sizeselex %>%
  filter(Fleet > 9) %>%
  filter(Yr == 2014) %>%
  select(-c(Factor, Fleet, Yr, Sex, Label))



#Get age-length key and flip bc it is inverted
ALK <- as.data.frame(report.$ALK)
ALK <- ALK[, c(1:15)]
ALK <- apply(ALK, 2, rev)
trans.prob <- c(.011, 4, 2)
ageerror <- report.$age_error_sd[-1, 2]


# F for competition is going to be based on the RS abundance index
#for(iteration in 13:15){

N <- data.frame(
  Year = seq(2014, 2064, by = 0.5),
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

catch.se <- c(report.$catch_error[1:3], .1)

f.list <- list()

z.list <- list()

#Shrimp bycatch forecast
smp <- data.frame(Year = seq(2014, 2116), Seas = rep(1, 103), Fleet = rep(4, 103), f = rep(0.07356127, 103))
full.forecast.f <- smp

##Catch dataframes
catch <- matrix(NA, nrow = 4, ncol = 15)
###Maybe try reintroducing this but change the shrimp bycatch value from -2
catch.err <- report.$catch_error[1:4]
catch.by.year <- matrix(data = NA, nrow = Nyrs, ncol = Nages)
catch.by.fleet <- matrix(data = NA, nrow = 4, ncol = Nages)
catch.fleet.year <- matrix(data = NA, nrow = Nyrs, ncol = 4)
.datcatch <- matrix(NA, nrow = 100, ncol = 4)
init.catch <- 1555

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
q_order <- c(7,8,3,4,10,11)
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
  filter(Fleet == 7 |
           Fleet == 8 |
           Fleet == 3 | 
           Fleet == 4 | 
           Fleet == 10 |
           Fleet == 11) %>%
  summarise_all(mean)
CPUE.se <- CPUE.se[c(3, 4, 1, 2, 5, 6), ]

fecund <- report.$ageselex %>%
  filter(Factor == "Fecund") %>%
  filter(Yr == 2014) %>%
  select(-c(1:7))
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
  year_seq = seq(2014, 2064, by = 0.5),
  N_areas = 1,
  N_totalfleet = 4,
  catch_proportions = c(0.35, 0.17, 0.48), 
  full_forecast = full.forecast.f,
  files.to.copy = list("forecast.ss",
                       "starter.ss",
                       "VS.dat",
                       "VS.ctl",
                       "Report.sso",
                       "Forecast-report.sso",
                       "SS3.par",
                       "CompReport.sso",
                       "covar.sso",
                       "wtatage.ss_new")
  
)


rebuild.mat <- c()
## Reference points
ref.points <- list()
iteration = as.numeric(args[1])
T.target <- 0
recommend_catch = F
x <- 1
scenario = args[2]
rs.scen = args[3]
set.seed(seed[iteration])

##Prep files
dir.it <- paste0(dir., "/Assessments/", scenario, "/", iteration)
dir.create(file.path(dir., "Assessments", scenario))
dir.create(file.path(dir., "Assessments", scenario, paste(iteration)))
files. <- list("VS.dat", "VS.ctl", "forecast.ss", "starter.ss")
file.copy(file.path(dir.,"/initial_pop/initial_no_comp", files.), file.path(dir.it), overwrite = T)
ss3.file <- list("SS3")
file.copy(file.path(dir., ss3.file), file.path(dir.it))
#unlink(here("one_plus", "Assessments"), recursive = T)
writeLines(c(rs.scen, scenario), paste0("one_plus/Assessments/", scenario, "/", iteration, "/README.txt"))

for (year in Year.vec[1:50]) {
  
  if (recommend_catch == F) {
    .datcatch[year,1] <- rtruncnorm(1, a = 0, b = Inf, 633.09, 341.489)
    .datcatch[year,2] <- rtruncnorm(1, a = 0, b = Inf, 335.9493, 60.46)
    .datcatch[year,3] <- rtruncnorm(1, a = 0, b = Inf, 674.75, 126.1625)
    .datcatch[year,4] <- runif(1, c_4[1], c_4[2])
    
    f.by.fleet[1] <-
      H_rate(.datcatch[year,1], 1, N, year, dat.list, F)
    f.by.fleet[2] <-
      H_rate(.datcatch[year,2], 2, N, year, dat.list, F)
    f.by.fleet[3] <-
      H_rate(.datcatch[year,3], 3, N, year, dat.list, F)
    f.by.fleet[4] <-
      H_rate(.datcatch[year,4], 4, N, year, dat.list, F, F)
    
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
    
  }
  
  ##SSB caluclated at beginning of each year
  ssb_tot[year - 1] <- sum(N[year - 1, -1] * fecund)
  
  f.list[[year]] <- f.by.fleet
  
  #Catches
  catch.by.fleet[,] <- exploit_catch(f.by.fleet, N, year, dat.list)
  ncatch.by.fleet <- rbind(num_to_bio(catch.by.fleet[c(1,2),], dat.list, Natage = F, age0 = T), catch.by.fleet[c(3:4),])
  catch.by.year[year,] <- colSums(ncatch.by.fleet)  #catch from 2014.0 N
  catch.fleet.year[year,] <- rowSums(catch.by.fleet)
  
  hist.catch.list[[year]] <- catch.by.fleet
  
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
               .datcatch = catch.fleet.year,
               dir. = dir.it,
               write = T)
    
    ss <- run_ss(dir.it, lin = TRUE)
    
    converge.params <- check_convergence(dir.it)
    
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
    
    if (MSST_rel < 1) {
      rebuild = T
      rebuild.mat[year / 10] <- 1
    } else{
      rebuild = F
    }
    
    if (rebuild == T & year > x+T.target) {
      T.target <-
        rebuild_ttarg(dir.it, dat.list)
      OFL[[year]] <-
        rebuild_f(year, dir.it, dat.list, T.target)
      P.star <- p_star()
      ABC <- OFL[[year]]$catch * P.star
      
      c_1[x:(x+T.target-1)] <- ABC * dat.list$catch_proportions[1] 
      c_2[x:(x+T.target-1)] <- ABC * dat.list$catch_proportions[2]
      c_3[x:(x+T.target-1)] <- ABC * dat.list$catch_proportions[3]
      
      recommend_catch <- T
      reduced_catch <- T
    }
    
    if (rebuild == F) {
      rebuild.mat[year / 10] <- 0
      OFL[[year]] <- 1555
      ABC <- OFL[[year]] * .75
      c_1[x:(x+5)] <- ABC * dat.list$catch_proportions[1]
      c_2[x:(x+5)] <- ABC * dat.list$catch_proportions[2]
      c_3[x:(x+5)] <- ABC * dat.list$catch_proportions[3]
      
      recommend_catch <- T
      reduced_catch <- F
    }
    
    copy_files(dat.list, dir.it, paste0(dir.it, "/", dat.list$year_seq[year]))
    
    save.list <- list(N = N, 
                      I= I, 
                      catch.data = catch.fleet.year, 
                      F. = do.call(rbind, compact(f.list)), 
                      OFL = do.call(rbind, compact(OFL)), 
                      ref.points = do.call(rbind, compact(ref.points)),
                      rebuild.matrix = rebuild.mat)
    
    for(i in names(save.list)){
      write.csv(save.list[[i]], paste0(dir.it, "/", i, ".csv"))
    }
    start <- SS_readstarter(file = paste0(dir.it, "/starter.ss"))
    start$init_values_src <- 0
    SS_writestarter(start, dir = dir.it, overwrite = T)
    rm(rep.file, dat.)
    jit <- FALSE
    final_try <- FALSE
  }
}
#}
