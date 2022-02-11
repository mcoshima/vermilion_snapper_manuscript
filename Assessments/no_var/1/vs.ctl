#V3.30.xx.yy;_safe;_compile_date:_Dec  7 2021;_Stock_Synthesis_by_Richard_Methot_(NOAA)_using_ADMB_12.3
#_Stock_Synthesis_is_a_work_of_the_U.S._Government_and_is_not_subject_to_copyright_protection_in_the_United_States.
#_Foreign_copyrights_may_apply._See_copyright.txt_for_more_information.
#_User_support_available_at:NMFS.Stock.Synthesis@noaa.gov
#_User_info_available_at:https://vlab.noaa.gov/group/stock-synthesis
#_Source_code_at:_https://github.com/nmfs-stock-synthesis/stock-synthesis

#_data_and_control_files: vs_envM2.dat // vs2.ctl
0  # 0 means do not read wtatage.ss; 1 means read and use wtatage.ss and also read and use growth parameters
1  #_N_Growth_Patterns (Growth Patterns, Morphs, Bio Patterns, GP are terms used interchangeably in SS3)
1 #_N_platoons_Within_GrowthPattern 
#_Cond 1 #_Platoon_within/between_stdev_ratio (no read if N_platoons=1)
#_Cond  1 #vector_platoon_dist_(-1_in_first_val_gives_normal_approx)
#
4 # recr_dist_method for parameters:  2=main effects for GP, Area, Settle timing; 3=each Settle entity; 4=none (only when N_GP*Nsettle*pop==1)
1 # not yet implemented; Future usage: Spawner-Recruitment: 1=global; 2=by area
1 #  number of recruitment settlement assignments 
0 # unused option
#GPattern month  area  age (for each settlement assignment)
 1 1 1 0
#
#_Cond 0 # N_movement_definitions goes here if Nareas > 1
#_Cond 1.0 # first age that moves (real age at begin of season, not integer) also cond on do_migration>0
#_Cond 1 1 1 2 4 10 # example move definition for seas=1, morph=1, source=1 dest=2, age1=4, age2=10
#
3 #_Nblock_Patterns
 3 4 1 #_blocks_per_pattern 
# begin and end years of blocks
 1990 2004 2005 2007 2008 2017
 1990 1996 1997 2004 2005 2007 2008 2017
 2008 2016
#
# controls for all timevary parameters 
1 #_time-vary parm bound check (1=warn relative to base parm bounds; 3=no bound check); Also see env (3) and dev (5) options to constrain with base bounds
#
# AUTOGEN
 1 1 1 1 1 # autogen: 1st element for biology, 2nd for SR, 3rd for Q, 4th reserved, 5th for selex
# where: 0 = autogen time-varying parms of this category; 1 = read each time-varying parm line; 2 = read then autogen if parm min==-12345
#
#_Available timevary codes
#_Block types: 0: P_block=P_base*exp(TVP); 1: P_block=P_base+TVP; 2: P_block=TVP; 3: P_block=P_block(-1) + TVP
#_Block_trends: -1: trend bounded by base parm min-max and parms in transformed units (beware); -2: endtrend and infl_year direct values; -3: end and infl as fraction of base range
#_EnvLinks:  1: P(y)=P_base*exp(TVP*env(y));  2: P(y)=P_base+TVP*env(y);  3: P(y)=f(TVP,env_Zscore) w/ logit to stay in min-max;  4: P(y)=2.0/(1.0+exp(-TVP1*env(y) - TVP2))
#_DevLinks:  1: P(y)*=exp(dev(y)*dev_se;  2: P(y)+=dev(y)*dev_se;  3: random walk;  4: zero-reverting random walk with rho;  5: like 4 with logit transform to stay in base min-max
#_DevLinks(more):  21-25 keep last dev for rest of years
#
#_Prior_codes:  0=none; 6=normal; 1=symmetric beta; 2=CASAL's beta; 3=lognormal; 4=lognormal with biascorr; 5=gamma
#
# setup for M, growth, wt-len, maturity, fecundity, (hermaphro), recr_distr, cohort_grow, (movement), (age error), (catch_mult), sex ratio 
#_NATMORT
1 #_natM_type:_0=1Parm; 1=N_breakpoints;_2=Lorenzen;_3=agespecific;_4=agespec_withseasinterpolate;_5=BETA:_Maunder_link_to_maturity
15 #_N_breakpoints
 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 # age(real) at M breakpoints
#
1 # GrowthModel: 1=vonBert with L1&L2; 2=Richards with L1&L2; 3=age_specific_K_incr; 4=age_specific_K_decr; 5=age_specific_K_each; 6=NA; 7=NA; 8=growth cessation
0.5 #_Age(post-settlement)_for_L1;linear growth below this
999 #_Growth_Age_for_L2 (999 to use as Linf)
-999 #_exponential decay for growth above maxage (value should approx initial Z; -999 replicates 3.24; -998 to not allow growth above maxage)
0  #_placeholder for future growth feature
#
0 #_SD_add_to_LAA (set to 0.1 for SS2 V1.x compatibility)
1 #_CV_Growth_Pattern:  0 CV=f(LAA); 1 CV=F(A); 2 SD=F(LAA); 3 SD=F(A); 4 logSD=F(A)
#
1 #_maturity_option:  1=length logistic; 2=age logistic; 3=read age-maturity matrix by growth_pattern; 4=read age-fecundity; 5=disabled; 6=read length-maturity
1 #_First_Mature_Age
2 #_fecundity option:(1)eggs=Wt*(a+b*Wt);(2)eggs=a*L^b;(3)eggs=a*Wt^b; (4)eggs=a+b*L; (5)eggs=a+b*W
0 #_hermaphroditism option:  0=none; 1=female-to-male age-specific fxn; -1=male-to-female age-specific fxn
1 #_parameter_offset_approach for M, G, CV_G:  1- direct, no offset**; 2- male=fem_parm*exp(male_parm); 3: male=female*exp(parm) then old=young*exp(parm)
#_** in option 1, any male parameter with value = 0.0 and phase <0 is set equal to female parameter
#
#_growth_parms
#_ LO HI INIT PRIOR PR_SD PR_type PHASE env_var&link dev_link dev_minyr dev_maxyr dev_PH Block Block_Fxn
# Sex: 1  BioPattern: 1  NatMort
 0.0001 1e+06 0.234 0.2 0 0 -1 0 0 0 0 0 0 0 # NatM_break_1_Fem_GP_1
 0.0001 1e+06 0.342 0.2 0 0 -1 0 0 0 0 0 0 0 # NatM_break_2_Fem_GP_1
 0.0001 1e+06 0.287 0.2 0 0 -1 0 0 0 0 0 0 0 # NatM_break_3_Fem_GP_1
 0.0001 1e+06 0.257 0.2 0 0 -1 0 0 0 0 0 0 0 # NatM_break_4_Fem_GP_1
 0.0001 1e+06 0.239 0.2 0 0 -1 0 0 0 0 0 0 0 # NatM_break_5_Fem_GP_1
 0.0001 1e+06 0.228 0.2 0 0 -1 0 0 0 0 0 0 0 # NatM_break_6_Fem_GP_1
 0.0001 1e+06 0.22 0.2 0 0 -1 0 0 0 0 0 0 0 # NatM_break_7_Fem_GP_1
 0.0001 1e+06 0.215 0.2 0 0 -1 0 0 0 0 0 0 0 # NatM_break_8_Fem_GP_1
 0.0001 1e+06 0.212 0.2 0 0 -1 0 0 0 0 0 0 0 # NatM_break_9_Fem_GP_1
 0.0001 1e+06 0.209 0.2 0 0 -1 0 0 0 0 0 0 0 # NatM_break_10_Fem_GP_1
 0.0001 1e+06 0.207 0.2 0 0 -1 0 0 0 0 0 0 0 # NatM_break_11_Fem_GP_1
 0.0001 1e+06 0.206 0.2 0 0 -1 0 0 0 0 0 0 0 # NatM_break_12_Fem_GP_1
 0.0001 1e+06 0.205 0.2 0 0 -1 0 0 0 0 0 0 0 # NatM_break_13_Fem_GP_1
 0.0001 1e+06 0.204 0.2 0 0 -1 0 0 0 0 0 0 0 # NatM_break_14_Fem_GP_1
 0.0001 1e+06 0.204 0.2 0 0 -1 0 0 0 0 0 0 0 # NatM_break_15_Fem_GP_1
# Sex: 1  BioPattern: 1  Growth
 0.0001 1e+06 11.83 11.83 0 0 -1 0 0 0 0 0 0 0 # L_at_Amin_Fem_GP_1
 0.0001 1e+06 34.4 34.4 0 0 -1 0 0 0 0 0 0 0 # L_at_Amax_Fem_GP_1
 0 1e+06 0.3254 0.3254 0 0 -1 0 0 0 0 0 0 0 # VonBert_K_Fem_GP_1
 0 1e+06 0.2535 0.0001 0 0 -1 0 0 0 0 0 0 0 # CV_young_Fem_GP_1
 0 1e+06 0.2535 0.0001 0 0 -1 0 0 0 0 0 0 0 # CV_old_Fem_GP_1
# Sex: 1  BioPattern: 1  WtLen
 0 1e+06 2.19e-05 2.19e-05 0 0 -1 0 0 0 0 0 0 0 # Wtlen_1_Fem_GP_1
 0 1e+06 2.916 2.916 0 0 -1 0 0 0 0 0 0 0 # Wtlen_2_Fem_GP_1
# Sex: 1  BioPattern: 1  Maturity&Fecundity
 0 1e+06 14.087 14.087 0 0 -1 0 0 0 0 0 0 0 # Mat50%_Fem_GP_1
 -1 1e+06 -0.574 -0.574 0 0 -1 0 0 0 0 0 0 0 # Mat_slope_Fem_GP_1
 0 1e+06 278.715 278.715 0 0 -1 0 0 0 0 0 0 0 # Eggs_scalar_Fem_GP_1
 0 1e+06 3.042 3.042 0 0 -1 0 0 0 0 0 0 0 # Eggs_exp_len_Fem_GP_1
# Hermaphroditism
#  Recruitment Distribution  
#  Cohort growth dev base
 0.1 10 1 1 1 0 -1 0 0 0 0 0 0 0 # CohortGrowDev
#  Movement
#  Age Error from parameters
#  catch multiplier
#  fraction female, by GP
 1e-06 0.999999 0.5 0.5 0.5 0 -1 0 0 0 0 0 0 0 # FracFemale_GP_1
#  M2 parameter for each predator fleet
 -0.2 0.5 0.0565212 0 0 0 3 201 0 0 0 0 0 0 # M2_pred1
#
# timevary MG parameters 
#_ LO HI INIT PRIOR PR_SD PR_type  PHASE
 0 1 0.0510646 0 0 0 3 # M2_pred1_ENV_add
# info on dev vectors created for MGparms are reported with other devs after tag parameter section 
#
#_seasonal_effects_on_biology_parms
 0 0 0 0 0 0 0 0 0 0 #_femwtlen1,femwtlen2,mat1,mat2,fec1,fec2,Malewtlen1,malewtlen2,L1,K
#_ LO HI INIT PRIOR PR_SD PR_type PHASE
#_Cond -2 2 0 0 -1 99 -2 #_placeholder when no seasonal MG parameters
#
3 #_Spawner-Recruitment; Options: 1=NA; 2=Ricker; 3=std_B-H; 4=SCAA; 5=Hockey; 6=B-H_flattop; 7=survival_3Parm; 8=Shepherd_3Parm; 9=RickerPower_3parm
1  # 0/1 to use steepness in initial equ recruitment calculation
0  #  future feature:  0/1 to make realized sigmaR a function of SR curvature
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type      PHASE    env-var    use_dev   dev_mnyr   dev_mxyr     dev_PH      Block    Blk_Fxn #  parm_name
             0         13.82       10.6723          6.91             0             0          1          0          0          0          0          0          0          0 # SR_LN(R0)
          0.22          0.96      0.473371           0.6          0.74             0          2          0          0          0          0          0          0          0 # SR_BH_steep
             0             2           0.3           0.2             0             0         -3          0          0          0          0          0          0          0 # SR_sigmaR
            -5             5             0             0             0             0         -3          0          0          0          0          0          0          0 # SR_regime
             0           0.5             0             0             0             0         -2          0          0          0          0          0          0          0 # SR_autocorr
#_no timevary SR parameters
1 #do_recdev:  0=none; 1=devvector (R=F(SSB)+dev); 2=deviations (R=F(SSB)+dev); 3=deviations (R=R0*dev; dev2=R-f(SSB)); 4=like 3 with sum(dev2) adding penalty
1994 # first year of main recr_devs; early devs can preceed this era
2015 # last year of main recr_devs; forecast devs start in following year
3 #_recdev phase 
1 # (0/1) to read 13 advanced options
 0 #_recdev_early_start (0=none; neg value makes relative to recdev_start)
 -4 #_recdev_early_phase
 5 #_forecast_recruitment phase (incl. late recr) (0 value resets to maxphase+1)
 1 #_lambda for Fcast_recr_like occurring before endyr+1
 1966.0 #_last_yr_nobias_adj_in_MPD; begin of ramp
 1999.3 #_first_yr_fullbias_adj_in_MPD; begin of plateau
 2014.8 #_last_yr_fullbias_adj_in_MPD
 2018.2 #_end_yr_for_ramp_in_MPD (can be in forecast to shape ramp, but SS3 sets bias_adj to 0.0 for fcast yrs)
 0.9307 #_max_bias_adj_in_MPD (typical ~0.8; -3 sets all years to 0.0; -2 sets all non-forecast yrs w/ estimated recdevs to 1.0; -1 sets biasadj=1.0 for all yrs w/ recdevs)
 0 #_period of cycles in recruitment (N parms read below)
 -5 #min rec_dev
 5 #max rec_dev
 0 #_read_recdevs
#_end of advanced SR options
#
#_placeholder for full parameter lines for recruitment cycles
# read specified recr devs
#_Yr Input_value
#
# all recruitment deviations
#  1994R 1995R 1996R 1997R 1998R 1999R 2000R 2001R 2002R 2003R 2004R 2005R 2006R 2007R 2008R 2009R 2010R 2011R 2012R 2013R 2014R 2015R 2016F 2017F 2018F 2019F 2020F 2021F 2022F
#  -0.640542 -0.360788 -0.352397 -0.288281 -0.323908 0.264475 0.172854 0.133309 0.12218 0.131315 -0.127889 0.0298662 0.264976 -0.17175 -0.266162 -0.41824 -0.108542 0.230987 0.339283 0.188272 0.315108 0.865874 0.462062 -0.102891 0 0 0 0 0
#
#Fishing Mortality info 
0.5 # F ballpark value in units of annual_F
-2001 # F ballpark year (neg value to disable)
4 # F_Method:  1=Pope midseason rate; 2=F as parameter; 3=F as hybrid; 4=fleet-specific parm/hybrid (#4 is superset of #2 and #3 and is recommended)
4 # max F (methods 2-4) or harvest fraction (method 1)
# read list of fleets that do F as parameter; unlisted fleets stay hybrid, bycatch fleets must be included with start_PH=1, high F fleets should switch early
# (A) fleet, (B) F_initial_value (used if start_PH=1), (C) start_PH for parms (99 to stay in hybrid)
# (A) (B) (C)  (terminate list with -9999 for fleet)
 1 0.05 3 # CM_E
 2 0.05 3 # CM_W
 3 0.05 3 # REC
 4 0.05 1 # SMP_BYC
-9999 1 1 # end of list
3 #_number of loops for hybrid tuning; 4 good; 3 faster; 2 enough if switching to parms is enabled
#
#_initial_F_parms; for each fleet x season that has init_catch; nest season in fleet; count = 0
#_for unconstrained init_F, use an arbitrary initial catch and set lambda=0 for its logL
#_ LO HI INIT PRIOR PR_SD  PR_type  PHASE
#
# F rates by fleet x season
# Yr:  1950 1951 1952 1953 1954 1955 1956 1957 1958 1959 1960 1961 1962 1963 1964 1965 1966 1967 1968 1969 1970 1971 1972 1973 1974 1975 1976 1977 1978 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012 2013 2014 2015 2016 2017 2018 2019 2020 2021 2022
# seas:  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
# CM_E 4.06692e-05 8.12487e-05 0.000123009 0.000165134 0.000208307 0.000252204 0.000296426 0.000341974 0.000388843 0.000436834 0.000486591 0.000535787 0.000582513 0.000624927 0.000680791 0.000674051 0.000350564 0.000708211 0.00140516 0.00178856 0.00166901 0.00182244 0.00160525 0.00269718 0.00255004 0.00553563 0.00485619 0.00658258 0.00567432 0.00434216 0.0031863 0.00464 0.00482974 0.00765885 0.0108952 0.0140841 0.014902 0.0119711 0.0114216 0.0114168 0.0284555 0.0238017 0.0315431 0.0455904 0.0459332 0.0465185 0.0386887 0.0384719 0.0325462 0.0402045 0.0323036 0.0367929 0.0428539 0.0477813 0.0374353 0.0550844 0.059315 0.0604352 0.0681364 0.101215 0.0487156 0.0944859 0.0672166 0.0383592 0.0424965 0.0229636 0.0253698 0.0259983 0.183016 0.183016 0.183016 0.183016 0.183016
# CM_W 3.96647e-05 7.94393e-05 0.000119563 0.00016038 0.000201981 0.000244256 0.000287101 0.000330613 0.000374923 0.000420334 0.000466939 0.000514387 0.00056123 0.000605254 0.00063026 0.000552934 0.000176366 0.000415434 0.00132103 0.000709898 0.00116168 0.00125587 0.00121357 0.00143041 0.00173357 0.00282912 0.00156255 0.00503789 0.00422138 0.00570972 0.00385313 0.00301061 0.00384147 0.00427754 0.0227413 0.0202668 0.0266488 0.0296423 0.0305988 0.0321096 0.0322103 0.027919 0.0378045 0.0386111 0.0412783 0.0273165 0.0285889 0.0504861 0.0462027 0.0623958 0.0467047 0.0586207 0.0659321 0.0800823 0.0740423 0.0675952 0.0455189 0.082418 0.0563383 0.0523333 0.0406097 0.0379749 0.0476067 0.0290966 0.0378957 0.0370506 0.0369371 0.0294302 0.221373 0.221373 0.221373 0.221373 0.221373
# REC 0.000137692 0.00037285 0.000614719 0.000862374 0.00111426 0.0013714 0.00163372 0.00190345 0.00218436 0.00247827 0.00278102 0.00307611 0.00333794 0.00357266 0.0038125 0.0040634 0.00432292 0.00458498 0.00485339 0.0051291 0.00541484 0.00569247 0.00594194 0.006178 0.00641228 0.00665639 0.00690547 0.00717255 0.00746965 0.00778505 0.00810671 0.0084296 0.0186484 0.00719021 0.0111278 0.0215974 0.0306685 0.0390245 0.058829 0.0326624 0.0452418 0.0590443 0.0748159 0.05772 0.0498776 0.0699082 0.0296677 0.054319 0.0302044 0.0631304 0.0380079 0.102505 0.0811773 0.073215 0.0723704 0.0882173 0.0822731 0.0807313 0.0387354 0.0563704 0.0411824 0.0972136 0.0721943 0.118047 0.102503 0.0827288 0.0791097 0.0940245 0.561394 0.561394 0.561394 0.561394 0.561394
# SMP_BYC 0.0314791 0.0429227 0.0506946 0.0533068 0.0691052 0.072035 0.0920951 0.105445 0.12914 0.139221 0.139183 0.10542 0.101515 0.115736 0.122263 0.135724 0.133588 0.145561 0.147939 0.168179 0.158479 0.151204 0.147277 0.149905 0.149244 0.150182 0.155803 0.170594 0.180216 0.188644 0.193 0.17605 0.168706 0.171831 0.196304 0.190179 0.2108 0.180689 0.166644 0.1845 0.16572 0.175991 0.223865 0.238662 0.362602 0.258098 0.25259 0.260053 0.316508 0.203932 0.167769 0.186684 0.22733 0.193043 0.182754 0.150531 0.0977749 0.0688866 0.0490056 0.0769732 0.0567064 0.0670413 0.0599576 0.0672588 0.0543364 0.0453144 0.0466827 0.050227 0.0621485 0.0621485 0.0621485 0.0621485 0.0621485
#
#_Q_setup for fleets with cpue or survey data
#_1:  fleet number
#_2:  link type: (1=simple q, 1 parm; 2=mirror simple q, 1 mirrored parm; 3=q and power, 2 parm; 4=mirror with offset, 2 parm)
#_3:  extra input for link, i.e. mirror fleet# or dev index number
#_4:  0/1 to select extra sd parameter
#_5:  0/1 for biasadj or not
#_6:  0/1 to float
#_   fleet      link link_info  extra_se   biasadj     float  #  fleetname
         1         1         0         0         0         1  #  CM_E
         2         1         0         0         0         1  #  CM_W
         3         1         0         0         0         1  #  REC
         4         1         0         0         0         0  #  SMP_BYC
         5         1         0         0         0         1  #  HB_E
         6         1         0         0         0         1  #  HB_W
         7         1         0         0         0         1  #  LARVAL
         8         1         0         0         0         1  #  VIDEO
         9         1         0         0         0         1  #  SEAMAP
-9999 0 0 0 0 0
#
#_Q_parms(if_any);Qunits_are_ln(q)
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type      PHASE    env-var    use_dev   dev_mnyr   dev_mxyr     dev_PH      Block    Blk_Fxn  #  parm_name
           -25            25      -9.34322             0             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_CM_E(1)
           -25            25      -8.99573             0             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_CM_W(2)
           -25            25      -9.94041             0             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_REC(3)
           -10            20       1.84347             0             0             0          2          0          0          0          0          0          0          0  #  LnQ_base_SMP_BYC(4)
           -25            25      -10.3281             0             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_HB_E(5)
           -25            25      -10.2511             0             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_HB_W(6)
           -30            25      -26.6364             0             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_LARVAL(7)
           -25            25      -10.9251             0             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_VIDEO(8)
           -25            25      -11.0472             0             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_SEAMAP(9)
#_no timevary Q parameters
#
#_size_selex_patterns
#Pattern:_0;  parm=0; selex=1.0 for all sizes
#Pattern:_1;  parm=2; logistic; with 95% width specification
#Pattern:_2;  parm=6; modification of pattern 24 with improved sex-specific offset
#Pattern:_5;  parm=2; mirror another size selex; PARMS pick the min-max bin to mirror
#Pattern:_11; parm=2; selex=1.0  for specified min-max population length bin range
#Pattern:_15; parm=0; mirror another age or length selex
#Pattern:_6;  parm=2+special; non-parm len selex
#Pattern:_43; parm=2+special+2;  like 6, with 2 additional param for scaling (average over bin range)
#Pattern:_8;  parm=8; double_logistic with smooth transitions and constant above Linf option
#Pattern:_9;  parm=6; simple 4-parm double logistic with starting length; parm 5 is first length; parm 6=1 does desc as offset
#Pattern:_21; parm=2+special; non-parm len selex, read as pairs of size, then selex
#Pattern:_22; parm=4; double_normal as in CASAL
#Pattern:_23; parm=6; double_normal where final value is directly equal to sp(6) so can be >1.0
#Pattern:_24; parm=6; double_normal with sel(minL) and sel(maxL), using joiners
#Pattern:_25; parm=3; exponential-logistic in length
#Pattern:_27; parm=special+3; cubic spline in length; parm1==1 resets knots; parm1==2 resets all 
#Pattern:_42; parm=special+3+2; cubic spline; like 27, with 2 additional param for scaling (average over bin range)
#_discard_options:_0=none;_1=define_retention;_2=retention&mortality;_3=all_discarded_dead;_4=define_dome-shaped_retention
#_Pattern Discard Male Special
 0 2 0 0 # 1 CM_E
 0 2 0 0 # 2 CM_W
 0 2 0 0 # 3 REC
 0 3 0 0 # 4 SMP_BYC
 0 0 0 0 # 5 HB_E
 0 0 0 0 # 6 HB_W
 0 0 0 0 # 7 LARVAL
 24 0 0 0 # 8 VIDEO
 24 0 0 0 # 9 SEAMAP
 0 0 0 0 # 10 COMP
#
#_age_selex_patterns
#Pattern:_0; parm=0; selex=1.0 for ages 0 to maxage
#Pattern:_10; parm=0; selex=1.0 for ages 1 to maxage
#Pattern:_11; parm=2; selex=1.0  for specified min-max age
#Pattern:_12; parm=2; age logistic
#Pattern:_13; parm=8; age double logistic
#Pattern:_14; parm=nages+1; age empirical
#Pattern:_15; parm=0; mirror another age or length selex
#Pattern:_16; parm=2; Coleraine - Gaussian
#Pattern:_17; parm=nages+1; empirical as random walk  N parameters to read can be overridden by setting special to non-zero
#Pattern:_41; parm=2+nages+1; // like 17, with 2 additional param for scaling (average over bin range)
#Pattern:_18; parm=8; double logistic - smooth transition
#Pattern:_19; parm=6; simple 4-parm double logistic with starting age
#Pattern:_20; parm=6; double_normal,using joiners
#Pattern:_26; parm=3; exponential-logistic in age
#Pattern:_27; parm=3+special; cubic spline in age; parm1==1 resets knots; parm1==2 resets all 
#Pattern:_42; parm=2+special+3; // cubic spline; with 2 additional param for scaling (average over bin range)
#Age patterns entered with value >100 create Min_selage from first digit and pattern from remainder
#_Pattern Discard Male Special
 12 0 0 0 # 1 CM_E
 12 0 0 0 # 2 CM_W
 20 0 0 0 # 3 REC
 19 0 0 0 # 4 SMP_BYC
 15 0 0 3 # 5 HB_E
 15 0 0 3 # 6 HB_W
 0 0 0 0 # 7 LARVAL
 0 0 0 0 # 8 VIDEO
 0 0 0 0 # 9 SEAMAP
 11 0 0 0 # 10 COMP
#
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type      PHASE    env-var    use_dev   dev_mnyr   dev_mxyr     dev_PH      Block    Blk_Fxn  #  parm_name
# 1   CM_E LenSelex
            10           100         10.16         10.16            -1             0         -3          0          0          0          0          0          1          2  #  Retain_L_infl_CM_E(1)
            -1            20         1e-06             1            -1             0         -3          0          0          0          0          0          0          0  #  Retain_L_width_CM_E(1)
           -10            10            10            10            -1             0         -2          0          0          0          0          0          1          2  #  Retain_L_asymptote_logit_CM_E(1)
            -1             2             0             0            -1             0         -4          0          0          0          0          0          0          0  #  Retain_L_maleoffset_CM_E(1)
           -10            10            -5            -5            -1             0         -2          0          0          0          0          0          0          0  #  DiscMort_L_infl_CM_E(1)
            -1             2         1e-06             1            -1             0         -4          0          0          0          0          0          0          0  #  DiscMort_L_width_CM_E(1)
            -1             2          0.15          0.15            -1             0         -2          0          0          0          0          0          3          2  #  DiscMort_L_level_old_CM_E(1)
            -1             2             0             0            -1             0         -4          0          0          0          0          0          0          0  #  DiscMort_L_male_offset_CM_E(1)
# 2   CM_W LenSelex
            10           100         10.16         10.16            -1             0         -3          0          0          0          0          0          1          2  #  Retain_L_infl_CM_W(2)
            -1            20         1e-06             1            -1             0         -3          0          0          0          0          0          0          0  #  Retain_L_width_CM_W(2)
           -10            10            10            10            -1             0         -2          0          0          0          0          0          1          2  #  Retain_L_asymptote_logit_CM_W(2)
            -1             2             0             0            -1             0         -4          0          0          0          0          0          0          0  #  Retain_L_maleoffset_CM_W(2)
           -10            10            -5            -5            -1             0         -2          0          0          0          0          0          0          0  #  DiscMort_L_infl_CM_W(2)
            -1             2         1e-06             1            -1             0         -4          0          0          0          0          0          0          0  #  DiscMort_L_width_CM_W(2)
            -1             2          0.15          0.15            -1             0         -2          0          0          0          0          0          3          2  #  DiscMort_L_level_old_CM_W(2)
            -1             2             0             0            -1             0         -4          0          0          0          0          0          0          0  #  DiscMort_L_male_offset_CM_W(2)
# 3   REC LenSelex
            10           100         10.16         10.16            -1             0         -3          0          0          0          0          0          2          2  #  Retain_L_infl_REC(3)
            -1            20         1e-06             1            -1             0         -3          0          0          0          0          0          0          0  #  Retain_L_width_REC(3)
           -10            10            10            10            -1             0         -2          0          0          0          0          0          2          2  #  Retain_L_asymptote_logit_REC(3)
            -1             2             0             0            -1             0         -4          0          0          0          0          0          0          0  #  Retain_L_maleoffset_REC(3)
           -10            10            -5            -5            -1             0         -2          0          0          0          0          0          0          0  #  DiscMort_L_infl_REC(3)
            -1             2         1e-06             1            -1             0         -4          0          0          0          0          0          0          0  #  DiscMort_L_width_REC(3)
            -1             2          0.15          0.15            -1             0         -2          0          0          0          0          0          3          2  #  DiscMort_L_level_old_REC(3)
            -1             2             0             0            -1             0         -4          0          0          0          0          0          0          0  #  DiscMort_L_male_offset_REC(3)
# 4   SMP_BYC LenSelex
# 5   HB_E LenSelex
# 6   HB_W LenSelex
# 7   LARVAL LenSelex
# 8   VIDEO LenSelex
           7.5          52.5       19.2557          42.7          0.05             0          2          0          0          0          0        0.5          0          0  #  Size_DblN_peak_VIDEO(8)
           -10             3      -1.37276          -0.4          0.05             0          3          0          0          0          0        0.5          0          0  #  Size_DblN_top_logit_VIDEO(8)
            -6            12       1.10772           5.5          0.05             0          3          0          0          0          0        0.5          0          0  #  Size_DblN_ascend_se_VIDEO(8)
            -6             6      0.313286           5.1          0.05             0          3          0          0          0          0        0.5          0          0  #  Size_DblN_descend_se_VIDEO(8)
           -15             5      -1.53307          -4.2          0.05             0          2          0          0          0          0        0.5          0          0  #  Size_DblN_start_logit_VIDEO(8)
            -8             5      0.693019           0.4          0.05             0          2          0          0          0          0        0.5          0          0  #  Size_DblN_end_logit_VIDEO(8)
# 9   SEAMAP LenSelex
           7.5          52.5       14.8115            13          0.05             0          2          0          0          0          0        0.5          0          0  #  Size_DblN_peak_SEAMAP(9)
           -10             3      -3.93027          -1.1          0.05             0          3          0          0          0          0        0.5          0          0  #  Size_DblN_top_logit_SEAMAP(9)
            -6            12       1.29723           3.1          0.05             0          3          0          0          0          0        0.5          0          0  #  Size_DblN_ascend_se_SEAMAP(9)
            -4             6       3.11953             5          0.05             0          3          0          0          0          0        0.5          0          0  #  Size_DblN_descend_se_SEAMAP(9)
           -15             5      -1.24781          -4.5          0.05             0          2          0          0          0          0        0.5          0          0  #  Size_DblN_start_logit_SEAMAP(9)
            -8             5      -5.22542           0.1          0.05             0          2          0          0          0          0        0.5          0          0  #  Size_DblN_end_logit_SEAMAP(9)
# 10   COMP LenSelex
# 1   CM_E AgeSelex
           0.5            14       2.13977          2.66             0             0          3          0          0          0          0          0          0          0  #  Age_inflection_CM_E(1)
           0.5            14      0.920405        7.2774             0             0          1          0          0          0          0          0          0          0  #  Age_95%width_CM_E(1)
# 2   CM_W AgeSelex
           0.5            14       3.66543          2.66             0             0          3          0          0          0          0          0          0          0  #  Age_inflection_CM_W(2)
           0.5            14       2.04896        7.2774             0             0          1          0          0          0          0          0          0          0  #  Age_95%width_CM_W(2)
# 3   REC AgeSelex
             1            10       3.30967           4.3          0.05             0          2          0          0          0          0        0.5          0          0  #  Age_DblN_peak_REC(3)
           -10             3      -9.25829          -4.6          0.05             0          3          0          0          0          0        0.5          0          0  #  Age_DblN_top_logit_REC(3)
            -6            12      0.527117           0.7          0.05             0          3          0          0          0          0        0.5          0          0  #  Age_DblN_ascend_se_REC(3)
            -5             6        3.0289           2.7          0.05             0          3          0          0          0          0        0.5          0          0  #  Age_DblN_descend_se_REC(3)
           -15             5       -12.055         -11.2          0.05             0          2          0          0          0          0        0.5          0          0  #  Age_DblN_start_logit_REC(3)
            -8             5      -2.17918          -3.3          0.05             0          2          0          0          0          0        0.5          0          0  #  Age_DblN_end_logit_REC(3)
# 4   SMP_BYC AgeSelex
         1e-07             2           0.5           0.5             0             0         -4          0          0          0          0          0          0          0  #  AgeSel_P1_SMP_BYC(4)
           0.5         1e+07           100           100             0             0         -4          0          0          0          0          0          0          0  #  AgeSel_P2_SMP_BYC(4)
           0.3             3           1.5           1.5             0             0         -4          0          0          0          0          0          0          0  #  AgeSel_P3_SMP_BYC(4)
           0.5         1e+07        2.4096        2.4096             0             0         -4          0          0          0          0          0          0          0  #  AgeSel_P4_SMP_BYC(4)
            -1             1             0             0             0             0         -4          0          0          0          0          0          0          0  #  AgeSel_P5_SMP_BYC(4)
            -1             1             0             0             0             0         -4          0          0          0          0          0          0          0  #  AgeSel_P6_SMP_BYC(4)
# 5   HB_E AgeSelex
# 6   HB_W AgeSelex
# 7   LARVAL AgeSelex
# 8   VIDEO AgeSelex
# 9   SEAMAP AgeSelex
# 10   COMP AgeSelex
            -1             3             1             0             0             0        -99          0          0          0          0          0          0          0  #  minage@sel=1_COMP(10)
             1            20            14             0             0             0        -99          0          0          0          0          0          0          0  #  maxage@sel=1_COMP(10)
#_No_Dirichlet parameters
# timevary selex parameters 
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type    PHASE  #  parm_name
            10           100         20.32         20.32            -1             0      -4  # Retain_L_infl_CM_E(1)_BLK1repl_1990
            10           100         27.94         27.94            -1             0      -4  # Retain_L_infl_CM_E(1)_BLK1repl_2005
            10           100          25.4          25.4            -1             0      -4  # Retain_L_infl_CM_E(1)_BLK1repl_2008
           -10            10            10            10            -1             0      -4  # Retain_L_asymptote_logit_CM_E(1)_BLK1repl_1990
           -10            10            10            10            -1             0      -4  # Retain_L_asymptote_logit_CM_E(1)_BLK1repl_2005
           -10            10            10            10            -1             0      -4  # Retain_L_asymptote_logit_CM_E(1)_BLK1repl_2008
            -1             2          0.15          0.15            -1             0      -4  # DiscMort_L_level_old_CM_E(1)_BLK3repl_2008
            10           100         20.32         20.32            -1             0      -4  # Retain_L_infl_CM_W(2)_BLK1repl_1990
            10           100         27.94         27.94            -1             0      -4  # Retain_L_infl_CM_W(2)_BLK1repl_2005
            10           100          25.4          25.4            -1             0      -4  # Retain_L_infl_CM_W(2)_BLK1repl_2008
           -10            10            10            10            -1             0      -4  # Retain_L_asymptote_logit_CM_W(2)_BLK1repl_1990
           -10            10            10            10            -1             0      -4  # Retain_L_asymptote_logit_CM_W(2)_BLK1repl_2005
           -10            10            10            10            -1             0      -4  # Retain_L_asymptote_logit_CM_W(2)_BLK1repl_2008
            -1             2          0.15          0.15            -1             0      -4  # DiscMort_L_level_old_CM_W(2)_BLK3repl_2008
            10           100         20.32         20.32            -1             0      -4  # Retain_L_infl_REC(3)_BLK2repl_1990
            10           100          25.4          25.4            -1             0      -4  # Retain_L_infl_REC(3)_BLK2repl_1997
            10           100         27.94         27.94            -1             0      -4  # Retain_L_infl_REC(3)_BLK2repl_2005
            10           100          25.4          25.4            -1             0      -4  # Retain_L_infl_REC(3)_BLK2repl_2008
           -10            10            10            10            -1             0      -4  # Retain_L_asymptote_logit_REC(3)_BLK2repl_1990
           -10            10            10            10            -1             0      -4  # Retain_L_asymptote_logit_REC(3)_BLK2repl_1997
           -10            10            10            10            -1             0      -4  # Retain_L_asymptote_logit_REC(3)_BLK2repl_2005
           -10            10            10            10            -1             0      -4  # Retain_L_asymptote_logit_REC(3)_BLK2repl_2008
            -1             2          0.15          0.15            -1             0      -4  # DiscMort_L_level_old_REC(3)_BLK3repl_2008
# info on dev vectors created for selex parms are reported with other devs after tag parameter section 
#
0   #  use 2D_AR1 selectivity(0/1)
#_no 2D_AR1 selex offset used
#
# Tag loss and Tag reporting parameters go next
0  # TG_custom:  0=no read and autogen if tag data exist; 1=read
#_Cond -6 6 1 1 2 0.01 -4 0 0 0 0 0 0 0  #_placeholder if no parameters
#
# deviation vectors for timevary parameters
#  base   base first block   block  env  env   dev   dev   dev   dev   dev
#  type  index  parm trend pattern link  var  vectr link _mnyr  mxyr phase  dev_vector
#      1    29     1     0     0     2     1     0     0     0     0     0
#      5     1     2     1     2     0     0     0     0     0     0     0
#      5     3     5     1     2     0     0     0     0     0     0     0
#      5     7     8     3     2     0     0     0     0     0     0     0
#      5     9     9     1     2     0     0     0     0     0     0     0
#      5    11    12     1     2     0     0     0     0     0     0     0
#      5    15    15     3     2     0     0     0     0     0     0     0
#      5    17    16     2     2     0     0     0     0     0     0     0
#      5    19    20     2     2     0     0     0     0     0     0     0
#      5    23    24     3     2     0     0     0     0     0     0     0
     #
# Input variance adjustments factors: 
 #_1=add_to_survey_CV
 #_2=add_to_discard_stddev
 #_3=add_to_bodywt_CV
 #_4=mult_by_lencomp_N
 #_5=mult_by_agecomp_N
 #_6=mult_by_size-at-age_N
 #_7=mult_by_generalized_sizecomp
#_Factor  Fleet  Value
 -9999   1    0  # terminator
#
10 #_maxlambdaphase
1 #_sd_offset; must be 1 if any growthCV, sigmaR, or survey extraSD is an estimated parameter
# read 5 changes to default Lambdas (default value is 1.0)
# Like_comp codes:  1=surv; 2=disc; 3=mnwt; 4=length; 5=age; 6=SizeFreq; 7=sizeage; 8=catch; 9=init_equ_catch; 
# 10=recrdev; 11=parm_prior; 12=parm_dev; 13=CrashPen; 14=Morphcomp; 15=Tag-comp; 16=Tag-negbin; 17=F_ballpark; 18=initEQregime
#like_comp fleet  phase  value  sizefreq_method
 2 1 1 0 1
 2 2 1 0 1
 2 3 1 0 1
 1 10 1 0 1
 8 10 1 0 1
-9999  1  1  1  1  #  terminator
#
# lambdas (for info only; columns are phases)
#  1 1 1 1 1 1 1 1 1 1 #_CPUE/survey:_1
#  1 1 1 1 1 1 1 1 1 1 #_CPUE/survey:_2
#  1 1 1 1 1 1 1 1 1 1 #_CPUE/survey:_3
#  1 1 1 1 1 1 1 1 1 1 #_CPUE/survey:_4
#  1 1 1 1 1 1 1 1 1 1 #_CPUE/survey:_5
#  1 1 1 1 1 1 1 1 1 1 #_CPUE/survey:_6
#  1 1 1 1 1 1 1 1 1 1 #_CPUE/survey:_7
#  1 1 1 1 1 1 1 1 1 1 #_CPUE/survey:_8
#  1 1 1 1 1 1 1 1 1 1 #_CPUE/survey:_9
#  0 0 0 0 0 0 0 0 0 0 #_CPUE/survey:_10
#  0 0 0 0 0 0 0 0 0 0 #_discard:_1
#  0 0 0 0 0 0 0 0 0 0 #_discard:_2
#  0 0 0 0 0 0 0 0 0 0 #_discard:_3
#  1 1 1 1 1 1 1 1 1 1 #_discard:_4
#  0 0 0 0 0 0 0 0 0 0 #_discard:_5
#  0 0 0 0 0 0 0 0 0 0 #_discard:_6
#  0 0 0 0 0 0 0 0 0 0 #_discard:_7
#  0 0 0 0 0 0 0 0 0 0 #_discard:_8
#  0 0 0 0 0 0 0 0 0 0 #_discard:_9
#  0 0 0 0 0 0 0 0 0 0 #_discard:_10
#  0 0 0 0 0 0 0 0 0 0 #_lencomp:_1
#  0 0 0 0 0 0 0 0 0 0 #_lencomp:_2
#  0 0 0 0 0 0 0 0 0 0 #_lencomp:_3
#  0 0 0 0 0 0 0 0 0 0 #_lencomp:_4
#  0 0 0 0 0 0 0 0 0 0 #_lencomp:_5
#  0 0 0 0 0 0 0 0 0 0 #_lencomp:_6
#  0 0 0 0 0 0 0 0 0 0 #_lencomp:_7
#  1 1 1 1 1 1 1 1 1 1 #_lencomp:_8
#  1 1 1 1 1 1 1 1 1 1 #_lencomp:_9
#  0 0 0 0 0 0 0 0 0 0 #_lencomp:_10
#  1 1 1 1 1 1 1 1 1 1 #_agecomp:_1
#  1 1 1 1 1 1 1 1 1 1 #_agecomp:_2
#  1 1 1 1 1 1 1 1 1 1 #_agecomp:_3
#  0 0 0 0 0 0 0 0 0 0 #_agecomp:_4
#  0 0 0 0 0 0 0 0 0 0 #_agecomp:_5
#  0 0 0 0 0 0 0 0 0 0 #_agecomp:_6
#  0 0 0 0 0 0 0 0 0 0 #_agecomp:_7
#  0 0 0 0 0 0 0 0 0 0 #_agecomp:_8
#  0 0 0 0 0 0 0 0 0 0 #_agecomp:_9
#  0 0 0 0 0 0 0 0 0 0 #_agecomp:_10
#  1 1 1 1 1 1 1 1 1 1 #_init_equ_catch1
#  1 1 1 1 1 1 1 1 1 1 #_init_equ_catch2
#  1 1 1 1 1 1 1 1 1 1 #_init_equ_catch3
#  1 1 1 1 1 1 1 1 1 1 #_init_equ_catch4
#  1 1 1 1 1 1 1 1 1 1 #_init_equ_catch5
#  1 1 1 1 1 1 1 1 1 1 #_init_equ_catch6
#  1 1 1 1 1 1 1 1 1 1 #_init_equ_catch7
#  1 1 1 1 1 1 1 1 1 1 #_init_equ_catch8
#  1 1 1 1 1 1 1 1 1 1 #_init_equ_catch9
#  1 1 1 1 1 1 1 1 1 1 #_init_equ_catch10
#  1 1 1 1 1 1 1 1 1 1 #_recruitments
#  1 1 1 1 1 1 1 1 1 1 #_parameter-priors
#  1 1 1 1 1 1 1 1 1 1 #_parameter-dev-vectors
#  1 1 1 1 1 1 1 1 1 1 #_crashPenLambda
#  0 0 0 0 0 0 0 0 0 0 # F_ballpark_lambda
0 # (0/1/2) read specs for more stddev reporting: 0 = skip, 1 = read specs for reporting stdev for selectivity, size, and numbers, 2 = add options for M,Dyn. Bzero, SmryBio
 # 0 2 0 0 # Selectivity: (1) fleet, (2) 1=len/2=age/3=both, (3) year, (4) N selex bins
 # 0 0 # Growth: (1) growth pattern, (2) growth ages
 # 0 0 0 # Numbers-at-age: (1) area(-1 for all), (2) year, (3) N ages
 # -1 # list of bin #'s for selex std (-1 in first bin to self-generate)
 # -1 # list of ages for growth std (-1 in first bin to self-generate)
 # -1 # list of ages for NatAge std (-1 in first bin to self-generate)
999

