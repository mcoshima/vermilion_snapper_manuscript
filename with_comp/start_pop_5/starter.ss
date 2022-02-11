#C
vs_envM2.dat #file name of data file
vs.ctl #file name of control file
0 #inital parameter values 0 = use values in control file; 1 = use ss.par after reading setup in control file
1 #run display detail 0 = none other than ADMB outputs; 1 = one brief line of display for each iteration; 2 = fuller display per iteration 
1 #detailed age-structure report: 0 = minimal output for data-limited methods; 1 = include all output; 2 = brief output, no growth; 3 = custom output
0 #write 1st iteration details: 0 = omit; 1 = write detailed intermediate calculations to echoinput.sso 
0 #parameter trace: 0 = omit; 1 = write good iteration and active parameters; 2  = write good iterations and all params; 3 = write every iteration and all params; 4 = write every iteration and active params
1 #cumulative report: 0 = omit; 1 = brief; 2 = full 
0 #full priors: 0 = only calculate priors for active parameters; 1 = calculate priors for all paramters that have a defined prior
1 #soft bounds: 0 = omit; 1 = use 
1 #data file output: 0 = none, 1 = output an annotated replicate of the input data file; 2 = add a second data file containing the model's expected values with no added error; 3+ = add N-2 parametric bootstrap data files 
10 #turn off estimation: -1 = exit after reading input files; 0 = exit after one call to the calculation routines and production of sso and ss_new files; <positive value> = exit after completing this phase
100 #MCMC burn interval
1 #MCMC thin interval
0.0 #jitter
-1 #SD report start: -1 = begin annual SD report in start year; <year> = begin SD report this year 
-2 #SD report end: -1 = end annual SD report in end year; -2 = end annual SD report in last forecast year; <value> = end SD report in this year 
0 #extra SD report years: 0 = none; <value> = number of years to read 
#1940 1950 if extra SD report years > 0 add vector of years for additional SD reporting
0.0001 #final convergence 
0 #retrospective year 
1 #summary biomass min age 
1 #depletion basis: 0 = skip; 1 = X*SB0; 2 = X*SBmsy; 3 = X*SBstyr; 4 = X*SBendyr
1 #fraction for depletion denominator
4 #SPR report basis: 0 = skip; 1 = use 1-SPRtarget; 2 = use 1-SPR at MSY; 3 = use 1-SPR at Btarget; 4 = no denominator, report actual 1-SPR values 
2 #Annual F units: 0 = skip; 1 = exploitation rate in biomass; 2 = exploitation rate in numbers; 3 = sum(apical Fs by fleet); 4 = population F for range of ages; 5 = unweighted average F for range of ages 
#3 7 #if Fstd reporting >= 4 add age range 
0 #F report basis: 0 = not relative, report raw values; 1 = use F std values relative to SPRtarget; 2 = use F std value relative to FMSY; 3 = use Fstd value relative to Fbtarget 
0 #MCMC output detail: 0 = default; 1 =  output likelihood components and assocaited lambda values; 2 = expanded output; 3 = make output subdirectory for each MCMC vector
0 
#ALK tolerance level 
3.30 #model version 