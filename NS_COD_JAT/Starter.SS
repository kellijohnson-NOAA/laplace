# darkblotched starter file for SS v3.x
SimTest.dat 	# Data file
SimTest.ctl # Control file

1 # Read initial values from .par file: 0=no,1=yes
1	# DOS display detail: 0,1,2
0 	# Report file detail: 0,1,2 
0 	# Detailed checkup.sso file (0,1) 
0	# Write parameter iteration trace file during minimization
0	# Write cumulative report: 0=skip,1=short,2=full
1	# Include prior likelihood for non-estimated parameters
1 	# Use Soft Boundaries to aid convergence (0,1) (recommended)
0 	# N bootstrap datafiles to create
25 	# Last phase for estimation
0 	# MCMC burn-in
1 	# MCMC thinning interval
0 	# Jitter initial parameter values by this fraction
-1	# Min year for spbio sd_report (neg val = styr-2, virgin state)
-2	# Max year for spbio sd_report (-1=endyr+1, -2=entire forecast)
0 	# N individual SD years
0.0001 	# Ending convergence criteria
0 	# Retrospective year relative to end year (i.e. -4)
1 	# Min age for summary biomass
1 	# Depletion basis: denom is: 0=skip; 1=rel X*B0; 2=rel X*Bmsy; 3=rel X*B_styr
1 	# Fraction (X) for Depletion denominator (e.g. 0.4)
1 	# (1-SPR)_reporting:  0=skip; 1=rel(1-SPR); 2=rel(1-SPR_MSY); 3=rel(1-SPR_Btarget); 4=notrel
1 	# F_std reporting: 0=skip; 1=exploit(Bio); 2=exploit(Num); 3=sum(frates)
#0 45    #_min and max age over which average F will be calculated
0 	# F_report_basis: 0=raw; 1=rel Fspr; 2=rel Fmsy ; 3=rel Fbtgt

999 # end of file marker
