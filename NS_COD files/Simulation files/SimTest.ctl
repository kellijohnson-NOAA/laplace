#V3.24L
#_data_and_control_files: NS_COD.dat // NS_COD_Randwalk.ctl
#_SS-V3.24L-safe;_12/5/2012;_Stock_Synthesis_by_Richard_Methot_(NOAA)_using_ADMB_10.1
1  #_N_Growth_Patterns
1 #_N_Morphs_Within_GrowthPattern 
#_Cond 1 #_Morph_between/within_stdev_ratio (no read if N_morphs=1)
#_Cond  1 #vector_Morphdist_(-1_in_first_val_gives_normal_approx)
#
#_Cond 0  #  N recruitment designs goes here if N_GP*nseas*area>1
#_Cond 0  #  placeholder for recruitment interaction request
#_Cond 1 1 1  # example recruitment design element for GP=1, seas=1, area=1
#
#_Cond 0 # N_movement_definitions goes here if N_areas > 1
#_Cond 1.0 # first age that moves (real age at begin of season, not integer) also cond on do_migration>0
#_Cond 1 1 1 2 4 10 # example move definition for seas=1, morph=1, source=1 dest=2, age1=4, age2=10
#
1 #_Nblock_Patterns
 19 #_blocks_per_pattern 
# begin and end years of blocks
 1993 1993 1994 1994 1995 1995 1996 1996 1997 1997 1998 1998 1999 1999 2000 2000 2001 2001 2002 2002 2003 2003 2004 2004 2005 2005 2006 2006 2007 2007 2008 2008 2009 2009 2010 2010 2011 2011
#
0.5 #_fracfemale 
1 #_natM_type:_0=1Parm; 1=N_breakpoints;_2=Lorenzen;_3=agespecific;_4=agespec_withseasinterpolate
5 #_N_breakpoints
 1 2 3 4 5 # age(real) at M breakpoints
3 # GrowthModel: 1=vonBert with L1&L2; 2=Richards with L1&L2; 3=age_speciific_K; 4=not implemented
2 #_Growth_Age_for_L1
999 #_Growth_Age_for_L2 (999 to use as Linf)
14 # number of K multipliers to read
 2 3 4 5 6 7 8 9 10 11 12 13 14 15 # ages for K multiplier
0 #_SD_add_to_LAA (set to 0.1 for SS2 V1.x compatibility)
0 #_CV_Growth_Pattern:  0 CV=f(LAA); 1 CV=F(A); 2 SD=F(LAA); 3 SD=F(A); 4 logSD=F(A)
3 #_maturity_option:  1=length logistic; 2=age logistic; 3=read age-maturity matrix by growth_pattern; 4=read age-fecundity; 5=read fec and wt from wtatage.ss
#_Age_Maturity by growth pattern
 0 0.01 0.05 0.23 0.62 0.86 1 1 1 1 1 1 1 1 1 1
1 #_First_Mature_Age
1 #_fecundity option:(1)eggs=Wt*(a+b*Wt);(2)eggs=a*L^b;(3)eggs=a*Wt^b; (4)eggs=a+b*L; (5)eggs=a+b*W
0 #_hermaphroditism option:  0=none; 1=age-specific fxn
1 #_parameter_offset_approach (1=none, 2= M, G, CV_G as offset from female-GP1, 3=like SS2 V1.x)
1 #_env/block/dev_adjust_method (1=standard; 2=logistic transform keeps in base parm bounds; 3=standard w/ no bound check)
#
#_growth_parms
#_LO HI INIT PRIOR PR_type SD PHASE env-var use_dev dev_minyr dev_maxyr dev_stddev Block Block_Fxn
 0.01 1.8 0 0 -1 0.8 -3 -1 0 0 0 0.5 0 0 # NatM_p_1_Fem_GP_1
 0.01 1.8 0 0 -1 0.8 -3 -2 0 0 0 0.5 0 0 # NatM_p_2_Fem_GP_1
 0.01 1.8 0 0 -1 0.8 -3 -3 0 0 0 0.5 0 0 # NatM_p_3_Fem_GP_1
 0.01 1.8 0 0 -1 0.8 -3 -4 0 0 0 0.5 0 0 # NatM_p_4_Fem_GP_1
 0.01 1.8 0.2 0.2 -1 0.8 -3 0 0 0 0 0.5 0 0 # NatM_p_5_Fem_GP_1
 10 80 35.0851 30.8 -1 0.2 2 0 0 0 0 0.5 0 0 # L_at_Amin_Fem_GP_1
 25 250 106.204 200.1 -1 0.2 5 0 0 0 0 0.5 0 0 # L_at_Amax_Fem_GP_1
 0.01 2 0.135374 0.251 -1 0.8 2 0 0 0 0 0.5 0 0 # VonBert_K_Fem_GP_1
 0.25 4 1 1 0 0.8 -2 0 0 0 0 0.5 0 0 # Age_K_Fem_GP_1_a_2
 0.25 4 3.19611 1 0 0.8 2 0 0 0 0 0.5 0 0 # Age_K_Fem_GP_1_a_3
 0.25 4 0.61557 1 0 0.8 2 0 0 0 0 0.5 0 0 # Age_K_Fem_GP_1_a_4
 0.25 4 1.91847 1 0 0.8 2 0 0 0 0 0.5 0 0 # Age_K_Fem_GP_1_a_5
 0.25 4 0.734884 1 0 0.8 2 0 0 0 0 0.5 0 0 # Age_K_Fem_GP_1_a_6
 0.25 4 1.02755 1 0 0.8 2 0 0 0 0 0.5 0 0 # Age_K_Fem_GP_1_a_7
 0.25 4 1.06879 1 0 0.8 2 0 0 0 0 0.5 0 0 # Age_K_Fem_GP_1_a_8
 0.25 4 1.32637 1 0 0.8 2 0 0 0 0 0.5 0 0 # Age_K_Fem_GP_1_a_9
 0.25 4 1.1256 1 0 0.8 2 0 0 0 0 0.5 0 0 # Age_K_Fem_GP_1_a_10
 0.25 4 1.09201 1 0 0.8 2 0 0 0 0 0.5 0 0 # Age_K_Fem_GP_1_a_11
 0.25 4 1.58639 1 0 1.8 2 0 0 0 0 0.5 0 0 # Age_K_Fem_GP_1_a_12
 0.25 4 1.35917 1 0 1.8 2 0 0 0 0 0.5 0 0 # Age_K_Fem_GP_1_a_13
 0.25 4 1.03905 1 0 1.8 2 0 0 0 0 0.5 0 0 # Age_K_Fem_GP_1_a_14
 0.25 4 1.00191 1 0 1.8 2 0 0 0 0 0.5 0 0 # Age_K_Fem_GP_1_a_15
 0.01 0.5 0.2 0.42 -1 0.8 -3 0 0 0 0 0.5 0 0 # CV_young_Fem_GP_1
 0.01 0.5 0.2 0.2 -1 0.8 -5 0 0 0 0 0.5 0 0 # CV_old_Fem_GP_1
 0 3 2.0475e-005 2.0475e-005 -1 0.2 -3 0 0 0 0 0.5 0 0 # Wtlen_1_Fem
 2.5 3.5 2.857 2.98 -1 0.2 -3 0 0 0 0 0.5 0 0 # Wtlen_2_Fem
 25 35 28.9 0.879 -1 0.8 -3 0 0 0 0 0.5 0 0 # Mat50%_Fem
 -3 3 -0.42 -1.14 -1 0.8 -3 0 0 0 0 0.5 0 0 # Mat_slope_Fem
 -3 3 1 1 -1 0.8 -3 0 0 0 0 0.5 0 0 # Eggs/kg_inter_Fem
 -3 4 0 0 -1 0.8 -3 0 0 0 0 0.5 0 0 # Eggs/kg_slope_wt_Fem
 -4 4 0 0 -1 0 -4 0 0 0 0 0 0 0 # RecrDist_GP_1
 -4 4 0 0 -1 0 -4 0 0 0 0 0 0 0 # RecrDist_Area_1
 -4 4 0 0 -1 0 -4 0 0 0 0 0 0 0 # RecrDist_Seas_1
 -4 4 1 1 -1 0 -4 0 3 1964 2010 0.5 0 0 # CohortGrowDev
#
1 #_custom_MG-env_setup (0/1)
 -2 2 1 0 -1 99 -2 # NatM_p_1_Fem_GP_1_ENV_add
 -2 2 1 0 -1 99 -2 # NatM_p_2_Fem_GP_1_ENV_add
 -2 2 1 0 -1 99 -2 # NatM_p_3_Fem_GP_1_ENV_add
 -2 2 1 0 -1 99 -2 # NatM_p_4_Fem_GP_1_ENV_add
#
#_Cond 0  #custom_MG-block_setup (0/1)
#_Cond -2 2 0 0 -1 99 -2 #_placeholder when no MG-block parameters
#_Cond No MG parm trends 
#
#_seasonal_effects_on_biology_parms
 0 0 0 0 0 0 0 0 0 0 #_femwtlen1,femwtlen2,mat1,mat2,fec1,fec2,Malewtlen1,malewtlen2,L1,K
#_Cond -2 2 0 0 -1 99 -2 #_placeholder when no seasonal MG parameters
#
#DisplayOnly -0.0283546 # CohortGrowDev_DEVrwalk_1964
#DisplayOnly -0.027921 # CohortGrowDev_DEVrwalk_1965
#DisplayOnly -0.00769853 # CohortGrowDev_DEVrwalk_1966
#DisplayOnly -0.0348144 # CohortGrowDev_DEVrwalk_1967
#DisplayOnly 0.024816 # CohortGrowDev_DEVrwalk_1968
#DisplayOnly -0.0115681 # CohortGrowDev_DEVrwalk_1969
#DisplayOnly -0.0266095 # CohortGrowDev_DEVrwalk_1970
#DisplayOnly 0.0469508 # CohortGrowDev_DEVrwalk_1971
#DisplayOnly 0.108721 # CohortGrowDev_DEVrwalk_1972
#DisplayOnly -0.0327253 # CohortGrowDev_DEVrwalk_1973
#DisplayOnly -0.0875729 # CohortGrowDev_DEVrwalk_1974
#DisplayOnly -0.0269338 # CohortGrowDev_DEVrwalk_1975
#DisplayOnly 0.103011 # CohortGrowDev_DEVrwalk_1976
#DisplayOnly -0.0192614 # CohortGrowDev_DEVrwalk_1977
#DisplayOnly 0.00966106 # CohortGrowDev_DEVrwalk_1978
#DisplayOnly -0.133257 # CohortGrowDev_DEVrwalk_1979
#DisplayOnly 0.0475624 # CohortGrowDev_DEVrwalk_1980
#DisplayOnly 0.0606031 # CohortGrowDev_DEVrwalk_1981
#DisplayOnly -0.0348526 # CohortGrowDev_DEVrwalk_1982
#DisplayOnly -0.0930036 # CohortGrowDev_DEVrwalk_1983
#DisplayOnly -0.0285435 # CohortGrowDev_DEVrwalk_1984
#DisplayOnly 0.0419926 # CohortGrowDev_DEVrwalk_1985
#DisplayOnly -0.0160221 # CohortGrowDev_DEVrwalk_1986
#DisplayOnly 0.132825 # CohortGrowDev_DEVrwalk_1987
#DisplayOnly -0.0808965 # CohortGrowDev_DEVrwalk_1988
#DisplayOnly 0.193959 # CohortGrowDev_DEVrwalk_1989
#DisplayOnly 0.0201439 # CohortGrowDev_DEVrwalk_1990
#DisplayOnly -0.13269 # CohortGrowDev_DEVrwalk_1991
#DisplayOnly -0.0467171 # CohortGrowDev_DEVrwalk_1992
#DisplayOnly -0.0556397 # CohortGrowDev_DEVrwalk_1993
#DisplayOnly 0.0122748 # CohortGrowDev_DEVrwalk_1994
#DisplayOnly -0.138701 # CohortGrowDev_DEVrwalk_1995
#DisplayOnly -0.0629802 # CohortGrowDev_DEVrwalk_1996
#DisplayOnly 0.0991759 # CohortGrowDev_DEVrwalk_1997
#DisplayOnly 0.135335 # CohortGrowDev_DEVrwalk_1998
#DisplayOnly -0.189558 # CohortGrowDev_DEVrwalk_1999
#DisplayOnly 0.00975739 # CohortGrowDev_DEVrwalk_2000
#DisplayOnly 0.0671067 # CohortGrowDev_DEVrwalk_2001
#DisplayOnly 0.00505334 # CohortGrowDev_DEVrwalk_2002
#DisplayOnly 0.0885133 # CohortGrowDev_DEVrwalk_2003
#DisplayOnly 0.00424363 # CohortGrowDev_DEVrwalk_2004
#DisplayOnly -0.00617196 # CohortGrowDev_DEVrwalk_2005
#DisplayOnly 0.112747 # CohortGrowDev_DEVrwalk_2006
#DisplayOnly 0.121361 # CohortGrowDev_DEVrwalk_2007
#DisplayOnly 0.0665081 # CohortGrowDev_DEVrwalk_2008
#DisplayOnly -0.200588 # CohortGrowDev_DEVrwalk_2009
#DisplayOnly 0 # CohortGrowDev_DEVrwalk_2010
6 #_MGparm_Dev_Phase
#
#_Spawner-Recruitment
3 #_SR_function: 2=Ricker; 3=std_B-H; 4=SCAA; 5=Hockey; 6=B-H_flattop; 7=survival_3Parm
#_LO HI INIT PRIOR PR_type SD PHASE
 13 20 17 10.3 -1 10 1 # SR_LN(R0)
 0.2 1 0.964496 0.7 1 0.05 4 # SR_BH_steep
 0 2 0.7 0.8 -1 0.8 -4 # SR_sigmaR
 -5 5 0.1 0 -1 1 -3 # SR_envlink
 -5 5 -0.169901 0 -1 1 4 # SR_R1_offset
 0 0 0 0 -1 0 -99 # SR_autocorr
0 #_SR_env_link
0 #_SR_env_target_0=none;1=devs;_2=R0;_3=steepness
1 #do_recdev:  0=none; 1=devvector; 2=simple deviations
1955 # first year of main recr_devs; early devs can preceed this era
2011 # last year of main recr_devs; forecast devs start in following year
2 #_recdev phase 
1 # (0/1) to read 13 advanced options
 0 #_recdev_early_start (0=none; neg value makes relative to recdev_start)
 -4 #_recdev_early_phase
 0 #_forecast_recruitment phase (incl. late recr) (0 value resets to maxphase+1)
 1 #_lambda for Fcast_recr_like occurring before endyr+1
 1900 #_last_early_yr_nobias_adj_in_MPD
 1900 #_first_yr_fullbias_adj_in_MPD
 2010 #_last_yr_fullbias_adj_in_MPD
 2011 #_first_recent_yr_nobias_adj_in_MPD
 1 #_max_bias_adj_in_MPD (-1 to override ramp and set biasadj=1.0 for all estimated recdevs)
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
#DisplayOnly -0.0280539 # Main_InitAge_8
#DisplayOnly -0.00599562 # Main_InitAge_7
#DisplayOnly 0.049106 # Main_InitAge_6
#DisplayOnly 0.368923 # Main_InitAge_5
#DisplayOnly -0.212201 # Main_InitAge_4
#DisplayOnly -0.560151 # Main_InitAge_3
#DisplayOnly -0.0390276 # Main_InitAge_2
#DisplayOnly -0.864373 # Main_InitAge_1
#DisplayOnly -0.466825 # Main_RecrDev_1963
#DisplayOnly 0.255243 # Main_RecrDev_1964
#DisplayOnly 0.45336 # Main_RecrDev_1965
#DisplayOnly 0.277965 # Main_RecrDev_1966
#DisplayOnly -0.533705 # Main_RecrDev_1967
#DisplayOnly -1.11031 # Main_RecrDev_1968
#DisplayOnly 0.38459 # Main_RecrDev_1969
#DisplayOnly 0.817532 # Main_RecrDev_1970
#DisplayOnly -0.814579 # Main_RecrDev_1971
#DisplayOnly -0.159684 # Main_RecrDev_1972
#DisplayOnly 0.0264713 # Main_RecrDev_1973
#DisplayOnly 0.542827 # Main_RecrDev_1974
#DisplayOnly 0.384896 # Main_RecrDev_1975
#DisplayOnly 1.43877 # Main_RecrDev_1976
#DisplayOnly -0.0358882 # Main_RecrDev_1977
#DisplayOnly 1.21693 # Main_RecrDev_1978
#DisplayOnly 1.66261 # Main_RecrDev_1979
#DisplayOnly 0.530016 # Main_RecrDev_1980
#DisplayOnly 1.21537 # Main_RecrDev_1981
#DisplayOnly 0.348522 # Main_RecrDev_1982
#DisplayOnly 1.09188 # Main_RecrDev_1983
#DisplayOnly -0.604048 # Main_RecrDev_1984
#DisplayOnly 1.12521 # Main_RecrDev_1985
#DisplayOnly -0.0864913 # Main_RecrDev_1986
#DisplayOnly -0.446777 # Main_RecrDev_1987
#DisplayOnly 0.456535 # Main_RecrDev_1988
#DisplayOnly -0.705419 # Main_RecrDev_1989
#DisplayOnly -0.468324 # Main_RecrDev_1990
#DisplayOnly 0.585302 # Main_RecrDev_1991
#DisplayOnly -0.133122 # Main_RecrDev_1992
#DisplayOnly 1.18736 # Main_RecrDev_1993
#DisplayOnly 0.417146 # Main_RecrDev_1994
#DisplayOnly -0.31472 # Main_RecrDev_1995
#DisplayOnly 1.0601 # Main_RecrDev_1996
#DisplayOnly -1.15712 # Main_RecrDev_1997
#DisplayOnly -0.552619 # Main_RecrDev_1998
#DisplayOnly -0.0066168 # Main_RecrDev_1999
#DisplayOnly -1.0671 # Main_RecrDev_2000
#DisplayOnly -0.167679 # Main_RecrDev_2001
#DisplayOnly -1.171 # Main_RecrDev_2002
#DisplayOnly -0.455854 # Main_RecrDev_2003
#DisplayOnly -0.587544 # Main_RecrDev_2004
#DisplayOnly 0.323425 # Main_RecrDev_2005
#DisplayOnly -0.31163 # Main_RecrDev_2006
#DisplayOnly -0.450355 # Main_RecrDev_2007
#DisplayOnly -0.668736 # Main_RecrDev_2008
#DisplayOnly -0.445659 # Main_RecrDev_2009
#DisplayOnly -1.54351 # Main_RecrDev_2010
#DisplayOnly -0.0449621 # Main_RecrDev_2011
#
#Fishing Mortality info 
0.3 # F ballpark for tuning early phases
-2001 # F ballpark year (neg value to disable)
2 # F_Method:  1=Pope; 2=instan. F; 3=hybrid (hybrid is recommended)
4 # max F or harvest rate, depends on F_Method
# no additional F input needed for Fmethod 1
# if Fmethod=2; read overall start F value; overall phase; N detailed inputs to read
# if Fmethod=3; read N iterations for tuning for Fmethod 3
 0.4 1 0 # overall start F value; overall phase; N detailed inputs to read
#Fleet Year Seas F_value se phase (for detailed setup of F_Method=2)

#
#_initial_F_parms
#_LO HI INIT PRIOR PR_type SD PHASE
 0 2 0.865773 0.01 0 99 1 # InitF_1Fishery

# F rates for Fmethod=2
# 0.373519 F_fleet_1_YR_1963_s_1
# 0.32791 F_fleet_1_YR_1964_s_1
# 0.447589 F_fleet_1_YR_1965_s_1
# 0.473531 F_fleet_1_YR_1966_s_1
# 0.447915 F_fleet_1_YR_1967_s_1
# 0.494334 F_fleet_1_YR_1968_s_1
# 0.37001 F_fleet_1_YR_1969_s_1
# 0.452813 F_fleet_1_YR_1970_s_1
# 0.725824 F_fleet_1_YR_1971_s_1
# 0.787645 F_fleet_1_YR_1972_s_1
# 0.743955 F_fleet_1_YR_1973_s_1
# 0.816384 F_fleet_1_YR_1974_s_1
# 0.788309 F_fleet_1_YR_1975_s_1
# 0.996486 F_fleet_1_YR_1976_s_1
# 0.946062 F_fleet_1_YR_1977_s_1
# 0.7807 F_fleet_1_YR_1978_s_1
# 1.144 F_fleet_1_YR_1979_s_1
# 1.3818 F_fleet_1_YR_1980_s_1
# 0.975835 F_fleet_1_YR_1981_s_1
# 0.993713 F_fleet_1_YR_1982_s_1
# 0.790432 F_fleet_1_YR_1983_s_1
# 1.08462 F_fleet_1_YR_1984_s_1
# 0.799465 F_fleet_1_YR_1985_s_1
# 1.23347 F_fleet_1_YR_1986_s_1
# 0.873005 F_fleet_1_YR_1987_s_1
# 0.889417 F_fleet_1_YR_1988_s_1
# 1.07084 F_fleet_1_YR_1989_s_1
# 0.859945 F_fleet_1_YR_1990_s_1
# 0.81367 F_fleet_1_YR_1991_s_1
# 0.86907 F_fleet_1_YR_1992_s_1
# 0.708502 F_fleet_1_YR_1993_s_1
# 0.950457 F_fleet_1_YR_1994_s_1
# 1.05301 F_fleet_1_YR_1995_s_1
# 1.14269 F_fleet_1_YR_1996_s_1
# 0.987121 F_fleet_1_YR_1997_s_1
# 1.27911 F_fleet_1_YR_1998_s_1
# 1.3629 F_fleet_1_YR_1999_s_1
# 1.42925 F_fleet_1_YR_2000_s_1
# 0.95581 F_fleet_1_YR_2001_s_1
# 0.948153 F_fleet_1_YR_2002_s_1
# 1.14815 F_fleet_1_YR_2003_s_1
# 0.894557 F_fleet_1_YR_2004_s_1
# 1.14345 F_fleet_1_YR_2005_s_1
# 0.843872 F_fleet_1_YR_2006_s_1
# 0.941545 F_fleet_1_YR_2007_s_1
# 0.988352 F_fleet_1_YR_2008_s_1
# 0.91511 F_fleet_1_YR_2009_s_1
# 0.855724 F_fleet_1_YR_2010_s_1
# 0.795419 F_fleet_1_YR_2011_s_1
#
#_Q_setup
 # Q_type options:  <0=mirror, 0=float_nobiasadj, 1=float_biasadj, 2=parm_nobiasadj, 3=parm_w_random_dev, 4=parm_w_randwalk, 5=mean_unbiased_float_assign_to_parm
#_for_env-var:_enter_index_of_the_env-var_to_be_linked
#_Den-dep  env-var  extra_se  Q_type
 0 0 0 0 # 1 Fishery
 0 0 1 2 # 2 SURVEY
 0 0 1 2 # 3 RecruitSvy
#
#_Cond 0 #_If q has random component, then 0=read one parm for each fleet with random q; 1=read a parm for each year of index
#_Q_parms(if_any);Qunits_are_ln(q)
# LO HI INIT PRIOR PR_type SD PHASE
 0 0.9 0.0768743 0.01 -1 99 3 # Q_extraSD_2_SURVEY
 0 0.9 0.01 0.01 -1 99 -3 # Q_extraSD_3_RecruitSvy
 -25 5 -8.81463 -10 -1 10 1 # LnQ_base_2_SURVEY
 -25 5 -6.70339 -10 -1 10 -1 # LnQ_base_3_RecruitSvy
#
#_size_selex_types
#discard_options:_0=none;_1=define_retention;_2=retention&mortality;_3=all_discarded_dead
#_Pattern Discard Male Special
 1 1 0 0 # 1 Fishery
 0 0 0 0 # 2 SURVEY
 0 0 0 0 # 3 RecruitSvy
#
#_age_selex_types
#_Pattern ___ Male Special
 17 0 0 0 # 1 Fishery
 17 0 0 0 # 2 SURVEY
 17 0 0 0 # 3 RecruitSvy
#_LO HI INIT PRIOR PR_type SD PHASE env-var use_dev dev_minyr dev_maxyr dev_stddev Block Block_Fxn
 4 20 5.75131 6 0 0.5 2 0 0 0 0 0 0 0 # SizeSel_1P_1_Fishery
 0.1 10 8 1 0 99 -3 0 0 0 0 0 0 0 # SizeSel_1P_2_Fishery
 4 70 4 40 0 99 -3 0 0 0 0 0 0 0 # Retain_1P_1_Fishery
 0.1 10 0.1 1 0 99 -3 0 0 0 0 0 0 0 # Retain_1P_2_Fishery
 0.001 1 0.99 1 0 99 -3 0 0 0 0 0.5 0 0 # Retain_1P_3_Fishery
 -10 10 0 0 0 99 -3 0 0 0 0 0 0 0 # Retain_1P_4_Fishery
 -1000 10 -1000 -1000 -1 99 -2 0 0 0 0 0.5 0 0 # AgeSel_1P_1_Fishery
 -10 10 0 0 -1 99 -2 0 0 0 0 0.5 0 0 # AgeSel_1P_2_Fishery
 -10 10 1.25188 0 0 0.5 2 0 0 0 0 0.5 0 0 # AgeSel_1P_3_Fishery
 -10 10 -0.0169411 0 0 0.2 3 0 0 0 0 0.5 -1 0 # AgeSel_1P_4_Fishery
 -10 10 -0.115959 0 0 0.2 5 0 0 0 0 0.5 0 0 # AgeSel_1P_5_Fishery
 -10 10 -0.167069 0 0 0.2 5 0 0 0 0 0.5 0 0 # AgeSel_1P_6_Fishery
 -10 10 -0.181004 0 0 0.2 5 0 0 0 0 0.5 0 0 # AgeSel_1P_7_Fishery
 -10 10 0 0 0 0.5 -2 0 0 0 0 0.5 0 0 # AgeSel_1P_8_Fishery
 -10 10 0 0 0 0.5 -2 0 0 0 0 0.5 0 0 # AgeSel_1P_9_Fishery
 -10 10 0 0 0 0.5 -2 0 0 0 0 0.5 0 0 # AgeSel_1P_10_Fishery
 -10 10 0 0 0 0.5 -2 0 0 0 0 0.5 0 0 # AgeSel_1P_11_Fishery
 -10 10 0 0 0 0.5 -2 0 0 0 0 0.5 0 0 # AgeSel_1P_12_Fishery
 -10 10 0 0 0 0.5 -2 0 0 0 0 0.5 0 0 # AgeSel_1P_13_Fishery
 -10 10 0 0 0 0.5 -2 0 0 0 0 0.5 0 0 # AgeSel_1P_14_Fishery
 -10 10 0 0 0 0.5 -2 0 0 0 0 0.5 0 0 # AgeSel_1P_15_Fishery
 -10 10 0 0 0 0.5 -2 0 0 0 0 0.5 0 0 # AgeSel_1P_16_Fishery
 -1000 10 -1000 -1000 -1 99 -2 0 0 0 0 0.5 0 0 # AgeSel_2P_1_SURVEY
 -10 10 0 0 -1 99 -2 0 0 0 0 0.5 0 0 # AgeSel_2P_2_SURVEY
 -10 10 1.355 0 -1 99 2 0 0 0 0 0.5 0 0 # AgeSel_2P_3_SURVEY
 -10 10 0.46393 0 -1 99 2 0 0 0 0 0.5 0 0 # AgeSel_2P_4_SURVEY
 -10 10 0.328391 0 -1 99 2 0 0 0 0 0.5 0 0 # AgeSel_2P_5_SURVEY
 -10 10 0.285548 0 -1 99 2 0 0 0 0 0.5 0 0 # AgeSel_2P_6_SURVEY
 -10 10 0.268835 0 -1 99 2 0 0 0 0 0.5 0 0 # AgeSel_2P_7_SURVEY
 -10 10 0 0 -1 99 -2 0 0 0 0 0.5 0 0 # AgeSel_2P_8_SURVEY
 -10 10 0 0 -1 99 -2 0 0 0 0 0.5 0 0 # AgeSel_2P_9_SURVEY
 -10 10 0 0 -1 99 -2 0 0 0 0 0.5 0 0 # AgeSel_2P_10_SURVEY
 -10 10 0 0 -1 99 -2 0 0 0 0 0.5 0 0 # AgeSel_2P_11_SURVEY
 -10 10 0 0 -1 99 -2 0 0 0 0 0.5 0 0 # AgeSel_2P_12_SURVEY
 -10 10 0 0 -1 99 -2 0 0 0 0 0.5 0 0 # AgeSel_2P_13_SURVEY
 -10 10 0 0 -1 99 -2 0 0 0 0 0.5 0 0 # AgeSel_2P_14_SURVEY
 -10 10 0 0 -1 99 -2 0 0 0 0 0.5 0 0 # AgeSel_2P_15_SURVEY
 -10 10 0 0 -1 99 -2 0 0 0 0 0.5 0 0 # AgeSel_2P_16_SURVEY
 -1000 10 -1000 -1000 -1 99 -2 0 0 0 0 0.5 0 0 # AgeSel_3P_1_RecruitSvy
 -10 10 0 0 -1 99 -2 0 0 0 0 0.5 0 0 # AgeSel_3P_2_RecruitSvy
 -10 10 0.1 0 -1 99 -2 0 0 0 0 0.5 0 0 # AgeSel_3P_3_RecruitSvy
 -10 10 0.1 0 -1 99 -2 0 0 0 0 0.5 0 0 # AgeSel_3P_4_RecruitSvy
 -10 10 0.1 0 -1 99 -2 0 0 0 0 0.5 0 0 # AgeSel_3P_5_RecruitSvy
 -10 10 0.1 0 -1 99 -2 0 0 0 0 0.5 0 0 # AgeSel_3P_6_RecruitSvy
 -10 10 0.1 0 -1 99 -2 0 0 0 0 0.5 0 0 # AgeSel_3P_7_RecruitSvy
 -10 10 0 0 -1 99 -2 0 0 0 0 0.5 0 0 # AgeSel_3P_8_RecruitSvy
 -10 10 0 0 -1 99 -2 0 0 0 0 0.5 0 0 # AgeSel_3P_9_RecruitSvy
 -10 10 0 0 -1 99 -2 0 0 0 0 0.5 0 0 # AgeSel_3P_10_RecruitSvy
 -10 10 0 0 -1 99 -2 0 0 0 0 0.5 0 0 # AgeSel_3P_11_RecruitSvy
 -10 10 0 0 -1 99 -2 0 0 0 0 0.5 0 0 # AgeSel_3P_12_RecruitSvy
 -10 10 0 0 -1 99 -2 0 0 0 0 0.5 0 0 # AgeSel_3P_13_RecruitSvy
 -10 10 0 0 -1 99 -2 0 0 0 0 0.5 0 0 # AgeSel_3P_14_RecruitSvy
 -10 10 0 0 -1 99 -2 0 0 0 0 0.5 0 0 # AgeSel_3P_15_RecruitSvy
 -10 10 0 0 -1 99 -2 0 0 0 0 0.5 0 0 # AgeSel_3P_16_RecruitSvy
#_Cond 0 #_custom_sel-env_setup (0/1) 
#_Cond -2 2 0 0 -1 99 -2 #_placeholder when no enviro fxns
#_Cond 0 #_custom_sel-blk_setup (0/1) 
#_Cond -2 2 0 0 -1 99 -2 #_placeholder when no block usage
#_seltrend_parms 
 -10 10 0.0473021 0 -1 99 4 # AgeSel_1P_4_Fishery_TrendFinal_Offset
 1970 2009 1995.34 2000 0 5 4 # AgeSel_1P_4_Fishery_TrendInfl_
 0.5 20 2.75721 4 0 2 4 # AgeSel_1P_4_Fishery_TrendWidth_
# -0.820629 # Retain_1P_3_Fishery_DEVrwalk_1993
# -0.738523 # Retain_1P_3_Fishery_DEVrwalk_1994
# -0.707061 # Retain_1P_3_Fishery_DEVrwalk_1995
# -0.0668886 # Retain_1P_3_Fishery_DEVrwalk_1996
# 0.309739 # Retain_1P_3_Fishery_DEVrwalk_1997
# -0.076297 # Retain_1P_3_Fishery_DEVrwalk_1998
# -0.137264 # Retain_1P_3_Fishery_DEVrwalk_1999
# 0.0850417 # Retain_1P_3_Fishery_DEVrwalk_2000
# 0.1549 # Retain_1P_3_Fishery_DEVrwalk_2001
# 0.0025791 # Retain_1P_3_Fishery_DEVrwalk_2002
# -0.483674 # Retain_1P_3_Fishery_DEVrwalk_2003
# 0.260686 # Retain_1P_3_Fishery_DEVrwalk_2004
# -0.160058 # Retain_1P_3_Fishery_DEVrwalk_2005
# 0.0972438 # Retain_1P_3_Fishery_DEVrwalk_2006
# -0.00880021 # Retain_1P_3_Fishery_DEVrwalk_2007
# -0.0489344 # Retain_1P_3_Fishery_DEVrwalk_2008
# 0.120645 # Retain_1P_3_Fishery_DEVrwalk_2009
# 0.106888 # Retain_1P_3_Fishery_DEVrwalk_2010
# -0.0564513 # Retain_1P_3_Fishery_DEVrwalk_2011
#4 #_selparmdev-phase
#2 #_env/block/dev_adjust_method (1=standard; 2=logistic trans to keep in base parm bounds; 3=standard w/ no bound check)
#
# Tag loss and Tag reporting parameters go next
0  # TG_custom:  0=no read; 1=read if tags exist
#_Cond -6 6 1 1 2 0.01 -4 0 0 0 0 0 0 0  #_placeholder if no parameters
#
0 #_Variance_adjustments_to_input_values
#_fleet: 1 2 3 
#_Cond  0 0 0 #_add_to_survey_CV
#_Cond  0 0 0 #_add_to_discard_stddev
#_Cond  0 0 0 #_add_to_bodywt_CV
#_Cond  1 1 1 #_mult_by_lencomp_N
#_Cond  1 1 1 #_mult_by_agecomp_N
#_Cond  1 1 1 #_mult_by_size-at-age_N
#
4 #_maxlambdaphase
1 #_sd_offset
#
2 # number of changes to make to default Lambdas (default value is 1.0)
# Like_comp codes:  1=surv; 2=disc; 3=mnwt; 4=length; 5=age; 6=SizeFreq; 7=sizeage; 8=catch; 
# 9=init_equ_catch; 10=recrdev; 11=parm_prior; 12=parm_dev; 13=CrashPen; 14=Morphcomp; 15=Tag-comp; 16=Tag-negbin
#like_comp fleet/survey  phase  value  sizefreq_method
 1 3 1 0 1
 9 1 4 0.01 1
#
# lambdas (for info only; columns are phases)
#  0 0 0 0 #_CPUE/survey:_1
#  1 1 1 1 #_CPUE/survey:_2
#  0 0 0 0 #_CPUE/survey:_3
#  1 1 1 1 #_discard:_1
#  0 0 0 0 #_discard:_2
#  0 0 0 0 #_discard:_3
#  1 1 1 1 #_agecomp:_1
#  1 1 1 1 #_agecomp:_2
#  0 0 0 0 #_agecomp:_3
#  1 1 1 1 #_size-age:_1
#  0 0 0 0 #_size-age:_2
#  0 0 0 0 #_size-age:_3
#  1 1 1 0.01 #_init_equ_catch
#  1 1 1 1 #_recruitments
#  1 1 1 1 #_parameter-priors
#  1 1 1 1 #_parameter-dev-vectors
#  1 1 1 1 #_crashPenLambda
0 # (0/1) read specs for more stddev reporting 
 # 0 1 -1 5 1 5 1 -1 5 # placeholder for selex type, len/age, year, N selex bins, Growth pattern, N growth ages, NatAge_area(-1 for all), NatAge_yr, N Natages
 # placeholder for vector of selex bins to be reported
 # placeholder for vector of growth ages to be reported
 # placeholder for vector of NatAges ages to be reported
999

