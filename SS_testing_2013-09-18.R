

# File structure
#PlatformFile = "C:/Users/James.Thorson/"     # Laptop
PlatformFile = "c:/"                         # KFJ
if(grepl("Users", PlatformFile)) {
  #SourceFile = MiscFile = PlatformFile = "C:/Users/thorsonja/Dropbox/Laplace approximation/"   # NWFSC desktop
    #SourceFile = MiscFile = PlatformFile = "C:/Users/James.Thorson/Dropbox/Laplace approximation/"  # Laptop dropbox
      #SourceFile = MiscFile = PlatformFile = "C:/Users/James Thorson/Dropbox/Laplace approximation/"  # Laptop dropbox
MiscFile = paste(PlatformFile,"Desktop/UW Hideaway/Software/R/Thorson/",sep="")
SourceFile = paste(PlatformFile,"Desktop/UW Hideaway/Laplace approximation/",sep="")
SsFile = paste(SourceFile,"Simulation files/",sep="")
r4ss.name <- paste0(MiscFile,"R4SS/")
  } else{
      MiscFile = paste(PlatformFile, "~/ss3/", sep="")
      SourceFile = paste(PlatformFile, "laplace/", sep="")
      SsFile = paste(SourceFile,"Simulation files/",sep="") 
      r4ss.name <- NULL      
  }

# Load packages
library(mvtnorm)
library(pso)                                              
library(corpcor) # pseudoinverse
library(snowfall) # pseudoinverse
library(r4ss) # pseudoinverse

update_r4ss_files(local = r4ss.name, save=FALSE)
# integrate(f=function(x){DtruncnormFn(x,truncmean=0.6,truncsd=0.1,a=0.2,b=1)},lower=0.2,upper=1)

# Load scripts
source(paste0(SourceFile,"Fn_Laplace_Approx_2013-08-15b.R"))
source(paste0(SourceFile,"Fn_SS_testing_2013-08-26.R"))

# Define species
Species = c("NS_COD", "Hake_33", "Hake_2012", "Hake_2012_Explicit_F")[1]
Version = 1 # 1=Determinant of whole Hessian; 5=Determinant of just a subset of Hessian
ReRun = FALSE # if TRUE, then always runs; if FALSE, then only runs if it hasn't already
RandomSeed = ceiling(runif(1, min=1, max=1e6))

# File structure
Date = Sys.Date()
  #Date = "2013-08-12_200_rep_SigmaR=1.5"
  #Date = "2013-08-13_SigmaR=0.5"
  #Date = "2013-08-14_SigmaK=0.02"
  DateFile = paste(SourceFile,"SS_testing_",Species,"--",Date,"/",sep="")
  dir.create(DateFile)

# Save stuff
RecordList = list(RandomSeed=RandomSeed, Version=Version, ReRun=ReRun)
capture.output(RecordList, file=paste(DateFile,"RecordList.txt",sep=""))

# Settings 
Nreps = 200
Ncores = 2

# Copy files to DateFile
if(Species=="NS_COD"){
  SsFile = paste(SourceFile,"NS_COD files/Simulation files/",sep="")
  file.copy(from=paste(SsFile,c("ss3.exe","Starter.SS","ss3_0.par","SimTest.ctl","SimTest.dat","Forecast.ss"),sep=""), to=paste(DateFile,c("ss3.exe","Starter.SS","ss3_0.par","SimTest.ctl","SimTest.dat","Forecast.ss"),sep=""), overwrite=TRUE)  
  FileSet_From=c("ss3.exe","Starter.SS","ss3_0.par","SimTest.ctl","SimTest.dat","Forecast.ss")
  FileSet_To=c("ss3.exe","Starter.SS","ss3_init.par","SimTest.ctl","SimTest.dat","Forecast.ss")
  K_set_in_par = 40:86
  SigmaK = 0.02               # SHOULD BE 0.1
  RecDev_set_in_par = 93:149
  SigmaR = 0.7              # SHOULD BE 0.7
  F_set_in_par = 150:199
  SigmaF = 0.1              # SHOULD BE 0.1
  Fmin = 0.1
  Fmax = 0.6
  Flim = 0.8
  Input_SD = 0.5
  Interval = c(0.01, 0.5)
  CTL_linenum_List = list(77)
  ESTPAR_num_List = list(17:63)
  PAR_num_Vec = NA
  ReDoBiasRamp = TRUE
  BiasRamp_linenum_Vec = 162:166
  CTL_linenum_Type = NULL 
  StartFromPar = TRUE 
  NLL_threshold = 0  # Should be zero
  Gradient_threshold = 1  # Should be one
  VarianceAdjustment_in_CTL = NA # Variance adjustment is turned-off
}
if(Species=="Hake_2012" | Species=="Hake_2012_Explicit_F"){
  Type = c("With_Age","Without_Age")[1]
  if(Species=="Hake_2012" & Type=="With_Age") SsFile = paste(SourceFile,"Hake/81_2012_revised_acoustic_base_mle/",sep="")      # With marginal-age data
  if(Species=="Hake_2012" & Type=="Without_Age") SsFile = paste(SourceFile,"Hake/81_2012_revised_acoustic_base_no_age/",sep="")     # Without marginal-age data
  if(Species=="Hake_2012_Explicit_F") SsFile = paste(SourceFile,"Hake/81_2012_revised_acoustic_base_mle_explicit_F/",sep="")      # With marginal-age data
  file.copy(from=paste(SsFile,c("forecast.ss","2012_hake_control.ss","2012_hake_data.ss","SS3.exe","ss3.par","starter.ss","wtatage.ss"),sep=""), to=paste(DateFile,c("forecast.ss","2012_hake_control.ss","2012_hake_data.ss","ss3.exe","ss3_0.par","starter.ss","wtatage.ss"),sep=""), overwrite=TRUE)
  FileSet_From=c("forecast.ss","2012_hake_control.ss","2012_hake_data.ss","ss3.exe","ss3_0.par","starter.ss","wtatage.ss")
  FileSet_To=c("forecast.ss","2012_hake_control.ss","2012_hake_data.ss","ss3.exe","ss3_init.par","starter.ss","wtatage.ss")
  K_set_in_par = vector(length=0)
  SigmaK = 0.0
  RecDev_set_in_par = 24:89
  SigmaR = 0.5         ########## SHOULD BE 1.5
  #Interval = c(0.5, 3.0)
  Interval = c(0.1, 1.5)
  if(Species=="Hake_2012"){
    F_set_in_par = vector(length=0)
    SigmaF = 0.1
    Fmin = 0.02
    Fmax = 0.2
    Flim = 0.4
  }
  if(Species=="Hake_2012_Explicit_F"){
    F_set_in_par = 97:142
    SigmaF = 0.1
    Fmin = 0.1
    Fmax = 0.5
    Flim = 0.7
  }
  Input_SD = 0.5
  CTL_linenum_List = list(67)
  ESTPAR_num_List = list(4:72)
  PAR_num_Vec = c(20)
  ReDoBiasRamp = TRUE                      
  BiasRamp_linenum_Vec = 85:89
  CTL_linenum_Type = NULL
  StartFromPar = TRUE       ######### SHOULD BE TRUE 
  if(Type=="With_Age"){
    Gradient_threshold = 1  # Should be one
    NLL_threshold = 0 # Should be zero
  }
  if(Type=="Without_Age"){
    Gradient_threshold = 5  # Should be one
    NLL_threshold = -Inf # Should be zero
  }
  VarianceAdjustment_in_CTL = 188:193
}

###############
# Run replicated simulations
###############
if(Ncores>1){
  sfInit(parallel=TRUE, cpus=Ncores)
  sfCpus()
  sfExportAll()
  sfLibrary(mvtnorm)
  sfLibrary(corpcor)
  sfLibrary(r4ss)
  sfSapply(x=sample(1:Nreps), fun=SS_testing, Version=Version)
}else{
  for(RepI in 1:Nreps) SS_testing(RepI, Version=Version)
}



  