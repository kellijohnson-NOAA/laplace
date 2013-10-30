
SS_testing = function(RepI, Version){
  # Make directory for replicate
  RepFile = paste0(DateFile, "Rep=", RepI, "/")
  dir.create(RepFile)

  # Only re-generate the data if necessary
  file.copy(from = paste0(DateFile, FileSet_From), 
            to   = paste0(RepFile,  FileSet_To), overwrite = TRUE)
  if(ReRun == TRUE | 
     !all(c("BootData.ss", "ss3_0.par") %in% list.files(RepFile))) {
    # Copy files
    STARTER = SS_readstarter(paste0(DateFile,"Starter.SS"), verbose = FALSE)
    
    # Deal with randomization 
    # (i.e. its necessary to de-synchronize generating random values
    # using ADMB and parallelized code)
    set.seed(RandomSeed+RepI)
    Sleep <- runif(1, min = 1e-6, max = 3)
    Sys.sleep(Sleep)
                                                                         
    # Loop until generate valid data set
    Converged = FALSE
    while(Converged==FALSE){
      # Change starter file
      Starter = SS_readstarter(paste0(RepFile,"Starter.SS"), verbose = FALSE)
        Starter$init_values_src = 1
        Starter$last_estimation_phase = 0
        Starter$N_bootstraps = 3
        file.remove(paste(RepFile,"Starter.SS",sep=""))
      SS_writestarter(Starter, dir=RepFile, verbose=FALSE, overwrite=TRUE)
    
      # Change PAR file
      Par = scan(paste(RepFile,"ss3_init.par",sep=""), comment.char="#")
        # Change K
        Par[K_set_in_par] = rnorm(length(K_set_in_par), mean=0, sd=SigmaK)
        # Change RecDev
        Par[RecDev_set_in_par] = rnorm(length(RecDev_set_in_par), mean=0, sd=SigmaR)
      # Write file
      write(Par, file=paste(RepFile, "ss3.par", sep=""))
      Sys.sleep(1)
      file.copy(from=paste(RepFile, "ss3.par", sep=""), to=paste(RepFile,"ss3_true.par",sep=""), overwrite=TRUE)
    
      # Change variance adjustments to default
      if( !is.na(VarianceAdjustment_in_CTL[1]) ){
        CTL = readLines(paste(RepFile,STARTER$ctlfile,sep=""))
        VarAdj_Lines = CTL[VarianceAdjustment_in_CTL]
        for(i in seq_along(VarAdj_Lines)){
          Temp = as.vector(unlist(sapply(CTL[VarianceAdjustment_in_CTL[i]], FUN=function(Char){strsplit(Char," ")[[1]]})))
          Temp = as.vector(unlist(sapply(Temp, FUN=function(Char){strsplit(Char,"\t")[[1]]})))      
          Temp = Temp[which(Temp!="")]
          if(length(grep("#",Temp))>=1) Temp = Temp[-(grep("#",Temp):length(Temp))]
          Temp = as.numeric(Temp)
          if(i<=3) CTL[VarianceAdjustment_in_CTL[i]] = paste(rep(0, length(Temp)), sep=" ", collapse=" ")
          if(i>=4) CTL[VarianceAdjustment_in_CTL[i]] = paste(rep(1, length(Temp)), sep=" ", collapse=" ")
        }
        writeLines(CTL, paste(RepFile,STARTER$ctlfile,sep=""))    
      }
      
      # Run SS (to generate bootstrap replicate)
      setwd(RepFile)
      shell("ss3")
      Sys.sleep(1)
    
      # Revert variance adjustments to their original values
      if( !is.na(VarianceAdjustment_in_CTL[1]) ){
        CTL = readLines(paste(RepFile,STARTER$ctlfile,sep=""))
        CTL[VarianceAdjustment_in_CTL] = VarAdj_Lines
        writeLines(CTL, paste(RepFile,STARTER$ctlfile,sep=""))    
      }
  
      # write to new DAT file
      SS_splitdat(inpath=RepFile, outpath=RepFile, MLE=FALSE)
  
      # Change starter file
      Starter = SS_readstarter(paste0(RepFile,"Starter.SS"), verbose = FALSE)
        Starter$init_values_src = 1
        Starter$last_estimation_phase = 25
        Starter$N_bootstraps = 0
        Starter$datfile = "BootData.ss"
        file.remove(paste(RepFile,"Starter.SS",sep=""))
      SS_writestarter(Starter, dir=RepFile, verbose=FALSE, overwrite=FALSE)
    
      # Run SS (to estimate model for first time without hessian)
      setwd(RepFile)
      shell("ss3 -nohess")
      Sys.sleep(1)
            
      # Check for (1) convergence, and (2) a realistic level of log-likelihood given the quantity of data
      if("ss3.par" %in% list.files(RepFile)){
        PAR = scan(paste(RepFile,"ss3.par",sep=""), what="character", quiet=TRUE)
        if(!is.na(as.numeric(PAR[11])) && as.numeric(PAR[11])>NLL_threshold && as.numeric(PAR[16])<Gradient_threshold) Converged = TRUE
      }
    } # End While-Loop
    file.copy(from=paste(RepFile,"ss3.par",sep=""), to=paste(RepFile,"ss3_0.par",sep=""), overwrite=TRUE)
  }
  
  # Only re-estimate the model if necessary
  if(ReRun == TRUE | !("Save_Rep.RData" %in% list.files(RepFile))){
    # Explore Methot and Taylor (2011) algorithm for SigmaR
    if(FALSE){
      # Run SS
      shell("ss3")
      Sys.sleep(1)
      SsOutput = SS_output(RepFile, covar = TRUE, forecast = FALSE)
      # Analyze
      Param = SsOutput$parameters
      RecrDev = Param[grep("RecrDev", Param$Label), ]
      RMSE = sqrt(mean(RecrDev[, 'Value']^2))
      c('Fixed'=Param[grep("SR_sigmaR",Param$Label),'Value'], 'RMSE'=RMSE)
      # Methot and Taylor Eq. 11 -- sigmaR^2 = SE(r)^2 + SD(r)^2
      sqrt(sd(RecrDev[, 'Value'])^2 + mean(RecrDev[, 'Parm_StDev'])^2)
    }
    
    # Run Laplace Approximation function
    if(TRUE){
      Iteration = 0
      assign("Iteration", Iteration, envir=.GlobalEnv)
        # Input_SD_Group_Vec=Input_SD; File=RepFile; CTL_linenum_List=CTL_linenum_List; ESTPAR_num_List=ESTPAR_num_List; PAR_num_Vec=PAR_num_Vec; Int_Group_List=list(1); Version=Version; StartFromPar=StartFromPar; Intern=TRUE; ReDoBiasRamp=ReDoBiasRamp; BiasRamp_linenum_Vec=BiasRamp_linenum_Vec; CTL_linenum_Type=CTL_linenum_Type
      NegLogInt_Fn(Input_SD_Group_Vec=Input_SD, File=RepFile, CTL_linenum_List=CTL_linenum_List, ESTPAR_num_List=ESTPAR_num_List, PAR_num_Vec=PAR_num_Vec, Int_Group_List=list(1), Version=Version, StartFromPar=StartFromPar, Intern=TRUE, ReDoBiasRamp=ReDoBiasRamp, BiasRamp_linenum_Vec=BiasRamp_linenum_Vec, CTL_linenum_Type=CTL_linenum_Type)
    }
  
    # Optimize
    Iteration = 0
    assign("Iteration", Iteration, envir=.GlobalEnv)
    if("Optimization_record.txt" %in% list.files(RepFile)) file.remove(paste(RepFile,"Optimization_record.txt",sep=""))
      # f=NegLogInt_Fn; interval=Interval; maximum=FALSE; CTL_linenum_List=CTL_linenum_List; ESTPAR_num_List=ESTPAR_num_List; Int_Group_List=list(1); PAR_num_Vec=PAR_num_Vec; Version=Version; Intern=FALSE
      Opt = try( optimize(f=NegLogInt_Fn, interval=Interval, File=RepFile, maximum=FALSE, CTL_linenum_List=CTL_linenum_List, ESTPAR_num_List=ESTPAR_num_List, Int_Group_List=list(1), PAR_num_Vec=PAR_num_Vec, Version=Version, StartFromPar=StartFromPar, Intern=FALSE, ReDoBiasRamp=ReDoBiasRamp, BiasRamp_linenum_Vec=BiasRamp_linenum_Vec, CTL_linenum_Type=CTL_linenum_Type) ) 
    
    # If successful
    if( class(Opt)!='try-error' ){
      # Save stuff
      capture.output(Opt, file=paste(RepFile,"Optimization_solution.txt",sep=""))                   
      Save_Rep = list(Opt=Opt)
      save(Save_Rep, file=paste(RepFile,"Save_Rep.RData",sep=""))
      # Read into R4SS and make plots
      #SsOutput = SS_output(RepFile, covar=TRUE, forecast=FALSE)
      #SS_plots(SsOutput, png=TRUE, uncertainty=TRUE, aalresids=TRUE, datplot=FALSE) #aalyear, aalbin
    }
  }
}
