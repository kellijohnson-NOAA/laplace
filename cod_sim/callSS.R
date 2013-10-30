###############################################################################
###############################################################################
#-----------------------------------------------------------------------------#
####Author      : Johnson, Kelli Faye, and Thorson, James A.
####Contact     : kellifayejohnson@gmail.com
####Lastupdate  : 2013-10-25
####Purpose     : CJFAS simulation - cod M w/ Laplace approximation using SS3
####Packages    : corpcor; mvtnorm; pso; r4ss; snowfall;
####Inputs      : "Fn_Laplace_Approx_2013-08-15b.R"
####Inputs      : "Fn_SS_testing_2013-08-26.R"
####Outputs     : 
####Remarks     : Character width == 80
#-----------------------------------------------------------------------------#
###############################################################################
###############################################################################

###############################################################################
## Step 01
## Set the variable inputs
###############################################################################
 
# Simulation specs
  Nreps   <- 2     # Default was 200
  Version <- 1     # 1 == Determinant of whole Hessian;
                   # 5 == Determinant of just a subset of Hessian
  ReRun   <- FALSE # TRUE == always runs; FALSE == only runs if it hasn't already

  Species <- "coda"
  design <- c("comp30.mr", "comp15.mr", "comp30.mt", "comp15.mt")
  ss.param <- "NatM_p_1_Fem_GP_1"
  
###############################################################################
## Step 02
## Set the PlatformFile
###############################################################################  
# File Structure
  PlatformFile <- "c:/laplace/"                   # KFJ
  sim.fold <- paste0(PlatformFile, "cod_sim/")
  out.fold <- paste0(sim.fold,"output/")
  om.fold <- paste0(PlatformFile, Species, "_om/")
  DateFile <- paste0(out.fold, design, "/")
  
###############################################################################
## Step 03
## Load the necessary packages
###############################################################################
  
# Check to make sure all packages are installed
  needed.packages <- c("corpcor", "mvtnorm", "pso", "r4ss", "snowfall",
                       "ss3sim")
  is.installed <- is.element(needed.packages, installed.packages()[,1])
  for ( i in seq_along(needed.packages)) {
    if(is.installed[i]) {
      require(needed.packages[i], quietly = TRUE, character.only = TRUE)
    } else {
        install.packages(needed.packages[i])
        require(needed.packages[i], quietly = TRUE, character.only = TRUE)
      }
  }
    
###############################################################################
## Step 04
## Load the necessary scripts
###############################################################################
  
# Load scripts
  source(paste0(sim.fold, "Fn_Laplace_Approx_2013-08-15b.R"))
  source(paste0(sim.fold, "Fn_SS_testing_2013-08-26.R"))  
  
###############################################################################
## Step 05
## Random tasks prior to the start of the simulation
###############################################################################
  Date <- Sys.Date()
  Ncores  <- as.numeric(Sys.getenv('NUMBER_OF_PROCESSORS'))  
  RandomSeed <- ceiling(runif(1, min=1, max=1e6))
  
  sapply(DateFile, dir.create, showWarnings = FALSE, recursive = TRUE)

# Save stuff
  RecordList <- list(RandomSeed = RandomSeed,
                     Version = Version,
                     ReRun = ReRun,
                     Machine = Sys.info()["nodename"])
  sapply(paste0(DateFile, "RecordList.txt"), function(x) {
    capture.output(RecordList, file = x)
    }
  )
  
###############################################################################
## Step 06
## Copy files and set up the simulation based on the species set in Step 01
###############################################################################

# Determine the number of years for the simulation
 years <- scan(paste0(om.fold, "ss3.dat"), 
               comment.char = "#", nmax=2, quiet= TRUE)

# Set up the data files for the given sampling schemes
# 30 years
  curr.folds <- DateFile[grep("30", DateFile)]
  for(q in seq(curr.folds)) {
    change_agecomp(SS_readdat(paste0(om.fold, "ss3.dat"),
                              verbose = FALSE),
                   paste0(curr.folds[q], "ss3.dat"),
                   fleets = c(1, 2), 
                   Nsamp = list(rep(40, 20), 
                                c(rep(10, 5), rep(40, 20), rep(100, 5))),
                   years = list(1906:1925, 1901:1930), cpar = c(1, 1), 
                   agebin_vector = NULL, write_file = TRUE)
  }
# 15 years
  curr.folds <- DateFile[grep("15", DateFile)]
  for(q in seq(curr.folds)) {
    change_agecomp(SS_readdat(paste0(om.fold, "ss3.dat"),
                              verbose = FALSE),
                   file.path(curr.folds[q], "ss3.dat"),
                   fleets = c(1, 2), 
                   Nsamp = list(rep(80, 10), 
                                rep(c(20, 50, 200), each = 5)),
                   years = list(1916:1925, 1916:1930), cpar = c(1, 1), 
                   agebin_vector = NULL, write_file = TRUE)
  }
  
# Configure the control file for the given parameter
    ss3sim.ctl <- readLines(paste0(om.fold, "om.ctl"))
      param.line <- grep(ss.param, ss3sim.ctl)
      CTL_linenum_List = list(param.line)

      BiasRamp_linenum_Vec <- seq(
        grep("last_early_yr_nobias_adj_in_MPD",ss3sim.ctl),
        grep("#_max_bias_adj_in_MPD", ss3sim.ctl))
      param <- unlist(strsplit(ss3sim.ctl[param.line], " "))
      param <- param[-nchar(param) < 0]
      param[9] <- 3
      param[10] <- years[1]
      param[11] <- years[2]
      param[12] <- 0.5
      ss3sim.ctl[param.line] <- paste(param, collapse = " ")
      ss3sim.ctl[grep("MGparm_Dev_Phase", ss3sim.ctl)] <- 
        "4 #_MGparm_Dev_Phase"
      sapply(paste0(DateFile, "om.ctl"), function(x) {
        writeLines(ss3sim.ctl, x)
        }
      )
      
# Configure the par file for the given parameter
    ss3sim.par <- readLines(paste0(om.fold, "ss3.par"))
      sapply(paste0(DateFile, "ss3_0.par"), function(x) {
        writeLines(ss3sim.par, x)
        }
      )
      ## Find the specific elements of par file, where .par is read in
      ## scan("ss3.par", comment.char="#")
      # '# MGparm_dev'
      num.mgparm <- length(grep("MGparm", ss3sim.par))
      K_set_in_par <- seq(num.mgparm + 2,
                          num.mgparm + diff(years) + 2)
      # '# recdev1:'
      num.srparm <- length(grep("SR\\_parm", ss3sim.par)) 
      RecDev_set_in_par <- seq(num.srparm + num.mgparm + diff(years) + 2,
                               num.srparm + num.mgparm + diff(years)*2 + 3)
      # '# init_F[x]:'
      F_set_in_par <- seq(rev(RecDev_set_in_par)[1] + 1,
                          rev(RecDev_set_in_par)[1] + diff(years) + 1)

                
    
    files.from <- c(file.path(om.fold, c("starter.ss", 
                                        "forecast.ss")))
    sapply(DateFile, function(x) {
        files.to   <- paste0(x, c("Starter.SS",
                                  "forecast.ss"))
        file.copy(from = files.from, to = files.to, overwrite = TRUE)
      }
    )           
  
    FileSet_From <- c("Starter.SS", "ss3_0.par",
                      "om.ctl", "ss3.dat", "Forecast.ss")
    FileSet_To   <- c("Starter.SS", "ss3_init.par",
                      "om.ctl", "ss3.dat", "Forecast.ss")
    SigmaK = 0.1                # SHOULD BE 0.1
    SigmaR = 0.4                # SHOULD BE 0.4
    SigmaF = 0.1                # SHOULD BE 0.1
    Fmin = 0.1
    Fmax = 0.6
    Flim = 0.8
    Input_SD = 0.5
    Interval = c(0.01, 0.5)

    PAR_num_Vec = NA
    ReDoBiasRamp = TRUE
    CTL_linenum_Type = NULL 
    StartFromPar = TRUE 
    NLL_threshold = 0              # Should be zero
    Gradient_threshold = 1         # Should be one
    VarianceAdjustment_in_CTL = NA # Variance adjustment is turned-off

    # Determine where the parameter and the mg deviation vector lie in 
    # the .std file.
    std.fold <- paste0(om.fold, "getstd/")
    std.from <- dir(om.fold, full.names = TRUE)
    dir.create(std.fold)
    curr.wd <- getwd()
    setwd(om.fold)
    file.copy(std.from, std.fold)
    setwd(std.fold)
    system("SS3 -nohess", invisible = TRUE, ignore.stdout = TRUE, 
           show.output.on.console = FALSE)
    start.ss <- SS_readstarter(paste0(std.fold,"Starter.SS"), verbose = FALSE)
        start.ss$last_estimation_phase <- 20
        start.ss$N_bootstraps <- 0
      SS_writestarter(start.ss, dir=std.fold, verbose=FALSE, overwrite=TRUE)
    SS_splitdat(inpath=std.fold, outpath=std.fold, 
                outpattern = "ss3.dat", verbose = FALSE)
    file.copy("ss3.dat.ss", "ss3.dat", overwrite = TRUE)
    system("SS3", invisible = TRUE, ignore.stdout = TRUE, 
           show.output.on.console = FALSE)
    std.file <- readLines("ss3.std")[-1]
    results <- SS_output(getwd(), verbose = FALSE, 
                         warn = FALSE, printstats = FALSE)
    ESTPAR_num_List <- list(grep(ss.param,
                                 results[["stdtable"]][,"name"], 
                                 fixed = TRUE))
    
    setwd(curr.wd)
    unlink(std.fold, recursive = TRUE)
###############################################################################
## Step 07
## Run replicated simulations
###############################################################################
DateFile <- DateFile[1]
# Run the simulation based on the number of cores (Ncores) set in step 01
  if(Ncores > 1) {
    sfInit(parallel = TRUE, cpus = Ncores)
    sfCpus()
    sfExportAll()
    sfLibrary(mvtnorm)
    sfLibrary(corpcor)
    sfLibrary(r4ss)
    sfSapply(x = sample(1 : Nreps), fun = SS_testing, Version = Version)
  }else{
    for(RepI in 1 : Nreps) SS_testing(RepI, Version = Version)
  }

###############################################################################
## Step 08
## End of the simulation
###############################################################################

  