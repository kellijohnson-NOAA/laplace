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
## Variable inputs
###############################################################################
Species <- "coda"


###############################################################################
## Step 02
## Set the PlatformFile
###############################################################################  
# File Structure
  PlatformFile <- "c:/laplace/"                   # KFJ
  sim.fold <- paste0(PlatformFile, "cod_sim/")
  out.fold <- paste0(sim.fold,"output/")
  om.fold <- paste0(PlatformFile, Species, "_om/")
  DateFile <- paste0(out.fold, dir(out.fold, pattern = "comp"), "/")
  
###############################################################################
## Step 03
## Read in the data
###############################################################################

solution <- do.call("rbind", 
  lapply(DateFile, function(x) {
      old.wd <- getwd()
      f.list <- dir(x, pattern = "Rep")
      fold   <- paste0(x, f.list, "/")
      temp <- 
      sapply(paste0(fold, "Optimization_solution.txt"), function(p) {
          use.me <- readLines(p)
          use.me <- strsplit(grep("\\[", use.me, value = TRUE), " ")
          use.me <- as.numeric(do.call("rbind", use.me)[,2])
        }
      )
      print.me <- cbind(temp[1,], temp[2,])
      rownames(print.me) <- rep(x, nrow(print.me))
      colnames(print.me) <- c("par", "negll")
      print.me
    }
  )
)

sol.mean <- tapply(solution[, "par"], rownames(solution), mean)
sol.sd   <- tapply(solution[, "par"], rownames(solution), sd)

n.samp <- lapply(paste0(DateFile, "ss3.dat"), function(x) {
    temp <- SS_readdat(x, verbose = FALSE)$agecomp
    temp <- subset(temp, select=c("Yr", "FltSvy", "Nsamp"))
    temp
  }
)

png("c:/laplace/data.png")
plot(n.samp[[3]]$Yr, n.samp[[3]]$Nsamp, col = "black", pch = n.samp[[3]]$FltSvy+20, cex = 2, 
     xlab = "year", ylab = "number of age comp samples", las = 1, ylim = c(0, max(n.samp[[1]]$Nsamp)))
points(n.samp[[1]]$Yr, n.samp[[1]]$Nsamp, col = "red", pch = n.samp[[1]]$FltSvy+20, cex = 2)
lines(1906:1925, rep(-1,20), lwd = 4)
legend("topleft", c("fishery", "survey", "fishing"), pch = c(21, 22, 32),
       lty = c(0,0,1), bty = "n")
dev.off()

png("c:/laplace/results.png", width = 4, height = 4, units = "in", res = 400)
plotCI(c(1,2), sol.mean[c(1,3)], sol.sd[c(1,3)], xaxt = "n", las = 1,
       xlab = "fewer versus more sampling years", 
       ylab = "mean sigma for natural mortality", 
       cex = 4, lwd = 4, col = c("red", "black"), xlim = c(0,3))
abline(h = 0.10, col = "blue", lwd = 4)
legend("topright", "true sigma", lty = 1, col = "blue", bty = "n")
dev.off()
