##---------------------------------------------------------------------------------------------
## CMSY and BSM analysis ----
## Developed by Rainer Froese, Gianpaolo Coro and Henning Winker in 2016, version of November 2019
## PDF creation added by Gordon Tsui and Gianpaolo Coro
## Default rules for biomass windows and range of q were improved for better batch processing
## Time series within 1950-2020 are stored in csv file 
## Correction for effort creep added by RF
## Multivariate normal r-k priors added to CMSY by HW, RF and GP in October 2019
## Multivariate normal plus observation error on catch added to BSM by HW in November 2019
## Retrospective analysis added by GP in November 2019
##---------------------------------------------------------------------------------------------

# Automatic installation of missing packages
list.of.packages <- c("R2jags","coda","parallel","foreach","doParallel","gplots","mvtnorm",'snpar')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(R2jags)  # Interface with JAGS
library(coda) 
library(parallel)
library(foreach)
library(doParallel)
library(gplots)
library(mvtnorm)
library(snpar)

#-----------------------------------------
# Some general settings ----
#-----------------------------------------
# set.seed(999) # use for comparing results between runs
rm(list=ls(all=FALSE)) # clear previous variables etc
options(digits=3) # displays all numbers with three significant digits as default
graphics.off() # close graphics windows from previous sessions
FullSchaefer <- F    # initialize variable; automatically set to TRUE if enough abundance data are available
n.chains     <- ifelse(detectCores() > 2,3,2) # set 3 chains in JAGS if more than 2 cores are available
ncores_for_computation=detectCores() # cores to be used for parallel processing of CMSY
cl           <- makeCluster(ncores_for_computation)
registerDoParallel(cl, cores = ncores_for_computation)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set working directory to source file location

#-----------------------------------------
# Required settings, File names ----
#-----------------------------------------
catch_file  <- "Stocks_Catch_EU_1.csv" #  name of file containing "Stock", "yr", "ct", and optional "bt"
id_file     <- "Stocks_ID_EU_8.csv" #  name of file containing stock-specific info and settings for the analysis
outfile     <- paste("Out_",format(Sys.Date(),format="%B%d%Y_"),id_file,sep="") # default name for output file

#----------------------------------------
# Select stock to be analyzed ----
#----------------------------------------
stocks      <-NA
# If the input files contain more than one stock, specify below the stock to be analyzed
# If the line below is commented out (#), all stocks in the input file will be analyzed
 stocks <- "fle-2732" # 

#-----------------------------------------
# General settings for the analysis ----
#-----------------------------------------
dataUncert   <- 0.3  # set observation error as uncertainty in catch - default is SD=0.3
sigmaR       <- 0.1 # overall process error for CMSY; SD=0.1 is the default
cor.log.rk   <- -0.607 # empirical value of log r-k correlation in 140 stocks analyzed with BSM (without r-k correlation)
rk.cor.beta  <- c(2.52,3.37) # beta.prior for rk cor+1
nbk          <- 3 # Number of B/k priors to be used by BSM, with options 1 (first year), 2 (first & intermediate), 3 (first, intermediate & final bk priors)  
q.biomass.pr <- c(0.9,1.1) # if btype=="biomass" this is the prior range for q
n            <- 20000 #20000 # initial number of r-k pairs #reduced for testing
# n.new        <- n # initialize n.new
ni           <- 3 # iterations for r-k-startbiomass combinations, to test different variability patterns; no improvement seen above 3
nab          <- 3 # recommended=5; minimum number of years with abundance data to run BSM
bw           <- 3 # default bandwidth to be used by ksmooth() for catch data
mgraphs      <- T # set to TRUE to produce additional graphs for management
e.creep.line <- T # set to TRUE to display uncorrected CPUE in biomass graph
kobe.plot    <- T # set to TRUE to produce additional kobe status plot; management graph needs to be TRUE for Kobe to work 
BSMfits.plot <- T # set to TRUE to plot fit diagnostics for BSM
pp.plot      <- T # set to TRUE to plot Posterior and Prior distributions 
retros       <- F # set to TRUE to enable retrospective analysis (1-3 years less in the time series)
save.plots   <- F # set to TRUE to save graphs to JPEG files
close.plots  <- F # set to TRUE to close on-screen plots after they are saved, to avoid "too many open devices" error in batch-processing
write.output <- T # set to TRUE if table with results in output file is wanted; expects years 2004-2014 to be available
write.pdf    <- F # set to TRUE if PDF output of results is wanted. See more instructions at end of code.
select.yr    <- NA # option to display F, B, F/Fmsy and B/Bmsy for a certain year; default NA



if(write.pdf == FALSE) save.plots=TRUE

#----------------------------------------------
#  FUNCTIONS ----
#----------------------------------------------
# Monte Carlo filtering with Schaefer Function ----
#----------------------------------------------
SchaeferParallelSearch<-function(ni, nyr,sigR,duncert,ct,int.yr,intbio, startbt, ki,i, ri,int.yr.i, nstartbt, yr, end.yr, endbio, npoints, pt){
  ptm<-proc.time()
  # create vectors for viable r, k and bt
  inmemorytable <- vector()
  # parallelised for the points in the r-k space
  inmemorytable <- foreach (i = 1 : npoints, .combine='rbind', .packages='foreach', .inorder=TRUE) %dopar%{
    nsbt = length(startbt)
    VP   <- FALSE
    for(nj in 1:nsbt) {  
      # create empty vector for annual biomasses
      bt <- vector()
      j<-startbt[nj]
      # set initial biomass, including 0.1 process error to stay within bounds
      bt[1]=j*ki[i]*exp(rnorm(1,0, 0.1*sigR))  ## set biomass in first year
      # repeat test of r-k-startbt combination to allow for different random error
      for(re in 1:ni)   {
        #loop through years in catch time series
        for (t in 1:nyr)  {  # for all years in the time series
          xt=rnorm(1,0, sigR) # set new process error for every year  
          zlog.sd = sqrt(log(1+(duncert)^2))
          zt=rlnorm(1,meanlog = 0, sdlog = zlog.sd) # model the catch error as a log normal distribution.
          # calculate biomass as function of previous year's biomass plus surplus production minus catch
          bt[t+1] <- ifelse(bt[t]/ki[i] >= 0.25,
                            bt[t]+ri[i]*bt[t]*(1-bt[t]/ki[i])*exp(xt)-ct[t]*zt,
                            bt[t]+(4*bt[t]/ki[i])*ri[i]*bt[t]*(1-bt[t]/ki[i])*exp(xt)-ct[t]*zt) # assuming reduced r at B/k < 0.25
          
          # if biomass < 0.01 k, discard r-k-startbt combination
          if(bt[t+1] < 0.01*ki[i]) { 
            break
          } # stop looping through years, go to next upper level
          # intermediate year check
          if ((t+1)==int.yr.i && (bt[t+1]>(intbio[2]*ki[i]) || bt[t+1]<(intbio[1]*ki[i]))) { 
            break 
          }  
        } # end of loop of years
        # if loop was broken or last biomass falls outside of expected ranges 
        # do not store results, go directly to next startbt
        if(t < nyr || bt[yr==end.yr] > (endbio[2]*ki[i]) || bt[yr==end.yr] < (endbio[1]*ki[i]) ) { 
          next 
        } else {
          #each vector will be finally appended to the others found by the threads - this is done by the .combine='rbind' option
          inmemorytablerow<-c(i,j,ri[i],ki[i],bt[1:(nyr+1)]/ki[i])
          if (length(inmemorytablerow)==(4+nyr+1)){
            if (VP==FALSE)
            {
              inmemorytable<-inmemorytablerow
            }
            else
            {
              inmemorytable<-rbind(inmemorytable,inmemorytablerow)
            }
            VP<-TRUE
          }
        }
      } # end of repetition for random error
    } # end of j-loop of initial biomasses 
    # instruction necessary to make the foreach loop see the variable:
    if (length(inmemorytable)==0)
    {inmemorytable<-vector(length=4+nyr+1)*NA}
    else
    {inmemorytable}
  }#end loop on points
  
  #create the output matrix
  mdat        <- matrix(data=NA, nrow = npoints*nstartbt, ncol = 2+nyr+1) 
  npointsmem = dim(inmemorytable)[1]
  npointscols = dim(inmemorytable)[2]
  #reconstruction of the processing matrix after the parallel search
  if (npointsmem>0 && npointscols>0){
    for (idxr in 1:npointsmem){
      i = inmemorytable[idxr,1]
      if (!is.na(i)){
        j = inmemorytable[idxr,2]
        mdatindex<-((i-1)*nstartbt)+which(startbt==j)  
        mdat[mdatindex,1]           <- inmemorytable[idxr,3]
        mdat[mdatindex,2]           <- inmemorytable[idxr,4]
        mdat[mdatindex,3:(2+nyr+1)] <- inmemorytable[idxr,5:(4+nyr+1)]
        if(pt==T) points(x=ri[i], y=ki[i], pch=".", cex=4, col="gray")
      }
    }
  }
  ptm<-proc.time()-ptm
  mdat <- na.omit(mdat)
  return(mdat)
}

SchaeferMC <- function(ri, ki, startbio, int.yr, intbio, endbio, sigR, pt, duncert, startbins, ni) {
  
  # create vector for initial biomasses
  startbt     <- seq(from =startbio[1], to=startbio[2], by = (startbio[2]-startbio[1])/startbins)
  nstartbt    <- length(startbt)
  npoints     <- length(ri)
  # get index of intermediate year
  int.yr.i     <- which(yr==int.yr) 
  
  #loop through r-k pairs with parallel search
  mdat<-SchaeferParallelSearch(ni, nyr,sigR,duncert,ct,int.yr,intbio, startbt, ki, i, ri, int.yr.i, nstartbt, yr, end.yr, endbio, npoints,pt)
  
  cat("\n")
  return(list(mdat))
} # end of SchaeferMC function

#-------------------------------------------------------------
# Function to create multivariate-normal distribution for r-k 
#-------------------------------------------------------------
mvn   <- function(n,mean.log.r,sd.log.r,mean.log.k,sd.log.k) {
  cov.log.rk <- cor.log.rk*sd.log.r*sd.log.k # covariance with empirical correlation and prior variances  covar.log.rk = matrix(NA, ncol=2,nrow=2)   # contract covariance matrix
  covar.log.rk      <- matrix(NA, ncol=2,nrow=2) # covariance matrix
  covar.log.rk[1,1] <- sd.log.r^2                # position [1,1] is variance of log.r
  covar.log.rk[2,2] <- sd.log.k^2               # position [2,2] is variance of log.k
  covar.log.rk[1,2] = covar.log.rk[2,1] = cov.log.rk     # positions [1,2] and [2,1] are correlations
  mu.log.rk  <- (c(mean.log.r,mean.log.k))      # vector of log.means
  mvn.log.rk <- rmvnorm(n,mean=mu.log.rk,sigma=covar.log.rk,method="svd") 
  return(mvn.log.rk)
}

#---------------------------------------------
# END OF FUNCTIONS
#---------------------------------------------

#-----------------------------------------
# Start output to screen
#-----------------------------------------
cat("-------------------------------------------\n")
cat("CMSY Analysis,", date(),"\n")
cat("-------------------------------------------\n")

#------------------------------------------
# Read data and assign to vectors
#------------------------------------------
# create headers for data table file
if(write.output==T){
  outheaders = data.frame("Group","Region", "Subregion","Name","SciName","Stock","start.yr","end.yr","btype",
                          "MaxCatch","LastCatch","MSY_BSM","lcl","ucl","r_BSM","lcl","ucl","log_var",
                          "k_BSM","lcl","ucl","log_var","log_rk_cor","log_rk_cov","q_BSM","lcl","ucl","rel_B_BSM","lcl","ucl","rel_F_BSM",
                          "r_CMSY","lcl","ucl","k_CMSY","lcl","ucl","MSY_CMSY","lcl","ucl",
                          "rel_B_CMSY","2.5th","97.5th","rel_F_CMSY",
                          "F_msy","lcl","ucl","curF_msy","lcl","ucl",
                          "MSY","lcl","ucl","Bmsy","lcl","ucl",
                          "B","lcl","ucl","B_Bmsy","lcl","ucl",
                          "F","lcl","ucl","F_Fmsy","lcl","ucl",
                          "sel_B","sel_B_Bmsy","sel_F","sel_F_Fmsy",
                          # create columns for catch, F/Fmsy and Biomass for 1970 to 2020
                          "c50","c51","c52","c53","c54","c55","c56","c57","c58","c59",
                          "c60","c61","c62","c63","c64","c65","c66","c67","c68","c69",
                          "c70","c71","c72","c73","c74","c75","c76","c77","c78","c79",
                          "c80","c81","c82","c83","c84","c85","c86","c87","c88","c89",
                          "c90","c91","c92","c93","c94","c95","c96","c97","c98","c99",
                          "c00","c01","c02","c03","c04","c05","c06","c07","c08","c09",
                          "c10","c11","c12","c13","c14","c15","c16","c17","c18","c19","c20",
                          "F.Fmsy50","F.Fmsy51","F.Fmsy52","F.Fmsy53","F.Fmsy54","F.Fmsy55","F.Fmsy56","F.Fmsy57","F.Fmsy58","F.Fmsy59",
                          "F.Fmsy60","F.Fmsy61","F.Fmsy62","F.Fmsy63","F.Fmsy64","F.Fmsy65","F.Fmsy66","F.Fmsy67","F.Fmsy68","F.Fmsy69",
                          "F.Fmsy70","F.Fmsy71","F.Fmsy72","F.Fmsy73","F.Fmsy74","F.Fmsy75","F.Fmsy76","F.Fmsy77","F.Fmsy78","F.Fmsy79",
                          "F.Fmsy80","F.Fmsy81","F.Fmsy82","F.Fmsy83","F.Fmsy84","F.Fmsy85","F.Fmsy86","F.Fmsy87","F.Fmsy88","F.Fmsy89",
                          "F.Fmsy90","F.Fmsy91","F.Fmsy92","F.Fmsy93","F.Fmsy94","F.Fmsy95","F.Fmsy96","F.Fmsy97","F.Fmsy98","F.Fmsy99",
                          "F.Fmsy00","F.Fmsy01","F.Fmsy02","F.Fmsy03","F.Fmsy04","F.Fmsy05","F.Fmsy06","F.Fmsy07","F.Fmsy08","F.Fmsy09",
                          "F.Fmsy10","F.Fmsy11","F.Fmsy12","F.Fmsy13","F.Fmsy14","F.Fmsy15","F.Fmsy16","F.Fmsy17","F.Fmsy18","F.Fmsy19","F.Fmsy20",
                          "B50","B51","B52","B53","B54","B55","B56","B57","B58","B59",
                          "B60","B61","B62","B63","B64","B65","B66","B67","B68","B69",
                          "B70","B71","B72","B73","B74","B75","B76","B77","B78","B79",
                          "B80","B81","B82","B83","B84","B85","B86","B87","B88","B89",
                          "B90","B91","B92","B93","B94","B95","B96","B97","B98","B99",
                          "B00","B01","B02","B03","B04","B05","B06","B07","B08","B09",
                          "B10","B11","B12","B13","B14","B15","B16","B17","B18","B19","B20")
  
  write.table(outheaders,file=outfile, append = T, sep=",",row.names=F,col.names=F)
}

cat("Parallel processing will use",ncores_for_computation,"cores\n")
# Read data
cdat         <- read.csv(catch_file, header=T, dec=".", stringsAsFactors = FALSE)
cinfo        <- read.csv(id_file, header=T, dec=".", stringsAsFactors = FALSE)
cat("Files", catch_file, ",", id_file, "read successfully","\n")

#---------------------------------
# Analyze stock(s)
#---------------------------------
if(is.na(stocks[1])==TRUE){
  # stocks         <- as.character(cinfo$Stock) # Analyze stocks in sequence of ID file
  # stocks         <- sort(as.character(cinfo$Stock[cinfo$Stock>="ghl-arct"])) # Analyze in alphabetic order after a certain stock
   stocks         <- sort(as.character(cinfo$Stock)) # Analyze stocks in alphabetic order
  # stocks         <- as.character(cinfo$Stock[cinfo$Subregion=="Sardinia"]) # Analyze stocks in Region
}

# analyze one stock after the other
for(stock in stocks) {

	cat("Processing",stock,",", as.character(cinfo$ScientificName[cinfo$Stock==stock]),"\n")

	#retrospective analysis
	retros.nyears<-ifelse(retros==T,3,0) #retrospective analysis
	FFmsy.retrospective<-list() #retrospective analysis
	BBmsy.retrospective<-list() #retrospective analysis
	years.retrospective<-list() #retrospective analysis


  for (retrosp.step in 0:retros.nyears){ #retrospective analysis loop
		
  
  # assign data from cinfo to vectors
  res          <- as.character(cinfo$Resilience[cinfo$Stock==stock])
  start.yr     <- as.numeric(cinfo$StartYear[cinfo$Stock==stock])
  end.yr       <- as.numeric(cinfo$EndYear[cinfo$Stock==stock])
  end.yr.orig  <- end.yr
  end.yr 	     <- end.yr-retrosp.step #retrospective analysis
  r.low        <- as.numeric(cinfo$r.low[cinfo$Stock==stock])
  r.hi         <- as.numeric(cinfo$r.hi[cinfo$Stock==stock])
  user.log.r   <- ifelse(is.na(r.low)==F & is.na(r.hi)==F,TRUE,FALSE)     
  stb.low      <- as.numeric(cinfo$stb.low[cinfo$Stock==stock])
  stb.hi       <- as.numeric(cinfo$stb.hi[cinfo$Stock==stock])
  int.yr       <- as.numeric(cinfo$int.yr[cinfo$Stock==stock])
  intb.low     <- as.numeric(cinfo$intb.low[cinfo$Stock==stock])
  intb.hi      <- as.numeric(cinfo$intb.hi[cinfo$Stock==stock])
  endb.low     <- as.numeric(cinfo$endb.low[cinfo$Stock==stock])
  endb.hi      <- as.numeric(cinfo$endb.hi[cinfo$Stock==stock])
  q.start      <- cinfo$q.start[cinfo$Stock==stock]
  q.end        <- cinfo$q.end[cinfo$Stock==stock]
  e.creep      <- as.numeric(cinfo$e.creep[cinfo$Stock==stock])
  btype        <- as.character(cinfo$btype[cinfo$Stock==stock])
  force.cmsy   <- cinfo$force.cmsy[cinfo$Stock==stock]
  comment      <- as.character(cinfo$Comment[cinfo$Stock==stock])
  source       <- as.character(cinfo$Source[cinfo$Stock==stock])
  # set global defaults for uncertainty
  duncert      <- dataUncert
  sigR         <- sigmaR
  
	if (retros==T && retrosp.step==0){
		cat("* ",ifelse(btype!="None","BSM","CMSY")," retrospective analysis for ",stock," has been enabled\n",sep="") #retrospective analysis
	}

	
  if (retros==T){
		cat("* Retrospective analysis: step n. ",(retrosp.step+1),"/",(retros.nyears+1),". Range of years: [",start.yr ," - ",end.yr,"]\n",sep="") #retrospective analysis
	}
  # check for common errors
  if(length(btype)==0){
    cat("ERROR: Could not find the stock in the ID input file - check that the stock names match in ID and Catch files and that commas are used (not semi-colon)")
    return (NA) }
  if(length(cdat$yr[cdat$Stock==stock])==0){
    cat("ERROR: Could not find the stock in the Catch file - check that the stock names match in ID and Catch files and that commas are used (not semi-colon)")
    return (NA) }
  if(start.yr < cdat$yr[cdat$Stock==stock][1]){
    cat("ERROR: start year in ID file before first year in catch file\n")
    return (NA)}
  
  # extract data on stock
  yr           <- as.numeric(cdat$yr[cdat$Stock==stock & cdat$yr >= start.yr & cdat$yr <= end.yr])
  
  if(length(yr)==0){
    cat("ERROR: Could not find the stock in the Catch input files - Please check that the code is written correctly")
    return (NA) }
  if(btype %in% c("None","CPUE","biomass")==FALSE){
    cat("ERROR: In ID file, btype must be None, CPUE, or biomass.")
    return (NA) }
  if(length(yr) != (end.yr-start.yr+1)) {
    cat("ERROR: indicated year range is of different length than years in catch file\n")
    return (NA)}
  
  ct.raw   <- as.numeric(cdat$ct[cdat$Stock==stock & cdat$yr >= start.yr & cdat$yr <= end.yr])/1000  ## assumes that catch is given in tonnes, transforms to '000 tonnes
  if(btype=="biomass" | btype=="CPUE" ) {
    bt.raw <- as.numeric(cdat$bt[cdat$Stock==stock & cdat$yr >= start.yr & cdat$yr <= end.yr])/1000  ## assumes that biomass is in tonnes, transforms to '000 tonnes
    bt     <- bt.raw #ksmooth(x=yr,y=bt.raw,kernel="normal",n.points=length(yr),bandwidth=3)$y
    if(length(bt[is.na(bt)==F])==0) {
      cat("ERROR: No CPUE or biomass data in the Catch input file")
      return (NA) }
  } else {bt <- NA}
 
  # apply correction for effort-creep to commercial(!) CPUE 
  if(btype=="CPUE" && is.na(e.creep)==FALSE) {
    cpue.first  <- min(which(is.na(bt)==F)) # ifelse(is.na(q.start)==T,min(which(is.na(bt)==F)),which(yr==q.start))
    cpue.last   <- max(which(is.na(bt)==F)) # ifelse(is.na(q.end)==T,max(which(is.na(bt)==F)),which(yr==q.end))
#	if (is.na(cpue.last) && retros == T){ # modification for retrospective analysis
#		cat("Error: q.end should be set to maximum three years before last year when running retrospective analysis")
#	}
		
    cpue.length <- cpue.last - cpue.first
    bt.cor      <- bt
    for(i in 1:(cpue.length)) {
      bt.cor[cpue.first+i]  <- bt[cpue.first+i]*(1-e.creep/100)^i # equation for decay in %
    }
    bt <- bt.cor
  }

  # set q.start and q.end both to NA if q.end falls within the last three years
  if(is.na(q.end)==F && force.cmsy==F && btype!="None" && retros==T ) { #adjustment for retrospective analysis
	if ( (end.yr.orig-q.end) <= 3){
			q.start<-NA
			q.end<-NA
			if (retrosp.step == 0){
				cat("Warning: User range for q.start-q.end overlaps with last three years, changed to NA for retrospective analysis\n")
			}
	}
  }
  
  if(retros==T && force.cmsy == F && (btype !="None" & length(bt[is.na(bt)==F])<nab) ) { #stop retrospective analysis if cpue is < nab
	cat("Warning: Cannot run retrospective analysis for ",end.yr,", number of remaining ",btype," values is too low (<",nab,")\n",sep="")
	#retrosp.step<-retros.nyears
	break
  }
  
  if(is.na(mean(ct.raw))){
    cat("ERROR: Missing value in Catch data; fill or interpolate\n")  
  }
  
  
  
  
  nyr          <- length(yr) # number of years in the time series
  
  # apply kernel smoothing with a bandwidth of bw
  ct           <- ksmooth(x=yr,y=ct.raw,kernel="normal",n.points=length(yr),bandwidth=bw)$y
  
  # initialize vectors for viable r, k, bt, and all in a matrix
  mdat.all    <- matrix(data=vector(),ncol=2+nyr+1)
  
  # initialize other vectors anew for each stock
  current.attempts <- NA
  
  # use start.yr if larger than select year
  if(is.na(select.yr)==F) {
    sel.yr <- ifelse(start.yr > select.yr,start.yr,select.yr)
  } else sel.yr <- NA
  
  #----------------------------------------------------
  # Determine initial ranges for parameters and biomass
  #----------------------------------------------------
  # initial range of r from input file
  if(is.na(r.low)==F & is.na(r.hi)==F) {
    prior.r <- c(r.low,r.hi)
  } else {
    # initial range of r based on resilience
    if(res == "High") {
      prior.r <- c(0.6,1.5)} else if(res == "Medium") {
        prior.r <- c(0.2,0.8)}    else if(res == "Low") {
          prior.r <- c(0.05,0.5)}  else { # i.e. res== "Very low"
            prior.r <- c(0.015,0.1)} 
  }
  
  # get index of years with lowest and highest catch between start+3 and end-3 years
  min.yr.i     <- which.min(ct[4:(length(ct)-3)])+3
  max.yr.i     <- which.max(ct[4:(length(ct)-3)])+3
  min.ct       <- ct[min.yr.i]
  max.ct       <- ct[max.yr.i]
  
  # use initial biomass range from input file if stated
  if(is.na(stb.low)==F & is.na(stb.hi)==F) {
    startbio <- c(stb.low,stb.hi)
  } else {
    # if catch < 0.1 max catch, assume nearly unexploited biomass
    if(ct[1] < 0.1*max.ct) { startbio <- c(0.9,1)
    # if catch < 0.25 max catch, assume high biomass
    } else if(ct[1] < 0.25*max.ct) { startbio <- c(0.8,1)
    # if catch < 0.33 max catch, assume high biomass
    } else if(ct[1] < 0.33*max.ct) { startbio <- c(0.6,1)
    # if catch < 0.66 max catch, assume medium to high biomass
    } else if(ct[1] < 0.66*max.ct | start.yr <=1960) { startbio <- c(0.4,0.8)
    # otherwise assume low to medium biomass
    } else startbio <- c(0.2,0.6) 
  }
  
  # use year and biomass range for intermediate biomass from input file
  if(is.na(intb.low)==F & is.na(intb.hi)==F) {
    int.yr   <- int.yr
    intbio   <- c(intb.low,intb.hi)
    
    # if contrast in catch is low, use initial range again in mid-year
  } else if(min(ct)/max(ct) > 0.6) {
    int.yr    <- as.integer(mean(c(start.yr, end.yr)))
    intbio    <- startbio 
    
    # else if year of minimum catch is after max catch then use min catch
  } else if(min.yr.i > max.yr.i) {
    int.yr    <- yr[min.yr.i-1]
    if(startbio[1]>=0.5 &  (int.yr-start.yr) < (end.yr-int.yr) & 
       (min.ct/max.ct) > 0.3) intbio <- c(0.2,0.6) else intbio <- c(0.01,0.4)
       
       # else use max catch  
  } else {
    # assume that biomass range in year before maximum catch was high or medium
    int.yr    <- yr[max.yr.i-1]
    intbio    <- if((startbio[1]>=0.5 & (int.yr-start.yr) < (end.yr-int.yr))| # if initial biomass is high, assume same for intermediate
                    (((max.ct-min.ct)/max.ct)/(max.yr.i-min.yr.i) > 0.04)) c(0.5,0.9) else 
                      if(ct[which(yr==int.yr)-5]/max.ct<0.33) {c(0.4,0.8)} else c(0.2,0.6) } # if incease is steep, assume high, else medium
  # end of intbio setting
  
  # final biomass range from input file
  if(is.na(endb.low)==F & is.na(endb.hi)==F) {
    endbio   <- c(endb.low,endb.hi)
  } else {
    rawct.ratio=ct.raw[nyr]/max(ct)
    # else use mean final catch/max catch to estimate final biomass
    endbio  <- if(ct[nyr]/max(ct) > 0.8) {c(0.4,0.8)} else if(rawct.ratio < 0.5) {c(0.01,0.4)} else {c(0.2,0.6)}
    
    # if final catch is <=  catch at int.yr, endbio may not exceed intbio
    if(endbio[1] > intbio[1] & ct[nyr] <= 1.1*ct[which(yr==int.yr)]) {endbio <- intbio}
    
    # if endbio is less than 5 years after intbio endbio may not exceed intbio 
    if(endbio[1] > intbio[1]) {endbio <- intbio}
    
    # if default endbio is low (0.01-0.4), check whether the upper bound should be lower than 0.4 for depleted stocks
    if(endbio[2]==0.4){
      if(rawct.ratio< 0.05) {endbio[2] <- 0.1} else
        if(rawct.ratio< 0.15) {endbio[2] <- 0.2} else
          if(rawct.ratio< 0.35) {endbio[2] <- 0.3} else {endbio[2] <- 0.4}
    }
  } # end of final biomass setting

  if(mean(endbio) <= 0.5) {
   prior.k <- c(2*max(ct)/mean(prior.r),6*max(ct)/mean(prior.r))} else {
    prior.k <- c(2*max(ct)/mean(prior.r),10*max(ct)/mean(prior.r))}

  cat("startbio=",startbio,ifelse(is.na(stb.low)==T,"default","expert"),
      ", intbio=",int.yr,intbio,ifelse(is.na(intb.low)==T,"default","expert"),
      ", endbio=",endbio,ifelse(is.na(endb.low)==T,"default","expert"),"\n")
  
  #----------------------------------------------------------------
  # Multivariate normal sampling of r-k log space
  #----------------------------------------------------------------
  # turn numerical ranges into log-normal distributions 
  
  mean.log.r=mean(log(prior.r))
  sd.log.r=(log(prior.r[2])-log(prior.r[1]))/4  # assume range covers 4 SD
  
  mean.log.k <- mean(log(prior.k))
  sd.log.k   <- (log(prior.k[2])-log(prior.k[1]))/4 # assume range covers 4 SD
  
  mvn.log.rk <- mvn(n=n,mean.log.r=mean.log.r,sd.log.r=sd.log.r,mean.log.k=mean.log.k,sd.log.k=sd.log.k)
  ri1    <- exp(mvn.log.rk[,1])
  ki1    <- exp(mvn.log.rk[,2]) 
 
  #-----------------------------------------------------------------
  #Plot data and progress -----
  #-----------------------------------------------------------------
  # check for operating system, open separate window for graphs if Windows
  if(grepl("win",tolower(Sys.info()['sysname']))) {windows(14,9)}
  par(mfrow=c(2,3),mar=c(5.1,4.5,4.1,2.1))
  # plot catch ----
  plot(x=yr, y=ct.raw, 
       ylim=c(0,max(ifelse(substr(id_file,1,3)=="Sim",
                           1.1*true.MSY,0),1.2*max(ct.raw))),
       type ="l", bty="l", main=paste("A: Catch",stock), xlab="", ylab="Catch (1000 tonnes/year)", lwd=2, cex.main = 1.8, cex.lab = 1.55, cex.axis = 1.5)
  lines(x=yr,y=ct,col="blue", lwd=1)
  points(x=yr[max.yr.i], y=max.ct, col="red", lwd=2)
  points(x=yr[min.yr.i], y=min.ct, col="red", lwd=2)

  # (b): plot r-k graph 
  plot(x=ri1, y=ki1, xlim = c(0.95*quantile(ri1,0.001),1.2*quantile(ri1,0.999)), 
       ylim = c(0.95*quantile(ki1,0.001),1.2*quantile(ki1,0.999)),
       log="xy", xlab="r", ylab="k (1000 tonnes)", main="B: Finding viable r-k", pch=".", cex=3, bty="l", 
       col="grey95", cex.main = 1.8, cex.lab = 1.55, cex.axis = 1.5)
  lines(x=c(prior.r[1],prior.r[2],prior.r[2],prior.r[1],prior.r[1]), # plot original prior range
        y=c(prior.k[1],prior.k[1],prior.k[2],prior.k[2],prior.k[1]),
        lty="dotted") 
  
  #---------------------------------------------------------------------
  # 1 - Call CMSY-SchaeferMC function to preliminary explore the r-k space ----
  #---------------------------------------------------------------------
  cat("First Monte Carlo filtering of r-k space with ",n," points...\n")
  MCA <-  SchaeferMC(ri=ri1, ki=ki1, startbio=startbio, int.yr=int.yr, intbio=intbio, endbio=endbio, sigR=sigR, 
                     pt=T, duncert=dataUncert, startbins=10, ni=ni)
  mdat.all <- rbind(mdat.all,MCA[[1]])
  rv.all   <- mdat.all[,1]
  kv.all   <- mdat.all[,2]
  btv.all  <- mdat.all[,3:(2+nyr+1)]
  # count viable trajectories and r-k pairs ----
  n.viable.b   <- length(mdat.all[,1])
  n.viable.pt <- length(unique(mdat.all[,1]))
  cat("Found ",n.viable.b," viable trajectories for", n.viable.pt," r-k pairs\n")
  
  #----------------------------------------------------------------------- 
  # 2 - if the lower bound of k is too high, reduce it by half and rerun ----
  #-----------------------------------------------------------------------
  # if overall points are few and mostly found in the lower-left prior space, then reduce lower bound of k
  if(length(kv.all < 200 && kv.all[kv.all < 1.1*prior.k[1] & rv.all < mean(prior.r)]) > 20) { 
    cat("Reducing lower bound of k, resampling area with",n,"additional points...\n")
    prior.k    <- c(0.5*prior.k[1],prior.k[2])
    mean.log.k <- mean(log(prior.k))
    sd.log.k   <- (log(prior.k[2])-log(prior.k[1]))/4 # assume range covers 4 SD
    
    mvn.log.rk <- mvn(n=n,mean.log.r=mean.log.r,sd.log.r=sd.log.r,mean.log.k=mean.log.k,sd.log.k=sd.log.k)
    ri1        <- exp(mvn.log.rk[,1])
    ki1        <- exp(mvn.log.rk[,2])
   
    MCA <-  SchaeferMC(ri=ri1, ki=ki1, startbio=startbio, int.yr=int.yr, intbio=intbio, endbio=endbio, sigR=sigR, 
                       pt=T, duncert=dataUncert, startbins=10, ni=ni)
    mdat.all <- rbind(mdat.all,MCA[[1]])
    rv.all   <- mdat.all[,1]
    kv.all   <- mdat.all[,2]
    btv.all  <- mdat.all[,3:(2+nyr+1)]
    n.viable.b   <- length(mdat.all[,1])
    n.viable.pt <- length(unique(mdat.all[,1]))
    cat("Found altogether",n.viable.b," viable trajectories for", n.viable.pt," r-k pairs\n")
  }
  
  #-------------------------------------------------------------------
  # 3 - if few points were found then resample with more points
  #-------------------------------------------------------------------
  if (n.viable.b <= 400){
    log.prior.k.new  <- log(prior.k) 
    max.attempts     <- 3
    current.attempts <- 1
    startbins        <- 10  
    while (n.viable.b <= 400 && current.attempts <= max.attempts){
      n.new      <- n*current.attempts #add more points
      mvn.log.rk <- mvn(n=n.new,mean.log.r=mean.log.r,sd.log.r=sd.log.r,
                        mean.log.k=mean.log.k,sd.log.k=sd.log.k)
      ri1        <- exp(mvn.log.rk[,1])
      ki1        <- exp(mvn.log.rk[,2])
      
      cat("Repeating analysis with more points...\n")
      cat("Attempt ",current.attempts," of ",max.attempts," with ",n.new," additional points...","\n")
      if(current.attempts==2 & n.viable.b < 50){
        duncert   <- 2*dataUncert
        sigR      <- 2*sigmaR
        startbins <- 20
        bw        <- 4 # increase bandwidth of smoothing
        ct        <- ksmooth(x=yr,y=ct.raw,kernel="normal",n.points=length(yr),bandwidth=bw)$y
        cat("Increasing startbins, smoothing, catch and process error, and number of variability patterns \n")   
      }
      MCA <-  SchaeferMC(ri=ri1, ki=ki1, startbio=startbio, int.yr=int.yr, intbio=intbio, endbio=endbio, sigR=sigR, 
                         pt=T, duncert=duncert, startbins=startbins, ni=2*ni)
      mdat.all <- rbind(mdat.all,MCA[[1]])
      rv.all   <- mdat.all[,1]
      kv.all   <- mdat.all[,2]
      btv.all  <- mdat.all[,3:(2+nyr+1)]
      n.viable.b   <- length(mdat.all[,1])
      n.viable.pt <- length(unique(mdat.all[,1]))
      cat("Found altogether",n.viable.b," viable trajectories for", n.viable.pt," r-k pairs\n")
      current.attempts=current.attempts+1 #increment the number of attempts
    } # end of 3 attempts loop 
  } # end of loop 3
  
  # --------------------------------------------------------------------------------
  # 4 - if too few points are found, remove intermediate filter by setting it to 0-1 ----
  # --------------------------------------------------------------------------------
  if(n.viable.b < 5) {
    cat("Setting intermediate biomass to 0-1... \n")
    int.yr=yr[as.integer(nyr/2)] 
    intbio=c(0,1)
    MCA <-  SchaeferMC(ri=ri1, ki=ki1, startbio=startbio, int.yr=int.yr, intbio=intbio, endbio=endbio, sigR=sigR, 
                       pt=T, duncert=dataUncert, startbins=10, ni=ni)
    mdat.all <- rbind(mdat.all,MCA[[1]])
    rv.all   <- mdat.all[,1]
    kv.all   <- mdat.all[,2]
    btv.all  <- mdat.all[,3:(2+nyr+1)]
    # count viable trajectories and r-k pairs ----
    n.viable.b   <- length(mdat.all[,1])
    n.viable.pt <- length(unique(mdat.all[,1]))
    cat("Found ",n.viable.b," viable trajectories for", n.viable.pt," r-k pairs\n")
  }
  if(n.viable.b < 5) {
    cat("Only",n.viable.pt,"viable r-k pairs found, check data and settings \n")
    if(retros==TRUE) {break} else {next}
  }  

  #-----------------------------------------------
  # For CMSY, get best r and k as 75th percentiles
  #-----------------------------------------------  
  unique.rk         <- unique(mdat.all[,1:2]) # get unique r-k pairs
  rs                <- unique.rk[,1]
  ks                <- unique.rk[,2]
  r.est             <- quantile(rs,0.75)
  ucl.r.est         <- quantile(rs,0.9875)
  # use symmetrical confidence limits in log space
  ucl.dist.log.r    <- log(ucl.r.est)-log(r.est) 
  lcl.r.est         <- exp(log(r.est)-ucl.dist.log.r)
  k.est             <- quantile(ks,0.25)
  lcl.k.est         <- quantile(ks,0.0125)
  lcl.dist.log.k    <- log(k.est)-log(lcl.k.est) 
  ucl.k.est         <- exp(log(k.est)+lcl.dist.log.k)
  # get MSY from r-k pairs within the approximate confidence limits
  MSYs              <- rs[rs>=lcl.r.est & rs<=ucl.r.est & ks>=lcl.k.est & ks<=ucl.k.est]*
    ks[rs>=lcl.r.est & rs<=ucl.r.est & ks>=lcl.k.est & ks<=ucl.k.est]/4
  MSY.est           <- median(MSYs)
  lcl.MSY.est       <- quantile(MSYs,0.025)
  ucl.MSY.est       <- quantile(MSYs,0.975)
  
    #-----------------------------------------  
  # get predicted biomass vectors as median and quantiles 
  # only use biomass trajectories from r-k pairs within the confidence limits
  
  #><> HW select rk region for BRP calculations 
  rem = which(rs>=lcl.r.est & rs<=ucl.r.est & ks>=lcl.k.est & ks<=ucl.k.est)
  # Get B/Bmsy CMSY posterior
  rem.btv.all    <- mdat.all[rem,3:(2+nyr+1)]
  #><> HW get select r-k CMSY posteriors (from previous version)
  rem.log.r      <- log(unique.rk[,1][rem])
  rem.log.k      <- log(unique.rk[,2][rem])
  
  median.btv <- apply(rem.btv.all,2, median)
  median.btv.lastyr  <- median.btv[length(median.btv)-1]
  nextyr.bt  <- median.btv[length(median.btv)]
  lcl.btv    <- apply(rem.btv.all,2, quantile, probs=0.025)
  q.btv      <- apply(rem.btv.all,2, quantile, probs=0.25)
  ucl.btv    <- apply(rem.btv.all,2, quantile, probs=0.975)
  lcl.median.btv.lastyr <- lcl.btv[length(lcl.btv)-1]
  ucl.median.btv.lastyr <- ucl.btv[length(lcl.btv)-1]  
  lcl.nextyr.bt <- lcl.btv[length(lcl.btv)]
  ucl.nextyr.bt <- ucl.btv[length(lcl.btv)]
  
  # get F derived from predicted CMSY biomass
  F.CMSY      <- ct.raw/(median.btv[1:nyr]*k.est)
  lcl.F.CMSY  <- ct.raw/(ucl.btv[1:nyr]*k.est)
  ucl.F.CMSY  <- ct.raw/(lcl.btv[1:nyr]*k.est)
  Fmsy.CMSY   <- r.est/2;lcl.Fmsy.CMSY<-lcl.r.est/2;ucl.Fmsy.CMSY<-ucl.r.est/2 # Fmsy from CMSY   
  F_Fmsy.CMSY <- vector();lcl.F_Fmsy.CMSY<-vector();ucl.F_Fmsy.CMSY<-vector() # initialize vectors
  for(z in 1: length(F.CMSY)) {
    F_Fmsy.CMSY[z]     <- F.CMSY[z]/ifelse(median.btv[z]<0.25,Fmsy.CMSY*4*median.btv[z],Fmsy.CMSY)
    lcl.F_Fmsy.CMSY[z] <- lcl.F.CMSY[z]/ifelse(median.btv[z]<0.25,Fmsy.CMSY*4*median.btv[z],Fmsy.CMSY)
    ucl.F_Fmsy.CMSY[z] <- ucl.F.CMSY[z]/ifelse(median.btv[z]<0.25,Fmsy.CMSY*4*median.btv[z],Fmsy.CMSY)
  }
  
  # show CMSY estimate in prior space of graph B
  points(x=r.est, y=k.est, pch=19, col="blue")
  lines(x=c(lcl.r.est, ucl.r.est),y=c(k.est,k.est), col="blue")
  lines(x=c(r.est,r.est),y=c(lcl.k.est, ucl.k.est), col="blue")
  lines(x=c(prior.r[1],prior.r[2],prior.r[2],prior.r[1],prior.r[1]), # plot original prior range
        y=c(prior.k[1],prior.k[1],prior.k[2],prior.k[2],prior.k[1]),lty="dotted")
  

  # ------------------------------------------------------------------
  # Bayesian analysis of catch & biomass (or CPUE) with Schaefer model ----
  # ------------------------------------------------------------------
  FullSchaefer <- F
  if(btype != "None" & length(bt[is.na(bt)==F])>=nab) {
  
    FullSchaefer <- T
        
    # set inits for r-k in lower right corner of log r-k space
    init.r      <- prior.r[1]+0.8*(prior.r[2]-prior.r[1])
    init.k      <- prior.k[1]+0.1*(prior.k[2]-prior.k[1])
    
    # vector with no penalty (=0) if predicted biomass is within viable range, else a penalty of 10 is set 
    pen.bk = pen.F = rep(0,length(ct))
    
    # Add biomass priors
    b.yrs = c(1,length(start.yr:int.yr),length(start.yr:end.yr))
    b.prior = rbind(matrix(c(startbio[1],startbio[2],intbio[1],intbio[2],endbio[1],endbio[2]),2,3),rep(0,3)) # last row includes the 0 pen
    
    #><> Add Catch CV
    CV.C = dataUncert/2 
    #><> Add minimum realistic cpue CV 
    CV.cpue = dataUncert/2
    
    cat("Running MCMC analysis....\n")
    
      # ---------------------------------------------------------------------
      # Schaefer model for Catch & CPUE or biomass
      # ---------------------------------------------------------------------

      # get prior for q from stable catch/biomass period, min 5 years; get range of years from input file
      if(is.na(q.start)==F & is.na(q.end)==F) {
        mean.last.ct      <-mean(ct[yr >= q.start & yr <= q.end], na.rm=T) # get mean catch of indicated years
        mean.last.cpue    <-mean(bt[yr >= q.start & yr <= q.end], na.rm=T) # get mean of CPUE of indicated years
      } else { 
        # get prior range for q from mean catch and mean CPUE in recent years 
        lyr               <- ifelse(mean(prior.r)>=0.5,5,10)  # determine number of last years to use, 5 for normal and 10 for slow growing fish    
        # determine last year with CPUE data
        lbt <- max(which(bt>0))
        mean.last.ct      <-mean(ct[(lbt-lyr):lbt],na.rm=T) # get mean catch of last years
        mean.last.cpue    <-mean(bt[(lbt-lyr):lbt],na.rm=T) # get mean of CPUE of last years
        
      }
      gm.prior.r      <- exp(mean(log(prior.r))) # get geometric mean of prior r range
      if(btype=="biomass") {
        q.prior <- q.biomass.pr
        init.q  <- mean(q.prior)
      } else {
        if(mean(endbio) >= 0.5) {  # if biomass is high  
          q.1           <- mean.last.cpue*0.25*gm.prior.r/mean.last.ct
          q.2           <- mean.last.cpue*0.5*prior.r[2]/mean.last.ct
        } else {
          q.1           <- mean.last.cpue*0.5*gm.prior.r/mean.last.ct
          q.2           <- mean.last.cpue*prior.r[2]/mean.last.ct
        }
        
        q.prior         <- c(q.1,q.2)
        init.q          <- mean(q.prior)
      }
      # Data to be passed on to JAGS
      jags.data        <- c('ct','bt','nyr', 'prior.r', 'prior.k', 'startbio', 'q.prior',
                            'init.q','init.r','init.k','pen.bk','pen.F','b.yrs','b.prior','CV.C','CV.cpue','nbk','rk.cor.beta')
      # Parameters to be returned by JAGS
      jags.save.params <- c('r','k','q', 'P','ct.jags','cpuem','proc.logB') 
      
      # JAGS model ----
      Model = "model{
    # to reduce chance of non-convergence, Pmean[t] values are forced >= eps
    eps<-0.01
    #><> Add Catch.CV
    for(t in 1:nyr){    
      ct.jags[t] ~ dlnorm(log(ct[t]),pow(CV.C,-2))
    }
      

    penm[1] <- 0 # no penalty for first biomass
    Pmean[1] <- log(alpha)
    P[1] ~ dlnorm(Pmean[1],itau2)

    for (t in 2:nyr) {
      Pmean[t] <- ifelse(P[t-1] > 0.25,
      log(max(P[t-1] + r*P[t-1]*(1-P[t-1]) - ct.jags[t-1]/k,eps)),  # Process equation
      log(max(P[t-1] + 4*P[t-1]*r*P[t-1]*(1-P[t-1]) - ct.jags[t-1]/k,eps))) # assuming reduced r at B/k < 0.25
      P[t] ~ dlnorm(Pmean[t],itau2) # Introduce process error
      penm[t]  <- ifelse(P[t]<(eps+0.001),log(q*k*P[t])-log(q*k*(eps+0.001)),ifelse(P[t]>1,log(q*k*P[t])-log(q*k*(0.99)),0)) # penalty if Pmean is outside viable biomass
    }
    
    # Get Process error deviation 
    for(t in 1:nyr){
      proc.logB[t] <- log(P[t]*k)-log(exp(Pmean[t])*k)} 
      

    # Biomass priors/penalties are enforced as follows
    for (i in 1:nbk) {
      penb[i]  <- ifelse(P[b.yrs[i]]<b.prior[1,i],log(q*k*P[b.yrs[i]])-log(q*k*b.prior[1,i]),ifelse(P[b.yrs[i]]>b.prior[2,i],log(q*k*P[b.yrs[i]])-log(q*k*b.prior[2,i]),0)) 
      b.prior[3,i] ~ dnorm(penb[i],100)
    }
    
    for (t in 1:nyr){
      Fpen[t]   <- ifelse(ct[t]>(0.9*k*P[t]),ct[t]-(0.9*k*P[t]),0) # Penalty term on F > 1, i.e. ct>B
      pen.F[t]  ~ dnorm(Fpen[t],1000)
      pen.bk[t] ~ dnorm(penm[t],10000) 
      cpuem[t]  <- log(q*P[t]*k);
      bt[t]     ~ dlnorm(cpuem[t],pow(sigma2,-1));
    }
    
  # priors
  log.alpha               <- log((startbio[1]+startbio[2])/2) # needed for fit of first biomass
  sd.log.alpha            <- (log.alpha-log(startbio[1]))/4
  tau.log.alpha           <- pow(sd.log.alpha,-2)
  alpha                   ~  dlnorm(log.alpha,tau.log.alpha)

  # set realistic prior for q
  log.qm              <- mean(log(q.prior))
  sd.log.q            <- (log.qm-log(q.prior[1]))/2 # previous 4
  tau.log.q           <- pow(sd.log.q,-2)
  q                   ~  dlnorm(log.qm,tau.log.q)
    
  # define process (tau) and observation (sigma) variances as inversegamma prios
  itau2 ~ dgamma(4,0.01)
  tau2  <- 1/itau2
  tau   <- pow(tau2,0.5)
    
  isigma2 ~ dgamma(2,0.01)
  sigma2 <- 1/isigma2+pow(CV.cpue,2) # Add minimum realistic CPUE CV
  sigma <- pow(sigma2,0.5)

  
  log.rm              <- mean(log(prior.r))
  sd.log.r         <- abs(log.rm - log(prior.r[1]))/2
  tau.log.r           <- pow(sd.log.r,-2)
  #r                   ~  dlnorm(log.rm-0.5*pow(sd.log.r,2),tau.log.r) # bias-corrected   
  
  # bias-correct lognormal for k
  log.km              <- mean(log(prior.k))
  sd.log.k            <- abs(log.km-log(prior.k[1]))/2
  tau.log.k           <- pow(sd.log.k,-2)
  #k                   ~  dlnorm(log.km-0.5*pow(sd.log.k,2),tau.log.k) # bias-correct
  
  # Construct Multivariate lognormal (MVLN) prior
  mu.rk[1] <- log.rm
  mu.rk[2] <- log.km
  
  # Prior for correlation log(r) vs log(k)
  rho1 ~ dbeta(rk.cor.beta[1],rk.cor.beta[2])
  rho <- rho1-1

  # Construct Covariance matrix
  cov.rk[1,1] <- sd.log.r * sd.log.r
  cov.rk[1,2] <- sd.log.r * sd.log.k* rho
  cov.rk[2,1] <- sd.log.r * sd.log.k* rho
  cov.rk[2,2] <- sd.log.k * sd.log.k  
  
  # MVLN prior for r-k
  log.rk[1:2] ~ dmnorm(mu.rk[],inverse(cov.rk[,]))
  r <- exp(log.rk[1])
  k <- exp(log.rk[2])

} "    # end of JAGS model                  
      
    # Write JAGS model to file ----
    cat(Model, file="r2jags.bug")  
    
    
    j.inits <- function(){list("log.rk"=c(log(rnorm(1,mean=init.r,sd=0.2*init.r)),log(rnorm(1,mean=init.k,sd=0.1*init.k))),
                  "q"=rnorm(1,mean=init.q,sd=0.2*init.q),"itau2"=1000,"isigma2"=1000)}
    # run model ----
    jags_outputs <- jags.parallel(data=jags.data, 
                                  working.directory=NULL, inits=j.inits, 
                                  parameters.to.save=jags.save.params, 
                                  model.file="r2jags.bug", n.chains = n.chains, 
                                  n.burnin = 30000, n.thin = 10, 
                                  n.iter = 60000)
    
    # ------------------------------------------------------
    # Results from JAGS Schaefer ----
    # ------------------------------------------------------
    r_raw            <- as.numeric(mcmc(jags_outputs$BUGSoutput$sims.list$r))
    k_raw            <- as.numeric(mcmc(jags_outputs$BUGSoutput$sims.list$k))
    ct.jags          <- jags_outputs$BUGSoutput$sims.list$ct.jags
    cpue.jags        <- exp(jags_outputs$BUGSoutput$sims.list$cpuem)
    pe.logbt.jags   <- (jags_outputs$BUGSoutput$sims.list$proc.logB)
    
    # get catch predicted=adapted within error range by jags
    predC            <- apply(ct.jags,2,quantile,c(0.5,0.025,0.975)) 
    ct.jags          <- predC[1,]
    lcl.ct.jags      <- predC[2,]
    ucl.ct.jags      <- predC[3,]
    
    # get cpue predicted
    pred.cpue            <- apply(cpue.jags,2,quantile,c(0.5,0.025,0.975)) 
    cpue.jags          <- pred.cpue[1,]
    lcl.cpue.jags      <- pred.cpue[2,]
    ucl.cpue.jags      <- pred.cpue[3,]
    
    # get process error on log(biomass)   pred.cpue            <- apply(cpue.jags,2,quantile,c(0.5,0.025,0.975)) 
    pred.pe         <- apply(pe.logbt.jags,2,quantile,c(0.5,0.025,0.975)) 
    pe.jags         <- pred.pe[1,]
    lcl.pe.jags     <- pred.pe[2,]
    ucl.pe.jags     <- pred.pe[3,]

    #------------------------------------------------------------------
    mean.log.r.jags  <- mean(log(r_raw))
    sd.log.r.jags    <- sd(log(r_raw))
    r.jags           <- exp(mean.log.r.jags)
    lcl.r.jags       <- exp(mean.log.r.jags - 1.96*sd.log.r.jags)
    ucl.r.jags       <- exp(mean.log.r.jags + 1.96*sd.log.r.jags)
    mean.log.k.jags  <- mean(log(k_raw))
    sd.log.k.jags    <- sd(log(k_raw))
    k.jags           <- exp(mean.log.k.jags)
    lcl.k.jags       <- exp(mean.log.k.jags - 1.96*sd.log.k.jags)
    ucl.k.jags       <- exp(mean.log.k.jags + 1.96*sd.log.k.jags)
    MSY.posterior     <- r_raw*k_raw/4 
    mean.log.MSY.jags <- mean(log(MSY.posterior))
    sd.log.MSY.jags   <- sd(log(MSY.posterior))
    MSY.jags          <- exp(mean.log.MSY.jags)
    lcl.MSY.jags      <- exp(mean.log.MSY.jags - 1.96*sd.log.MSY.jags)
    ucl.MSY.jags      <- exp(mean.log.MSY.jags + 1.96*sd.log.MSY.jags)

      q_out           <- as.numeric(mcmc(jags_outputs$BUGSoutput$sims.list$q))
      mean.log.q      <- mean(log(q_out))
      sd.log.q        <- sd(log(q_out))
      mean.q          <- exp(mean.log.q)
      lcl.q           <- exp(mean.log.q-1.96*sd.log.q)
      ucl.q           <- exp(mean.log.q+1.96*sd.log.q)
      F.bt.jags       <- mean.q*ct.raw/bt # F from raw data
      Fmsy.jags       <- r.jags/2 
      F.bt_Fmsy.jags  <- vector() # initialize vector
      for(z in 1: length(F.bt.jags)) {
        F.bt_Fmsy.jags[z] <- ifelse(is.na(bt[z])==T,NA,F.bt.jags[z]/
                               ifelse(((bt[z]/mean.q)/k.jags)<0.25,Fmsy.jags*4*(bt[z]/mean.q)/k.jags,Fmsy.jags))}

    # get relative biomass P=B/k as predicted by BSM, including predictions for years with NA abundance
    all.P    <- jags_outputs$BUGSoutput$sims.list$P # matrix with P distribution by year
    quant.P  <- apply(all.P,2,quantile,c(0.025,0.5,0.975),na.rm=T)
    
    # get k, r posterior
    all.k  <- jags_outputs$BUGSoutput$sims.list$k # matrix with k distribution by year
    all.r  <- jags_outputs$BUGSoutput$sims.list$r # matrix with r distribution by year
    
    # get B/Bmys posterior 
    all.b_bmsy=NULL
    for(t in 1:ncol(all.P)){
      all.b_bmsy  <- cbind(all.b_bmsy,all.P[,t]*2)}
    
    # get F/Fmsy posterior 
    all.F_Fmsy=NULL
    for(t in 1:ncol(all.P)){
      all.F_Fmsy      <- cbind(all.F_Fmsy,(ct.jags[t]/(all.P[,t]*all.k))/
                            ifelse(all.P[,t]>0.25,all.r/2,all.r/2*4*all.P[,t]))}
    quant.all.F_Fmsy  <- apply(all.F_Fmsy,2,quantile,c(0.025,0.5,0.975),na.rm=T)
    F_Fmsy.jags       <- quant.all.F_Fmsy[2,]
    lcl.F_Fmsy.jags   <- quant.all.F_Fmsy[1,]
    ucl.F_Fmsy.jags   <- quant.all.F_Fmsy[3,]
    
    # get variance and correlation between log(r) and log(k)
    log.r.var    <- var(x=log(r_raw))
    log.k.var    <- var(x=log(k_raw))
    log.rk.cor   <- cor(x=log(r_raw),y=log(k_raw))
    log.rk.cov   <- cov(x=log(r_raw),y=log(k_raw))
    
  } # end of MCMC Schaefer loop 
  
  # --------------------------------------------
  # Get results for management ----
  # --------------------------------------------
  if(FullSchaefer==F | force.cmsy==T) { # if only CMSY is available or shall be used
    MSY   <-MSY.est; lcl.MSY<-lcl.MSY.est; ucl.MSY<-ucl.MSY.est 
    Bmsy  <-k.est/2; lcl.Bmsy<-lcl.k.est/2; ucl.Bmsy<-ucl.k.est/2
    Fmsy  <-r.est/2; lcl.Fmsy<-lcl.r.est/2; ucl.Fmsy<-ucl.r.est/2
    F.Fmsy<-F_Fmsy.CMSY;lcl.F.Fmsy<-lcl.F_Fmsy.CMSY; ucl.F.Fmsy<-ucl.F_Fmsy.CMSY
    B.Bmsy<-2*median.btv[1:nyr];lcl.B.Bmsy<-2*lcl.btv[1:nyr];ucl.B.Bmsy<-2*ucl.btv[1:nyr]
    if(is.na(sel.yr)==F){B.Bmsy.sel<-2*median.btv[yr==sel.yr]}
    
  } else { # if FullSchaefer is TRUE
    MSY   <-MSY.jags; lcl.MSY<-lcl.MSY.jags; ucl.MSY<-ucl.MSY.jags 
    Bmsy  <-k.jags/2; lcl.Bmsy<-lcl.k.jags/2; ucl.Bmsy<-ucl.k.jags/2
    Fmsy  <-r.jags/2; lcl.Fmsy<-lcl.r.jags/2; ucl.Fmsy<-ucl.r.jags/2
    F.Fmsy<-F_Fmsy.jags; lcl.F.Fmsy<-lcl.F_Fmsy.jags; ucl.F.Fmsy<-ucl.F_Fmsy.jags
    B.Bmsy<-2*quant.P[2,];lcl.B.Bmsy<-2*quant.P[1,];ucl.B.Bmsy<-2*quant.P[3,]
    if(is.na(sel.yr)==F) {B.Bmsy.sel<-2*quant.P[2,][yr==sel.yr]}
  }
  # the following code works for CMSY and for BSM
  B          <-B.Bmsy*Bmsy;lcl.B<-lcl.B.Bmsy*Bmsy;ucl.B<-ucl.B.Bmsy*Bmsy
  B.last     <-B[nyr];lcl.B.last<-lcl.B[nyr];ucl.B.last<-ucl.B[nyr]
  B.Bmsy.last<-B.Bmsy[nyr];lcl.B.Bmsy.last<-lcl.B.Bmsy[nyr];ucl.B.Bmsy.last<-ucl.B.Bmsy[nyr]
  
  if(FullSchaefer==T & force.cmsy==F){cm=ct.jags} else{cm=ct.raw}
  
  Fm           <- cm/B;lcl.F<-cm/ucl.B;ucl.F<-cm/lcl.B
  Fmsy.vec     <- ifelse(B.Bmsy>0.5,Fmsy,Fmsy*2*B.Bmsy)
  lcl.Fmsy.vec <- ifelse(B.Bmsy>0.5,lcl.Fmsy,lcl.Fmsy*2*B.Bmsy)
  ucl.Fmsy.vec <- ifelse(B.Bmsy>0.5,ucl.Fmsy,ucl.Fmsy*2*B.Bmsy)
  
  F.last     <-Fm[nyr];lcl.F.last<-lcl.F[nyr];ucl.F.last<-ucl.F[nyr]
  Fmsy.last  <-Fmsy.vec[nyr];lcl.Fmsy.last<-lcl.Fmsy.vec[nyr];ucl.Fmsy.last<-ucl.Fmsy.vec[nyr]
  F.Fmsy.last<-F.Fmsy[nyr];lcl.F.Fmsy.last<-lcl.F.Fmsy[nyr];ucl.F.Fmsy.last<-ucl.F.Fmsy[nyr]
  
  if(is.na(sel.yr)==F){
    B.sel<-B.Bmsy.sel*Bmsy
    F.sel<-ct.raw[yr==sel.yr]/B.sel
    F.Fmsy.sel<-F.sel/Fmsy.vec[yr==sel.yr]
  }
  
  # ------------------------------------------
  # print input and results to screen ----
  #-------------------------------------------
  cat("---------------------------------------\n")
  cat("Species:", cinfo$ScientificName[cinfo$Stock==stock], ", stock:",stock,"\n")
  cat(cinfo$Name[cinfo$Stock==stock], "\n")
  cat("Region:",cinfo$Region[cinfo$Stock==stock],",",cinfo$Subregion[cinfo$Stock==stock],"\n")
  cat("Catch data used from years", min(yr),"-", max(yr),", abundance =", btype, "\n")
  cat("Prior initial relative biomass =", startbio[1], "-", startbio[2],ifelse(is.na(stb.low)==T,"default","expert"), "\n")
  cat("Prior intermediate rel. biomass=", intbio[1], "-", intbio[2], "in year", int.yr,ifelse(is.na(intb.low)==T,"default","expert"), "\n")
  cat("Prior final relative biomass   =", endbio[1], "-", endbio[2],ifelse(is.na(endb.low)==T,"default","expert"), "\n")
  cat("Prior range for r =", format(prior.r[1],digits=2), "-", format(prior.r[2],digits=2),ifelse(is.na(r.low)==T,"default","expert,"),  
      ", prior range for k =", prior.k[1], "-", prior.k[2],"\n")
  # if Schaefer and CPUE, print prior range of q
  if(FullSchaefer==T) {
    cat("Prior range of q =",q.prior[1],"-",q.prior[2],", assumed effort creep",e.creep,"%\n")
  }
  # results of CMSY analysis
  cat("\nResults of CMSY analysis \n")
  cat("-------------------------\n")
  cat("Altogether", n.viable.b, "viable trajectories for", n.viable.pt," r-k pairs were found \n")
  cat("r   =", r.est,", 95% CL =", lcl.r.est, "-", ucl.r.est,", k =", k.est,", 95% CL =", lcl.k.est, "-", ucl.k.est,"\n")
  cat("MSY =", MSY.est,", 95% CL =", lcl.MSY.est, "-", ucl.MSY.est,"\n")
  cat("Relative biomass in last year =", median.btv.lastyr, "k, 2.5th perc =", lcl.median.btv.lastyr, 
      ", 97.5th perc =", ucl.median.btv.lastyr,"\n")
  cat("Exploitation F/(r/2) in last year =", F_Fmsy.CMSY[nyr],", 2.5th perc =",lcl.F_Fmsy.CMSY[nyr],
      ", 97.5th perc =",ucl.F_Fmsy.CMSY[nyr],"\n\n")
  
  # print results from full Schaefer if available
  if(FullSchaefer==T) {
    cat("Results from Bayesian Schaefer model (BSM) using catch &",btype,"\n")
    cat("------------------------------------------------------------\n")
    cat("q   =", mean.q,", lcl =", lcl.q, ", ucl =", ucl.q,"\n")
    cat("r   =", r.jags,", 95% CL =", lcl.r.jags, "-", ucl.r.jags,", k =", k.jags,", 95% CL =", lcl.k.jags, "-", ucl.k.jags,", r-k log correlation =", log.rk.cor,"\n")
    cat("MSY =", MSY.jags,", 95% CL =", lcl.MSY.jags, "-", ucl.MSY.jags,"\n")
    cat("Relative biomass in last year =", quant.P[2,][nyr], "k, 2.5th perc =",quant.P[1,][nyr], 
        ", 97.5th perc =", quant.P[3,][nyr],"\n")
    cat("Exploitation F/(r/2) in last year =", F_Fmsy.jags[nyr],", 2.5th perc =",lcl.F_Fmsy.jags[nyr],
        ", 97.5th perc =",ucl.F_Fmsy.jags[nyr],"\n\n")
  }
  
  # print results to be used in management
  cat("Results for Management (based on",ifelse(FullSchaefer==F | force.cmsy==T,"CMSY","BSM"),"analysis) \n")
  cat("-------------------------------------------------------------\n")
  if(force.cmsy==T) cat("Mangement results based on CMSY because abundance data seem unrealistic\n")
  cat("Fmsy =",Fmsy,", 95% CL =",lcl.Fmsy,"-",ucl.Fmsy,"(if B > 1/2 Bmsy then Fmsy = 0.5 r)\n") 
  cat("Fmsy =",Fmsy.last,", 95% CL =",lcl.Fmsy.last,"-",ucl.Fmsy.last,"(r and Fmsy are linearly reduced if B < 1/2 Bmsy)\n")
  cat("MSY  =",MSY,", 95% CL =",lcl.MSY,"-",ucl.MSY,"\n") 
  cat("Bmsy =",Bmsy,", 95% CL =",lcl.Bmsy,"-",ucl.Bmsy,"\n") 
  cat("Biomass in last year =",B.last,", 2.5th perc =", lcl.B.last, ", 97.5 perc =",ucl.B.last,"\n")
  cat("B/Bmsy in last year  =",B.Bmsy.last,", 2.5th perc =", lcl.B.Bmsy.last, ", 97.5 perc =",ucl.B.Bmsy.last,"\n")
  cat("Fishing mortality in last year =",F.last,", 2.5th perc =", lcl.F.last, ", 97.5 perc =",ucl.F.last,"\n")
  cat("Exploitation F/Fmsy  =",F.Fmsy.last,", 2.5th perc =", lcl.F.Fmsy.last, ", 97.5 perc =",ucl.F.Fmsy.last,"\n")
  
  # show stock status and exploitation for optional selected year
  if(is.na(sel.yr)==F) {
    cat("\nStock status and exploitation in",sel.yr,"\n")
    cat("Biomass =",B.sel, ", B/Bmsy =",B.Bmsy.sel,", F =",F.sel,", F/Fmsy =",F.Fmsy.sel,"\n") }
  
  # indicate if less than 5 years of biomass or CPUE are available
  if(btype !="None" & length(bt[is.na(bt)==F])<nab) {
    cat("Less than",nab,"years with abundance data available, shown on second axis\n")
  }
  cat("Comment:", comment,"\n")
  cat("----------------------------------------------------------\n")
  
  
  # -----------------------------------------
  # Plot results ----
  # -----------------------------------------
  # plot best r-k from full Schaefer analysis in prior space of graph B
  if(FullSchaefer==T) {
    points(x=r.jags, y=k.jags, pch=19, col="red")  
    lines(x=c(lcl.r.jags, ucl.r.jags),y=c(k.jags,k.jags), col="red")
    lines(x=c(r.jags,r.jags),y=c(lcl.k.jags, ucl.k.jags), col="red")
  }
  lines(x=c(prior.r[1],prior.r[2],prior.r[2],prior.r[1],prior.r[1]), # plot original prior range
        y=c(prior.k[1],prior.k[1],prior.k[2],prior.k[2],prior.k[1]),lty="dotted")
  
  
  # (c) Analysis of viable r-k plot -----
  # ----------------------------
  max.y    <- max(c(ifelse(FullSchaefer==T,max(k_raw,ucl.k.jags),NA), max(kv.all)), 
                  ifelse(substr(id_file,1,3)=="Sim",1.2*true.k,max(kv.all)),na.rm=T)
  min.y    <- min(c(ifelse(FullSchaefer==T,min(k_raw),NA), 0.9*min(kv.all)), 
                  ifelse(substr(id_file,1,3)=="Sim",0.8*true.k,0.9*min(kv.all)),na.rm=T)
  max.x    <- max(c(ifelse(FullSchaefer==T,max(r_raw),NA),max(rv.all)),na.rm=T)
  min.x    <- min(c(ifelse(FullSchaefer==T,min(r_raw),NA),0.9*lcl.r.est,prior.r[1]),na.rm=T)
    
  plot(x=rv.all, y=kv.all, xlim=c(min.x,max.x), 
       ylim=c(min.y,max.y), 
       pch=16, col="gray",log="xy", bty="l",
       xlab="", ylab="k (1000 tonnes)", main="C: Analysis of viable r-k",  cex.main = 1.8, cex.lab = 1.55, cex.axis = 1.5)
  title(xlab = "r", line = 2.25, cex.lab = 1.55)
 
  # plot r-k pairs from MCMC
  if(FullSchaefer==T) {points(x=r_raw, y=k_raw, pch=16,cex=0.5)}
  
  # plot best r-k from full Schaefer analysis
  if(FullSchaefer==T) {
    points(x=r.jags, y=k.jags, pch=19, col="red")  
    lines(x=c(lcl.r.jags, ucl.r.jags),y=c(k.jags,k.jags), col="red")
    lines(x=c(r.jags,r.jags),y=c(lcl.k.jags, ucl.k.jags), col="red")
  }
  
  # plot blue dot for CMSY r-k, with 95% CL lines 
  points(x=r.est, y=k.est, pch=19, col="blue")
  lines(x=c(lcl.r.est, ucl.r.est),y=c(k.est,k.est), col="blue")
  lines(x=c(r.est,r.est),y=c(lcl.k.est, ucl.k.est), col="blue")
  
  # (d) Pred. biomass plot ----
  #--------------------
  # determine k to use for red line in b/k plot 
  if(FullSchaefer==T)  {k2use <- k.jags} else {k2use <- k.est}
  # determine hight of y-axis in plot
  max.y  <- max(c(ifelse(btype=="biomass",max(bt/k2use,na.rm=T),NA),
                  ifelse(btype=="CPUE" & length(bt[is.na(bt)==F])>=nab,max(bt/(mean.q*k2use),na.rm=T),NA),
                  max(ucl.btv),0.6,startbio[2], endbio[2]),
                ifelse(FullSchaefer==T & btype=="biomass",max(bt[is.na(bt)==F]/lcl.k.jags,na.rm=T),NA),
                ifelse(FullSchaefer==T & btype=="CPUE",1.1*max(bt/(mean.q*lcl.k.jags),na.rm=T),NA), na.rm=T)
  
  # Main plot of relative CMSY biomass
  plot(x=yr,y=median.btv[1:nyr], lwd=1.5, xlab="", ylab="Relative biomass B/k", type="l",
       ylim=c(0,max.y), bty="l", main="D: Stock size",col="blue",  cex.main = 1.8, cex.lab = 1.55, cex.axis = 1.5)
  lines(x=yr, y=lcl.btv[1:nyr],type="l",lty="dotted",col="blue")
  lines(x=yr, y=ucl.btv[1:nyr],type="l",lty="dotted",col="blue")
  # plot lines for 0.5 and 0.25 biomass
  abline(h=0.5, lty="dashed")
  abline(h=0.25, lty="dotted")
  # Add BSM
  if(FullSchaefer==T){
  bk.bsm = apply(all.P,2,quantile,c(0.025,0.5,0.975))
  lines(x=yr, y=bk.bsm[2,1:nyr],type="l",lty=1,col="red")
  lines(x=yr, y=bk.bsm[1,1:nyr],type="l",lty="dotted",col="red")
  lines(x=yr, y=bk.bsm[3,1:nyr],type="l",lty="dotted",col="red")
  points(x=yr,y=bt/(mean.q*k.jags),pch=21,bg="grey")  
  }
  # plot biomass windows
  lines(x=c(yr[1],yr[1]), y=startbio, col="blue")
  lines(x=c(int.yr,int.yr), y=intbio, col="blue")  
  lines(x=c(max(yr),max(yr)), y=endbio, col="blue")  
  
  # if CPUE has been corrected for effort creep, display uncorrected CPUE
  if(btype=="CPUE" & FullSchaefer==T & e.creep.line==T & is.na(e.creep)==FALSE) {
    lines(x=yr,y=bt.raw/(mean.q*k.jags),type="l", col="green", lwd=1)
  }

  # (e) Exploitation rate plot ----
  # -------------------------
  # if CPUE data are available but fewer than nab years, plot on second axis
  if(btype == "CPUE" | btype=="biomass") {
    q=1/(max(median.btv[1:nyr][is.na(bt)==F],na.rm=T)*k.est/max(bt,na.rm=T))
    u.cpue      <- q*ct/bt
  }
  # determine upper bound of Y-axis
  max.y <- max(c(1.5,ucl.F_Fmsy.CMSY,ifelse(FullSchaefer==T,max(c(F.bt.jags/Fmsy.jags,ucl.F_Fmsy.jags),na.rm=T),NA),na.rm=T),na.rm=T)
  max.y <- ifelse(max.y>10,10,max.y)
  # plot F from CMSY
  plot(x=yr,y=F_Fmsy.CMSY, type="l", bty="l", lwd=1.5, ylim=c(0,max.y), xlab="", 
       ylab="F / Fmsy", main="E: Exploitation rate", col="blue",  cex.main = 1.8, cex.lab = 1.55, cex.axis = 1.5)
  lines(x=yr,y=lcl.F_Fmsy.CMSY,lty="dotted",col="blue")
  lines(x=yr,y=ucl.F_Fmsy.CMSY,lty="dotted",col="blue")
  abline(h=1, lty="dashed")
  
  # plot F/Fmsy as points from observed catch and CPUE and as red curves from BSM predicted catch and biomass
  if(FullSchaefer==T){
    points(x=yr, y=F.bt_Fmsy.jags, pch=21,bg="grey")
    lines(x=yr,y=F_Fmsy.jags, col="red")
    lines(x=yr,y=lcl.F_Fmsy.jags, col="red",lty="dotted")
    lines(x=yr,y=ucl.F_Fmsy.jags, col="red",lty="dotted")
  }
 
  # (f) Parabola plot ----
  #-------------------------
  max.y <- max(c(ct/MSY.est,ifelse(FullSchaefer==T,max(ct/MSY.jags),NA),1.2),na.rm=T)
  # plot parabola
  x=seq(from=0,to=2,by=0.001)
  y.c  <- ifelse(x>0.25,1,ifelse(x>0.125,4*x,exp(-10*(0.125-x))*4*x)) # correction for low recruitment below half and below quarter of Bmsy
  y=(4*x-(2*x)^2)*y.c
  plot(x=x, y=y, xlim=c(0,1), ylim=c(0,max.y), type="l", bty="l",xlab="", 
       ylab="Catch / MSY", main="F: Equilibrium curve",  cex.main = 1.8, cex.lab = 1.55, cex.axis = 1.5)
  title(xlab= "Relative biomass B/k", line = 2.25, cex.lab = 1.55)
  
  # plot catch against CMSY estimates of relative biomass
  lines(x=median.btv[1:nyr], y=ct/MSY.est, pch=16, col="blue", lwd=1)
  points(x=median.btv[1], y=ct[1]/MSY.est[1], pch=0, cex=2, col="blue")
  points(x=median.btv[nyr], y=ct[length(ct)]/MSY.est[length(MSY.est)],cex=2,pch=2,col="blue")
  
  # for CPUE, plot catch scaled by BSM MSY against observed biomass derived as q * CPUE scaled by BSM k
  if(FullSchaefer==T) {
    points(x=bt/(mean.q*k.jags), y=ct/MSY.jags, pch=21,bg="grey")
    lines(x=median.btv[1:nyr], y=predC[1,]/MSY.jags, pch=16, col="red",lwd=1)
    points(x=median.btv[1], y=predC[1,][1]/MSY.jags, pch=0, cex=2, col="red")
    points(x=median.btv[nyr], y=predC[1,][length(ct)]/MSY.jags[length(MSY.jags)], pch=2, cex=2,col="red")
  }
  
  #analysis.plot <- recordPlot()
  
  #save analytic chart to JPEG file
  if (save.plots==TRUE) 
  {
    jpgfile<-paste(stock,"_AN.jpg",sep="")
	
	if (retrosp.step>0) jpgfile<-gsub(".jpg", paste0("_retrostep_",retrosp.step,".jpg"), jpgfile) #modification added to save all steps in retrospective analysis
	
    dev.copy(jpeg,jpgfile,
             width = 1024, 
             height = 768, 
             units = "px", 
             pointsize = 18,
             quality = 95,
             res=80,
             antialias="cleartype")
    dev.off()
  }
  
  #---------------------------------------------
  # Plot Management-Graphs if desired ----
  #---------------------------------------------
  if(mgraphs==T) {
    # open window for plot of four panels
    if(grepl("win",tolower(Sys.info()['sysname']))) {windows(14,12)}
    par(mfrow=c(2,2))  
    # make margins narrower
    par(mar=c(3.1,4.2,2.1,2.1))
    
    #---------------------
    # plot catch with MSY ----
    #---------------------
    max.y <- max(c(1.1*max(cm),ucl.MSY),na.rm=T)
    plot(x=yr,rep(0,nyr),type="n",ylim=c(0,max.y), bty="l", main=paste("Catch",stock), 
         xlab="",ylab="Catch (1000 tonnes/year)",  cex.main = 1.6, cex.lab = 1.35, cex.axis = 1.35)
    rect(yr[1],lcl.MSY,yr[nyr],ucl.MSY,col="lightgray", border=NA)
    lines(x=c(yr[1],yr[nyr]),y=c(MSY,MSY),lty="dashed", col="black", lwd=1.5)
    lines(x=yr, y=cm, lwd=2) # 
    text("MSY",x=end.yr-1.5, y=MSY+MSY*0.1, cex = .75)
    
    #----------------------------------------
    # Plot of estimated biomass relative to Bmsy 
    #----------------------------------------
    # plot empty frame
    plot(yr, rep(0,nyr),type="n", ylim=c(0,max(c(2, max(ucl.B.Bmsy)))), ylab="B / Bmsy",xlab="", main="Stock size", bty="l",  cex.main = 1.6, cex.lab = 1.35, cex.axis = 1.35) 
    # plot gray area of uncertainty in predicted biomass
    polygon(c(yr,rev(yr)), c(lcl.B.Bmsy,rev(ucl.B.Bmsy)),col="lightgray", border=NA)
    # plot median biomass
    lines(yr,B.Bmsy,lwd=2)
    # plot lines for Bmsy and 0.5 Bmsy
    lines(x=c(yr[1],yr[nyr]),y=c(1,1), lty="dashed", lwd=1.5)
    lines(x=c(yr[1],yr[nyr]),y=c(0.5,0.5), lty="dotted", lwd=1.5)
    
    # -------------------------------------
    ## Plot of exploitation rate
    # -------------------------------------
    # plot empty frame
    plot(yr, rep(0,nyr),type="n", ylim=c(0,max(c(2,ucl.F.Fmsy))), 
         ylab="F / Fmsy",xlab="", main="Exploitation", bty="l",  cex.main = 1.6, cex.lab = 1.35, cex.axis = 1.35)  
    # plot gray area of uncertainty in predicted exploitation
    polygon(c(yr,rev(yr)), c(lcl.F.Fmsy,rev(ucl.F.Fmsy)),col="lightgray", border=NA)
    # plot median exploitation rate
    lines(x=yr,y=F.Fmsy,lwd=2)
    # plot line for u.msy
    lines(x=c(yr[1],yr[nyr]),y=c(1,1), lty="dashed", lwd=1.5)
    
    # -------------------------------------
    ## plot stock-status graph
    # -------------------------------------
    
    if(FullSchaefer==T & force.cmsy==F) {x.F_Fmsy = all.F_Fmsy[,nyr]
    y.b_bmsy = all.b_bmsy[,nyr]} else {
      log.rk = cbind(rem.log.r,rem.log.k)
      rem.log.btv.lastyr = log(mdat.all[rem,nyr])
      log.bbmsy = rem.log.btv.lastyr+log(2)
      log.ffmsy = (log(ct.raw[nyr])-(rem.log.btv.lastyr+rem.log.k))-(rem.log.r-log(2))
      # get mean after all the CMSY subsetting (can't match with biomass sbmsetting)
      mu.kobe = log(c(F.Fmsy.last,B.Bmsy.last))
      # Get covariance of the 2 vectors
      cov.kobe = cov(cbind(log.ffmsy,log.bbmsy)) 
      # Generate 10000 new random deviates from a MVN
      log.kobe.mvn = rmvnorm(10000 ,mean = mu.kobe,sigma = cov.kobe)
      kobe.mvn = exp(log.kobe.mvn)
      # Generate 10000 new random deviates from a MVN
      x.F_Fmsy =exp(log.kobe.mvn[,1])
      y.b_bmsy =exp(log.kobe.mvn[,2])
    }
    
    kernelF <- ci2d(x.F_Fmsy,y.b_bmsy,nbins=201,factor=2.2,ci.levels=c(0.50,0.80,0.75,0.90,0.95),show="none")
    c1 <- c(-1,100)
    c2 <- c(1,1)
    
    max.x1   <- max(c(2, max(kernelF$contours$"0.95"$x,F.Fmsy),na.rm =T))
    max.x    <- ifelse(max.x1 > 5,min(max(5,F.Fmsy*2),8),max.x1)
    max.y    <- max(max(2,quantile(y.b_bmsy,0.96)))
    
    plot(1000,1000,type="b", xlim=c(0,max.x), ylim=c(0,max.y),lty=3,xlab="",ylab="B / Bmsy", bty="l",  cex.main = 1.6, cex.lab = 1.35, cex.axis = 1.35)
    mtext("F / Fmsy",side=1, line=2, cex=1.05)
    #mtext("B / Bmsy",side=2, line=2.2,cex=1.15)
    
    # extract interval information from ci2d object
    # and fill areas using the polygon function
    polygon(kernelF$contours$"0.95",lty=2,border=NA,col="cornsilk4")
    polygon(kernelF$contours$"0.8",border=NA,lty=2,col="grey")
    polygon(kernelF$contours$"0.5",border=NA,lty=2,col="cornsilk2")
    
    ## Add points and trajectory lines
    lines(c1,c2,lty=3,lwd=0.7)
    lines(c2,c1,lty=3,lwd=0.7)
    lines(F.Fmsy,B.Bmsy, lty=1,lwd=1.)
    
    # points(F.Fmsy,B.Bmsy,cex=0.8,pch=4)
    points(F.Fmsy[1],B.Bmsy[1],col=1,pch=22,bg="white",cex=1.5)
    points(F.Fmsy[which(yr==int.yr)],B.Bmsy[which(yr==int.yr)],col=1,pch=21,bg="white",cex=1.5)
    points(F.Fmsy[nyr],B.Bmsy[nyr],col=1,pch=24,bg="white",cex=1.5)
    
    ## Add legend
    legend('topright', inset = .03, c(paste(start.yr),paste(int.yr),paste(end.yr),"50% C.I.","80% C.I.","95% C.I."), 
           lty=c(1,1,1,-1,-1,-1),pch=c(22,21,24,22,22,22),pt.bg=c(rep("white",3),"cornsilk2","grey","cornsilk4"), 
           col=1,lwd=.8,cex=0.85,pt.cex=c(rep(1.1,3),1.5,1.5,1.5),bty="n",y.intersp = 1.1)
    #End of Biplot
    
  } # end of management graphs
  
  #management.plot <- recordPlot()
  
  # save management chart to JPEG file
  if (save.plots==TRUE & mgraphs==TRUE) 
  {
    jpgfile<-paste(stock,"_MAN.jpg",sep="")
	if (retrosp.step>0) jpgfile<-gsub(".jpg", paste0("_retrostep_",retrosp.step,".jpg"), jpgfile) #modification added to save all steps in retrospective analysis
    dev.copy(jpeg,jpgfile,
             width = 1024, 
             height = 768, 
             units = "px", 
             pointsize = 18,
             quality = 95,
             res=80,
             antialias="cleartype")
    dev.off()
  }
  
  #----------------------------------------------------------
  #><> Optional prior - posterior plots 
  #---------------------------------------------------------
  if(pp.plot==T) {
    # open window for plot of four panels
    if(grepl("win",tolower(Sys.info()['sysname']))) {windows(17,12)} 
    # make margins narrower
    par(mfrow=c(2,3),mar=c(4.5,4.5,2,0.5))
    greycol = c(grey(0.7,0.5),grey(0.3,0.5)) # changed 0.6 to 0.7
  
    # plot PP-diagnostics for CMSY
    # r
    rk <- exp(mvn(n=10000,mean.log.r=mean.log.r,sd.log.r=sd.log.r,mean.log.k=mean.log.k,sd.log.k=sd.log.k))
    post.cmsy = exp(rem.log.r)
    nmc = length(post.cmsy) 
    rpr = rk[,1]
    pdf.cmsy = stats::density(post.cmsy,adjust=2)  
	prior.mean.log.r=mean(log(prior.r))
	prior.sd.log.r=(log(prior.r[2])-log(prior.r[1]))/4
    prior.samples.r<-rlnorm(3000, meanlog = prior.mean.log.r, sdlog = prior.sd.log.r)
	prior = stats::density(prior.samples.r,adjust=2) # stats::density(rk[,1],adjust=2)   # modification by GP 03/12/2019
	prior.r<-prior
    plot(pdf.cmsy,type="l",ylim=range(prior$y,pdf.cmsy$y*1.1),xlim=range(c(pdf.cmsy$x,rpr,max(pdf.cmsy$x,rpr)*1.1)),
         yaxt="n",xlab="r",ylab="",xaxs="i",yaxs="i",main="",bty="l",cex.lab = 1.55, cex.axis = 1.5)
    polygon(c((prior$x),rev(prior$x)),c(prior$y,rep(0,length(sort(prior$y)))),col=greycol[1])
    polygon(c(pdf.cmsy$x,rev(pdf.cmsy$x)),c(pdf.cmsy$y,rep(0,length(pdf.cmsy$y))),col=greycol[2])
    PPVR.cmsy = round((sd(post.cmsy)/mean(post.cmsy))^2/(sd(rpr)/mean(rpr))^2,2)  
    PPVM.cmsy = round(mean(post.cmsy)/mean(rpr),2)
    pp = c(paste("PPVR =",PPVR.cmsy))
    legend('right',c("Prior","Posterior"),pch=22,pt.cex=1.5,pt.bg = greycol,bty="n",cex=1.5)
    legend("topright",pp,cex=1.4,bty="n")  
    # k
	post.cmsy = exp(rem.log.k)
    nmc = length(post.cmsy) 
    rpr = rk[,2]
    pdf.cmsy = stats::density(post.cmsy,adjust=2)  
    prior.mean.log.k <- mean(log(prior.k))
	prior.sd.log.k   <- (log(prior.k[2])-log(prior.k[1]))/4
	prior.samples.k<-rlnorm(3000, meanlog = prior.mean.log.k, sdlog = prior.sd.log.k)
	prior = stats::density(prior.samples.k,adjust=2) # stats::density(rk[,2],adjust=2)   # modification by GP 03/12/2019
	prior.k<-prior
	plot(pdf.cmsy,type="l",ylim=range(prior$y,pdf.cmsy$y*1.1),xlim=range(c(pdf.cmsy$x,rpr,max(pdf.cmsy$x,rpr)*1.1)),yaxt="n",xlab="k (1000 tonnes)",ylab="",xaxs="i",yaxs="i",main="",bty="l",cex.lab = 1.55, cex.axis = 1.5)
    polygon(c((prior$x),rev(prior$x)),c(prior$y,rep(0,length(sort(prior$y)))),col=greycol[1])
    polygon(c(pdf.cmsy$x,rev(pdf.cmsy$x)),c(pdf.cmsy$y,rep(0,length(pdf.cmsy$y))),col=greycol[2])
    PPVR.cmsy = round((sd(post.cmsy)/mean(post.cmsy))^2/(sd(rpr)/mean(rpr))^2,2)  
    PPVM.cmsy = round(mean(post.cmsy)/mean(rpr),2)
    pp = c(paste("PPVR =",PPVR.cmsy))
    legend("topright",pp,cex=1.4,bty="n")  
    mtext(paste0("CMSY prior & posterior distributions for ",stock),  side=3,cex=1.5)
    # MSY
    post.cmsy = exp(rem.log.k)*exp(rem.log.r)/4
    rpr = rk[,1]*rk[,2]/4
    pdf.cmsy = stats::density(post.cmsy,adjust=2)  
	prior.cmsy.calc = (prior.r$y*prior.k$y)/4 # modification by GP 03/12/2019
	#prior.cmsy.sd = sd(prior.r$y*prior.k$y)/16
	prior = stats::density(rpr,adjust=2)
	
    plot(prior.cmsy.calc,type="l",ylim=range(prior$y,pdf.cmsy$y*1.1),xlim=range(c(pdf.cmsy$x,rpr,max(pdf.cmsy$x,rpr)*1.1)),yaxt="n",xlab="MSY (1000 tonnes/year)",ylab="",xaxs="i",yaxs="i",main="",bty="l",cex.lab = 1.55, cex.axis = 1.5)
    polygon(c((prior$x),rev(prior$x)),c(prior$y,rep(0,length(sort(prior$y)))),col=greycol[1])
    polygon(c(pdf.cmsy$x,rev(pdf.cmsy$x)),c(pdf.cmsy$y,rep(0,length(pdf.cmsy$y))),col=greycol[2])
    PPVR.cmsy = round((sd(post.cmsy)/mean(post.cmsy))^2/(sd(rpr)/mean(rpr))^2,2)  
    PPVM.cmsy = round(mean(post.cmsy)/mean(rpr),2)
    pp = c(paste("PPVR =",PPVR.cmsy))
    legend("topright",pp,cex=1.4,bty="n")  
    
    # bk1
    post.cmsy = rem.btv.all[,1]
    nmc = length(post.cmsy) 
    rpr = startbio
    pdf.cmsy = stats::density(post.cmsy,adjust=2)  
    prior = rpr
	prior.height<-1/(prior[2]-prior[1])	# modification by GP 03/12/2019
	    plot(pdf.cmsy,type="l",ylim=range(pdf.cmsy$y),xlim=range(c(pdf.cmsy$x,0.3*rpr,min(1.7*rpr[2],1.05),max(pdf.cmsy$x,rpr)*1.1)),yaxt="n",xlab=paste0("B/k ",yr[1]),ylab="",xaxs="i",yaxs="i",main="",bty="l",cex.lab = 1.55, cex.axis = 1.5)
    rect(prior[1],0,prior[2],prior.height,col=greycol[1])
    polygon(c(pdf.cmsy$x,rev(pdf.cmsy$x)),c(pdf.cmsy$y,rep(0,length(pdf.cmsy$y))),col=greycol[2])
    
    # bk2
    post.cmsy = rem.btv.all[,which(int.yr==yr)]
    rpr = intbio
    pdf.cmsy = stats::density(post.cmsy,adjust=2)  
    prior = rpr   
	prior.height<-1/(intbio[2]-intbio[1])	# modification by GP 03/12/2019
    plot(pdf.cmsy,type="l",ylim=range(pdf.cmsy$y),xlim=range(c(pdf.cmsy$x,0.3*rpr,min(1.7*rpr[2],1.1),max(pdf.cmsy$x,rpr)*1.2)),yaxt="n",xlab=paste0("B/k ", int.yr),ylab="",xaxs="i",yaxs="i",main="",bty="l",cex.lab = 1.55, cex.axis = 1.5)
    rect(prior[1],0,prior[2],prior.height,col=greycol[1])
    polygon(c(pdf.cmsy$x,rev(pdf.cmsy$x)),c(pdf.cmsy$y,rep(0,length(pdf.cmsy$y))),col=greycol[2])
    
    # bk3
    post.cmsy = rem.btv.all[,length(yr)]
    rpr = endbio
    pdf.cmsy = stats::density(post.cmsy,adjust=2)  
    prior = rpr   
	prior.height<-1/(endbio[2]-endbio[1])	# modification by GP 03/12/2019
    plot(pdf.cmsy,type="l",ylim=range(pdf.cmsy$y),xlim=range(c(pdf.cmsy$x,0.3*rpr,min(1.7*rpr[2],1.1),max(pdf.cmsy$x,rpr)*1.2)),yaxt="n",xlab=paste0("B/k ",max(yr)),ylab="",xaxs="i",yaxs="i",main="",bty="l",cex.lab = 1.55, cex.axis = 1.5)
    rect(prior[1],0,prior[2],prior.height,col=greycol[1])
    polygon(c(pdf.cmsy$x,rev(pdf.cmsy$x)),c(pdf.cmsy$y,rep(0,length(pdf.cmsy$y))),col=greycol[2])
    mtext(paste("Density"), side=2, outer=T, at=0.5,line=-2,cex=1.5)
    
    #save analytic chart to JPEG file
    if (save.plots==TRUE) 
    {
      jpgfile<-paste(stock,"_PP_CMSY.jpg",sep="")
	  if (retrosp.step>0) jpgfile<-gsub(".jpg", paste0("_retrostep_",retrosp.step,".jpg"), jpgfile) #modification added to save all steps in retrospective analysis
      dev.copy(jpeg,jpgfile,
               width = 1024, 
               height = 768, 
               units = "px", 
               pointsize = 18,
               quality = 95,
               res=80,
               antialias="cleartype")
      dev.off()
    } 
    
    # plot PP diagnostics for BSM if available
    if(FullSchaefer==T & force.cmsy==F){ # BSM PLOT
    # open window for plot of four panels
    if(grepl("win",tolower(Sys.info()['sysname']))) {windows(17,12)}
    # make margins narrower
    par(mfrow=c(2,3),mar=c(4.5,4.5,2,0.5))
    greycol = c(grey(0.7,0.5),grey(0.3,0.5)) 
     # r
     rk <- exp(mvn(n=10000,mean.log.r=mean.log.r,sd.log.r=sd.log.r,mean.log.k=mean.log.k,sd.log.k=sd.log.k))
     rpr = rk[,1]
     post.bsm = r_raw 
     pdf.bsm = stats::density(post.bsm,adjust=2)  
     pdf.cmsy = stats::density(post.cmsy,adjust=2)  
	 prior<-prior.r     # modification by GP 04/12/2019
	 plot(prior,type="l",ylim=range(prior$y,pdf.bsm$y*1.1),xlim=range(c(pdf.bsm$x,rpr,max(pdf.bsm$x,rpr)*1.1)),yaxt="n",xlab="r",ylab="",xaxs="i",yaxs="i",main="",bty="l",cex.lab = 1.55, cex.axis = 1.5)
     polygon(c((prior$x),rev(prior$x)),c(prior$y,rep(0,length(sort(prior$y)))),col=greycol[1])
     polygon(c(pdf.bsm$x,rev(pdf.bsm$x)),c(pdf.bsm$y,rep(0,length(pdf.bsm$y))),col=greycol[2])
     PPVR.bsm = round((sd(post.bsm)/mean(post.bsm))^2/(sd(rpr)/mean(rpr))^2,2)  
     PPVM.bsm = round(mean(post.bsm)/mean(rpr),2)
     pp = c(paste("PPVR =",PPVR.bsm))
     legend('right',c("Prior","Posterior"),pch=22,pt.cex=1.5,pt.bg = greycol,bty="n",cex=1.5)
     legend("topright",pp,cex=1.4,bty="n")  
     # k
     rpr = rk[,2]
     post.bsm = k_raw 
     pdf.bsm = stats::density(post.bsm,adjust=2)  
	 prior<-prior.k     # modification by GP 04/12/2019
     plot(prior,type="l",ylim=range(prior$y,pdf.bsm$y*1.1),xlim=range(c(pdf.bsm$x,rpr,max(pdf.bsm$x,rpr)*1.1)),yaxt="n",xlab="k (1000 tonnes)",ylab="",xaxs="i",yaxs="i",main="",bty="l",cex.lab = 1.55, cex.axis = 1.5)
     polygon(c((prior$x),rev(prior$x)),c(prior$y,rep(0,length(sort(prior$y)))),col=greycol[1])
     polygon(c(pdf.bsm$x,rev(pdf.bsm$x)),c(pdf.bsm$y,rep(0,length(pdf.bsm$y))),col=greycol[2])
     PPVR.bsm = round((sd(post.bsm)/mean(post.bsm))^2/(sd(rpr)/mean(rpr))^2,2)  
     PPVM.bsm = round(mean(post.bsm)/mean(rpr),2)
     pp = c(paste("PPVR =",PPVR.bsm))
     legend("topright",pp,cex=1.4,bty="n")  
     mtext(paste0("BSM prior & posterior distributions for ",stock),  side=3,cex=1.5)
    
     # MSY
     rpr = rk[,1]*rk[,2]/4
     post.bsm = k_raw*r_raw/4
     pdf.bsm = stats::density(post.bsm,adjust=2)  
     prior = stats::density(rpr,adjust=2)  
	 prior.cmsy.calc = (prior.r$y*prior.k$y)/4 # modification by GP 04/12/2019
	 
     plot(prior.cmsy.calc,type="l",ylim=range(prior$y,pdf.bsm$y*1.1),xlim=range(c(pdf.bsm$x,rpr,max(pdf.bsm$x,rpr)*1.1)),yaxt="n",xlab="MSY (1000 tonnes/year)",ylab="",xaxs="i",yaxs="i",main="",bty="l",cex.lab = 1.55, cex.axis = 1.5)
     polygon(c((prior$x),rev(prior$x)),c(prior$y,rep(0,length(sort(prior$y)))),col=greycol[1])
     polygon(c(pdf.bsm$x,rev(pdf.bsm$x)),c(pdf.bsm$y,rep(0,length(pdf.bsm$y))),col=greycol[2])
     PPVR.bsm = round((sd(post.bsm)/mean(post.bsm))^2/(sd(rpr)/mean(rpr))^2,2)  
     PPVM.bsm = round(mean(post.bsm)/mean(rpr),2)
     pp = c(paste("PPVR =",PPVR.bsm))
     legend("topright",pp,cex=1.4,bty="n")  
    
     # bk1
     post.bsm = all.P[,1]  
     rpr = startbio
     pdf.bsm = stats::density(post.bsm,adjust=2)  
     prior = rpr  
	 prior.height<-1/(prior[2]-prior[1])	# modification by GP 03/12/2019	 
     plot(pdf.bsm,type="l",ylim=range(pdf.bsm$y),xlim=range(c(prior,pdf.bsm$x,0.3*rpr,min(1.7*rpr[2],1.05),post.bsm,max(pdf.bsm$x,rpr)*1.1)),yaxt="n",xlab=paste0("B/k ",yr[1]),ylab="",xaxs="i",yaxs="i",main="",bty="l",cex.lab = 1.55, cex.axis = 1.5)
     rect(prior[1],0,prior[2],prior.height,col=greycol[1])
     polygon(c(pdf.bsm$x,rev(pdf.bsm$x)),c(pdf.bsm$y,rep(0,length(pdf.bsm$y))),col=greycol[2])
     # bk2
     post.bsm = all.P[,which(int.yr==yr)]
     rpr = intbio
     pdf.bsm = stats::density(post.bsm,adjust=2)  
     prior = rpr   
	 prior.height<-1/(prior[2]-prior[1])	# modification by GP 03/12/2019	 
     plot(pdf.bsm,type="l",ylim=range(pdf.bsm$y),xlim=range(c(prior,pdf.bsm$x,0.3*rpr,min(1.7*rpr[2],1.05),post.bsm,max(pdf.bsm$x,rpr)*1.1)),yaxt="n",xlab=paste0("B/k ",int.yr),ylab="",xaxs="i",yaxs="i",main="",bty="l",cex.lab = 1.55, cex.axis = 1.5)
     if(nbk>1){
      rect(prior[1],0,prior[2],prior.height,col=greycol[1])
     }else{
      abline(v=prior[1],lty=2)  
      abline(v=prior[2],lty=2)
     }
     polygon(c(pdf.bsm$x,rev(pdf.bsm$x)),c(pdf.bsm$y,rep(0,length(pdf.bsm$y))),col=greycol[2])
    
     # bk3
     post.bsm = all.P[,length(yr)]  
     pdf.bsm = stats::density(post.bsm,adjust=2)  
     rpr = endbio
     prior = rpr   
	 prior.height<-1/(prior[2]-prior[1])	# modification by GP 03/12/2019	 
     plot(pdf.bsm,type="l",ylim=range(pdf.bsm$y),xlim=range(c(prior,pdf.bsm$x,0.3*rpr,min(1.7*rpr[2],1.05),
        post.bsm,max(pdf.bsm$x,rpr)*1.1)),yaxt="n",xlab=paste0("B/k ",max(yr)),ylab="",xaxs="i",yaxs="i",main="",bty="l",cex.lab = 1.55, cex.axis = 1.5)
     if(nbk>2){
      rect(prior[1],0,prior[2],prior.height,col=greycol[1])
     }else{
      abline(v=prior[1],lty=2)  
      abline(v=prior[2],lty=2)
     }
     polygon(c(pdf.bsm$x,rev(pdf.bsm$x)),c(pdf.bsm$y,rep(0,length(pdf.bsm$y))),col=greycol[2])
     mtext(paste("Density"), side=2, outer=T, at=0.5,line=-2,cex=1.5)
    
     #save analytic chart to JPEG file
     if (save.plots==TRUE) 
     {
       jpgfile<-paste(stock,"_PP_BSM.jpg",sep="")
	   if (retrosp.step>0) jpgfile<-gsub(".jpg", paste0("_retrostep_",retrosp.step,".jpg"), jpgfile) #modification added to save all steps in retrospective analysis
       dev.copy(jpeg,jpgfile,
                width = 1024, 
                height = 768, 
                units = "px", 
                pointsize = 18,
                quality = 95,
                res=80,
                antialias="cleartype")
       dev.off()
     } 
    } # end of BSM plot
  } # End of posterior/prior plot
 
  #----------------------------------------------------------
  #><> Optional BSM diagnostic plot
  #---------------------------------------------------------
  if(BSMfits.plot==T & FullSchaefer==T & force.cmsy==F){
    #---------------------------------------------
    # open window for plot of four panels
    if(grepl("win",tolower(Sys.info()['sysname']))) {windows(9,6)}
    # make margins narrower
    par(mfrow=c(2,2),mar=c(3.1,4.1,2.1,2.1),cex=1)
    cord.x <- c(yr,rev(yr))
    # Observed vs Predicted Catch
    cord.y<-c(lcl.ct.jags,rev(ucl.ct.jags))
    plot(yr,ct,type="n",ylim=c(0,max(predC,na.rm=T)),lty=1,lwd=1.3,xlab="Year",
         ylab=paste0("Catch (1000 tonnes)"),main=paste("Catch fit",stock),bty="l")
    polygon(cord.x,cord.y,col="gray",border=0,lty=1)
    lines(yr,ct.jags,lwd=2,col=1)
    points(yr,(ct),pch=21,bg="white",cex=1.)
    legend("topright",c("Observed","Predicted","95%CIs"),pch=c(21,-1,22),pt.cex = c(1,1,1.5),
           pt.bg=c("white",-1,"grey"),lwd=c(-1,2,-1),col=c(1,1,"grey"),bty="n",y.intersp = 0.9)
    
    # Observed vs Predicted CPUE
    cord.y<-c(lcl.cpue.jags,rev(ucl.cpue.jags))
    plot(yr,bt,type="n",ylim=c(0,max(c(pred.cpue,bt),na.rm=T)),lty=1,lwd=1.3,xlab="Year",ylab=paste0("cpue"),
         main="cpue fit",bty="l")
    polygon(cord.x,cord.y,col="gray",border=0,lty=1)
    lines(yr,cpue.jags,lwd=2,col=1)
    points(yr,(bt),pch=21,bg="white",cex=1.)
    legend("topright",c("Observed","Predicted","95%CIs"),pch=c(21,-1,22),pt.cex = c(1,1,1.5),pt.bg=c("white",-1,"grey"),lwd=c(-1,2,-1),col=c(1,1,"grey"),bty="n",y.intersp = 0.9)
    
    # Process error log-biomass
    cord.y<-c(lcl.pe.jags,rev(ucl.pe.jags))
    plot(yr,rep(0,length(yr)),type="n",ylim=c(-max(c(abs(pred.pe),0.2),na.rm=T),max(c(abs(pred.pe),0.2),na.rm=T)),lty=1,lwd=1.3,xlab="Year",ylab=paste0("Deviation log(B)"),main="Process variation",bty="l")
    polygon(cord.x,cord.y,col="gray",border=0,lty=1)
    abline(h=0,lty=2)
    lines(yr,pe.jags,lwd=2)
    
    
    #-------------------------------------------------
    # Function to do runs.test and 3 x sigma limits  
    #------------------------------------------------
    runs.sig3 <- function(x,type="resid") {
      if(type=="resid"){mu = 0}else{mu = mean(x, na.rm = TRUE)} 
      # Average moving range
      mr  <- abs(diff(x - mu))
      amr <- mean(mr, na.rm = TRUE)
      # Upper limit for moving ranges
      ulmr <- 3.267 * amr
      # Remove moving ranges greater than ulmr and recalculate amr, Nelson 1982
      mr  <- mr[mr < ulmr]
      amr <- mean(mr, na.rm = TRUE)
      # Calculate standard deviation, Montgomery, 6.33
      stdev <- amr / 1.128
      # Calculate control limits
      lcl <- mu - 3 * stdev
      ucl <- mu + 3 * stdev
      if(nlevels(factor(sign(x)))>1){ 
        runstest = snpar::runs.test(resid) 
        pvalue = round(runstest$p.value,3)} else {
        pvalue = 0.001  
      }
      
      return(list(sig3lim=c(lcl,ucl),p.runs= pvalue))
    }
    
    # get residuals 
    resid = (log(bt)-log(cpue.jags))[is.na(bt)==F]  
    res.yr = yr[is.na(bt)==F]
    runstest = runs.sig3(resid)
    
    # CPUE Residuals with runs test
    plot(yr,rep(0,length(yr)),type="n",ylim=c(min(-0.25,runstest$sig3lim[1]*1.1),max(0.25,runstest$sig3lim[2]*1.1)),lty=1,lwd=1.3,xlab="Year",ylab=expression(log(cpue[obs])-log(cpue[pred])),main="Residual diagnostics",bty="l")
    abline(h=0,lty=2)
    RMSE = sqrt(mean(resid^2)) # Residual mean sqrt error
    if(RMSE>0.1){lims = runstest$sig3lim} else {lims=c(-1,1)}
    cols = c(rgb(1,0,0,0.5),rgb(0,1,0,0.5))[ifelse(runstest$p.runs<0.05,1,2)]
    if(RMSE>=0.1) rect(min(yr),lims[1],max(yr),lims[2],col=cols,border=cols) # only show runs if RMSE >= 0.1
    for(i in 1:length(resid)){
      lines(c(res.yr[i],res.yr[i]),c(0,resid[i]))  
    }
    points(res.yr,resid,pch=21,bg=ifelse(resid < lims[1] | resid > lims[2],2,"white"),cex=1)
    
    # save management chart to JPEG file
    if (save.plots==TRUE & FullSchaefer == T & BSMfits.plot==TRUE) 
    {
      jpgfile<-paste(stock,"_bsmfits.jpg",sep="")
	  if (retrosp.step>0) jpgfile<-gsub(".jpg", paste0("_retrostep_",retrosp.step,".jpg"), jpgfile) #modification added to save all steps in retrospective analysis
      dev.copy(jpeg,jpgfile,
               width = 1024, 
               height = 768, 
               units = "px", 
               pointsize = 18,
               quality = 95,
               res=80,
               antialias="cleartype")
      dev.off()
    }
  }
  
  
  #-------------------------------------
  # HW Produce optional kobe plot 
  #-------------------------------------
  
  if(kobe.plot==T){
    # open window for plot of four panels
    if(grepl("win",tolower(Sys.info()['sysname']))) {windows(7,7)}
    par(mfrow=c(1,1))  
    # make margins narrower
    par(mar=c(5.1,5.1,2.1,2.1))
    
    if(FullSchaefer==T & force.cmsy==F) {x.F_Fmsy = all.F_Fmsy[,nyr]
    y.b_bmsy = all.b_bmsy[,nyr]} else {
      log.rk = cbind(rem.log.r,rem.log.k)
      rem.log.btv.lastyr = log(mdat.all[rem,nyr])
      log.bbmsy = rem.log.btv.lastyr+log(2)
      log.ffmsy = (log(ct.raw[nyr])-(rem.log.btv.lastyr+rem.log.k))-(rem.log.r-log(2))
      # get mean after all the CMSY subsetting (can't match with biomass sbmsetting)
      mu.kobe = log(c(F.Fmsy.last,B.Bmsy.last))
      # Get covariance of the 2 vectors
      cov.kobe = cov(cbind(log.ffmsy,log.bbmsy)) 
      # Generate 10000 new random deviates from a MVN
      log.kobe.mvn = rmvnorm(10000 ,mean = mu.kobe,sigma = cov.kobe)
      kobe.mvn = exp(log.kobe.mvn)
      # Generate 10000 new random deviates from a MVN
      x.F_Fmsy =exp(log.kobe.mvn[,1])
      y.b_bmsy =exp(log.kobe.mvn[,2])
    }
    
    kernelF <- ci2d(y.b_bmsy,x.F_Fmsy,nbins=151,factor=2.2,ci.levels=c(0.50,0.80,0.75,0.90,0.95),show="none",col=1,xlab= ifelse(harvest.label=="Fmsy",expression(paste(F/F[MSY])),expression(paste(H/H[MSY]))),ylab=expression(paste(B/B[MSY])))
    
    max.y1   <- max(c(2, max(kernelF$contours$"0.95"$x,F.Fmsy),na.rm =T))
    max.y    <- ifelse(max.x1 > 5,min(max(5,F.Fmsy*2),8),max.x1)
    max.x    <- max(max(2,quantile(y.b_bmsy,0.96)))
    
    # -------------------------------------
    ## KOBE plot building
    # -------------------------------------
    #Create plot
    plot(1000,1000,type="b", xlim=c(0,max.x), ylim=c(0,max.y),lty=3,xlab="",ylab=expression(F/F[MSY]), bty="l",  cex.main = 2, cex.lab = 1.35, cex.axis = 1.35,xaxs = "i",yaxs="i")
    mtext(expression(B/B[MSY]),side=1, line=3, cex=1.3)
    c1 <- c(-1,100)
    c2 <- c(1,1)
    
    # extract interval information from ci2d object
    # and fill areas using the polygon function
    zb2 = c(0,1)
    zf2  = c(1,100)
    zb1 = c(1,100)
    zf1  = c(0,1)
    polygon(c(zb1,rev(zb1)),c(0,0,1,1),col="green",border=0)
    polygon(c(zb2,rev(zb2)),c(0,0,1,1),col="yellow",border=0)
    polygon(c(1,100,100,1),c(1,1,100,100),col="orange",border=0)
    polygon(c(0,1,1,0),c(1,1,100,100),col="red",border=0)
    
    polygon(kernelF$contours$"0.95",lty=2,border=NA,col="cornsilk4")
    polygon(kernelF$contours$"0.8",border=NA,lty=2,col="grey")
    polygon(kernelF$contours$"0.5",border=NA,lty=2,col="cornsilk2")
    points(B.Bmsy,F.Fmsy,pch=16,cex=1)
    lines(c1,c2,lty=3,lwd=0.7)
    lines(c2,c1,lty=3,lwd=0.7)
    lines(B.Bmsy,F.Fmsy, lty=1,lwd=1.)
    points(B.Bmsy[1],F.Fmsy[1],col=1,pch=22,bg="white",cex=1.5)
    points(B.Bmsy[which(yr==int.yr)],F.Fmsy[which(yr==int.yr)],col=1,pch=21,bg="white",cex=1.5)
    points(B.Bmsy[nyr],F.Fmsy[nyr],col=1,pch=24,bg="white",cex=1.5)
    # Get Propability
    Pr.green = sum(ifelse(y.b_bmsy>1 & x.F_Fmsy<1,1,0))/length(y.b_bmsy)*100
    Pr.red = sum(ifelse(y.b_bmsy<1 & x.F_Fmsy>1,1,0))/length(y.b_bmsy)*100
    Pr.yellow = sum(ifelse(y.b_bmsy<1 & x.F_Fmsy<1,1,0))/length(y.b_bmsy)*100
    Pr.orange = sum(ifelse(y.b_bmsy>1 & x.F_Fmsy>1,1,0))/length(y.b_bmsy)*100
    
    sel.years = c(yr[sel.yr])
    
    legend('topright', 
           c(paste(start.yr),paste(int.yr),paste(end.yr),"50% C.I.","80% C.I.","95% C.I.",paste0(round(c(Pr.red,Pr.yellow,Pr.orange,Pr.green),1),"%")), 
           lty=c(1,1,1,rep(-1,8)),pch=c(22,21,24,rep(22,8)),pt.bg=c(rep("white",3),"cornsilk2","grey","cornsilk4","red","yellow","orange","green"), 
           col=1,lwd=1.1,cex=1.1,pt.cex=c(rep(1.3,3),rep(1.7,3),rep(2.2,4)),bty="n",y.intersp = 1.)  
    
    if (save.plots==TRUE & kobe.plot==TRUE) 
    {
      jpgfile<-paste(stock,"_KOBE.jpg",sep="")
	  if (retrosp.step>0) jpgfile<-gsub(".jpg", paste0("_retrostep_",retrosp.step,".jpg"), jpgfile) #modification added to save all steps in retrospective analysis
      dev.copy(jpeg,jpgfile,
               width = 1024*0.7, 
               height = 1024*0.7, 
               units = "px", 
               pointsize = 18,
               quality = 95,
               res=80,
               antialias="cleartype")
      dev.off()
    }  
  }
  
  #HW Kobe plot end
  
  # -------------------------------------
  ## Write results into csv outfile
  # -------------------------------------
  if(write.output == TRUE && retrosp.step==0) { #account for retrospective analysis - write only the last result
    
    # fill catches from 1970 to 2020
    # if leading catches are missing, set them to zero; if trailing catches are missing, set them to NA
    ct.out     <- vector()
    F.Fmsy.out <- vector()
    bt.out     <- vector()
    
    j <- 1
    for(i in 1950 : 2020) {
      if(yr[1]>i) {
        ct.out[j]     <-0 
        F.Fmsy.out[j] <-0
        bt.out[j]     <-2*Bmsy
      } else {
        if(i>yr[length(yr)]) {
          ct.out[j]     <-NA 
          F.Fmsy.out[j] <-NA
          bt.out[j]     <-NA } else {
            ct.out[j]     <- ct.raw[yr==i] 
            F.Fmsy.out[j] <- F.Fmsy[yr==i]
            bt.out[j]     <- B[yr==i]}
      }
      j=j+1
    }
    
    # write data into csv file
    output = data.frame(as.character(cinfo$Group[cinfo$Stock==stock]),
                        as.character(cinfo$Region[cinfo$Stock==stock]),
                        as.character(cinfo$Subregion[cinfo$Stock==stock]),
                        as.character(cinfo$Name[cinfo$Stock==stock]),
                        cinfo$ScientificName[cinfo$Stock==stock], 
                        stock, start.yr, end.yr, btype,
                        max(ct.raw),ct.raw[nyr],
                        ifelse(FullSchaefer==T,MSY.jags,NA), # full Schaefer
                        ifelse(FullSchaefer==T,lcl.MSY.jags,NA),
                        ifelse(FullSchaefer==T,ucl.MSY.jags,NA),
                        ifelse(FullSchaefer==T,r.jags,NA), 
                        ifelse(FullSchaefer==T,lcl.r.jags,NA),
                        ifelse(FullSchaefer==T,ucl.r.jags,NA),
                        ifelse(FullSchaefer==T,log.r.var,NA),
                        ifelse(FullSchaefer==T,k.jags,NA),                                          
                        ifelse(FullSchaefer==T,lcl.k.jags,NA),                    
                        ifelse(FullSchaefer==T,ucl.k.jags,NA),
                        ifelse(FullSchaefer==T,log.k.var,NA),
                        ifelse(FullSchaefer==T,log.rk.cor,NA),
                        ifelse(FullSchaefer==T,log.rk.cov,NA),
                        ifelse(FullSchaefer==T, mean.q,NA),
                        ifelse(FullSchaefer==T,lcl.q,NA),
                        ifelse(FullSchaefer==T,ucl.q,NA),
                        ifelse(FullSchaefer==T,quant.P[2,][nyr],NA), # last B/k JAGS
                        ifelse(FullSchaefer==T,quant.P[1,][nyr],NA),
                        ifelse(FullSchaefer==T,quant.P[3,][nyr],NA), 
                        ifelse(FullSchaefer==T,F_Fmsy.jags[nyr],NA), # last F/Fmsy JAGS
                        r.est, lcl.r.est, ucl.r.est, # CMSY r
                        k.est, lcl.k.est, ucl.k.est, # CMSY k     
                        MSY.est, lcl.MSY.est, ucl.MSY.est, # CMSY MSY
                        median.btv.lastyr, lcl.median.btv.lastyr,ucl.median.btv.lastyr, # CMSY B/k in last year with catch data
                        (F.CMSY/Fmsy.CMSY)[nyr],
                        Fmsy,lcl.Fmsy,ucl.Fmsy,Fmsy.last,lcl.Fmsy.last,ucl.Fmsy.last,
                        MSY,lcl.MSY,ucl.MSY,Bmsy,lcl.Bmsy,ucl.Bmsy,
                        B.last, lcl.B.last, ucl.B.last, B.Bmsy.last, lcl.B.Bmsy.last, ucl.B.Bmsy.last,
                        F.last, lcl.F.last, ucl.F.last, F.Fmsy.last, lcl.F.Fmsy.last, ucl.F.Fmsy.last,
                        ifelse(is.na(sel.yr)==F,B.sel,NA),
                        ifelse(is.na(sel.yr)==F,B.Bmsy.sel,NA),
                        ifelse(is.na(sel.yr)==F,F.sel,NA),
                        ifelse(is.na(sel.yr)==F,F.Fmsy.sel,NA),
                        ct.out[1],ct.out[2],ct.out[3],ct.out[4],ct.out[5],ct.out[6],ct.out[7],ct.out[8],ct.out[9],ct.out[10],          # 1950-1959
                        ct.out[11],ct.out[12],ct.out[13],ct.out[14],ct.out[15],ct.out[16],ct.out[17],ct.out[18],ct.out[19],ct.out[20], # 1960-1969
                        ct.out[21],ct.out[22],ct.out[23],ct.out[24],ct.out[25],ct.out[26],ct.out[27],ct.out[28],ct.out[29],ct.out[30], # 1970-1979
                        ct.out[31],ct.out[32],ct.out[33],ct.out[34],ct.out[35],ct.out[36],ct.out[37],ct.out[38],ct.out[39],ct.out[40], # 1980-1989
                        ct.out[41],ct.out[42],ct.out[43],ct.out[44],ct.out[45],ct.out[46],ct.out[47],ct.out[48],ct.out[49],ct.out[50], # 1990-1999
                        ct.out[51],ct.out[52],ct.out[53],ct.out[54],ct.out[55],ct.out[56],ct.out[57],ct.out[58],ct.out[59],ct.out[60], # 2000-2009
                        ct.out[61],ct.out[62],ct.out[63],ct.out[64],ct.out[65],ct.out[66],ct.out[67],ct.out[68],ct.out[69],ct.out[70],ct.out[71], # 2010-2020
                        F.Fmsy.out[1],F.Fmsy.out[2],F.Fmsy.out[3],F.Fmsy.out[4],F.Fmsy.out[5],F.Fmsy.out[6],F.Fmsy.out[7],F.Fmsy.out[8],F.Fmsy.out[9],F.Fmsy.out[10], # 1950-1959
                        F.Fmsy.out[11],F.Fmsy.out[12],F.Fmsy.out[13],F.Fmsy.out[14],F.Fmsy.out[15],F.Fmsy.out[16],F.Fmsy.out[17],F.Fmsy.out[18],F.Fmsy.out[19],F.Fmsy.out[20], # 1960-1969
                        F.Fmsy.out[21],F.Fmsy.out[22],F.Fmsy.out[23],F.Fmsy.out[24],F.Fmsy.out[25],F.Fmsy.out[26],F.Fmsy.out[27],F.Fmsy.out[28],F.Fmsy.out[29],F.Fmsy.out[30], # 1970-1979
                        F.Fmsy.out[31],F.Fmsy.out[32],F.Fmsy.out[33],F.Fmsy.out[34],F.Fmsy.out[35],F.Fmsy.out[36],F.Fmsy.out[37],F.Fmsy.out[38],F.Fmsy.out[39],F.Fmsy.out[40], # 1980-1989
                        F.Fmsy.out[41],F.Fmsy.out[42],F.Fmsy.out[43],F.Fmsy.out[44],F.Fmsy.out[45],F.Fmsy.out[46],F.Fmsy.out[47],F.Fmsy.out[48],F.Fmsy.out[49],F.Fmsy.out[50], # 1990-1999
                        F.Fmsy.out[51],F.Fmsy.out[52],F.Fmsy.out[53],F.Fmsy.out[54],F.Fmsy.out[55],F.Fmsy.out[56],F.Fmsy.out[57],F.Fmsy.out[58],F.Fmsy.out[59],F.Fmsy.out[60], # 2000-2009
                        F.Fmsy.out[61],F.Fmsy.out[62],F.Fmsy.out[63],F.Fmsy.out[64],F.Fmsy.out[65],F.Fmsy.out[66],F.Fmsy.out[67],F.Fmsy.out[68],F.Fmsy.out[69],F.Fmsy.out[70],F.Fmsy.out[71], # 2010-2020
                        bt.out[1],bt.out[2],bt.out[3],bt.out[4],bt.out[5],bt.out[6],bt.out[7],bt.out[8],bt.out[9],bt.out[10],           # 1950-1959
                        bt.out[11],bt.out[12],bt.out[13],bt.out[14],bt.out[15],bt.out[16],bt.out[17],bt.out[18],bt.out[19],bt.out[20],  # 1960-1969
                        bt.out[21],bt.out[22],bt.out[23],bt.out[24],bt.out[25],bt.out[26],bt.out[27],bt.out[28],bt.out[29],bt.out[30],  # 1970-1979
                        bt.out[31],bt.out[32],bt.out[33],bt.out[34],bt.out[35],bt.out[36],bt.out[37],bt.out[38],bt.out[39],bt.out[40],  # 1980-1989 
                        bt.out[41],bt.out[42],bt.out[43],bt.out[44],bt.out[45],bt.out[46],bt.out[47],bt.out[48],bt.out[49],bt.out[50],  # 1990-1999
                        bt.out[51],bt.out[52],bt.out[53],bt.out[54],bt.out[55],bt.out[56],bt.out[57],bt.out[58],bt.out[59],bt.out[60],  # 2000-2009
                        bt.out[61],bt.out[62],bt.out[63],bt.out[64],bt.out[65],bt.out[66],bt.out[67],bt.out[68],bt.out[69],bt.out[70],bt.out[71]) # 2010-2020
    
    write.table(output, file=outfile, append = T, sep = ",", 
                dec = ".", row.names = FALSE, col.names = FALSE)
  }  
  
  #----------------------------------------------------------------------------------
  # The code below creates a report in PDF format if write.pdf is TRUE ----
  #----------------------------------------------------------------------------------
  ## To generate reports in PDF format, install a LaTeX program. For Windows, you can use https://miktex.org/howto/install-miktex (restart after installation)
  ## Set write.pdf to 'TRUE' if you want pdf output.
  
  # Using MarkdownReports, this creates a markdown file for each stock then using rmarkdown to render each markdown file into a pdf file.
  if(write.pdf == TRUE) {
    library(knitr)
    library(tinytex)

    docTemplate <- "\\documentclass[12pt,a4paper]{article}
\\setlength\\parindent{0pt}
  \\usepackage{geometry}
  \\geometry{margin=0.5in}
  \\begin{document}
  
  \\section*{#TITLE#}
  
  
  #INTRO#
  
  \\begin{figure}[ht]
  \\centering
  \\includegraphics[width=1.00\\textwidth]{#IMAGE1#}
  \\end{figure}
  
  #MANAGEMENT#
  
  \\pagebreak
  
  \\begin{figure}[ht]
  \\centering
  \\includegraphics[width=1.00\\textwidth]{#IMAGE2#}
  \\end{figure}
  
  #ANALYSIS#
  
  \\end{document}"
    
    title = cinfo$Name[cinfo$Stock==stock]
    
    intro = (paste("Species: \\\\emph{",cinfo$ScientificName[cinfo$Stock==stock],"}, Stock code: ", stock,".", sep=""))
    intro = (paste(intro,"\n\n","Region: ",cinfo$Region[cinfo$Stock==stock],".", sep=""))
    intro = (paste(intro,"\n\n","Marine Ecoregion: ",cinfo$Subregion[cinfo$Stock==stock],".", sep="" ))
    intro = (paste(intro,"\n\n","Reconstructed catch data used from years ", min(yr)," - ", max(yr),sep=""))
    intro = (paste(intro,"\n\n","For figure captions and method see http://www.seaaroundus.org/cmsy-method"))
    
    
    docTemplate<-gsub("#TITLE#", title, docTemplate)
    docTemplate<-gsub("#INTRO#", intro, docTemplate)
    
    
    management_text<-paste("\\\\textbf{Results for management (based on",ifelse(FullSchaefer==F | force.cmsy==T,"CMSY","BSM"),"analysis)}\\\\\\\\")
    management_text<-(paste(management_text,"\n\n","Fmsy = ",format(Fmsy, digits =3),", 95% CL = ",format(lcl.Fmsy, digits =3)," - ",format(ucl.Fmsy, digits =3)," (if B > 1/2 Bmsy then Fmsy = 0.5 r)", sep=""))
    management_text<-(paste(management_text,"\n\n","Fmsy = ",format(Fmsy.last, digits =3),", 95% CL = ",format(lcl.Fmsy.last, digits =3)," - ",format(ucl.Fmsy.last, digits =3)," (r and Fmsy are linearly reduced if B < 1/2 Bmsy)",sep=""))
    management_text<-(paste(management_text,"\n\n","MSY = ",format(MSY, digits =3),",  95% CL = ",format(lcl.MSY, digits =3)," - ",format(ucl.MSY, digits =3),'; Bmsy = ',format(Bmsy, digits =3),",  95% CL = ",format(lcl.Bmsy, digits =3)," - ",format(ucl.Bmsy, digits =3)," (1000 tonnes)",sep="")) 
    management_text<-(paste(management_text,"\n\n","Biomass in last year = ",format(B.last, digits =3),", 95% CL = ", format(lcl.B.last, digits =3), " - ",format(ucl.B.last, digits =3)," (1000 tonnes)",sep=""))
    management_text<-(paste(management_text,"\n\n","B/Bmsy in last year = " ,format(B.Bmsy.last, digits =3),", 95% CL = ", format(lcl.B.Bmsy.last, digits =3), " - ",format(ucl.B.Bmsy.last, digits =3),sep=""))
    management_text<-(paste(management_text,"\n\n","Fishing mortality in last year = ",format(F.last, digits =3),", 95% CL =", format(lcl.F.last, digits =3), " - ",format(ucl.F.last, digits =3),sep=""))
    management_text<-(paste(management_text,"\n\n","F/Fmsy  = ",format(F.Fmsy.last, digits =3),", 95% CL = ", format(lcl.F.Fmsy.last, digits =3), " - ",format(ucl.F.Fmsy.last, digits =3),sep=""))
    management_text<-(paste(management_text,"\n\n","Comment:", comment, ""))
    docTemplate<-gsub("#MANAGEMENT#", management_text, docTemplate)
    
    analysis_text<-(paste("\\\\textbf{Results of CMSY analysis with altogether ",n.viable.b, " viable trajectories for ", n.viable.pt," r-k pairs}\\\\\\\\",sep=""))
    analysis_text<-(paste(analysis_text,"\n\n","r = ", format(r.est, digits =3),", 95% CL = ", format(lcl.r.est, digits =3), " - ", format(ucl.r.est, digits =3),"; k = ", format(k.est, digits =3),", 95% CL = ", format(lcl.k.est, digits =3), " - ", format(ucl.k.est, digits =3)," (1000 tonnes)",sep=""))
    analysis_text<-(paste(analysis_text,"\n\n","MSY = ", format(MSY.est, digits =3),", 95% CL = ", format(lcl.MSY.est, digits =3), " - ", format(ucl.MSY.est, digits =3)," (1000 tonnes/year)",sep=""))
    analysis_text<-(paste(analysis_text,"\n\n","Relative biomass last year = ", format(median.btv.lastyr, digits =3), " k, 95% CL = ", format(lcl.median.btv.lastyr, digits =3), " - ", format(ucl.median.btv.lastyr, digits =3),sep=""))
    analysis_text<-(paste(analysis_text,"\n\n","Exploitation F/(r/2) in last year = ", format((F.CMSY/Fmsy.CMSY)[length(median.btv)-1], digits =3),sep=""))
    
    if(FullSchaefer==T) {
      analysis_text <- paste(analysis_text,"\\\\\\\\")
      analysis_text<-(paste(analysis_text,"\n\n", "\\\\textbf{Results from Bayesian Schaefer model using catch and ",btype,"}\\\\\\\\",sep=""))
      analysis_text<-(paste(analysis_text,"\n\n","r = ", format(r.jags, digits =3),", 95% CL = ", format(lcl.r.jags, digits =3), " - ", format(ucl.r.jags, digits =3),"; k = ", format(k.jags, digits =3),", 95% CL = ", format(lcl.k.jags, digits =3), " - ", format(ucl.k.jags, digits =3),sep=""))
      analysis_text<-(paste(analysis_text,"\n\n","MSY = ", format(MSY.jags, digits =3),", 95% CL = ", format(lcl.MSY.jags, digits =3), " - ", format(ucl.MSY.jags, digits =3)," (1000 tonnes/year)",sep=""))
      analysis_text<-(paste(analysis_text,"\n\n","Relative biomass in last year = ", format(quant.P[2,][nyr], digits =3), " k, 95% CL = ",format(quant.P[1,][nyr], digits =3)," - ", format(quant.P[3,][nyr], digits =3),sep=""))
      analysis_text<-(paste(analysis_text,"\n\n","Exploitation F/(r/2) in last year = ", format((ct.raw[nyr]/(quant.P[2,][nyr]*k.jags))/(r.jags/2), digits =3),sep=""))
    }
      analysis_text<-(paste(analysis_text,"\n\n","q = ", format(mean.q, digits =3),", 95% CL = ", format(lcl.q, digits =3), " - ", format(ucl.q, digits =3),sep=""))
    
    
    if(FullSchaefer==T) {
      analysis_text<-(paste(analysis_text,"\n\n","Prior range of q = ",format(q.prior[1], digits =3)," - ",format(q.prior[2], digits =3),sep=""))
    }
    # show stock status and exploitation for optional selected year
    if(is.na(sel.yr)==F) {
      analysis_text<-(paste(analysis_text,"\n\n","Stock status and exploitation in ",sel.yr,sep=""))
      analysis_text<-(paste(analysis_text,"\n\n","Biomass = ",format(B.sel, digits =3), ", B/Bmsy = ",format(B.Bmsy.sel, digits =3),", fishing mortality F = ",format(F.sel, digits =3),", F/Fmsy = ",format(F.Fmsy.sel, digits =3),sep=""))
    }
    
    if(btype !="None" & length(bt[is.na(bt)==F])<nab) {
      analysis_text<-(paste(analysis_text,"\n\n","Less than ",nab," years with abundance data available, shown on second axis",sep="")) }
    
    
    analysis_text<-(paste(analysis_text,"\n\n","Relative abundance data type = ", format(btype, digits =3),sep=""))
    analysis_text<-(paste(analysis_text,"\n\n","Prior initial relative biomass = ", startbio[1], " - ", startbio[2],ifelse(is.na(stb.low)==T," default"," expert"),sep=""))
    analysis_text<-(paste(analysis_text,"\n\n","Prior intermediate relative biomass = ", intbio[1], " - ", intbio[2], " in year ", int.yr,ifelse(is.na(intb.low)==T," default"," expert"),sep=""))
    analysis_text<-(paste(analysis_text,"\n\n","Prior final relative biomass = ", endbio[1], " - ", endbio[2],ifelse(is.na(endb.low)==T,", default"," expert"),sep=""))
    analysis_text<-(paste(analysis_text,"\n\n","Prior range for r = ", format(prior.r[1],digits=2), " - ", format(prior.r[2],digits=2),ifelse(is.na(r.low)==T," default"," expert"),", prior range for k = " , format(prior.k[1], digits =3), " - ", format(prior.k[2], digits =3)," (1000 tonnes) default",sep=""))
    analysis_text<-(paste(analysis_text,"\n\n","Source for relative biomass: \n\n",source,"",sep=""))
    
    docTemplate<-gsub("#ANALYSIS#", analysis_text, docTemplate)
    
    docTemplate<-gsub("_", "\\\\_", docTemplate)
    docTemplate<-gsub("%", "\\\\%", docTemplate)
    
    
    analysischartfile<-paste(stock,"_AN.jpg",sep="")
    managementchartfile<-paste(stock,"_MAN.jpg",sep="")
    docTemplate<-gsub("#IMAGE1#", managementchartfile, docTemplate)
    docTemplate<-gsub("#IMAGE2#", analysischartfile, docTemplate)
    
    # unique filenames to prevent error if files exists from previous run
    documentfile<-paste(stock,substr(as.character(Sys.time()),1,10),"-",sub(":","",substr(as.character(Sys.time()),12,16)),".RnW",sep="") # concatenated hours and minutes added to file name
    cat(docTemplate,file=documentfile,append=F)
    
    knit(documentfile)
    knitr::knit2pdf(documentfile)
    
    cat("PDF document is ",gsub(".RnW",".pdf",documentfile))
    
  }  
  # end of loop to write text to file
  
  
  if(close.plots==T) graphics.off() # close on-screen graphics windows after files are saved
  
  FFmsy.retrospective[[retrosp.step+1]]<-F.Fmsy #retrospective analysis
  BBmsy.retrospective[[retrosp.step+1]]<-B.Bmsy #retrospective analysis
  years.retrospective[[retrosp.step+1]]<-yr #retrospective analysis
  
  } #retrospective analysis - end loop
	
	#retrospective analysis plots
	if (retros == T){
	  
	   if(grepl("win",tolower(Sys.info()['sysname']))) {windows(14,7)}
		par(mfrow=c(1,2), mar=c(4,5,4,5),  oma=c(2,2,2,2))  

	  allyears<-years.retrospective[[1]]
	  nyrtotal<-length(allyears)
	  legendyears<-c("All years")
	  #CHECK IF ALL YEARS HAVE BEEN COMPUTED
	  for (ll in 1:4){
	    if (ll>length(FFmsy.retrospective)){
	      FFmsy.retrospective[[ll]]<-c(0)
	      BBmsy.retrospective[[ll]]<-c(0)
	    }
	    else {
	      if(ll>1)
	        legendyears<-c(legendyears,allyears[nyrtotal-ll+1])
	    }
	  }
	  
	  #PLOT FFMSY RETROSPECTIVE ANALYSIS
	  plot(x=allyears[1:nyrtotal],y=FFmsy.retrospective[[1]], main="", ylim=c(0,max(max(FFmsy.retrospective[[1]],na.rm=T),
	                                                                                max(FFmsy.retrospective[[2]],na.rm=T),
	                                                                                max(FFmsy.retrospective[[3]],na.rm=T),
	                                                                                max(FFmsy.retrospective[[4]],na.rm=T))),
	       lwd=2, xlab="Year", ylab="F/Fmsy", type="l", bty="l",  cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5) #, xaxs="i",yaxs="i",xaxt="n",yaxt="n")
	  #PLOT ONLY THE TIME SERIES THAT ARE COMPLETE
	  if (length(FFmsy.retrospective[[2]])>1 || FFmsy.retrospective[[2]]!=0)
	    lines(x=allyears[1:(nyrtotal-1)],y=FFmsy.retrospective[[2]], type = "o", pch=15, col="red")
	  if (length(FFmsy.retrospective[[3]])>1 || FFmsy.retrospective[[3]]!=0)
	    lines(x=allyears[1:(nyrtotal-2)],y=FFmsy.retrospective[[3]], type = "o", pch=16, col="green") 
	  if (length(FFmsy.retrospective[[4]])>1 || FFmsy.retrospective[[4]]!=0)
	    lines(x=allyears[1:(nyrtotal-3)],y=FFmsy.retrospective[[4]], type = "o", pch=17, col="blue") 
	  legend("bottomleft", legend = legendyears, 
	         col=c("black","red", "green", "blue"), lty=1, pch=c(-1,15,16,17))
	  #PLOT BBMSY RETROSPECTIVE ANALYSIS
	  plot(x=allyears[1:(nyrtotal)],y=BBmsy.retrospective[[1]],main="", ylim=c(0,max(max(BBmsy.retrospective[[1]],na.rm=T),
	                                                                                 max(BBmsy.retrospective[[2]],na.rm=T),
	                                                                                 max(BBmsy.retrospective[[3]],na.rm=T),
	                                                                                 max(BBmsy.retrospective[[4]],na.rm=T))), 
	       lwd=2, xlab="Year", ylab="B/Bmsy", type="l", bty="l",cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5) #, xaxs="i",yaxs="i",xaxt="n",yaxt="n")
	  if (length(BBmsy.retrospective[[2]])>1 || BBmsy.retrospective[[2]]!=0)
	    lines(x=allyears[1:(nyrtotal-1)],y=BBmsy.retrospective[[2]], type = "o", pch=15, col="red")
	  if (length(BBmsy.retrospective[[3]])>1 || BBmsy.retrospective[[3]]!=0)
	    lines(x=allyears[1:(nyrtotal-2)],y=BBmsy.retrospective[[3]], type = "o", pch=16, col="green") 
	  if (length(BBmsy.retrospective[[4]])>1 || BBmsy.retrospective[[4]]!=0)
	    lines(x=allyears[1:(nyrtotal-3)],y=BBmsy.retrospective[[4]], type = "o", pch=17, col="blue") 
	  legend("bottomleft", legend = legendyears, 
	         col=c("black","red", "green", "blue"), lty=1, pch=c(-1,15,16,17))
	  
	  mtext(paste0("Retrospective analysis for ",stock),  outer = T , cex=1.5)
	  
	  #save analytic chart to JPEG file
    if (save.plots==TRUE) 
    {
      jpgfile<-paste(stock,"_RetrospectiveAnalysis.jpg",sep="")
      dev.copy(jpeg,jpgfile,
               width = 1024, 
               height = 576, 
               units = "px", 
               pointsize = 10,
               quality = 95,
               res=80,
               antialias="default")
      dev.off()
    } 
	
	  if(close.plots==T) graphics.off() # close on-screen graphics windows after files are saved
	} #retrospective analysis plots - end

} # end of stocks loop

#stop parallel processing clusters
stopCluster(cl)
stopImplicitCluster()


