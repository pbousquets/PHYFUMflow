if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load_gh("adamallo/rwty")
p_load(rstan, LaplacesDemon, ggplot2, cowplot, ggrepel, data.table, HDInterval)


#My functions
#######################################

#' Summary function for point estimates
#' 
#' @param x data
#' @param type character indicating if we want to use the median, mean, or mode
#' @value point estimate
sumFun=function(x,type=c("mean","median","mode")){
  type=match.arg(type)
  if(type=="mean"){
    return(mean(x))
  } else if (type=="median"){
    return(median(x))
  } else if (type=="mode") {
    return(mode(x))
  }
}

#' Method-of-moments estimator of AICM
#' Posterior simulation-based analogue of Akaike's information criterion (AIC)
#' through Markov chain Monte Carlo
#' Raftery et al. 2007
#' Translated from BEAST's logMarginalLikelihoodAICM function
#' 
#' @param loglikelihoods vector of posterior sample of loglikelihoods
#' @retuns AICM
AICM=function(loglikelihoods){
  return(2*var(loglikelihoods)-2*mean(loglikelihoods))
}

#' Detects constant parameters
#' @param theseData data.table
#' @param params list of columns to consider
#' @retuns vector of names of parameters with sd == 0
detectConstants <- function(theseData,params){
  params[t(theseData[,lapply(.SD,FUN = function(x){sd(x)==0}),.SDcols = params])]}

#' Generates a table with convergence statistics using Rstan
#'
#' @param stanTable MCMC data for one parameter, with as many columns as independent chains of same length
#' @param credMass Mass for the high posterior density interval
#' @return data.table with a row per parameter and colums for Rhat, ESSB and ESST
stanStats <- function(stanTable,credMass = 0.95) {
  stanMatrix <- as.matrix(stanTable)
  thisHDI <- hdi(stanMatrix,credMass = credMass)
  data.table("Rhat"=Rhat(stanMatrix),
             "essB"=ess_bulk(stanMatrix),
             "essT"=ess_tail(stanMatrix),
             "medianP"=median(stanMatrix), 
             "meanP"=mean(stanMatrix),
             "HDILower"=thisHDI[1],
             "HDIUpper"=thisHDI[2])
}

#' Generates a trace of tree distances from a focal tree selected at random to calculate a pseudo ESS with it (or plotting, etc.)
#'
#' @param tree.list list of trees
#' @param treedist Tree distance to be used to calculate the distance to the focal tree ("PD" or "RF")
#' @return vector of tree distances
treeDistanceTrace <- function (tree.list, treedist = "PD", iTree = NA) 
{
  if (is.na(iTree)){
    iTree <- sample(1:length(tree.list), 1)
  }
  distances <- rwty:::tree.distances(tree.list, iTree, treedist = treedist)
  return(distances)
}

#' Generates a data.table merging data from multiple MCMC chains
#' 
#' @param chains MCMC data. List of rwty.chains (RWTY's load.trees value)
#' @param burninP Proportion of MCMC samples to discard as burnin
#' @return data.table with same columns as chain$ptable + a chain factor
getDataTable=function(chains,burninP){
  returnTable=data.table(chains[[1]]$ptable,chain=1)
  if(length(chains)>1){
    for(ichain in 2:length(chains)){
      returnTable=rbind(returnTable,data.table(chains[[ichain]]$ptable,chain=ichain))
    }}
  returnTable[,`:=`(chain=as.factor(chain),burnin=ifelse(state<max(state)*burninP,T,F))]
  return(returnTable)
}

#' Regenerates multiple MCMC chains from a data.table created by getDataTable
#' 
#' @param chains MCMC data. List of rwty.chains (RWTY's load.trees value)
#' @param theData Data table with modified chain data
#' @return new chains
updateChainsFromDataTable=function(chains,theData){
  chainCols <- names(theData)
  chainCols <- chainCols[!chainCols %in% c("chain","burnin","treedistance")] 
  for(ichain in 1:length(chains)){
    chains[[ichain]]$ptable <- data.table(theData[chain==ichain,.SD,.SDcols = chainCols])
  }
  return(chains)
}

#' Merges a list of rwty.chains into one, after proper burnin
#' 
#' @param chains MCMC data. List of rwty.chains (RWTY's load.trees value)
#' @param burninP Proportion of MCMC samples to discard as burnin
#' @return rwty.chain with the chains concatenated with the burnin trees and all paramter values (ptable) removed.
mergeTreeTrace=function(chains,burninP){
  returnChain=chains[[1]]
  returnChain$ptable=NULL
  last=length(returnChain$trees)
  first=round(burninP*last)
  returnChain$trees=returnChain$trees[first:last]
  if(length(chains)>1){
  for(ichain in 2:length(chains)){
    returnChain$trees=c(returnChain$trees,chains[[ichain]]$trees[first:last])
  }}
  return(returnChain)
}

#' Generates a data.table with convergence statistics for all parameters
#' 
#' @param theData MCMC data in data.table format after getDataTable
#' @param params List of parameter names to analyze
#' @param credMass Mass for the high posterior density interval
#' @param byChain Calculates Rhat, and the two ESS parameters both using all chains and by chain (otherwise, using all of them only)
#' @return data.table with a row per parameter and colums for Rhat, ESSB and ESST
getConvergenceTableFromDataTable=function(theData,params,credMass = 0.95,byChain = T){
  returnTable <- rbindlist(lapply(1:length(params),FUN=function(iparam){
          thisParam <- params[iparam]
          stanStats(dcast(theData[,c("state",thisParam,"chain"),with=F],state~chain,value.var = thisParam)[,-1], credMass = credMass)[,`:=`(param = thisParam)]}))
  setkey(returnTable, param)
  returnTable[,`:=`(chain = "ALL")]
  if(byChain == T){
    returnTableByChain <- rbindlist(lapply(theData[,unique(chain)],FUN = function(thisChain){
      thisReturnTable <- rbindlist(lapply(1:length(params),FUN=function(iparam){
        thisParam <- params[iparam]
        stanStats(dcast(theData[chain == thisChain,c("state",thisParam,"chain"),with=F],state~chain,value.var = thisParam)[,-1], credMass = credMass)[,`:=`(param = thisParam)]}))
      thisReturnTable[,`:=`(chain = thisChain)]
      return(thisReturnTable)}))
    setkey(returnTableByChain,param)
    returnTable <- returnTable[returnTableByChain[,.(medianPCV=sd(medianP)/mean(medianP),meanPCV=sd(meanP)/mean(meanP)),by=param]]
    returnTable <- rbindlist(list(returnTable,returnTableByChain),fill = T)
  }
  return(returnTable)
}

#' Generates a table with tree convergence statistics using the pseudo ESS approach
#'
#' @param chains MCMC data. List of rwty.chains (RWTY's load.trees value)
#' @param burnin Number of MCMC samples to discard as burnin
#' @param n Number of focal trees to calculate tree distances, for which convergence parameters are calculated. Then, the mean is reported
#' @param treedist Tree distance to be used to calculate the distance to the focal tree
#' @param credMass Mass for the high posterior density interval
#' @param byChain Calculates Rhat, and the two ESS parameters both using all chains and by chain (otherwise, using all of them only)
#' @return data.table with a row per parameter and colums for Rhat, ESSB and ESST
treeConvergenceStats <- function(chains, burnin = 0, n = 20, treedist = "PD", credMass = 0.95, byChain = T){
  sdCols <- c("Rhat","essB","essT","medianP","meanP","HDILower","HDIUpper")
  chains <- check.chains(chains)
  chain <- chains[[1]]
  indices <- seq(from = burnin + 1, to = length(chain$trees), 
                 by = 1)
  trees <- lapply(chains, function(x) x[["trees"]][indices])
  replicatedStanStats <- rbindlist(lapply(1:n,FUN = function(x){stanStats(do.call(data.table,sapply(trees,treeDistanceTrace,treedist = treedist)),credMass = credMass)}))
  returnTable <- replicatedStanStats[,lapply(.SD,mean),.SDcols = sdCols]
  returnTable[,`:=`(chain = "ALL")]
  if(byChain == T) {
    returnTableByChain <- rbindlist(lapply(1:length(trees),FUN = function(thisChain){
      theseTrees <- list(trees[[thisChain]])
      theseChainstanStats <- rbindlist(lapply(1:n,FUN = function(x){stanStats(do.call(data.table,sapply(theseTrees,treeDistanceTrace,treedist = treedist)),credMass = credMass)}))
      theseChainsummaryStanStats <- theseChainstanStats[,lapply(.SD,mean),.SDcols = sdCols]
      theseChainsummaryStanStats[,`:=`(chain = thisChain)]}))
    returnTable <- cbind(returnTable,returnTableByChain[,.(medianPCV=sd(medianP)/mean(medianP),meanPCV=sd(meanP)/mean(meanP))])
    returnTable <- rbind(returnTable,returnTableByChain,fill=T)
  }
  returnTable[,`:=`(param = paste(sep="_",treedist,"treeTopology"))][]
  setkey(returnTable,param)
  return(returnTable)
}

#' Generates tree-distance traces for all MCMC chains for plotting
#' To calculate ESS using the pseudo ESS approach, it is much better to run several replicates and calculate the mean. See treeConvergenceStats
#'
#' @param chains MCMC data. List of rwty.chains (RWTY's load.trees value)
#' @param burnin Number of MCMC samples to discard as burnin
#' @param treedist Tree distance to be used to calculate the distance to the focal tree
#' @return data.table with a row per MCMC sample and columns for treedistance and chain
treeDistanceTraces <- function(chains,burnin,treedist = "PD"){
  chains <- check.chains(chains)
  chain <- chains[[1]]
  indices <- seq(from = burnin + 1, to = length(chain$trees), 
                 by = 1)
  trees <- lapply(chains, function(x) x[["trees"]][indices])
  returnTable <- rbindlist(lapply(1:length(trees),FUN = function(ichain,treedist){tree.list <- trees[[ichain]]
  data.table(chain = ichain, treedistance=treeDistanceTrace(tree.list,iTree = 1,treedist = treedist)$topological.distance)},treedist = treedist))
}

#' Plots posterior sample and estimated values for a parameter
#' 
#' @param theseData data.table, with burnin removed if needed
#' @param parameterNames vector with the names of the parameters to print
#' @param outDir output directory
#' @param outName output name
#' @param ndigits number of digits to round numbers printed in the plot
#' @param width width in inches
#' @param height height in inches
#' @param ndigits number of digits to round the estimated and simulated values labeled on the plot
#' @return void
writePlotContinuousParameters=function(theseData,parameterNames,outDir,outName,ndigits=3,width=7,height=7){
  pdf(file=paste0(outDir,"/",outName),width = width, height = height)
  for(param in parameterNames){
    theplot=ggplot(data=theseData,aes(x=.data[[param]],group=chain,fill=chain,color=chain)) +
      geom_density(alpha=0.3,linewidth=1) +
      stat_summary(aes(xintercept = after_stat(x),y=0,linetype="Estimated",fill=NULL),fun=sumFun,geom="vline",orientation="y", show.legend = F) +
      stat_summary(aes(xintercept = after_stat(x),y=0,linetype="Estimated",fill=NULL,group=NULL),fun=sumFun,geom="vline",orientation="y",color="black") +
      stat_summary(fun=sumFun,geom="text_repel",aes(x=.data[[param]],label=round(after_stat(x),digits=ndigits),y=0),orientation="y",hjust=-0.5,direction="y", show.legend = F)+
      stat_summary(fun=sumFun,geom="text_repel",aes(x=.data[[param]],label=round(after_stat(x),digits=ndigits),y=0),orientation="y",hjust=1.5,inherit.aes=F, show.legend = F)+
      scale_linetype(name="") +
      scale_color_brewer(name="Chain",type = "qual",palette = 2) +
      scale_fill_brewer(name="Chain", type = "qual",palette = 2) +
      scale_x_continuous(name="Parameter value") +
      scale_y_continuous(name="Density")+
      labs(title=param)+
      theme_cowplot()
    print(theplot)
    #save_plot(filename = paste0(outDir,"/",param,format),theplot,base_height = 6)
  }
  dev.off()
  
}
#######################################

#Hardcoded config
#######################################
#Convergence and mixing standards
minESS <- 200
maxRhat <- 1.1
#maxCV <- 0.01 #TODO what is a good maxCV?
credMass <- 0.95

#Output config
problematicOutputSuffix <- "problematicParams.csv"
allOutputSuffix <- "allParams.csv"
plotName.continuousParameters <- "plotPosteriorDensities.pdf"
plotName.rwty <- "plotsRWTY.pdf"
plotName.prefix.correlations <- "plotCorrelations."
MLEOutputSuffix <- "MLE.csv"
#######################################

#Parsing command-line arguments and checking input
#######################################
args <- commandArgs(trailingOnly = TRUE) #Comment for debugging, uncomment to run with Rscript

#DEBUG
#######################################
#baseDir="~/projects/flipFlop/infantCrypts"
#args=c(0,1,paste(sep="/",baseDir,"analyses/overS/AI_colon_2_overS"),paste(sep="/",baseDir,"runs/combined/overS/AI_colon_2_overS.trees"))

#baseDir="~/projects/flipFlop/simulationStudy/simStudy"
#args=c(0.1,8,paste(sep="/",baseDir,"analysis"),paste(sep="/",baseDir,"333/sim_3_0.1_0.001_0.001_100_9_3cells.trees"),paste(sep="/",baseDir,"666/sim_3_0.1_0.001_0.001_100_9_3cells.trees"),paste(sep="/",baseDir,"999/sim_3_0.1_0.001_0.001_100_9_3cells.trees"))
#args=c(0.1,8,paste(sep="/",baseDir,"analysis"),paste(sep="/",baseDir,"333/sim_10_1.0_0.01_0.01_100_9_10cells.trees"),paste(sep="/",baseDir,"666/sim_10_1.0_0.01_0.01_100_9_10cells.trees"),paste(sep="/",baseDir,"999/sim_10_1.0_0.01_0.01_100_9_10cells.trees"))

#baseDir <- "~/Downloads/short_worflow_res/169"
#args <- c(0.1, 8, baseDir, paste(sep="/",baseDir,"169.3cells.trees"))

#baseDir <- "~/Downloads/short_worflow_res/test"
#inDir <- "~/Documents/Ciencia/Postdoc/projects/flipFlop/beFirstRun/"
#args <- c(0.1, 8, baseDir, paste(sep="/",inDir,"666/169flipflopS3.trees"),paste(sep="/",inDir,"999/169flipflopS3.trees"))

#baseDir <- "~/Downloads/short_worflow_res/test"
#inDir <- "~/Documents/Ciencia/Postdoc/projects/flipFlop/beFirstRun/"
#args <- c(0.1, 8, baseDir, paste(sep="/",inDir,"666/169flipflopS4.trees"),paste(sep="/",inDir,"999/169flipflopS4.trees"))

#baseDir <- "~/Downloads/short_worflow_res/test"
#inDir <- "~/Documents/Ciencia/Postdoc/projects/flipFlop/beFirstRun/"
#args <- c(0.1, 8, baseDir, paste(sep="/",inDir,"666/169flipflopS5.trees"),paste(sep="/",inDir,"999/169flipflopS5.trees"))

#baseDir <- "~/Downloads/short_worflow_res/test"
#inDir <- "~/Documents/Ciencia/Postdoc/projects/flipFlop/beFirstRun/"
#args <- c(0.1, 8, baseDir, paste(sep="/",inDir,"666/169flipflopS6.trees"),paste(sep="/",inDir,"999/169flipflopS6.trees"))

#baseDir <- "~/Downloads/short_worflow_res/test"
#inDir <- "~/Documents/Ciencia/Postdoc/projects/flipFlop/beFirstRun/"
#args <- c(0.1, 8, baseDir, paste(sep="/",inDir,"666/169flipflopS7.trees"),paste(sep="/",inDir,"999/169flipflopS7.trees"))

#baseDir <- "~/Downloads/short_worflow_res/test"
#inDir <- "~/Documents/Ciencia/Postdoc/projects/flipFlop/beFirstRun/"
#args <- c(0.1, 8, baseDir, paste(sep="/",inDir,"666/169flipflopS8.trees"),paste(sep="/",inDir,"999/169flipflopS8.trees"))

#######################################

if (length(args) <=3) {
  stop("Usage: script burnin_proportion n_cores baseOutdir sample1.tree ... sampleN.tree")
} else {
  ##ForPablo: Remove if we don't want this to be verbose
  print(paste0("Running the script in the trace files (",paste(collapse =", ",args[-c(1:3)]),") with a burnin of ",args[1]," proportion of the trees and ",args[2]," processors"))
}

##ForPablo: some sections of the code can run in parallel, but I assume we will not use this in the pipeline
##We could just modify this section and not get a n_cores argument, fixing rwty.processors to 1, or have the 1 as default in the call to this script. Up to you.
rwty.processors=as.numeric(args[2])
if(rwty.processors > 1){
  RNGkind("L'Ecuyer-CMRG") #Different random number generator, so that mclapply is reproducible, only needed if rwty.processors > 1
}
setDTthreads(rwty.processors)
burninP=as.numeric(args[1]) #ForPablo: the default should be 0.1
baseDir=args[3]
baseName=args[4]
files=args[-c(1:4)]

if(!all(file.exists(files))){
  stop("ERROR: Not all input files exist. ",paste(files,collapse=", "))
}
#######################################

#Preparing output
#######################################
outDir <- baseDir #paste(sep="/",baseDir,baseName)
dir.create(outDir,recursive=T)
print(paste0("Saving outputs in ",outDir)) #ForPablo: remove?
#######################################

#Parsing input traces
#######################################
sink(nullfile()) #shut up!
chains=invisible(lapply(files,FUN=function(x){
  logfile <- gsub(x,pattern = ".trees",replacement = ".log")
  if(!file.exists(logfile)) stop(sprintf("ERROR: logfile %s for tree file %s not found",logfile,x)) ##TODO this will not be visible because I am trying to mute the rest of the code in this section
  chain <- load.trees(x,format="BEAST",logfile = logfile);
  chain}))
nSamples=unique(lapply(chains,FUN=function(x){length(x$trees)}))
if(length(nSamples) != 1){
  stop("ERROR: chains are not of the same length")
}
#######################################

#Data reorganization
#######################################
theseData=getDataTable(chains,burninP)
burninStates=theseData[,max(state)]*burninP
burninSamples=theseData[chain==1 & burnin==T,.N]
#######################################

#Data augmentation
#######################################
theseData[,`:=`(flipflop.Rmu=flipflop.mu/flipflop.lambda,flipflop.Rgamma=flipflop.gamma/flipflop.lambda,flipflop.MDbias=flipflop.mu/flipflop.gamma)]
chains <- updateChainsFromDataTable(chains = chains,theData = theseData)
#######################################

#Assess tree convergence 
#######################################
#WARNING: We are doing this separate to be able to run 20 different pseudoESS estimations and report their mean instead of using just one
topologyConvergence <- treeConvergenceStats(chains, burninSamples, credMass = credMass, byChain = ifelse(length(chains)>1,T,F))

#Make one tree-distance trace from the tree with index 1 for tree-trace plots
treeTopologyData <- treeDistanceTraces(chains,0)[,`:=`(burnin=ifelse(.I<.N*burninP,T,F))][]
theseData <- cbind(theseData,treeTopologyData[,.(treedistance)])
#######################################

#Assessing convergence of continuous variables
#######################################
paramsForConvergence <- names(theseData)
paramsForConvergence <- paramsForConvergence[!paramsForConvergence%in%c("state","chain","burnin","treedistance")]#Never use this treedistance data for convergence, only for plotting
constantParameters <- detectConstants(theseData,c(paramsForConvergence,"treedistance"))
warning(sprintf("The parameters (%s) are constant and not taken into account for mixing and convergence evaluation",paste(collapse=",",constantParameters))) #ForPablo: This information should be relayed to the user
continuousConvergence <- getConvergenceTableFromDataTable(theseData[burnin==F,], paramsForConvergence[!paramsForConvergence%in%constantParameters], credMass = credMass, byChain = ifelse(length(chains)>1,T,F))
convergence <- rbind(continuousConvergence,topologyConvergence)
setkey(convergence,param)
#######################################

#Plotting continuous parameters
#######################################
#Posterior sample for parameters without true values
writePlotContinuousParameters(theseData[burnin==F],paramsForConvergence[!paramsForConvergence%in%constantParameters],outDir,paste(sep=".",baseName,plotName.continuousParameters))
#######################################

#Get marginal likelihood estimates from the posterior sample 
#######################################
thisLML=LML(LL=theseData[burnin==F,likelihood],method="HME")$LML
thisAICM=AICM(theseData[burnin==F,likelihood])
#######################################

#RWTY plots including traces
#######################################
rwtyPlots <- analyze.rwty(chains, burnin=burninSamples, fill.color = 'posterior', params = paramsForConvergence[!paramsForConvergence%in%constantParameters])

#Print all but the correlations plots in pdf
pdf(file = paste0(outDir,"/",paste(sep=".",baseName,plotName.rwty)), width = 10, height = 7)
print(rwtyPlots[-grep(".correlations",names(rwtyPlots))])
dev.off()

#Correlations in jpeg to avoid computers crashing
for (plotname in names(rwtyPlots[grep(".correlations",names(rwtyPlots),value = T)])){
  jpeg(filename = paste0(outDir,"/",paste(sep=".",baseName,plotname),".jpeg"),width = 1200, height = 1200,type = "quartz")
  print(rwtyPlots[[plotname]])
  dev.off()
}
#######################################

#Write output
#######################################

#All convergence statistics
outColumns=c("param","chain","Rhat","essB","essT","medianP","meanP","HDILower","HDIUpper")
allOutputFileName=paste0(outDir,"/",paste(sep=".",baseName,allOutputSuffix))
write.csv(convergence[,.SD,.SDcols = outColumns],file = allOutputFileName,quote = F,row.names = F)

#Parameters that do not meet certain standards

#For now, we consider that a param has a problem if:
# 1. It has not achieved the minimum ess in all chains (including all combined)
# 2. It has a Rhat > than the limit with all chains combined
#TODO 3. The estimates from the different chains have a CV< ??? Unsure what this value should be
paramProblems <- cbind(convergence[chain=="ALL",.(Rhat=is.na(Rhat)|Rhat>maxRhat),keyby=param],
      convergence[,lapply(.SD,FUN = function(x){any(is.na(x)) | any(x<minESS)}),keyby=param,.SDcols=c("essB","essT")][,-"param"])#,
      #convergence[chain=="ALL",lapply(.SD,FUN = function(x){any(is.na(x)) | any(x>maxCV)}),keyby=param,.SDcols=c("medianPCV","meanPCV")][,-"param"])
paramProblems <- paramProblems[,.(param,problem=any(.SD)),.SDcols=colnames(paramProblems)[!colnames(paramProblems)%in%c("param")],by=.I][problem==T,param]
problematicOutputFileName=paste0(outDir,"/",paste(sep=".",baseName,problematicOutputSuffix))
write.csv(convergence[param %in% paramProblems],file = problematicOutputFileName,quote = F,row.names = F)

#ForPablo: I do not know if you prefer to do this here or in your python code etc. Also, you can parse stdout of this script or directly read the tables written above. Up to you :)
#Analyze and communicate the number of problematic parameters
#if only tree topology, warn that they should still be checked by eye and that it may be ok if the number of leaves is small. Show the number of leaves?
#if there are more, warn that the results should be checked by hand and something in the analysis should be tweaked to improve mixing-convergence before proceeding
if(length(paramProblems)==0){
  print("No insights of problematic parameters. Please, still review the plots to confirm the analysis worked properly before interpeting or publishing these results")
} else if(length(grep("treeTopology",paramProblems,value = T)) == length(paramProblems)){
  warning(sprintf("there were not any insights of problematic continuous parameters. However, the tree topology has mixing and/or convergence issues. This may not be an issue for cases with very small trees (less than a handful of samples). For reference, the trees had %d leaves in these analyses",length(chains[[1]]$trees[[1]]$tip.label))
)} else {
  warning("WARNING: convergence statistics show convergence problems. You must re-run the analysis after tweaking the MCMC parameters to solve these issues. These results must not be used for biological interpretation or publication")
}

#MLE information
mleTable <- data.table(cond=c(baseName),method=c("HME","AICm"),lML=c(thisLML,thisAICM))
MLEOutputFileName=paste0(outDir,"/",paste(sep=".",baseName,MLEOutputSuffix))
write.csv(mleTable,file = MLEOutputFileName,quote = F,row.names = F)
#######################################
sink()
