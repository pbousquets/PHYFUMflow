#!/usr/bin/env Rscript
if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load(optparse, LDATS, ape, LaplacesDemon, data.table, lubridate)#, posterior) #only for second algorithm not yet implemented

option_list = list(
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="outTree filename", metavar="outTree"),
  make_option(c("-m","--mle"), type="character", default=NULL, 
              help="MLE input file. If used, the priority is SS>PS>HME. If not used, HME calculated here is used instead", metavar="mleFile"),
  make_option(c("-b","--burnin"), type="double", default=0.1, 
              help="burnin proportion [default= %default]", metavar="burninProportion"),
  make_option(c("-s","--seed"), type="double", default=NULL, 
              help="random number generator seed", metavar="seed")
);

opt_parser <- OptionParser(usage = "%prog [options] treefile1 ... treefileN", option_list=option_list)

#DEBUG
#######################################
#Equivalent: Rscript integratingOverS.R -o ~/projects/flipFlop/infantCrypts/runs/combined/overS/AI_colon_2_overS.trees -b 0 -s 20 -m ~/projects/flipFlop/infantCrypts/analyses/AI_colon_2/AI_colon_2.allS.csv ~/projects/flipFlop/infantCrypts//runs/combined/AI_colon_2_s*.trees
# baseDir="~/projects/flipFlop/infantCrypts/"
# arguments=list(options=list(burnin=0,
#                             out=paste(sep="/",baseDir,"runs/combined/overS/AI_colon_2_overS.trees"),
#                             mle=paste(sep="/",baseDir,"analyses/AI_colon_2/AI_colon_2.allS.csv"),
#                             seed=20),
#                args=c(paste(sep="/",baseDir,"runs/combined/AI_colon_2_s10.trees"),paste(sep="/",baseDir,"runs/combined/AI_colon_2_s2.trees"),paste(sep="/",baseDir,"runs/combined/AI_colon_2_s3.trees"),paste(sep="/",baseDir,"runs/combined/AI_colon_2_s4.trees"),paste(sep="/",baseDir,"runs/combined/AI_colon_2_s5.trees"),paste(sep="/",baseDir,"runs/combined/AI_colon_2_s6.trees"),paste(sep="/",baseDir,"runs/combined/AI_colon_2_s7.trees"),paste(sep="/",baseDir,"runs/combined/AI_colon_2_s8.trees"),paste(sep="/",baseDir,"runs/combined/AI_colon_2_s9.trees")))
#######################################
arguments <- parse_args(opt_parser, positional_arguments = c(2,Inf))
opt <- arguments$options
treeFiles <- arguments$args

if (is.null(opt$out)){
  print_help(opt_parser)
  stop("ERROR: output name and at least 2 input traces are required", call.=FALSE)
}

#nSamples <- NULL
burnin <- opt$burnin
outTreeFile <- opt$out
mleFile <- opt$mle
if(!is.null(opt$seed))
  set.seed(opt$seed)
outLogFile <- gsub(pattern = ".trees", replacement = ".log", outTreeFile)
logFiles <- gsub(pattern = ".trees", replacement = ".log", treeFiles)

if(!all(file.exists(logFiles)))
  stop("Logfiles not found")
if(!is.null(mleFile))
  if(!file.exists(mleFile))
    stop("The MLE file can't be found")

#Loading data and removing burnin
if(burnin >=1 || burnin < 0)
  stop("Burnin must be given as a proportion")
S <- as.numeric(gsub(pattern = ".*_s([0-9]*).trees",replacement = "\\1",basename(treeFiles)))
treeLists <- lapply(treeFiles,FUN = function(x){
  trees <- read.nexus(x)
  nRemove <- ceiling(length(trees)*burnin)
  return(trees[-seq(1,nRemove)])
})
logList <- lapply(logFiles,FUN = function(x){
  thisLog <- fread(x)
  nRemove <- ceiling(nrow(thisLog)*burnin)
  return(thisLog[-seq(1,nRemove)])
})

#Making sure posterior samples are complete
traceLengthMatrix <- matrix(c(sapply(logList,nrow),sapply(treeLists,length)),ncol = 2) #Lengths of log [,1] and tree [,2] traces for each chain
if(!all(apply(traceLengthMatrix,1,FUN = function(x){length(unique(x))==1})))
  stop("Posterior samples of trees and continuous parameters must be of the same length")

#Parsing or calculating MLEs and calculating sampling probabilities
MLEs <- NULL
if(!is.null(mleFile)) {
  mleData <- fread(mleFile)
  mleData <- merge(data.table(S=S),mleData,by="S",sort=F,all.x=T)
  
  #We make sure that we get all the lMLs using the same method in the pre-specified order
  if(nrow(mleData[method=="SS",]) == length(S)){
    MLEs <- mleData[method=="SS",lML]
  } else if (nrow(mleData[method=="PS",]) == length(S)) {
    MLEs <- mleData[method=="PS",lML]
  } else if (nrow(mleData[method=="HME",])== length(S)) {
    MLEs <- mleData[method=="HME",lML]
  } else {
    print("The MLE file does not contain all the lMLs needed. They will estimated using the HME method here")
  }
}

if(is.null(MLEs)){
  MLEs <- sapply(logList,FUN = function(x){LML(LL=x$likelihood,method="HME")$LML})
  if(!all(is.finite(MLEs)))
    stop("MLEs calculated here using the HME method are not finite. Make sure to use burnin if necessary")
}

probs <- softmax(MLEs)

#Info for the user
cat(paste0("Resampling ",length(S)," traces:\nS, P\n",paste(collapse="\n",sep=", ",S,round(probs,3))))

#Sampling
if(length(unique(c(sapply(logList,nrow),sapply(treeLists,length)))) == 1) {
  #All traces are the same length, so we will just sample what trace to get each sample from
  #We just copy the trace with the highest prob and replace the needed samples from others
  #Select default and use it
  defChain <- which(probs==max(probs))
  if(length(defChain)>1){
    defChain <- sample(x=defChain,size=1)
  }
  outTrees <- treeLists[[defChain]]
  outLog <- logList[[defChain]]
  theColnames <- colnames(outLog)
  #Replace when needed
  rSample <- apply(rmultinom(length(outTrees),1,probs),2,FUN=function(x){which(x==1)})
  for (iSample in 1:length(rSample)) {
    if(rSample[iSample] != defChain) {
    outTrees[iSample] <- treeLists[[rSample[iSample]]][iSample]
    outLog[iSample,(theColnames):=logList[[rSample[iSample]]][iSample,]]
    }
  }
} else {
  stop("Algorithm not implemented")
  #We will need a pre-specified number of samples
  #nSamples
  #We can just add weights to the traces and sample them based on them using for example
  #posterior:::.resample_stratified() #this should have less variance than a random sampling
  #To reduce autocorrelation, we could do something smarter but hackier: first we sample the number of samples per chain, then we thin the chains to the closest larger number of samples, resample to fit, and then mix them up (unsure what's the best way to do the last step)
  #samplesPerTrace <- as.numeric(rmultinom(1,nSamples,probs))
  #maybe use posterior:thin_draws 
}
dir.create(dirname(outTreeFile),showWarnings = F,recursive = T)
write.nexus(outTrees, file = outTreeFile)
writeLines(c("# Integrating over S after Phyfum",paste0("# ",now()),paste0(collapse=",",c("# arguments:",args))),outLogFile)
suppressWarnings(write.table(outLog, file = outLogFile,quote = F,row.names = F,sep = "\t",append = T))
