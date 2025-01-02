suppressPackageStartupMessages({ 
  library(pacman)
  p_load(data.table, ggplot2, cowplot, ggbeeswarm)
})
if(!interactive()) pdf(NULL) #Avoid generating empty Rplots.pdf when using save_plot

#Hardcoded config
#######################################
#Model comparison threshold
minBF <- 10

#Convergence and mixing standards
#minESS <- 200
#maxRhat <- 1.1

#Plotting
aspectRatio=4/2

#Input parsing
ncellsRegex <- ".*\\.([0-9]+)cells\\..*" #Phyfum version
basenameRegex <- "\\.[0-9]+cells\\..*" #Phyfum version

#Output config
plotName.allOutputSuffix <- "allS.pdf"
selectedOutputSuffix <- "selectedS.csv"
bestSOutputSuffix <- "bestFitS.csv"
allOutputSuffix <- "allS.csv"
#######################################

#Parsing command-line arguments and checking input
#######################################
args <- commandArgs(trailingOnly = TRUE) #Comment for debugging, uncomment to run with Rscript

#DEBUG
#######################################
# ncellsRegex <- ".*S([0-9]+).*" #DEBUG version
# basenameRegex <- "S[0-9]+.*" #DEBUG version
# baseDir <- "~/Downloads/short_worflow_res/test"
# inDir <- "~/Downloads/short_worflow_res/test"
# args <- c(baseDir, paste(sep="/",inDir,"169flipflopS3/169flipflopS3.MLE.csv"), 
#           paste(sep="/",inDir,"169flipflopS4/169flipflopS4.MLE.csv"), 
#           paste(sep="/",inDir,"169flipflopS5/169flipflopS5.MLE.csv"),
#           paste(sep="/",inDir,"169flipflopS6/169flipflopS6.MLE.csv"),
#           paste(sep="/",inDir,"169flipflopS7/169flipflopS7.MLE.csv"),
#           paste(sep="/",inDir,"169flipflopS8/169flipflopS8.MLE.csv"))
#######################################
if (length(args) <3) {
  stop("Usage: script outDir runNameSX.csv ... runNameSN.csv")
} else {
  ##ForPablo: Remove if we don't want this to be verbose
  print(paste0("Running the script in the MLE summaries (",paste(collapse =", ",args[-c(1)]),")"))
}

baseDir=args[1]
files=args[-c(1)]
baseName <- gsub(x = basename(files[1]),pattern = basenameRegex,replacement = "")

if(!all(file.exists(files))){
  stop("ERROR: Not all input files exist. ",paste(files,collapse=", "))  ##ForPablo: this is probably not necessary if the pipeline checks that files exist before using this script
}
#######################################

#Data collection and output
#######################################
theseDataRaw <- rbindlist(lapply(files,FUN = function(thisFile){
  thisS <- as.numeric(gsub(ncellsRegex,"\\1",basename(thisFile)))
  if(is.na(thisS)) stop(sprintf("The number of stem cells used for file %s was not detected properly",thisFile)) 
  thisTable <- fread(thisFile)
  thisTable[,`:=`(S=thisS)]
  }))

theseData <- theseDataRaw[method!="AICm",] #We are not using AIC for now
if(theseData[method=="PP",.N]) warning("Stepping stone and/or path sampling lML estimates not found. This script will use the less-accurate harmonic mean estimator, which is generally recommended against. In our simulation results, it has performed as well as PS and SS for the specific task of finding the optimum S parameter.")

#Data dump
allOutputSFileName <- paste0(baseDir,"/",paste(sep=".",baseName,allOutputSuffix))
write.csv(theseData,file = allOutputSFileName,quote = F,row.names = F)
#######################################

#Plotting
#######################################

#Plot all S values for LML, PP, and SS (or HME if it is the only present)
theseThresholds <- theseData[,.(minY=max(lML)-minBF),by=method]

thePlot <- ggplot(theseData,aes(x=S,y=lML,color=method))+
  geom_point(alpha=0.8) +
  geom_smooth(alpha=0.2,linewidth=0.3) +
  geom_hline(data=theseThresholds,aes(yintercept=minY,color=method,linetype="BF10")) +
  scale_y_continuous(name="Marginal Likelihood Estimation (logL)",n.breaks=10) +
  scale_x_continuous(name="Number of stem cells (S)") +
  scale_color_brewer(name="Estimation\nmethod",type = "qual",palette = 6) +
  scale_linetype_manual(name=NULL,values=2)+
  theme_cowplot()

save_plot(file=paste(sep="/",baseDir,paste(sep=".",baseName,plotName.allOutputSuffix)),thePlot,base_height = 6,base_asp = aspectRatio)
#######################################

#Model selection and final output
#######################################
bestFits <- theseData[order(-lML),.(lML=first(lML),S=first(S)),by=method]
bestSFileName <- paste0(baseDir,"/",paste(sep=".",baseName,bestSOutputSuffix))
write.csv(bestFits,file = bestSFileName,quote = F,row.names = F)

#Output the valid ones under each of them
validBFs <- theseData[order(-lML),.(valid=max(lML)-lML<=minBF,S,lML),by=method][valid==T,.(method,S,lML)]
selectedFileName <- paste0(baseDir,"/",paste(sep=".",baseName,selectedOutputSuffix))
write.csv(validBFs,file = selectedFileName,quote = F,row.names = F)

#Warn if the best fit value is in one of the extremes
extremeS <- theseData[,.(minS=min(S),maxS=max(S)),by=method]
for (thisMethod in theseData[,unique(method)]){ #ForPablo: these are important messages for the user
  if(bestFits[method==thisMethod,S] == extremeS[method==thisMethod,minS]) warning(sprintf("IMPORTANT: The best-fit S value %d is a extreme of the range explored using the %s method. The analysis should be re-run adding lower S values to find the global optimum",bestFits[method==thisMethod,S],thisMethod))
  if(bestFits[method==thisMethod,S] == extremeS[method==thisMethod,maxS]) warning(sprintf("IMPORTANT: The best-fit S value %d is a extreme of the range explored using the %s method. The analysis should be re-run adding higher S values to find the global optimum",bestFits[method==thisMethod,S],thisMethod))}

#Warn if any of the values within maxBF is one of the extremes
for (thisMethod in theseData[,unique(method)]){ #ForPablo: these are important messages for the user
  if(validBFs[method==thisMethod,min(S)] == extremeS[method==thisMethod,minS]) warning(sprintf("The smaller assayed S value %d cannot be rejected as a valid S value using the %s method. We recommend the analysis is re-run adding smaller S values to describe better the range of S values that may be valid",extremeS[method==thisMethod,minS],thisMethod))
  if(validBFs[method==thisMethod,max(S)] == extremeS[method==thisMethod,maxS]) warning(sprintf("The larget assayed S value %d cannot be rejected as a valid S value using the %s method. We recommend the analysis is re-run adding larger S values to describe better the range of S values that may be valid",extremeS[method==thisMethod,maxS],thisMethod))}
#######################################

#ForPablo: I would continue the analysis (trees etc.) with the best fit under the PS model, if present, otherwise HME
