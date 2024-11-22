#!/usr/bin/env Rscript

# Call this as follows:

## With logfile
# Rscript pisca-plot.R consensus.tree title,pdf_name,use_rate,burnin,age,my_pisca.log


suppressWarnings({
library(pacman)
p_load(data.table, ggplot2, cowplot, treeio, ggtree, HDInterval, lubridate, svglite)
pdf(NULL)

#Configuration
burnin=0.1
ppThreshold=0.8
credMass=0.95

folder = getwd()
args=commandArgs(trailingOnly = T)
# the args can be entered in 2 ways
# - from pisca-box we don't pass in the log file
# - from a script we might do so

mccTreeFile=args[1]
title=args[2]
outputfile=args[3]
burnin = as.numeric(args[4])
age=as.numeric(args[5])
#Parsing cenancestor information from the log file    
logDataFile=args[6]
logData=fread(logDataFile)
lucaBranch=mean(logData[,luca_branch][floor(nrow(logData)*burnin):nrow(logData)])
hpdLucaBranch=hdi(logData[,luca_branch][floor(nrow(logData)*burnin):nrow(logData)],credMass = credMass)

rm(logData)

print(paste("mccfile",mccTreeFile))
print(paste("title",title))
print(paste("outputfile",outputfile))
print(paste("burnin",burnin))
print(paste("age",age))
print(paste("lucaBranch",lucaBranch))
print(paste("hpdLucaBranch",hpdLucaBranch))


#rm(logData)

##Parsing the tree
mccTree=read.beast(mccTreeFile)
mccTreeDataFrame=fortify(mccTree)


#Adding luca with proper time
root=which(mccTreeDataFrame$parent==mccTreeDataFrame$node)
luca=nrow(mccTreeDataFrame)+1
mccTreeDataFrame[luca,]=mccTreeDataFrame[root,]
mccTreeDataFrame[root,]=mccTreeDataFrame[luca,]
mccTreeDataFrame[root,"parent"]=luca
mccTreeDataFrame[root,"branch.length"]=lucaBranch
mccTreeDataFrame[luca,"x"]=-lucaBranch
mccTreeDataFrame[luca,"node"]=luca

##Years to add - make x axis the age of sampling (tree file goes from 0-max height, not age)
yearsToAdd=age-max(mccTreeDataFrame[,"x"])

#Modifying times from time from the present to time from birth
mccTreeDataFrame$x=mccTreeDataFrame$x+yearsToAdd
mccTreeDataFrame$branch=mccTreeDataFrame$branch+yearsToAdd
mostRecentYearsOfAge=age
mccTreeDataFrame$height_0.95_HPD=lapply(mccTreeDataFrame$height_0.95_HPD,function(x){return(c(mostRecentYearsOfAge-x[2],mostRecentYearsOfAge-x[1]))})
mccTreeDataFrame$height_0.95_HPD[luca]=list(as.numeric(c(mostRecentYearsOfAge-hpdLucaBranch[2],mostRecentYearsOfAge-hpdLucaBranch[1])))

#Adding PP info
mccTreeDataFrame$posteriorAsterisk=as.character(ifelse(mccTreeDataFrame$posterior>=ppThreshold,"*",""))

#Removing PP data for LUCA
mccTreeDataFrame[luca,"posteriorAsterisk"]=NA
mccTreeDataFrame[luca,"posterior"]=NA

#Removing PP data for MRCA
mccTreeDataFrame[root,"posteriorAsterisk"]=NA
mccTreeDataFrame[root,"posterior"]=NA

#Plotting
#########

#Adjusting x limits to make sure tip names are shown completely
minX=min(mccTreeDataFrame[,"x"])
maxX=max(mccTreeDataFrame[,"x"])
extraX=0.2*(maxX-minX)
padMinX=2

minXError=min(as.vector(sapply(mccTreeDataFrame$height_0.95_HPD,function(x){return(c(x[1],x[2]))},simplify = TRUE)),na.rm = T)

#Getting y information
mccTreePlot=ggplot(mccTreeDataFrame,aes(x=x,y=y))+geom_tree(size=1.5)+
  theme_tree2(legend.position=c(0.1,0.7))+
  geom_nodelab(aes(label=round(posterior,2)),hjust=-.1,vjust=1.8,color="black") + ##PP values
  geom_range("height_0.95_HPD", color='black', size=5, alpha=.3) + #HPD height
  geom_tiplab(align = FALSE,color="black",size=5,hjust=-0.2) +
  labs(title=title) + 
  theme(plot.title= element_text(size = rel(2.5),hjust = 0.5),axis.text.x=element_text(size=rel(1.5)),axis.title.x=element_text(size=rel(1.5))) +
  scale_x_continuous(name="Patient age (years)",limits=c(min(minX,minXError)-padMinX,maxX+extraX))

save_plot(plot = mccTreePlot,filename = outputfile,base_height = 5,base_aspect_ratio = 1.6)

})
