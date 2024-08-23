#!/bin/bash

# Execute this script for each S independently


usage="Usage: $0 inputfile1 ... inputfileN outputFile.\n This script expects either a single (BEAST stdout) input file when only one chain per S was run, or multiple .MLE.log input files when multiple chains per S were run."

#ForPablo:
########## 
#this script should be run after mixingAndConvergence.R
#when only one chain per S is run, inputFile should be a log of BEAST's stdout (but alternatively, see comment within first if
#when multiple chains per S are run, we need to run beast to re-calculate the MLE merging the results from the multiple chains (done below). For this, inputFiles are the .MLE.log files that BEAST generates
#outputFile is always the MLE.csv table that mixingAndConvergence.R generates (this script concatenates to that file, that's why it must be run after it)

[ $# -eq 0 ] || [ "$1" == "-h" ] || [ "$1" == "--help" ] && echo -e "$usage" && exit 1
[ $# -lt 2 ] && echo -e "ERROR incorrect number of arguments: $usage" && exit 1

if [ $# -eq 2 ]
then
	#For cases with only one chain per condition, we can just get the information from the stdout log. As here. If we do not have it, we can remove this "if" and input the 
	resultContent=$(<$1)
	outputFile=$2	
else
	args=( "$@" )
	outputFile=${args[$#-1]}
	
	unset 'args[$#-1]'
	#inputString=$(echo $@ | sed "s/ [^ ][^ ]*$//g") #Alternative without using arrays, less safe?
	inputString=${args[@]}
	inputString=$(realpath $inputString)
	tempFile="${outputFile}_tmp.xml" #ForPablo: we could also keep this, but for now it is a temp file that is deleted below
	cat << EOF > $tempFile
<beast>
    <pathSamplingAnalysis fileName="$inputString">
        <likelihoodColumn name="pathLikelihood.delta"/>
        <thetaColumn name="pathLikelihood.theta"/>      
    </pathSamplingAnalysis>
    <steppingStoneSamplingAnalysis fileName="$inputString">
        <likelihoodColumn name="pathLikelihood.delta"/>
        <thetaColumn name="pathLikelihood.theta"/>      
    </steppingStoneSamplingAnalysis>
</beast>
EOF

	resultContent=$(beast -beagle_off $tempFile) #ForPablo: if beast is not in the PATH we will need to change this
	rm -f $tempFile	
fi

ps=$(echo "$resultContent" | grep 'path sampling)' | sed "s/^.*= \(.*\)$/\\1/g")
ss=$(echo "$resultContent" | grep 'stone sampling)' | sed "s/^.*= \(.*\)$/\\1/g")

echo -e "NA,SS,$ss\nNA,PS,$ps" >> $outputFile
