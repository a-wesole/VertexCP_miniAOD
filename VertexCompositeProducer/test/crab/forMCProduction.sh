#!/bin/bash

#######################
# Instructions
# How to run: ./forMCProduction.sh <username>
# will write to /store/user/<username>/D0MC_Run3

#Notes:
#Will submit MC crab jobs for all Dpt values below 

# make sure to establish voms proxy 

#will complain with ERROR: Invalid CRAB type 'awesolek' provided, valid values are prod, pre and dev.
#not a terminal issue - just ignore it 


########################

#Dpt=(0)
Dpt=(0 1 8 10 20 30 40 60 80) 

# Function to replace a string in a file using sed
replace_string() {
    local filename=$1
    local old_string=$2
    local new_string=$3

    if [ ! -f "$filename" ]; then
        echo "The file $filename does not exist."
        exit 1
    fi

    sed -i "s|$old_string|$new_string|g" "$filename"
}



# Navigate to the CMSSW directory and set up the environment
cd /home/awesole/VertexCP_clean/CMSSW_13_2_11/src || { echo "Directory not found"; exit 1; }
cmsenv
#scram b -j8
#source /cvmfs/cms.cern.ch/common/crab-setup.sh
#voms-proxy-init
#voms-proxy-init -voms cms

# Navigate to the crab directory
cd /home/awesole/VertexCP_clean/CMSSW_13_2_11/src/VertexCompositeAnalysis/VertexCompositeProducer/test/crab || { echo "Directory not found"; exit 1; }

for num in "${Dpt[@]}"; do


	# Replace strings in tmp_submit.py
	rm tmp_submit.py
	cp Lc_multicrab_1.py Dpt${num}_submit.py

	# Call replace_string function with the necessary arguments
	replace_string "Dpt${num}_submit.py" "awesolek" "wxie"
	replace_string "Dpt${num}_submit.py" "config.General.requestName = ''" "config.General.requestName = 'promptD0ToKPi_PT-$num'"
	if [ $num -ge 10 ]; then

		replace_string "Dpt${num}_submit.py" "config.Data.inputDataset = ''" "config.Data.inputDataset = '/promptD0ToKPi_PT-${num}_TuneCP5_5p36TeV_pythia8-evtgen/HINPbPbSpring23MiniAOD-132X_mcRun3_2023_realistic_HI_v9-v1/MINIAODSIM'"
	else 
		replace_string "Dpt${num}_submit.py" "config.Data.inputDataset = ''" "config.Data.inputDataset = '/promptD0ToKPi_PT-${num}_TuneCP5_5p36TeV_pythia8-evtgen/HINPbPbSpring23MiniAOD-132X_mcRun3_2023_realistic_HI_v9-v2/MINIAODSIM'"
	fi
	replace_string "Dpt${num}_submit.py" "config.Data.outLFNDirBase = '/store/user/awesolek/D0MC_Run3'" "config.Data.outLFNDirBase = '/store/user/wxie/D0MC_Run3'" 
	replace_string "Dpt${num}_submit.py" "config.JobType.psetName = '../run_VCProducer.py'" "config.JobType.psetName = '/home/awesole/VertexCP_clean/CMSSW_13_2_11/src/VertexCompositeAnalysis/VertexCompositeProducer/test/run_VCProducer.py'"

	#crab submit -c Dpt${num}_submit.py

	echo "   "
	echo " *****************************  "
	echo "Submitted CRAB jobs for Dpt> ${num} "
	echo " *****************************  "
	echo "   "

done


