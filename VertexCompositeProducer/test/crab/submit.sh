cd /home/awesole/VertexCP_clean/CMSSW_13_2_11/src || { echo "Directory not found"; exit 1; }
cmsenv
source /cvmfs/cms.cern.ch/common/crab-setup.sh

cd -

crab submit -c D0_mcTREE_submit_PT-0.py
crab submit -c D0_mcTREE_submit_PT-1.py
crab submit -c D0_mcTREE_submit_PT-8.py
crab submit -c D0_mcTREE_submit_PT-10.py
crab submit -c D0_mcTREE_submit_PT-20.py
crab submit -c D0_mcTREE_submit_PT-30.py
crab submit -c D0_mcTREE_submit_PT-40.py
crab submit -c D0_mcTREE_submit_PT-60.py
crab submit -c D0_mcTREE_submit_PT-80.py
#crab submit -c D0_dataTREE_submit.py

