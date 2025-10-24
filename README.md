---------------------------------------------------------------------------------------------------------------
<p align="center"> --------This code is written to process the 2023 PbPb Data-------------</p>
<p align="center"> ------------------The output is a updated TTree.root file----------------------</p>
<p align="center"> ------------------Authors: N. Saha, A. Wesolek -------------------------- </p>

 ---------------------------------------------------------------------------------------------------------------- 
<p align="center"> nihar.ranjan.saha@cern.ch </p>
<p align="center"> abigail.leigh.wesolek@cern.ch </p>


----------------------------------------------------------------------------------------------------------------
  

  <br>
 <br>
 <br>
 
This code was originally cloned from: <br>
https://github.com/NiharSaha/VertexCP_miniAOD/tree/D0analyzer_reproduction_TTree_skimmedEdm <br>
Commit: 133bfac <br>
Author: Nihar Saha <br>
 
-----------------------------------------------------------
**This branch of the code is used for reproducing the TTrees to include some missing variables. It was written by Nihar Saha. 
  It reads in the previously produced edm.root files and the correctly associated miniAOD file ('parent file') and produces the d0Analyzer tree, event info tree and zdc tree.**



## How to run

```bash 
#LXplus, bash, cmssw-el8 apptainer

cmsrel CMSSW_13_2_11

cd CMSSW_13_2_11/src
cmsenv

#clone the repo
git clone --branch TTreeReproduction_2023PbPb --single-branch https://github.com/a-wesole/VertexCP_miniAOD.git VertexCompositeAnalysis

cd VertexCompositeAnalysis
./setup.sh

#compile
cd ..
scram b -j12

#cd and run code that produces edm and TTree files 

cd VertexCompositeAnalysis/VertexCompositeProducer/test

cmsRun /run_edm_and_ttree_DATA_forD0_withParentFile_andZDC_andEP.py #for reproduction of D0 data TTrees, including EP, ZDC. reads edm.root files from Selector as input, instead of miniAOD 


```

