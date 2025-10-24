---------------------------------------------------------------------------------------------------------------
<p align="center"> --------This code is written to process the 2023 PbPb Data-------------</p>
<p align="center"> ------------------The output is a updated TTree.root file----------------------</p>
<p align="center"> ----------------The edm and root files are one-to-one------------------------</p>
<p align="center"> ------------------Authors: N. Saha, A. Wesolek -------------------------- </p>

 ---------------------------------------------------------------------------------------------------------------- 
<p align="center"> nihar.ranjan.saha@cern.ch </p>
<p align="center"> abigail.leigh.wesolek@cern.ch </p>


----------------------------------------------------------------------------------------------------------------
  

  <br>
 <br>
 <br>
 
This code was originally cloned from:
https://github.com/NiharSaha/VertexCP_miniAOD/tree/D0analyzer_reproduction_TTree_skimmedEdm
Commit: 133bfac
Author: Nihar Saha 
 
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
git clone https://github.com/a-wesole/VertexCP_miniAOD.git VertexCompositeAnalysis

#run the setup.sh script this will add the HeavyIonsAnalysis from CmsHI github that is needed for centrality 
cd VertexCompositeAnalysis
./setup.sh

#compile
cd ..
scram b -j12

#cd and run code that produces edm and TTree files 

cd VertexCompositeAnalysis/VertexCompositeProducer/test


cmsRun run_edm_and_ttree_MC_forD0.py #for D0 MC
cmsRun run_edm_and_ttree_DATA_forD0.py #for D0 data
cmsRun run_edm_and_ttree_MC_forLc.py #for Lc MC
cmsRun run_edm_and_ttree_DATA_forLc.py #for Lc data



```

