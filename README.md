---------------------------------------------------------------------------------------------------------------
<p align="center"> --------This code is written to process the 2023 PbPb Data-------------</p>
<p align="center"> ----The output is a sychronized edm.root file and TTree.root file----</p>
<p align="center"> ----------------The edm and root files are one-to-one------------------------</p>
<p align="center"> ------------------Authors: A.Wesolek, N. Saha, J. Lee-------------------------- </p>

 ---------------------------------------------------------------------------------------------------------------- 
<p align="center"> abigail.leigh.wesolek@cern.ch </p>
<p align="center"> nihar.ranjan.saha@cern.ch, junseok.lee@cern.ch </p>

----------------------------------------------------------------------------------------------------------------
  

  <br>
 <br>
 <br>
 
This code was originally forked from:
https://github.com/vince502/VertexCompositeAnalysis
 - branch = 13_2_X
- It has since been heavily edited and simplified.
- This code is updated for miniAOD input. i.e.) using ```packedPFCandidates``` and ```OfflineSlimmedPrimaryVertices``` instead of ```GeneralTracks``` and ```OfflinePrimaryVertices```
- see the notesOfAllChanges.txt for a mostly comprehensive list of the changes and updates made to this code from the original VCAnalyzer, some of them may be missing.
- For complete analysis, the HeavyIonsAnalysis package has been added. This is necessary for items such as Centrality.
    
- Resonace decay reconstruction algorithms with ```pat::CompositeCandiate``` collection in CMSSW 13_2_11. The package is fully orthogonal to the ```HiForest``` framework to be combined with other objects.

-----------------------------------------------------------

This branch of the code is used for reproducing the TTrees to include some missing variables. It was written by Nihar Saha. 
  It reads in the previously produced edm.root files and the correctly associated miniAOD file and produces the d0Analyzer tree, event info tree and zdc tree.  



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

