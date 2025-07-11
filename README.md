
Resonace decay reconstruction algorithms with ```VertexCompositeCandiate``` collection in cmssw. Compatible with 2023 PbPb datafomat. The package is fully orthogonal to the ```HiForest``` framework to be combined with other objects.

This branch was edited by Abigail Wesolek and currently will only process events of 

$D^{0} \to K+\pi$

But other branches can easily be added in the future.  You want to check the origianl branch for more information.

The code been edited to produce skimmedEDM and TTree simultaneously.  
The code is updated for miniAOD input. i.e.) using packedPFCandidates and OfflineSlimmedPrimaryVertices instead of GeneralTracks and OfflinePrimaryVertices.
VertexCompositeProducer/test/run_edm_TTree.py






## How to run

```bash 
#LXplus, bash, cmssw-el8 apptainer

cmsrel CMSSW_13_2_11

cd CMSSW_13_2_11/src
cmsenv

#clone the repo
git clone git@github.com:a-wesole/VertexCP_miniAOD.git VertexCompositeAnalysis

#run the setup.sh script this will add the HeavyIonsAnalysis from CmsHI github that is needed for centrality 
cd VertexCompositeAnalysis
./setup.sh

#compile
cd ..
scram b -j8

#cd and run code that produces edm and TTree files 
cd VertexCompositeAnalysis/VertexCompositeProducer/test
cmsRun run_edm_and_ttree_DATA.py 


```

