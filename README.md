* -----This code is written to process the 2023 PbPb Data----------
* --The output is a sychronized edm.root file and TTree.root file--
* ----------the edm and root files are one-to-one------------------
* -----------Authors: A.Wesolek, N. Saha, J. Lee-------------------

- abigail.leigh.wesolek@cern.ch
- nihar.ranjan.saha@cern.ch
- junseok.lee@cern.ch
 



This code was originally forked from:
https://github.com/vince502/VertexCompositeAnalysis
To reiterate, saving each of these is optional and can be modified in the configure file.
    -branch = 13_2_X
It has since been heavily edited and simplified.
For complete analysis, the HeavyIonsAnalysis package has been added. This is necessary for items such as Centrality.
Resonace decay reconstruction algorithms with ```pat::CompositeCandiate``` collection in CMSSW 13_2_11. The package is fully orthogonal to the ```HiForest``` framework to be combined with other objects.

-----------------------------------------------------------


At present, this branch supports 2 decay channels
(1.) D0 -> k,pi 
(2.) Lc -> p,k,pi 


-----------------------------------------------------------

The D0 channel was mainly edited by Abigail Wesolek and occurs in 3 steps.

- **Step 1 (VertexCompositeProducer/src/D0Fitter.cc) "Fitter"**

In short, this code first loops over all raw tracks inputted from the miniAOD data file and applies preliminary cuts to remove any tracks that are not of high quality.  Then, it loops over 2 oppositely charged tracks in a double loop, and assigns each one a PdgId and corresponding mass, these are now considered daughter tracks.

We perform 2 iterations of assigning PdgId and mass to the daughter candidiates so that we reconstruct both D0 and D0bar candidates, in turn this also creates the 'swap' component of the mass peak.
- First iteration:  (d1 = pi+, d2 = k-) (D0 candidate)
- Second iteration: (d1 = k+, d2 = pi-) (D0bar candidate)

We then reconstruct 'theD0' from the assigned daughter tracsk. The reconstructed D0 meson is of the type ```pat::CompositeCandidate```, this was a careful decision because of the capability to 'addUserFloat' to each candidate. These floats stay with each candidate throughout the process.  The daughters are of the type ```pat::PackedCandidates``` in accordance with the 2023 PbPb miniAOD data.

The output from D0Fitter.cc can be saved to an edm.root file, as written the output will be ``` "vector<pat::CompositeCandidate>      "generalD0CandidatesNew"   "D0" ```
It is important to note that the output from D0Fitter.cc is all D0Candidates, there are no cuts applies to the reconstrcted tracks in this step.
See more information in the config file (/VertexCompositeProducer/test/run_edm_and_ttree_DATA_forD0.py)

- **Step 2 (VertexCompositeAnalyzer/plugins/VCSelector_D02kpi.cc) "Selector"**

This code reads in the reconstructed tracks from Step 1 (mentioned above).  Although that can be adjusted in the configure file as well.
This code loops over all the D0 candidiates and applies a plethora of cuts. 
At present, most of the cuts are extremely loose and do not have an effect on the tracks.  If one wishes to apply a cut, they will need to apply a specific cut value in the config file.

We then apply 2 forms of BDT Traininng to the D0 candidates:
    - a.) TMVA BDT Training    -- this was implemented by Abigail Wesolek
    - b.) XGBoost BDT Training -- this was implemented by Junseok Lee

We then apply TMVA cuts and XGBoost cuts, these are implemented with an OR. Meaning, so long as the candidate passes one of the cuts, it is kept.

The output of this code consists of 3 collections all of which are stored in the edm.root file
  1. the selected D0 candidates after cuts ```  vector<pat::CompositeCandidate>      "d0selectorNewReduced"     "D0"              "ANASKIM"```
  2. the TMVA BDT values for each D0 candidate ```  vector<float>                        "d0selectorNewReduced"     "MVAValuesNewD0"   "ANASKIM"```
  3. the XGBoost BDT values for each D0 candidate```   vector<float>                        "d0selectorNewReduced"     "MVAValuesNewD02"   "ANASKIM"```

- To reiterate, saving each of these is optional and can be modified in the configure file.
The main reason to save the edm file is so that other colllaborators do not have to reproduce them on their own. 

- **Step 3 (VertexCompositeAnalyzer/plugins/VCTreeProducer_D02kpi.cc) "TreeProducer"**

Finally, this code reads in the D0Candidates that passed step 2.  It creates a TTree and populates it with many variables for each D0 candidiate.

** Each of these steps can be turned on/ff in the configure file, as well as can writing the output of each step **

-----------------------------------------------------------------


The Lc decay channel was built on top of the D0 branch by Nihar Saha.
- It consisits of the same 3 steps as above with the same general logic.  The details do differ; however this is mainlly because the Lc decay has 3 daughters and thus must begin with a triple loop.
- In due time, this will be updated with a link to more thorough notes written by Nihar Saha 
- Link will be posted here.



- This code is updated for miniAOD input. i.e.) using packedPFCandidates and OfflineSlimmedPrimaryVertices instead of GeneralTracks and OfflinePrimaryVertices.
- see the notesOfAllChanges.txt for a mostly comprehensive list of the changes and updates made to this code from the original VCAnalyzer, some of them may be missing.



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

