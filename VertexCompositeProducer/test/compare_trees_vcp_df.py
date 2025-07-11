import ROOT

ROOT.gROOT.SetBatch(False)  # interactive mode

# === Input files and tree paths ===
file1_path = "TTree_MC_test_1k.root"
file2_path = "/home/awesole/Dfinder/CMSSW_13_2_11/src/dfinder/dfinder_mc_newFile_numEvent1000.root"
tree_path_1 = "d0ana_newreduced/VertexCompositeNtuple"  # Tree inside file1
tree_path_2 = "Dfinder/ntDkpi"  # Tree inside file2

f1 = ROOT.TFile.Open(file1_path)
f2 = ROOT.TFile.Open(file2_path)

vcp = f1.Get(tree_path_1)
df = f2.Get(tree_path_2)
df.SetMarkerStyle(20)

# === Variable name mapping: ("tree2_variable", "tree1_variable")
matched_vars = [
    ("PvtxX", "PVx"),
    ("PvtxY", "PVy"),
    ("PvtxZ", "PVz"),
    ("BSx", "BSx"),
    ("BSy", "BSy"),
    ("PvtxXErr", "PVxE"),
    ("PvtxYErr", "PVyE"),
    ("PvtxZErr", "PVzE"),
    ("BSxErr", "BSxErr"),
    ("BSyErr", "BSyErr"),
    ("BSzErr", "BSzErr"),
    ("candSize", "Dsize"),
    ("pT", "Dpt"),
    ("y", "Dy"),
    ("phi", "Dphi"),
    ("mass", "Dmass"),
    ("eta", "Deta"),
    ("VtxProb", "Dchi2cl"),
    ("3DPointingAngle", "Dalpha"),
    ("Ddca", "Ddca"),
    ("3DDecayLengthSignificance", "Ddxyz / DdxyzErr"),
    ("3DDecayLength", "Ddxyz"),
    ("3DDecayLengthError", "DdxyzErr"),
    ("2DDecayLengthSignificance", "DsvpvDistance_2D / DsvpvDisErr_2D"),
    ("2DDecayLength", "DsvpvDistance_2D"),
    ("2DDecayLengthError", "DsvpvDisErr_2D"),
    ("DlxyBS", "DlxyBS"),
    ("DlxyBSErr", "DlxyBSErr"),
    ("pTD1 + pTD2", "Dtrk1Pt + Dtrk2Pt"),
    ("pTerrD1 + pTerrD2", "Dtrk1PtErr + Dtrk2PtErr"),
    ("EtaD1 + EtaD2", "Dtrk1Eta + Dtrk2Eta"),
    ("PhiD1 + PhiD2", "Dtrk1Phi + Dtrk2Phi"),
    ("Dtrk1Dz1 + Dtrk2Dz1", "Dtrk1Dz1 + Dtrk2Dz1"),
    ("Dtrk1DzError1 + Dtrk2DzError1", "Dtrk1DzError1 + Dtrk2DzError1"),
    ("Dtrk1Dxy1 + Dtrk2Dxy1", "Dtrk1Dxy1 + Dtrk2Dxy1"),
    ("Dtrk1DxyError1 + Dtrk2DxyError1", "Dtrk1DxyError1 + Dtrk2DxyError1")
    # Add more pairs if needed

]

canvas = ROOT.TCanvas("c", "", 800, 600)

for vcp_var, df_var in matched_vars:
    print(f"📊 Plotting vcp: {vcp_var} vs. df: {df_var}")

    canvas.Clear()
    # Draw tree2 first with automatic histogram
    vcp.Draw(vcp_var,"matchGEN==1 && isSwap==0")

    # Overlay df tree with marker style
    df.Draw(df_var, "Dgen==23333", "p same")

    canvas.Update()
    input("Press Enter to continue to next plot...")

# Cleanup
f1.Close()
f2.Close()

