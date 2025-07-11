import ROOT

# 🎯 Input files and tree path
file1_path = "TTree_wBDT_abbycheck.root"
file2_path = "/home/awesole/production_clone/CMSSW_13_2_11/src/VertexCompositeAnalysis/VertexCompositeProducer/test/TTree_wBDT_abbycheck.root"
tree_path = "d0ana_newreduced/VertexCompositeNtuple"  # Update with your actual nested path

# 📂 Open ROOT files
f1 = ROOT.TFile.Open(file1_path)
f2 = ROOT.TFile.Open(file2_path)

t1 = f1.Get(tree_path)
t2 = f2.Get(tree_path)

# 🧪 Get all branch names (variables)
branches = [b.GetName() for b in t1.GetListOfBranches()]

# 🖼️ Output canvas container
c = ROOT.TCanvas("c", "", 800, 600)

for var in branches:
    print(f"🔍 Plotting: {var}")
    
    # Draw histograms
    h1 = ROOT.TH1F("h1", f"{var}", 100, t1.GetMinimum(var), t1.GetMaximum(var))
    h2 = ROOT.TH1F("h2", f"{var}", 100, t2.GetMinimum(var), t2.GetMaximum(var))
    h1.SetMinimum(0.0)
    h2.SetMinimum(0.0)

    
    t1.Draw(f"{var} >> h1")
    t2.Draw(f"{var} >> h2")

    # Style tree2
    h2.SetMarkerStyle(20)
    

    # Plot together
    h1.Draw("hist")
    h2.Draw("hist p same")


    c.Update()
    input("Press Enter to continue to next plot...")

# ✅ Cleanup
f1.Close()
f2.Close()


