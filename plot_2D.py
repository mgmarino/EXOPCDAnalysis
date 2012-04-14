import cPickle as pick
import sys
import ROOT

def load_and_output(all_stuff, first):
    adict = {}
    for hists, num in all_stuff:
        adict[num] = hists
    
    keys = adict.keys()
    keys.sort()
    
    temp = adict[keys[0]][0]
    new_hist_ms = ROOT.TH2D("final_hist_ms", "MS", 4000, 0, 4000, len(keys), keys[0], keys[-1])  
    new_hist_ratio = ROOT.TH1D("final_hist_ratio", "RATIO", len(keys), keys[0], keys[-1])  
    
    new_hist_ss = new_hist_ms.Clone("final_hist_ss")
    new_hist_ss.SetTitle("SS")
    
    for key, hists in adict.items():
        ss, ms = hists
        for i in range(0,ss.GetNbinsX()+1):
            new_hist_ss.Fill(ss.GetBinCenter(i), key, ss.GetBinContent(i))
        for i in range(0,ms.GetNbinsX()+1):
            new_hist_ms.Fill(ms.GetBinCenter(i), key, ms.GetBinContent(i))
    new_file = ROOT.TFile(first + "out.root", "recreate")
    new_file.cd()
    new_hist_ss.Write()
    new_hist_ms.Write()
    new_file.Close()
    return 

    ROOT.gStyle.SetOptStat(0)
    canv = ROOT.TCanvas()
    canv.SetLogy()
    canv.SetGrid(1,1)
    ssName = first + "ss.ps"
    msName = first + "ms.ps"
    for num, aname, amax in [(0, ssName, 1e4), (1, msName, 1e5)]:
        canv.cd() 
        canv.Print(aname + "[")
        for key in keys:
            hist = adict[key][num] 
            hist.SetMaximum(amax)
            hist.SetMinimum(1)
            hist.SetLineColor(ROOT.kBlue)
            hist.SetLineWidth(2)
            hist.Rebin(8)
            hist.SetXTitle("Energy (keV)")
            hist.SetTitle("Cluster Radius: %g" % key)
            hist.Draw()
            canv.Update()
            canv.Print(aname)
        canv.Print(aname + "]")
    
    bothtogether = first + "both.ps"
    for aname, amax in [(bothtogether, 1e5)]:
        canv.cd() 
        canv.Print(aname + "[")
        for key in keys:
            ss, ms = adict[key]
            ss.SetLineColor(ROOT.kRed)
            ms.SetLineColor(ROOT.kBlue)
            leg = ROOT.TLegend(0.55, 0.7, 0.8, 0.85)
            leg.SetFillColor(0)
            leg.AddEntry(ss, "SS")
            leg.AddEntry(ms, "MS")
            same = ""
            new_hist_ratio.SetBinContent(new_hist_ratio.FindBin(key), ss.Integral()/float(ms.Integral()))
            for hist in (ms, ss):
                hist.SetMaximum(amax)
                hist.SetMinimum(1)
                hist.SetXTitle("Energy (keV)")
                hist.SetTitle("Cluster Radius: %g" % key)
                hist.Draw(same)
                same = "same"
             
            leg.Draw()
            canv.Update()
            canv.Print(aname)
        canv.Print(aname + "]")
    
    
    new_hist_ratio.SetLineColor(ROOT.kBlue)
    new_hist_ratio.SetXTitle("Cluster Radius (mm)")
    new_hist_ratio.SetYTitle("SS/MS ratio")
    new_hist_ratio.Draw()
    canv.Update()
    canv.Print(first + "ratio.ps")
    
if __name__ == "__main__":
    ROOT.gROOT.SetBatch()
    for afile in sys.argv[1:]:
        load_and_output(pick.load(afile), afile)
