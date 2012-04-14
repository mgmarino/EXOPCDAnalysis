import cPickle as pick
import sys
import ROOT

def load_and_output(first):
   
    open_file = ROOT.TFile(first)
    new_hist_ms = open_file.Get("final_hist_ms")
    new_hist_ss = open_file.Get("final_hist_ss")
    
    ROOT.gStyle.SetOptStat(0)
    canv = ROOT.TCanvas()
    canv.SetLogy()
    canv.SetGrid(1,1)

    for amax in [1e5]:
        canv.cd() 
        print new_hist_ms.GetNbinsY()
        for key in range(new_hist_ms.GetNbinsY()):
            ss, ms = new_hist_ss.ProjectionX("SS" + str(key), key, key + 1), new_hist_ms.ProjectionX("MS" + str(key), key, key + 1) 
            ss.Rebin(16)
            ms.Rebin(16)
            ss.SetLineColor(ROOT.kRed)
            ms.SetLineColor(ROOT.kBlue)
            leg = ROOT.TLegend(0.55, 0.7, 0.8, 0.85)
            leg.SetFillColor(0)
            leg.AddEntry(ss, "SS")
            leg.AddEntry(ms, "MS")
            same = ""
            for hist in (ms, ss):
                hist.SetMaximum(amax)
                hist.SetMinimum(1)
                hist.SetXTitle("Energy (keV)")
                hist.SetTitle("Cluster Radius: %g" % new_hist_ms.GetYaxis().GetBinCenter(key))
                hist.Draw(same)
                same = "same"
             
            leg.Draw()
            canv.Update()
            raw_input("E")
    
   
if __name__ == "__main__":
    for afile in sys.argv[1:]:
        load_and_output(afile)
