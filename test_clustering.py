import ROOT
from make_pickle_from_ROOT import grab_dict_from_afiles
from make_plots_for_pkl import GenClusters
import sys
import numpy


if __name__ == '__main__':
    ROOT.gSystem.Load("libEXOUtilities")
    ROOT.gSystem.Load("/Users/mgmarino/software/EXO_Fitting/EXO_Fitting/lib/libEXOFitting")
 
    culler = ROOT.EXOClusterCull()
    #adict = grab_dict_from_afiles(sys.argv[1:])
    afile = ROOT.TFile(sys.argv[1])  
    tree = afile.Get("tree")
    #cluster = GenClusters(adict.items(), 5.0, "temp")
    cluster = GenClusters(None, 5.0, "temp")
    
    for i in range(12, tree.GetEntries()):
        print "Entry: ", i
        tree.GetEntry(i)
        ed = tree.EventBranch
        mc = ed.fMonteCarloData
        if ed.GetNumScintillationClusters() != 1: continue
        temp_list = []
        for j in range(mc.GetNumPixelatedChargeDeposits()):
            pcd = mc.GetPixelatedChargeDeposit(j)
            cen = pcd.GetPixelCenter()
            temp_list.append((numpy.array([cen.GetX(), cen.GetY(), cen.GetZ()]), pcd.fTotalEnergy))
            
        for aval in temp_list: print aval 
        clustered = cluster.get_energy_in_ms_ss(temp_list, 3.0) 
        print "Clustered"
        for aval in temp_list: print aval 
        #print len(temp_list), temp_list 
        print clustered[0], sum(clustered[1])
        for j in range(ed.GetNumScintillationClusters()):
            energy = culler.ClusterCull(ed.GetScintillationCluster(j), True, True, ed)
            print "Energy: ", energy, "Sites: ", culler.GetNumSite(), "Missing Pos: ", culler.HasMissingPositions() 
            print "Fiducial: ", culler.IsFiducial()
        #for j in range(ed.GetNumChargeClusters()):
        #    ed.GetChargeCluster(j).Print()
        raw_input("E")
    

