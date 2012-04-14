import ROOT
import numpy
import cPickle as pick
import sys

ROOT.gSystem.Load("libEXOUtilities")


def grab_dict_from_afiles(file_list):
    chain = ROOT.TChain("tree")
    for afile in file_list: chain.Add(afile)
    
    chain.SetAlias("pkey", "fMonteCarloData.fPixelatedChargeDeposits.fCoordinateKey.GetCenter()")
    chain.SetAlias("penergy", "fMonteCarloData.fPixelatedChargeDeposits.fTotalEnergy")
    chain.SetEstimate(10*chain.GetEntries())
    print "Extracting from tree"
    entries = chain.Draw("pkey.GetX():pkey.GetY():pkey.GetZ():Entry$", "", "goff")
    x = chain.GetV1()
    y = chain.GetV2()
    z = chain.GetV3()
    entry = chain.GetV4()
    
    print "Putting in list"
    alist = [ (numpy.array([x[i], y[i], z[i]]), int(entry[i])) for i in range(entries) ]
    
    print "Grabbing energy"
    entries = chain.Draw("penergy", "", "goff")
    energy = chain.GetV1()
    energylist = [ energy[i] for i in range(entries) ]
    
    print "Grouping into events"
    aset = set([var[1] for var in alist])
    adict = dict([(i, []) for i in aset])
    map(lambda x: adict[x[0][1]].append((x[0][0], x[1])), zip(alist, energylist))
    return adict


if __name__ == '__main__':
    output_name = sys.argv[1]
    adict = grab_dict_from_afiles(sys.argv[2:])
    print "Saving to disk"
    pick.dump(adict, open(output_name, "w"))
