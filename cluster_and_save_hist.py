import cPickle as pick
import sys
import math
import numpy
import os
import multiprocessing, Queue
import ROOT
import gzip
import re

class GenClusters:
    def __init__(self, alist, rad, name):
        self.proc_list = alist
        self.radius = rad
        self.name = name
        self.energy_calibration = 0.6145

        self.limit =45./1000.
        self.width =15./1000.
        self.limit =50./1000.
        self.width =50./1000.
        self.v_wire_limit = 247.8/1000.
        self.v_wire_width = 67.12/1000.
        self.max_r = 163#ROOT.CATHODE_RADIUS/ROOT.CLHEP.mm 
        #self.max_r = ROOT.ACTIVE_XENON_RADIUS/ROOT.CLHEP.mm 
        self.min_z = 5.0 
        self.max_z = 182.0 
        self.culling_max_r = 173.
        self.volume_cut = ROOT.EXOClusterCull()

    def run(self):
        self.sshist, self.mshist = self.get_hists_for_list(self.proc_list, self.radius)

    def cluster(self, one, two):
        ener_one = one[1]
        ener_two = two[1]
        return ((one[0]*ener_one + two[0]*ener_two)/(ener_one + ener_two), ener_one + ener_two)

    def apply_active_volume_cut(self, event):
        # Returns True to keep
        arr = event[0] 
        return not self.volume_cut.HexagonalCulling(arr[0], arr[1], self.culling_max_r) 

    def apply_threshold(self, an_energy):
        # Returns false if we don't keep this
        #return True
        
        rand_energy = ROOT.gRandom.Rndm()
        return 0.5*(1 + ROOT.TMath.Erf((an_energy - self.limit)/(math.sqrt(2)*self.width))) > rand_energy
    
    def apply_v_wire_threshold(self, an_energy):
        # Returns false if we don't keep this
        rand_energy = ROOT.gRandom.Rndm()
        return 0.5*(1 + ROOT.TMath.Erf((an_energy - self.v_wire_limit)/(math.sqrt(2)*self.v_wire_width))) > rand_energy
        

    def limit_energy(self, events):
        # First apply the threshold that we actually see a cluster.
        remove_events = [] 
        i = 0
        for anev in events:
            if not self.apply_threshold(anev[1]): remove_events.append(i)
            i += 1
        remove_events.reverse()
        for ev in remove_events: events.pop(ev)
        return events

    def get_hists(self):
        return self.sshist, self.mshist
    
    def perform_fiducial_cuts(self, event):
        arr = event[0] 
        r = math.sqrt(arr[0]*arr[0] + arr[1]*arr[1])
        if r > self.max_r: return False
        return self.min_z < math.fabs(arr[2]) and math.fabs(arr[2]) < self.max_z

    def do_clustering(self, event, rad):
        adict = {}

	# We get here, this means everything was gotten rid of, or nothing was
	# in this event, return 0
	if len(event) == 0: return []

	# We get here, it means we have a single-site event, return
	if len(event) == 1: return event 
    
	# Else lets calculate the distance between all clusters
	for i in range(len(event)):
            for j in range(i+1, len(event)):
                arr = (event[i][0] - event[j][0])
                adict[math.sqrt(numpy.dot(arr,arr))] = (i,j) 
    
        rad_keys = adict.keys()
        rad_keys.sort() # should be now sorted least to biggest
        if rad_keys[0] < rad: 
            # we need to cluster this together    
            i, j = adict[rad_keys[0]]
            new_cluster = self.cluster(event[i], event[j])
            if i > j: 
                event.pop(i)
                event.pop(j)
            else:
                event.pop(j)
                event.pop(i)
            event.append(new_cluster)
            return self.do_clustering(event, rad)

	# if we get here it means we have no longer clustered and have a
	# multi=site, return it 
        return event 
 
    def remove_hex_events(self, event):
        # Now cull events outside the hexagonal 
        remove_events = [] 
        i = 0
        for anev in event:
            if not self.apply_active_volume_cut(anev): remove_events.append(i)
            i += 1
        remove_events.reverse()
        for ev in remove_events: event.pop(ev)
        return event


    def cull_and_cut(self, event):
	# This function returns a list of all events that passed the cuts 
	# it gets a list of the charge clusters

	# First limit energies, this removes clusters that are too small and
	# would have not been seen
	event = self.limit_energy(event)
        event = self.remove_hex_events(event)

	# Now we should have something that looks like data, we need to then
	# apply the v-wire thresholds and other cuts
        for anev in event:
            if not self.apply_v_wire_threshold(anev[1]):
                return 0, []


        for anev in event:
            if not self.perform_fiducial_cuts(anev): 
                return 0, []
        # Now return 
        if len(event) == 1: return event[0][1], []
        return 0, [x[1] for x in event]
    
        

    def get_energy_in_ms_ss(self, event, rad):
        event = self.do_clustering(event, rad)
        return self.cull_and_cut(event)
        
   
    def get_hists_for_list(self, alist, cut_radius):
        ss = ROOT.TH1D("ss" + self.name, "ss", 4000, 0, 4000)
        ms = ROOT.TH1D("ms" + self.name, "ms", 4000, 0, 4000)
    
        i = 0
        percentage = -1
        for total_entry in alist:
           if i*10/len(alist) != percentage: 
               percentage = i*10/len(alist)
               print " (%i) Percent done: %i " % (os.getpid(), percentage*10)
           #j = 0
           #remove_events = []
           #for an_entry in total_entry:
           #    if not self.apply_active_volume_cut(an_entry):
           #        remove_events.append(j)
           #    j += 1
           #remove_events.reverse()
           #for j in remove_events: total_entry.pop(j)
           ssenergy, msenergy = self.get_energy_in_ms_ss(total_entry, cut_radius) 
           sum_ms = sum(msenergy)
           if ssenergy != 0.0: ss.Fill(1000*ssenergy)
           if sum_ms != 0.0: ms.Fill(1000*sum_ms)
           #for ener in msenergy: ms.Fill(1000*ener)
           i += 1
        return ss, ms

def split_files(alist):
    new_dict = {}
    anre = re.compile(".*/P2_(.+)_[0-9]+.pkl.gz")
    max_size = 100*1024*1024
    for afile in alist:
        check = anre.match(afile)
        if not check: 
            print "Error, ", afile
            continue
        out_name = check.group(1)
        if out_name not in new_dict: new_dict[out_name] = [0, []]
        file_size = os.path.getsize(afile)
        if file_size + new_dict[out_name][0] > max_size: continue 
        new_dict[out_name][0] += file_size
        new_dict[out_name][1].append(afile)
    return new_dict
        
class FileLoader(multiprocessing.Process):
    def __init__(self, task_list, results_list):
        multiprocessing.Process.__init__(self)
        self.task_list = task_list
        self.results_list = results_list

    def run(self):
        while True: 
            tmp = self.task_list.get() 
            if tmp == None: break
            print "Doing job: ", tmp, os.getpid()
            try:
                res = pick.load(gzip.open(tmp))
            except IOError:
                res = pick.load(open(tmp))
            self.results_list.put(res)
            self.task_list.task_done()
        self.task_list.task_done()
        print "Done", os.getpid()


class ConsRunner(multiprocessing.Process):
    def __init__(self, task_list, mcd_list, results_list):
        multiprocessing.Process.__init__(self)
        self.task_list = task_list
        self.mcd_list = mcd_list
        self.results_list = results_list

    def run(self):
        while True: 
            tmp = self.task_list.get() 
            if tmp == None: break
            print "Doing job: ", tmp, os.getpid()
            new_clus = GenClusters(self.mcd_list, tmp, str(os.getpid()))
            new_clus.run()
            self.results_list.put(pick.dumps((new_clus.get_hists(), tmp)))
            self.task_list.task_done()
        self.task_list.task_done()
        print "Done", os.getpid()



if __name__ == '__main__':
    ROOT.gROOT.SetBatch()
    all_files = split_files(sys.argv[1:])
    for name, val in all_files.items():
        print "%s has %i items, with a size of %i MB" % \
          (name, len(val[1]), val[0]/(1024*1024))
    for name, val in all_files.items():
        mc_data = []
        output_file_name = name + ".pkl"
        if os.path.exists(output_file_name + "out.root"):
            print "Exists: ", output_file_name + "out.root"
            print "Skipping..."
            continue

        num_cpus = 16 
        q2 = multiprocessing.JoinableQueue()
        results = multiprocessing.Queue()
        children = []
        for cpu in range(num_cpus):
            p = FileLoader(q2, results)
            p.start()
            children.append(p)

        for afile in val[1]: 
            q2.put(afile) 

        for cpu in range(num_cpus): q2.put(None)

        q2.join()
        print "Consolidating data..."
        for _ in val[1]: 
            print "Pulling in:", _
            res = results.get()
            mc_data.extend(res) 
        
        print "Cleaning children..."
        for child in children:
            child.join()

        print "Beginning processing..."
        children = []
        radius_list = [0.1*i for i in range(10,100)]
        
        q = multiprocessing.JoinableQueue()
        results = multiprocessing.Queue()
        
        for cpu in range(num_cpus):
            p = ConsRunner(q, mc_data, results)
            p.start()
            children.append(p)
        
        for radius in radius_list: q.put(radius)
        for cpu in range(num_cpus): q.put(None)
        
        results_list = []
        q.join()
        for _ in radius_list: 
            res = results.get()
            results_list.append(pick.loads(res))
        
        for child in children:
            child.join()
        from plot_2D import load_and_output
        load_and_output(results_list, output_file_name)
        del results_list
        del mc_data
