import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from typing import List
from itertools import chain
import re
from ranges_radii_dense import get_ranges_radii
regex = re.compile(r'\d+')

particle_map = {0: "At211",
      1:  "Po211",
      2: "Po211*", 
      3: "Bi207",
      4: "Pb207",
      5: "Pb207*",
      6:"alphaAt211",
      7: "alphaPo211",
      8: "e-At211",
      9: "e-Bi207",
      10: "e-Pb207*",
      11: "gammaAt211",
      12: "gammaBi207",
      13: "gammaPb207",
      14: "gammaPb207*",
      15: "gammaPo211",
      16: "gammaPo211*",
   -1: "all"}


class ClusteringEvents:
    def __init__(self, infile):
        self.cl_file = infile
        self.header = pd.read_csv(self.cl_file, nrows=1) 
        self.data = pd.read_csv(self.cl_file, skiprows=2)
        self.data['ParticleName'] = list(map(lambda pid: particle_map[pid], self.data['ParticleID']))
    @property
    def events(self):
        return self.data['EventNo']
    @property
    def particle_ids(self):
        return self.data['ParticleID']
    @property
    def particle_avail_id(self):
        return np.unique(self.data['ParticleID'])
    @property
    def particle_avail(self):
        return np.unique(self.data['ParticleName'])
    @property
    def numBP(self):
        return float(self.header['numBP'])
    @property
    def numEvts(self):
        return int(self.header['numEvtIntersectingVolume'])
    
    def num_evts_pid(self, pid):
        # print(self.data.loc[self.data['ParticleID']==pid])
        num_evts = np.unique(self.data.loc[self.data['ParticleID']==pid]['NumEvts'])
        if len(num_evts):
            return num_evts[0]
        else:
            return np.nan
    def per_particle(self, pid, field_key):
        return self.data.loc[self.data["ParticleID"]==pid][field_key]
    
    def per_particlename(self, name, field_key):
        return self.data.loc[self.data["ParticleName"]==name][field_key]
        
    def get_damage(self, field_key):
        return self.data[field_key]
         

class ClusteringEventsRadii:
    def __init__(self, folder, fname_prefix, spacing=None, n_div_r=40, seed=None, keyword=None, particle=None, boxes_per_R = None, start_r = 10.5):
        self.datasets = {}
        self.n_r = n_div_r
        self.radii = range(n_div_r)        
        self.fname_base = f"{fname_prefix}"
        
        if spacing:
            self.fname_base += f"_{spacing}um"
        if seed: 
            self.fname_base+= f"_{seed}"
        if keyword:
            self.fname_base+=f'_{keyword}'
        if particle:
            self.fname_base+=f"_{particle}"
        if boxes_per_R:
            self.boxes_per_R = np.array(boxes_per_R)

        self.fnames_radii = {r: os.path.join(folder, self.fname_base+f"_{r}.csv") for r in self.radii if self.fname_base+f"_{r}.csv" in os.listdir(folder)}
        self.fnames = list(self.fnames_radii.values())
        
        if not len(self.fnames_radii):
            print(f"No data loaded for seed {seed}, particle {particle}")
            return None

        
        self.distance = np.zeros(self.n_r)
        for i in self.radii:
            if i in self.fnames_radii.keys():
                self.datasets[i] = ClusteringEvents(self.fnames_radii[i])
            
            if spacing:
                self.distance[i] = start_r+spacing* i
            else:
                self.distance[i] = start_r+2*i

    @property
    def particle_avail(self):
        return set(chain.from_iterable(dataset.particle_avail for dataset in self.datasets.values()))

    @property
    def particle_avail_id(self):
        return set(chain.from_iterable(dataset.particle_avail_id for dataset in self.datasets.values()))
    
    def get(self, damage_key):
        array = np.full(self.n_r, np.nan)
        
        for i, dataset in self.datasets.items():
            array[i] = np.array(dataset.get_damage(damage_key).sum())
        return array
    
    # def n_events(self, damage_key):
    #     nevents_r = np.full(self.n_r, np.nan)
    #     for i, dataset in self.datasets.items():
    #         _data = dataset.get_damage(damage_key)
    #         nevents_r[i] = len(_data[~np.isnan(_data)])
    #     return nevents_r

    # def n_events_pid(self, pid, damage_key):
    #     nevents_pid_r = np.full(self.n_r, np.nan)
    #     for i, dataset in self.datasets.items():
    #         _data = dataset.per_particlename(particle_map[pid], damage_key)
    #         nevents_pid_r[i] = len(_data[~np.isnan(_data)])
    #     return nevents_pid_r
    

    # def get_err(self, damage_key)

    def get_bydose(self, damage_key):
        array = np.full(self.n_r, np.nan)
        dose_arr = np.full(self.n_r, np.nan)
        for i, dataset in self.datasets.items():
            array[i] = np.array(dataset.get_damage(damage_key).sum())
            dose_arr[i] = np.array(dataset.get_damage("DoseGy").sum())
        if damage_key == 'DoseGy':
            return dose_arr
        else:
            arr_bydose = array/dose_arr
            #        np.divide(array, dose_arr, out=np.zeros_like(array), where=dose_arr)
            return arr_bydose
    
    def byparticle(self, pid, damage_key):
        array = np.full(self.n_r, np.nan)
        for i, dataset in self.datasets.items():
            array[i] = np.array(dataset.per_particlename(particle_map[pid], damage_key).sum())
        return array
    

    def byparticle_bydose(self, pid, damage_key):
        if damage_key=='DoseGy':
            return self.byparticle(pid, "DoseGy")
        else:
            damage = self.byparticle(pid, damage_key)
            dose = self.byparticle(pid, "DoseGy")
            return(damage/dose)
#            return np.divide(damage,dose, out=np.zeros_like(self.distance), where=dose)
    def num_evts_pid(self, pid):
        num_evts_pid_per_dataset = np.full(self.n_r, np.nan)
        for i, dataset in self.datasets.items():
            num_evts_pid_per_dataset[i] = dataset.num_evts_pid(pid)
        return num_evts_pid_per_dataset
    
    @property
    def numBP(self):
        numbp_per_dataset = np.unique([self.datasets[i].numBP for i in self.datasets])
        assert(len(numbp_per_dataset)==1)
        return float(numbp_per_dataset)
    
    @property
    def num_bpr(self):
        return self.boxes_per_R

    @property
    def numEvtsInters(self):
        num_evts_per_dataset = np.full(self.n_r, np.nan)
        for i, dataset in self.datasets.items():
            num_evts_per_dataset[i] = dataset.numEvts
        return num_evts_per_dataset
    
    def __str__(self) -> str:
        return f"fname: {self.fname_base}"

    def __repr__(self) -> str:
        return str(self)


particle_plotting_map = {
    -1: {"name": "Total",
       "colour": "mediumaquamarine",
       'marker': "s"},
    6: {"name": "alphaAt211",
        "colour": "darkblue",
        'marker': "*"},
    7: {"name": "alphaPo211",
        "colour": "mediumvioletred",
        'marker': "s"},
    8: {"name": "e-At211",
        "colour": "cornflowerblue",
        'marker': "x"},
    9: {"name": "e-Bi207",
        "colour": "cornflowerblue",
        'marker': "o"},
    10: {"name": "e-Pb207*",
       "colour": "cornflowerblue",
       'marker': "s"},
    11: {"name": "gammaAt211",
        "colour": "darkred",
        'marker': "x"},
    12: {"name": "gammaBi207",
        "colour": "darkred",
        'marker': "o"},
    14: {"name": "gammaPb207",
    "colour": "darkred",
    'marker': "s"},
    4: {"name": "Pb207",
       "colour": "gold",
       'marker': "s"},
    
    }
    
def plot_field(datadict, damage_type='DSBtotal'):
    # colors = ['red', "k", "mediumvioletred", "cornflowerblue", 'gold']
    damage = {}
    distance = {}
    for pmap in particle_plotting_map.keys():
        damage.setdefault(pmap, [])
        distance.setdefault(pmap, [])
    datasets_all = datadict['all']
    for p in datadict.keys():
            datasets = datadict[p]
            for seed_n, (dataset, dataset_all) in enumerate(zip(datasets, datasets_all)):
                for pid in dataset.particle_avail_id:
                    distance[pid].append(dataset.distance)
            
                    if damage_type=="DoseGy":
                        dam = dataset.byparticle(pid, damage_type)#*(dataset.n_events(damage_type)/1000)
                        
                    else:
                        # dam = dataset.get(damage_type)
                        damage_temp = dataset.byparticle(pid, damage_type)  
                        # damage_temp = np.where(dataset.n_events(damage_type)>2, np.nan, damage_temp)
                        dam = (damage_temp/dataset_all.get("DoseGy"))/(dataset.numBP*1e-9)
                            #np.divide(damage_temp, dataset_all.get("DoseGy"),#*(dataset_all.n_events('DoseGy')/1000), out=np.zeros_like(dataset.get(damage_type)), where=dataset_all.get("DoseGy")!=0))/(dataset.numBP*1e-9)
                        dam = np.where(dam==0., np.nan, dam)
                        
                    damage[pid].append(dam)
                    
                    
    for pid in damage.keys():
        
        std = np.nanstd(np.array(damage[pid]),axis=0)
        std = np.where(std==0., np.nan, std)

        dist = np.nanmean(np.asarray(distance[pid]), axis=0)
        mean = np.nanmean(np.array(damage[pid]), axis=0)
        if len(mean[~np.isnan(mean)])>1:
            plt.errorbar(dist, mean, yerr =std,
                    linestyle='',  alpha = 1 if pid!=-1 else 0.3,
                    capsize=1, marker=particle_plotting_map[pid]['marker'], 
                    markersize=3 if pid!=-1 else 8, linewidth=0.1, elinewidth=0.1,label=f"{particle_plotting_map[pid]['name']}", c=particle_plotting_map[pid]['colour'])
                
    return


def plot_rbe(datadict, damage_type='DSBtotal'):
    
    # colors = ['red', "k", "mediumvioletred", "cornflowerblue", 'gold']
    damage = {}
    distance = {}
    for pmap in particle_plotting_map.keys():
        damage.setdefault(pmap, [])
        distance.setdefault(pmap, [])
    datasets_all = datadict['all']
    for p in datadict.keys():
            datasets = datadict[p]
            for seed_n, (dataset, dataset_all) in enumerate(zip(datasets, datasets_all)):
                for pid in dataset.particle_avail_id:
                    distance[pid].append(dataset.distance)
            
                    if damage_type=="DoseGy":
                        dam = dataset.byparticle(pid, damage_type)#*(dataset.n_events(damage_type)/1000)
                        
                    else:
                        # dam = dataset.get(damage_type)
                        damage_temp = dataset.byparticle(pid, damage_type)  
                        # damage_temp = np.where(dataset.n_events(damage_type)>2, np.nan, damage_temp)
                        dam = (damage_temp/dataset_all.get("DoseGy"))/(dataset.numBP*1e-9)/6.9
                            #np.divide(damage_temp, dataset_all.get("DoseGy"),#*(dataset_all.n_events('DoseGy')/1000), out=np.zeros_like(dataset.get(damage_type)), where=dataset_all.get("DoseGy")!=0))/(dataset.numBP*1e-9)
                        dam = np.where(dam==0., np.nan, dam)
                        
                    damage[pid].append(dam)
                    
                    
    for pid in damage.keys():
        
        std = np.nanstd(np.array(damage[pid]),axis=0)
        std = np.where(std==0., np.nan, std)

        dist = np.nanmean(np.asarray(distance[pid]), axis=0)
        mean = np.nanmean(np.array(damage[pid]), axis=0)
        if len(mean[~np.isnan(mean)])>1:
            plt.errorbar(dist, mean, yerr =std,
                    linestyle='',  alpha = 1 if pid!=-1 else 0.3,
                    capsize=1, marker=particle_plotting_map[pid]['marker'], 
                    markersize=3 if pid!=-1 else 8, linewidth=0.1, elinewidth=0.1,label=f"{particle_plotting_map[pid]['name']}", c=particle_plotting_map[pid]['colour'])
                
    return
  
if __name__=="__main__": 

    datasets = []
    folder = "/home/cdesio/TAT/tat_shell_ps/output/test_At1k_shell_ps_10.5um_80R_continuous"
    data_folder = os.path.join(folder, "clustering_out")
    fname_prefix="out_AtDNA_1k_spacing"
    spacing= 1
    keyword = 'part'

    seeds = np.unique([int(regex.findall(fname.split(fname_prefix+f"_{spacing}um")[-1])[0]) for fname in os.listdir(data_folder) if keyword in fname])
    particles = np.unique([fname.split(fname_prefix+f"_{spacing}um")[-1].split("_")[-2] for fname in os.listdir(data_folder) if keyword in fname and "DSB" not in fname])

    print(seeds)
    print(particles)

    datadict = {}
    for particle in np.unique(particles):
        datadict.setdefault(particle, [])
    # boxes_per_R = get_ranges_radii(ndiv_R=100)[1]
    # print(boxes_per_R, len(boxes_per_R))
    for i, seed in enumerate(seeds):
        for particle in particles:
            
            datadict[particle].append(ClusteringEventsRadii(folder=data_folder, fname_prefix=fname_prefix, 
                                spacing=spacing, seed=int(seed), keyword=keyword, particle=particle, n_div_r=80, start_r = 0))  

    def plot_damage(damage_type, datadict):
        plt.figure(figsize=(8,5))
        plot_field(datadict, damage_type=damage_type)
        plt.legend(ncol=4, loc=(0,1.01))
        plt.xlabel("radial distance from the blood vessel (um)")
        if damage_type=="DoseGy":
            plt.ylabel(f"Dose ($Gy)$")
            plt.yscale('log')

        elif damage_type=="DSBtotal":
            plt.ylabel(f"total n. of DSB ($Gy^-1 Gbp^-1)$")      
            plt.ylim(-1, 16)
        elif damage_type=="TotalSBtotal":
            plt.ylabel(f"total n. of SB ($Gy^-1 Gbp^-1)$")  
            plt.ylim(-1, 200)

        plt.xlim(0, 80)
        #plt.legend(loc=(1,0.4))
        plt.savefig(f"{folder}/{damage_type}_1k_dense.png")

        return

    for damage in ["DSBtotal", "TotalSBtotal", "DoseGy"]:
        plot_damage(damage, datadict)  

    plt.figure(figsize=(8,5))
    plot_rbe(datadict, damage_type="DSBtotal")
    plt.legend()
    plt.xlabel("radial distance from blood vessel (um)")
    plt.ylabel("RBE")
    plt.xlim(0, 80)
    plt.ylim(0, 3)
    plt.savefig(f"{folder}/RBE_1k_dense.png")
    # damage_type="DSBtotal"

    # for i, seed in enumerate([6731, 7060]):
    #     testfolder = {}
    #     for particle in ['all', 'alpha', 'e-']:
    #         data = ClusteringEventsRadii(folder = "output/test_At1k_160k_boxes_10.5um/clustering_out", fname_prefix="out_AtDNA_1k_spacing",
    #                         spacing=2, seed=seed, keyword="part", particle=particle)
    #         if len(data.fnames_radii):
    #             testfolder[particle] = data
    #     plt.figure(figsize=(8,4))
        
    #     for p in ['all', 'alpha', 'e-']:
    #         if p in testfolder.keys():
    #             damage_arr = testfolder[p].get(damage_type)
    #             dose_arr = testfolder['all'].get("DoseGy")
    #             if np.any(damage_arr!=0):
    #                 plt.plot(testfolder[p].distance, np.divide(
    #                 damage_arr, dose_arr, out=np.zeros_like(damage_arr)/testfolder[p].numBP*1e-9, where=dose_arr != 0), marker='o', label=f'{p}', linewidth=0.3, markersize=4)
    #     #plt.plot(testfolder.distance, testfolder.get(damage_type),  marker='o',label=f'tot', linewidth=0.6, markersize=4)
    #     plt.xlabel("distance (um)")
    #     plt.ylabel(f"{damage_type} ($Gy^-1 Gbp^-1)$")
    #     plt.title(f"DSB : {seed}")
    #     plt.grid()
    #     plt.legend()
    #     plt.savefig(f"{damage_type}_1k_160k_{seed}.png")