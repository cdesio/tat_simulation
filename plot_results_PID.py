import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from typing import List
from itertools import chain
import re
from ranges_radii_dense import get_ranges_radii
regex = re.compile(r'\d+')

particle_map = {1: "alpha",
                2: "gamma",
                3: "e-",
                4: "nu_e",
                5: "At211",
                6: "Po211",
                7: "Bi207",
                8: "Pb207",
                9: "proton",
                10: "e+" ,
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

    def per_particle(self, pid, field_key):
        return self.data.loc[self.data["ParticleID"]==pid][field_key]
    
    def per_particlename(self, name, field_key):
        return self.data.loc[self.data["ParticleName"]==name][field_key]
        
    def get_damage(self, field_key):
        return self.data[field_key]
         

class ClusteringEventsRadii:
    def __init__(self, folder, fname_prefix, spacing=None, n_div_r=40, seed=None, keyword=None, particle=None, boxes_per_R = None):
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
                self.distance[i] = 10.5+spacing* i
            else:
                self.distance[i] = 10.5+2*i

    @property
    def particle_avail(self):
        return set(chain.from_iterable(dataset.particle_avail for dataset in self.datasets.values()))
    
    def get(self, damage_key):
        array = np.zeros(self.n_r)
        # std = np.zeros(self.n_r)
        for i, dataset in self.datasets.items():
            array[i] = np.array(dataset.get_damage(damage_key).sum())
            # std[i] = np.array(dataset.get_damage(damage_key).std())
        
        return array #, std

    # def get_err(self, damage_key)

    def get_bydose(self, damage_key):
        array = np.zeros(self.n_r)
        dose_arr = np.zeros(self.n_r)
        for i, dataset in self.datasets.items():
            array[i] = np.array(dataset.get_damage(damage_key).sum())
            dose_arr[i] = np.array(dataset.get_damage("DoseGy").sum())
            
        arr_bydose = np.divide(array, dose_arr, out=np.zeros_like(array), where=dose_arr != 0)
            
        return arr_bydose
    
    def byparticle(self, pid, damage_key):
        array = np.zeros(self.n_r)
        for i, dataset in self.datasets.items():
            array[i] = np.array(dataset.per_particlename(pid, damage_key).sum())
        if damage_key=='DoseGy':
            return array
        else:
            return array/self.numBP

    def byparticle_bydose(self, pid, damage_key):
        if damage_key=='DoseGy':
            return self.byparticle(pid, "DoseGy")
        else:
            damage = self.byparticle(pid, damage_key)
            dose = self.byparticle(pid, "DoseGy")
            return np.divide(damage,dose, out=np.zeros_like(self.distance), where=dose!=0)
    
    @property
    def numBP(self):
        numbp_per_dataset = np.unique([self.datasets[i].numBP for i in self.datasets])
        assert(len(numbp_per_dataset)==1)
        return float(numbp_per_dataset)
    
    @property
    def num_bpr(self):
        return self.boxes_per_R

    
    def __str__(self) -> str:
        return f"fname: {self.fname_base}"

    def __repr__(self) -> str:
        return str(self)

    
def plot_field(datadict, damage_type='DSBtotal'):
    colors = ['red', "k", "mediumvioletred", "cornflowerblue"]
    
    for i, p in enumerate(datadict.keys()):
        distance = []
        damage = []
        datasets = datadict[p]
        datasets_all = datadict['all']
        for _, (dataset, dataset_all) in enumerate(zip(datasets, datasets_all)):
            if len(dataset.particle_avail):
                distance.append(dataset.distance)
                if damage_type=="DoseGy":
                    dam = dataset.get(damage_type)/dataset.num_bpr
                else:
                    dam = (np.divide(dataset.get(damage_type), dataset_all.get("DoseGy"), 
                                     out=np.zeros_like(dataset.get(damage_type)), where=dataset_all.get("DoseGy")!=0))/(dataset.numBP*1e-9)
                dam = np.where(dam==0., np.nan, dam)
                damage.append(dam)
        
        if len(distance) and len(damage):
            plt.errorbar(np.mean(np.asarray(distance), axis=0), 
                        np.nanmean(np.asarray(damage) ,axis=0), 
                        yerr = np.nanstd(np.asarray(damage),axis=0), capsize=2, marker='o', 
                        markersize=4, linewidth=0.2, elinewidth=0.2,label=f"{p}", c=colors[i])
    return


def plot_rbe(datadict, damage_type='DSBtotal'):
    colors = ['red', "k", "mediumvioletred", "cornflowerblue"]
    
    for i, p in enumerate(datadict.keys()):
        distance = []
        damage = []
        datasets = datadict[p]
        datasets_all = datadict['all']
        for _, (dataset, dataset_all) in enumerate(zip(datasets, datasets_all)):
            if len(dataset.particle_avail):
                distance.append(dataset.distance)
            
                dam = (np.divide(dataset.get(damage_type), dataset_all.get("DoseGy"), 
                                out=np.zeros_like(dataset.get(damage_type)), where=dataset_all.get("DoseGy")!=0))/(dataset.numBP*1e-9)/6.9
                dam = np.where(dam==0., np.nan, dam)
                damage.append(dam)
        
        if len(distance) and len(damage):
            plt.errorbar(np.mean(np.asarray(distance), axis=0), 
                        np.nanmean(np.asarray(damage) ,axis=0), 
                        yerr = np.nanstd(np.asarray(damage),axis=0), capsize=2, marker='o', 
                        markersize=4, linewidth=0.2, elinewidth=0.2,label=f"{p}", c=colors[i])
    return
  
if __name__=="__main__": 

    datasets = []
    folder = "/home/cdesio/TAT/tat_dense/output/combined_all"
    #data_folder = os.path.join(folder, "clustering_out")
    fname_prefix="out_AtDNA_1k_spacing"
    spacing= 1
    keyword = 'part'

    seeds = np.unique([int(regex.findall(fname.split(fname_prefix+f"_{spacing}um")[-1])[0]) for fname in os.listdir(folder) if keyword in fname])
    particles = np.unique([fname.split(fname_prefix+f"_{spacing}um")[-1].split("_")[-2] for fname in os.listdir(folder) if keyword in fname and "DSB" not in fname])

    print(seeds)
    print(particles)

    datadict = {}
    for particle in np.unique(particles):
        datadict.setdefault(particle, [])
    boxes_per_R = get_ranges_radii(ndiv_R=120)[1]
    # print(boxes_per_R, len(boxes_per_R))
    for i, seed in enumerate(seeds):
        for particle in particles:
            
            datadict[particle].append(ClusteringEventsRadii(folder=folder, fname_prefix=fname_prefix, 
                                spacing=spacing, seed=int(seed), keyword=keyword, particle=particle, n_div_r=120, boxes_per_R=boxes_per_R))  

    def plot_damage(damage_type, datadict):
        plt.figure(figsize=(8,5))
        plot_field(datadict, damage_type=damage_type)
        plt.legend()
        plt.xlabel("distance (um)")
        if damage_type=="DoseGy":
            plt.ylabel(f"Dose ($Gy)$")
            
        elif damage_type=="DSBtotal":
            plt.ylabel(f"n. of {damage_type} ($Gy^-1 Gbp^-1)$")       
        
        elif damage_type=="TotalSBtotal":
            plt.ylim(0, 500)
        plt.xlim(0, 120)
        plt.savefig(f"{folder}/{damage_type}_1k_dense.png")

        return

    for damage in ["DSBtotal", "TotalSBtotal", "DoseGy"]:
        plot_damage(damage, datadict)  

    plt.figure(figsize=(8,5))
    plot_rbe(datadict, damage_type="DSBtotal")
    plt.legend()
    plt.xlabel("distance (um)")
    plt.ylabel("RBE")
    plt.xlim(0, 120)
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