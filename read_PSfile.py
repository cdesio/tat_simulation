import numpy as np


particle_map = {1: "alpha",
                2: "gamma",
                3: "e-",
                4: "nu_e",
                5: "At211",
                6: "Po211",
                7: "Bi207",
                8: "Pb207",
                9: "e+" }


from dataclasses import dataclass, InitVar, field
from typing import List

@dataclass
class Event:
    psdata: InitVar[List[float]]
    positions: List[float] = field(init=False)
    momentum: List[float] = field(init=False)
    energy: float = field(init=False)
    eventID: int = field(init=False)
    particleID: int = field(init=False)
    copyNo: int = field(init=False)
    time: float = field(init=False)
    parent: int = field(init=False)
    # excitation: float = field(init=False)

    def __post_init__(self, data_arr):
            if data_arr is not None:
                self.positions = data_arr[:3]
                self.momentum = data_arr[3:6]
                self.energy = data_arr[6]
                self.eventID = int(data_arr[7])
                self.particleID = int(data_arr[8])
                self.copyNo = int(data_arr[9])
                self.time = data_arr[10]
                self.parent = data_arr[11]
                # if len(data_arr)>11:
                #     self.excitation = data_arr[12]


class PSEvents:
    def __init__(self, psfile):
        self.events=[]
        self.psfile = psfile
        events_data = np.fromfile(self.psfile, dtype='float32').reshape(-1,12)
        for evt in events_data:
            self.events.append(Event(evt))
    @property
    def energies(self):
        if len(self.events):
            return np.array([evt.energy for evt in self.events])
    # @property
    # def exc_energies(self):
    #     if len(self.events):
    #         return np.array([evt.excitation for evt in self.events])
    @property
    def eventIDs(self):
        if len(self.events):
            return np.array([evt.eventID for evt in self.events])
    @property
    def particleIDs(self):
        if len(self.events):
            return np.array([evt.particleID for evt in self.events])
    @property
    def copyNos(self):
        if len(self.events):
            return np.array([evt.copyNo for evt in self.events])
    @property
    def parents(self):
        if len(self.events):
            return np.array([evt.parent for evt in self.events])
    @property
    def parent_names(self):
        return list(map(lambda pid: particle_map[pid], self.parents))
    @property
    def particle_names(self):
        return list(map(lambda pid: particle_map[pid], self.particleIDs))
    @property
    def times(self):
        if len(self.events):
            return np.array([evt.time for evt in self.events])

if __name__=="__main__":
    pass