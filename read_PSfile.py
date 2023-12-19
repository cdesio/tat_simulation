import numpy as np

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


from dataclasses import dataclass, InitVar, field
from typing import List

@dataclass
class Event:
    psdata: InitVar[List[float]]
    localpositions: List[float] = field(init=False)
    globalpositions: List[float] = field(init=False)
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
                if len(data_arr)==12:
                    self.localpositions = data_arr[:3]
                    self.momentum = data_arr[3:6]
                    self.energy = data_arr[6]
                    self.eventID = int(data_arr[7])
                    self.particleID = int(data_arr[8])
                    self.copyNo = int(data_arr[9])
                    self.time = data_arr[10]
                    self.parent = data_arr[11]
                elif len(data_arr)==15:
                    self.globalpositions = data_arr[:3]
                    self.localpositions = data_arr[3:6]
                    self.momentum = data_arr[6:9]
                    self.energy = data_arr[9]
                    self.eventID = int(data_arr[10])
                    self.particleID = int(data_arr[11])
                    self.copyNo = int(data_arr[12])
                    self.time = data_arr[13]
                    self.parent = data_arr[14]


class PSEvents:
    def __init__(self, psfile):
        self.events=[]
        self.psfile = psfile
        events_data = np.fromfile(self.psfile, dtype='float32').reshape(-1,15)
        for evt in events_data:
            self.events.append(Event(evt))
    @property
    def energies(self):
        if len(self.events):
            return np.array([evt.energy for evt in self.events])
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