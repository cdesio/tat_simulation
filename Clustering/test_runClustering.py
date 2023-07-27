# test case

from ROOT import TFile, TTree

from array import array
import unittest
from runClustering import *
import uproot

class Test(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """set up for test class creating root files used by many test mehtods
        """
        # create EventEdep TTree
        input_tree = TTree("EventEdep","EventEdep")
        input_tree_withExtras = TTree("EventEdep","EventEdep")

        eventNum = array('i', [0])
        input_tree.Branch('EventNum', eventNum, 'EventNum/I')
        input_tree_withExtras.Branch('EventNum', eventNum, 'EventNum/I')

        Edep_J = array('d', [0])
        input_tree.Branch('Edep_J', Edep_J, 'Edep_J/D')
        input_tree_withExtras.Branch('Edep_J', Edep_J, 'Edep_J/D')

        part1_EventNum = array('i', [0])
        input_tree.Branch('part1_EventNum', part1_EventNum, 'part1_EventNum/I')
        input_tree_withExtras.Branch('part1_EventNum', part1_EventNum, 'part1_EventNum/I')

        part1_CopyNum = array('i', [0])
        input_tree.Branch('part1_CopyNum', part1_CopyNum, 'part1_CopyNum/I')
        input_tree_withExtras.Branch('part1_CopyNum', part1_CopyNum, 'part1_CopyNum/I')

        PathLengthChromatin = array('d', [0])
        input_tree_withExtras.Branch('PathLengthChromatin', PathLengthChromatin, 'PathLengthChromatin/D')

        PrimaryKEEntrance = array('d', [0])
        input_tree_withExtras.Branch('PrimaryKEEntrance', PrimaryKEEntrance, 'PrimaryKEEntrance/D')

        PrimaryKEExit = array('d', [0])
        input_tree_withExtras.Branch('PrimaryKEExit', PrimaryKEExit, 'PrimaryKEExit/D')

        part1_particleSource = array('i', [0])
        input_tree.Branch('part1_particleSource', part1_particleSource, 'part1_particleSource/I')

        copies = [2, 2, 2, 4, 1, 0, 4, 4, 1, 0]
        particles = [13,14,20,8,4,5,6,13,2,3]
        part1_Evts = [10,21,33,34,35,36,37,38,39,40]

        for i in range(10):
            eventNum[0]=i
            part1_EventNum[0]= part1_Evts[i]
            part1_CopyNum[0]= copies[i]
            Edep_J[0] = float(i**2)+1
            PathLengthChromatin[0] = (i**0.5)
            PrimaryKEEntrance[0] = (i*10)
            PrimaryKEExit[0] = (i*9.5)
            part1_particleSource[0] = particles[i]

            input_tree.Fill()
            input_tree_withExtras.Fill()
        

        # Create info TTree

        finfo = TTree("Info","Info")
        finfoExtras = TTree("Info","Info")

        NumIntersecting = array('i', [0])
        finfo.Branch('NumIntersecting', NumIntersecting, 'NumIntersecting/I')
        finfoExtras.Branch('NumIntersecting', NumIntersecting, 'NumIntersecting/I')

        ChromatinVolume_m3 = array('d', [0])
        finfo.Branch('ChromatinVolume_m3', ChromatinVolume_m3, 'ChromatinVolume_m3/D')
        finfoExtras.Branch('ChromatinVolume_m3', ChromatinVolume_m3, 'ChromatinVolume_m3/D')

        NumBasepairs = array('i', [0])
        finfo.Branch('NumBasepairs', NumBasepairs, 'NumBasepairs/I')
        finfoExtras.Branch('NumBasepairs', NumBasepairs, 'NumBasepairs/I')

        GitHash = bytearray("abc123", "utf-8")
        finfo.Branch('GitHash', GitHash, 'GitHash/C')
        finfoExtras.Branch('GitHash', GitHash, 'GitHash/C')

        MeanLET = array('d', [0])
        finfo.Branch('MeanLET', MeanLET, 'MeanLET/D')
        finfoExtras.Branch('MeanLET', MeanLET, 'MeanLET/D')

        NumIntersecting[0] = 1849
        ChromatinVolume_m3[0] = 300e-9**3
        NumBasepairs[0] = 147
        MeanLET[0] = 100

        finfo.Fill()
        finfoExtras.Fill()

        fDirect = TTree("Direct","Direct")
        fDirectExtras = TTree("Direct","Direct")

        EventNum = array('i', [0])
        fDirect.Branch('EventNum', EventNum, 'EventNum/I')
        fDirectExtras.Branch('EventNum', EventNum, 'EventNum/I')

        eDep_eV = array('d', [0])
        fDirect.Branch('eDep_eV', eDep_eV, 'eDep_eV/D')
        fDirectExtras.Branch('eDep_eV', eDep_eV, 'eDep_eV/D')

        x = array('d', [0])
        fDirect.Branch('x', x, 'x/D')
        fDirectExtras.Branch('x', x, 'x/D')

        y = array('d', [0])
        fDirect.Branch('y', y, 'y/D')
        fDirectExtras.Branch('y', y, 'y/D')
    
        z = array('d', [0])
        fDirect.Branch('z', z, 'z/D')
        fDirectExtras.Branch('z', z, 'z/D')
    
        part1_particleSource = array('i', [0])
        fDirect.Branch('part1_particleSource', part1_particleSource, 'part1_particleSource/I')
        fDirectExtras.Branch('part1_particleSource', part1_particleSource, 'part1_particleSource/I')

        part1_CopyNum = array('i', [0])
        fDirect.Branch('part1_CopyNum', part1_CopyNum, 'part1_CopyNum/I')
        fDirectExtras.Branch('part1_CopyNum', part1_CopyNum, 'part1_CopyNum/I')


        edep = [50]*16
        evt = [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2]

        cls.mapping = {0: 10, 1: 21, 2: 33}

        xs = [10, 10, 10, 10, 10, 10, 20, 20, 20, 10, 20, 20, 10, 20, 20, 10]
        ys = [4, 3.6, -2, -8.8, -2.2, -10, 2.6, 7.8, 1, 6.2, 1.4, -8.4, 5.6, -5, -1.4, 3.4]
        zs = [10]*16
        part1_particleSources = [20, 8, 9, 13, 20, 12, 11, 18, 24, 10, 14, 15, 22, 12, 16, 16]
        part1copy = [0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3]

        for i in range(16):
            EventNum[0] = evt[i]
            eDep_eV[0] = edep[i]
            x[0] = xs[i]
            y[0] = ys[i]
            z[0] = zs[i]
            part1_particleSource[0] = part1_particleSources[i]
            part1_CopyNum[0] = part1copy[i]

            fDirect.Fill()
            fDirectExtras.Fill()

        fIndirect = TTree("Indirect","Indirect")
        fIndirectExtras = TTree("Indirect","Indirect")

        EventNum = array('i', [0])
        fIndirect.Branch('EventNum', EventNum, 'EventNum/I')
        fIndirectExtras.Branch('EventNum', EventNum, 'EventNum/I')

        x = array('d', [0])
        fIndirect.Branch('x', x, 'x/D')
        fIndirectExtras.Branch('x', x, 'x/D')

        y = array('d', [0])
        fIndirect.Branch('y', y, 'y/D')
        fIndirectExtras.Branch('y', y, 'y/D')
    
        z = array('d', [0])
        fIndirect.Branch('z', z, 'z/D')
        fIndirectExtras.Branch('z', z, 'z/D')
    
        part1_particleSource = array('i', [0])
        fIndirect.Branch('part1_particleSource', part1_particleSource, 'part1_particleSource/I')
        fIndirectExtras.Branch('part1_particleSource', part1_particleSource, 'part1_particleSource/I')

        part1_CopyNum = array('i', [0])
        fIndirect.Branch('part1_CopyNum', part1_CopyNum, 'part1_CopyNum/I')
        fIndirectExtras.Branch('part1_CopyNum', part1_CopyNum, 'part1_CopyNum/I')

        DNAmolecule = bytearray("Deoxyribose^0", "utf-8")
        fIndirect.Branch('DNAmolecule', DNAmolecule, 'DNAmolecule/C')
        fIndirectExtras.Branch('DNAmolecule', DNAmolecule, 'DNAmolecule/C')

        radical = bytearray("OH^0", "utf-8")
        fIndirect.Branch('radical', radical, 'radical/C')
        fIndirectExtras.Branch('radical', radical, 'radical/C')

        evt = [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2]

        xs = [10, 10, 10, 10, 10, 10, 20, 20, 20, 10, 20, 20, 10, 20, 20, 10]
        ys = [4.2, 3.8, -1.8, -8.6, -1.6, -9.8, 2.8, 8.0, 1.2, 6.4, 1.6, -8.2, 5.8, -4.8, -1.2, 3.6]
        zs = [10]*16
        part1_particleSources = [21, 9, 10, 13, 22, 11, 23, 15, 24, 12, 11, 19, 24, 9, 21, 13]
        part1copy = [0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3]

        for i in range(16):
            EventNum[0] = evt[i]
            x[0] = xs[i]
            y[0] = ys[i]
            z[0] = zs[i]
            part1_particleSource[0] = part1_particleSources[i]
            part1_CopyNum[0] = part1copy[i]

            fIndirect.Fill()
            fIndirectExtras.Fill()


        outHistFile = TFile.Open("test.root" ,"RECREATE")
        dir = outHistFile.mkdir("ntuple")
        dir.cd()
        # fDose.Write()
        fDirect.Write()
        fIndirect.Write()
        finfo.Write()
        input_tree.Write()
        outHistFile.Close()

        outHistFile = TFile.Open("testExtras.root" ,"RECREATE")
        dir = outHistFile.mkdir("ntuple")
        dir.cd()
        fDirectExtras.Write()
        fIndirectExtras.Write()
        finfoExtras.Write()
        input_tree_withExtras.Write()
        outHistFile.Close()

    @classmethod
    def tearDownClass(cls):
        """ Remove all results files and test root files
        """
        os.remove("test.csv")
        os.remove("testPhoton.csv")
        os.remove("testDecay.csv")
        os.remove("testDecayAlpha.csv")
        os.remove("DSBclusterSize_test.csv")
        os.remove("DSBclusterSize_testPhoton.csv")
        os.remove("DSBclusterSize_testDecay.csv")
        os.remove("DSBclusterSize_testDecayAlpha.csv")
        os.remove("test.root")
        os.remove("testExtras.root")
        os.remove("testPhoton.root")
        print()

    def test_doseCalculation_standalone(self):
        """ Test doseCalculation function with eventType = standalone 
        """
        pathLength = {}
        dosePerEvent = {}
        meanKEperEvent = {}
        eventType = "standalone"
        chromatinVolume = 1e-6

        fFile =  uproot.open("test.root")
        fDirectory = fFile['ntuple']
        input_tree = fDirectory["EventEdep"]
        energy, eventMapping, dosePerEvent, meanKEperEvent, pathLength = doseCalculation(input_tree, pathLength, dosePerEvent,meanKEperEvent, eventType, chromatinVolume, copyCheck = False) 

        self.assertEqual(dosePerEvent, dict([(i, (i**2+1)/(1000*chromatinVolume)) for i in range(10) if i**2 +1>0]), "Should be dictionary of events 1-9'")
        self.assertEqual(eventMapping, {}, "Should be empty dictionary")
        self.assertEqual(energy, "N/A", "Should be N/A")

    def test_doseCalculation_PSnocopy(self):
        """ Test doseCalculation function with eventType = PS and no copy number 
        """
        pathLength = {}
        dosePerEvent = {}
        meanKEperEvent = {}
        eventType = "PS"
        chromatinVolume = 1e-6

        fFile =  uproot.open("test.root")
        fDirectory = fFile['ntuple']
        input_tree = fDirectory["EventEdep"]

        energy, eventMapping, dosePerEvent, _, _ = doseCalculation(input_tree, pathLength, dosePerEvent,meanKEperEvent, eventType, chromatinVolume, copyCheck = False) 

        self.assertEqual(dosePerEvent, {10: 1000, 21: 2000, 33: 5000, 34: 10000, 35: 17000, 36: 26000, 37: 37000, 38: 50000, 39: 65000,40: 82000})

        part1_Evts = [10,21,33,34, 35, 36, 37, 38, 39, 40]

        self.assertEqual(eventMapping, dict([(i, part1_Evts[i]) for i in range(10)]), "Should be empty dictionary")
        self.assertEqual(energy, "N/A", "Should be N/A")

    def test_doseCalculation_PScopy(self):
        """ Test doseCalculation function with eventType = PS and copy number 
        """
        pathLength = {}
        dosePerEvent = {}
        meanKEperEvent = {}
        eventType = "PS"
        chromatinVolume = 1e-6

        fFile =  uproot.open("test.root")
        fDirectory = fFile['ntuple']
        input_tree = fDirectory["EventEdep"]

        energy, eventMapping, dosePerEvent, _, _ = doseCalculation(input_tree, pathLength, dosePerEvent,meanKEperEvent, eventType, chromatinVolume, copyCheck = True, copy = 2) 
        part1_Evts = [10,21,33,34, 35, 36, 37, 38, 39, 40]

        self.assertEqual(dosePerEvent, {10: (0**2+1)/(1000*chromatinVolume), 21: (1**2+1)/(1000*chromatinVolume), 33: (2**2+1)/(1000*chromatinVolume)}, "Should be dictionary of events 25 and 45")
        self.assertEqual(eventMapping, dict([(i, part1_Evts[i]) for i in range(10)]))
        self.assertEqual(energy, "N/A", "Should be N/A")

    def test_doseCalculation_standaloneWithExtraInfo(self):
        """ Test doseCalculation function with eventType = standalone, with a root file containing KE and pathlength per event 
        """
        pathLength = {}
        dosePerEvent = {}
        meanKEperEvent = {}
        eventType = "standalone"
        chromatinVolume = 1e-6

        fFile =  uproot.open("testExtras.root")
        fDirectory = fFile['ntuple']
        input_tree = fDirectory["EventEdep"]

        energy, eventMapping, dosePerEvent, meanKEperEvent , pathLength = doseCalculation(input_tree, pathLength, dosePerEvent,meanKEperEvent, eventType, chromatinVolume, copyCheck = False) 

        self.assertEqual(dosePerEvent, dict([(i, (i**2+1)/(1000*chromatinVolume)) for i in range(10) if i**2+1 >0]), "Should be dictionary of events 1-9'")
        self.assertEqual(eventMapping, {}, "Should be empty dictionary")
        self.assertEqual(energy, 43.875)
        self.assertEqual(pathLength, dict([(i, i**0.5) for i in range(10) if i**2 +1 >0]), "") # is original changed if not returned?
        self.assertEqual(meanKEperEvent, dict([(i, i*19.5/2) for i in range(10) if i**2 +1>0]), "") # is original changed if not returned?

    def test_doseCalculationParticle_noCopy(self):
        """ Test doseCalculation specifying particle function with no copy number 
        """
        dosePerEvent = {}
        pathLength = {}
        meanKEperEvent = {}

        chromatinVolume = 1e-6
        part1_particleSource = [13,14,15,16,17,28,19] # e-
        eventType = "PS"
        fFile =  uproot.open("test.root")
        fDirectory = fFile['ntuple']
        input_tree = fDirectory["EventEdep"]

        energy, eventMapping, dosePerEvent, meanKEperEvent, pathLength = doseCalculation(input_tree, pathLength, dosePerEvent, meanKEperEvent, eventType, chromatinVolume, part1_particleSource = part1_particleSource, copyCheck = False) 

        expected = {10:  1/(1000*chromatinVolume), 21: 2/(1000*chromatinVolume), 38: 50/(1000*chromatinVolume)}

        self.assertEqual(expected.keys(), dosePerEvent.keys())
        for k in expected:
            self.assertAlmostEqual(dosePerEvent[k], expected[k], places=7)

    def test_doseCalculationParticle_Copy(self):
        """ Test doseCalculation function with copy number and particle
        """
        dosePerEvent = {}
        pathLength = {}
        meanKEperEvent = {}

        chromatinVolume = 1e-6
        part1_particleSource = [13, 14, 15, 16, 17, 18, 19] # e-
        eventType = "PS"

        fFile =  uproot.open("test.root")
        fDirectory = fFile['ntuple']
        input_tree = fDirectory["EventEdep"]

        energy, eventMapping, dosePerEvent, meanKEperEvent, pathLength = doseCalculation(input_tree, pathLength, dosePerEvent, meanKEperEvent, eventType, chromatinVolume, part1_particleSource = part1_particleSource, copyCheck = True, copy = 2) 

        expected = {10:  1/(1000*chromatinVolume), 21: 2/(1000*chromatinVolume)}

        self.assertEqual(dosePerEvent, expected)

    def test_getNumIntersecting(self):
        """ Test getNumIntersecting
        """
        fFile =  uproot.open("test.root")
        fDirectory = fFile['ntuple']
        input_tree = fDirectory["Info"]
    
        NumIntersecting = getNumIntersecting(input_tree)
        self.assertEqual(NumIntersecting, 1849)

    def test_loadSugarFile(self):
        """ Test loadSugarFile
        """
        T0, T1 = loadSugarFile("sugarPosTest.bin")
        self.assertEqual(T0.n, 147)
        self.assertEqual(T1.n, 147)

    def test_AccumulateEdep_standalone(self):
        """ Test AccumulateEdep function with eventType = standalone
        """
        cumulatedEnergyDep = {}
        fFile =  uproot.open("test.root")
        fDirectory = fFile['ntuple']
        input_tree = fDirectory["Direct"]

        T0, T1 = loadSugarFile("sugarPosTest.bin")

        cumulatedEnergyDep = AccumulateEdep(input_tree, cumulatedEnergyDep, T0, T1, "standalone", eventMapping = False, copyCheck = False, part1_particleSource = False)

        expected = {
                    (0, 0, 70): 50,
                    (0, 0, 68): 50,
                    (0, 0, 40): 50,
                    (0, 0, 6): 50,
                    (0, 0, 39): 50,
                    (0, 0, 0): 50,
                    (1, 1, 63): 50,
                    (1, 1, 89): 50,
                    (1, 1, 55): 50,
                    (1, 0, 81): 50,
                    (1, 1, 57): 50,
                    (1, 1, 8): 50,
                    (1, 0, 78): 50,
                    (2, 1, 25): 50,
                    (2, 1, 43): 50,
                    (2, 0, 67): 50,
                }
        self.assertEqual(cumulatedEnergyDep, expected)

    def test_AccumulateEdep_EventMapping(self):
        """ Test AccumulateEdep function with eventType = PS and a mapping between event numbers and input PS event numbers
        """
        cumulatedEnergyDep = {}
        fFile =  uproot.open("test.root")
        fDirectory = fFile['ntuple']
        input_tree = fDirectory["Direct"]

        T0, T1 = loadSugarFile("sugarPosTest.bin")

        cumulatedEnergyDep = AccumulateEdep(input_tree, cumulatedEnergyDep, T0, T1, "PS", self.mapping, copyCheck = False, part1_particleSource = False)

        expected = {
                    (10, 0, 70): 50,
                    (10, 0, 68): 50,
                    (10, 0, 40): 50,
                    (10, 0, 6): 50,
                    (10, 0, 39): 50,
                    (10, 0, 0): 50,
                    (21, 1, 63): 50,
                    (21, 1, 89): 50,
                    (21, 1, 55): 50,
                    (21, 0, 81): 50,
                    (21, 1, 57): 50,
                    (21, 1, 8): 50,
                    (21, 0, 78): 50,
                    (33, 1, 25): 50,
                    (33, 1, 43): 50,
                    (33, 0, 67): 50,
                }

        self.assertEqual(cumulatedEnergyDep, expected)

    def test_AccumulateEdep_Copy(self):
        """ Test AccumulateEdep function with eventType = PS, a mapping between event numbers and input PS event numbers and specifying copy number
        """
        cumulatedEnergyDep = {}
        fFile =  uproot.open("test.root")
        fDirectory = fFile['ntuple']
        input_tree = fDirectory["Direct"]

        T0, T1 = loadSugarFile("sugarPosTest.bin")

        cumulatedEnergyDep = AccumulateEdep(input_tree, cumulatedEnergyDep, T0, T1, "PS", eventMapping = self.mapping, copyCheck = True, copy = 2, part1_particleSource = False)

        expected = {
                    (10, 0, 40): 50,
                    (21, 1, 63): 50,
                    (21, 1, 57): 50,
                    (33, 1, 43): 50,
                     }

        self.assertEqual(cumulatedEnergyDep, expected)

    def test_AccumulateEdep_part1_particleSource(self):
        """ Test AccumulateEdep function with eventType = PS, a mapping between event numbers and input PS event numbers and specifying particle from PS
        """
        cumulatedEnergyDep = {}
        fFile =  uproot.open("test.root")
        fDirectory = fFile['ntuple']
        input_tree = fDirectory["Direct"]

        T0, T1 = loadSugarFile("sugarPosTest.bin")

        cumulatedEnergyDep = AccumulateEdep(input_tree, cumulatedEnergyDep, T0, T1, "PS", eventMapping = self.mapping, copyCheck = False, part1_particleSource = [8,9,10,11,12])

        expected = {
                    (10, 0, 68): 50,
                    (10, 0, 40): 50,
                    (10, 0, 0): 50,
                    (21, 1, 63): 50,
                    (21, 0, 81): 50,
                    # (21, 1, 57): 50,
                    (33, 1, 25): 50,
                    # (33, 1, 43): 50,
                     }

        self.assertEqual(cumulatedEnergyDep, expected)

    def test_direct(self):
        """ Test direct function
        """
        cumulatedEnergyDep = {(0,0, 200): 0.5, (1, 0, 400): 40, (2, 1, 300): 45, (3, 0, 402): 0.1} # only test 0 or 1 probability test linear probability in IsEdepSufficient test
        fEMinDamage = 5
        fEMaxDamage = 37.5

        eventsListDirect, copyListDirect, strandListDirect = Direct(cumulatedEnergyDep,fEMinDamage, fEMaxDamage)

        self.assertEqual(eventsListDirect, [1,2])
        self.assertEqual(copyListDirect, [400,300])
        self.assertEqual(strandListDirect, [0,1])

    def test_indirect_standalone(self):
        """ Test indirect function, with eventType = standalone
        """
        fFile =  uproot.open("test.root")
        fDirectory = fFile['ntuple']
        input_tree = fDirectory["Indirect"]

        T0, T1 = loadSugarFile("sugarPosTest.bin")

        events = 0
        copies = 0
        strands = 0

        repeats = 10000
        for i in range(repeats):
            eventsListIndirect, copyListIndirect, strandListIndirect = Indirect(input_tree, "standalone", 0.5, T0, T1, eventMapping = False, copyCheck = False, part1_particleSource = False)

            events +=len(eventsListIndirect)
            copies +=len(copyListIndirect)
            strands +=len(strandListIndirect)

        # check length of each is half of number of reactions in root file (16)
        self.assertAlmostEqual(events/repeats, 8, delta = 0.01*8)
        self.assertAlmostEqual(copies/repeats, 8, delta = 0.01*8)
        self.assertAlmostEqual(strands/repeats, 8, delta = 0.01*8)

    def test_indirect_Mapping(self):
        """ Test indirect function, with eventType = PS with mapping between event number and PS event number
        """
        T0, T1 = loadSugarFile("sugarPosTest.bin")
        fFile =  uproot.open("test.root")
        fDirectory = fFile['ntuple']
        input_tree = fDirectory["Indirect"]

        eventsListIndirect, copyListIndirect, strandListIndirect = Indirect(input_tree, "PS", 1, T0, T1, eventMapping = self.mapping, copyCheck = False, part1_particleSource = False)

        self.assertEqual(eventsListIndirect, [10, 10, 10, 10, 10, 10, 21, 21, 21, 21, 21, 21, 21, 33, 33, 33])
        self.assertEqual(copyListIndirect, [71, 69, 41, 7, 42, 1, 64, 90, 56,82, 58, 9, 79, 26, 44, 68])
        self.assertEqual(strandListIndirect, [0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0])

    def test_indirect_Copy(self):
        """ Test indirect function, with eventType = PS with mapping between event number and PS event number and copy number
        """
        T0, T1 = loadSugarFile("sugarPosTest.bin")
        fFile =  uproot.open("test.root")
        fDirectory = fFile['ntuple']
        input_tree = fDirectory["Indirect"]

        eventsListIndirect, copyListIndirect, strandListIndirect = Indirect(input_tree, "PS", 1, T0, T1, eventMapping = self.mapping, copyCheck = True, copy = 0, part1_particleSource = False)

        self.assertEqual(eventsListIndirect, [10, 10, 21, 21])
        self.assertEqual(copyListIndirect, [71, 42, 56, 79])
        self.assertEqual(strandListIndirect, [0, 0, 1, 0])

    def test_indirect_Particle(self):
        """ Test indirect function, with eventType = PS with mapping between event number and PS event number and primary particle
        """
        T0, T1 = loadSugarFile("sugarPosTest.bin")
        fFile =  uproot.open("test.root")
        fDirectory = fFile['ntuple']
        input_tree = fDirectory["Indirect"]

        eventsListIndirect, copyListIndirect, strandListIndirect = Indirect(input_tree, "PS", 1, T0, T1, eventMapping = self.mapping, copyCheck = False, part1_particleSource = [13,14,15,16,17,18,19])

        self.assertEqual(eventsListIndirect, [10, 21, 21, 33])
        self.assertEqual(copyListIndirect, [7, 90, 9, 68])
        self.assertEqual(strandListIndirect, [0, 1, 1, 0])

    def test_IsEdepSufficient(self):
        """ Test IsEdepSufficient. This function generates a random number to implement the linear probability model, function is called many times to check the damage is returned with the expected probability.
        """
        fEMinDamage = 1
        fEMaxDamage = 2

        result = 0
        repeats = 100000

        # Test many times to check function as linear probability is as expected
        for i in range(repeats):
            result += IsEdepSufficient(1.5, fEMinDamage, fEMaxDamage)

        self.assertAlmostEqual(result, repeats/2, delta = 0.01*repeats) #1%
        self.assertEqual(IsEdepSufficient(0.5, fEMinDamage, fEMaxDamage), 0)
        self.assertEqual(IsEdepSufficient(2.5, fEMinDamage, fEMaxDamage), 1)

    def test_checkPoints(self):
        """ Test checkPoint returns the expected copy number and strand number
        """
        T0, T1 = loadSugarFile("sugarPosTest.bin")
        
        self.assertEqual(checkPoints([(50,50,50)], T0, T1), [[-1,-1]])
        self.assertEqual(checkPoints([(10, -9.8,10)], T0, T1), [[0, 1]])

    def test_getIndex(self):
        """ Test getIndex returns the expected copy number and strand number
        """
        T0, T1 = loadSugarFile("sugarPosTest.bin")
        
        self.assertEqual(getIndex([(10, -9.8,10)], T0, T1), [[1,0]])


    def test_clustering(self):
        """ Test clustering pybind module against test set
        """
        eventsListDirect = [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2]
        copyListDirect = [6, 7, 208, 6000, 820, 822, 600, 620, 650, 670, 6, 7, 211, 295, 250, 8546]
        strandListDirect = [0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0]
        eventsListIndirect = [0, 0, 0, 0, 0, 0, 1, 2, 2, 2]
        copyListIndirect= [6, 400, 410, 217, 823, 826, 200, 200, 300, 302]
        strandListIndirect = [1, 0, 1, 1, 1, 0, 1, 1, 0, 0]

        test = {
            # DoseGy,PathLength_nm,MeanKEMeV, TotalSBdirect,SSBdirect,cSSBdirect,DSBdirect,TotalSBindirect,SSBindirect,cSSBindirect,DSBindirect,TotalSBtotal,SSBtotal,cSSBtotal,DSBtotal
            0: ["N/A","N/A","N/A",6,2,1,1,6,2,0,2,12,1,1,3],
            1: ["N/A","N/A","N/A",7,5,0,1,1,1,0,0,8,6,0,1],
            2: ["N/A","N/A","N/A",3,3,0,0,3,1,1,0,6,3,0,1]
        }

        clusterSizetest = {
            0: ["N/A","N/A","N/A",
            0,0,0,0,0,0,0,0,0,0, #direct
            0,1,0,0,0,0,0,0,0,0, #indirect
            0,0,0,0,0,0,0,0,0,0, #hybrid
            0,0,1,1,0,0,0,0,0,0, #mixed
            0,1,1,1,0,0,0,0,0,0], #total

            1: ["N/A","N/A","N/A",
            0,1,0,0,0,0,0,0,0,0, #direct
            0,0,0,0,0,0,0,0,0,0, #indirect
            0,0,0,0,0,0,0,0,0,0, #hybrid
            0,0,0,0,0,0,0,0,0,0, #mixed
            0,1,0,0,0,0,0,0,0,0], #total

            2: ["N/A","N/A","N/A",
            0,0,0,0,0,0,0,0,0,0, #direct
            0,0,0,0,0,0,0,0,0,0, #indirect
            0,0,1,0,0,0,0,0,0,0, #hybrid
            0,0,0,0,0,0,0,0,0,0, #mixed
            0,0,1,0,0,0,0,0,0,0], #total
        }

        distancestest = {
            0: [399, 417], #6, 405, 821
            1: [],
            2: []
        }

        testResults = clustering([0,1,2],eventsListDirect,copyListDirect,strandListDirect,eventsListIndirect,copyListIndirect, strandListIndirect, True)

        clusteringResults = testResults[0] # strand break number results
        clusterSize = testResults[1] # DSB cluster size, cluster of size 1-10
        distances = testResults[2]
        
        for i in range(3):
            self.assertEqual (clusteringResults[i][1:],test[i][3:])
            self.assertEqual (clusterSize[i][1:],clusterSizetest[i][3:])
            self.assertEqual (distances[i],distancestest[i])

    def test_runClustering_standalone(self):
        """Test full runClustering function with example root files for a standalone simulation
        """
        runClustering(filename = "testExtras.root", 
                      outputFilename = "test.csv",
                      fEMinDamage = 5, 
                      fEMaxDamage = 37.5, 
                      probIndirect = 1, 
                      sugarPosFilename = "sugarPosTest.bin", 
                      simulationType = "standalone", 
                      filenamePhoton =False, 
                      continuous = True, 
                      part1_CopyNum = False, 
                      part1_particleSource = False)
        
        with open("test.csv", "r") as f:
            data = f.readlines()
        
        #  Expected output files
        expected = ["Filename,EnergyMeV,LET,numEvtIntersectingVolume,chromatinVolume,numBP,sugarPosFilename,GitHash,ClusteringGitHash",
        "testExtras.root,43.875,100.0,10,2.6999999999999997e-20,147,sugarPosTest.bin,abc123,M",
        "EventNum,DoseGy,PathLength_nm,MeanKEMeV,TotalSBdirect,SSBdirect,cSSBdirect,DSBdirect,TotalSBindirect,SSBindirect,cSSBindirect,DSBindirect,TotalSBtotal,SSBtotal,cSSBtotal,DSBtotal",
        "0,3.703703703703704e+16,0.0,0.0,6, 0, 3, 0, 6, 0, 3, 0, 12, 0, 3, 0",
        "1,7.407407407407408e+16,1.0,9.75,7, 1, 1, 1, 7, 1, 1, 1, 14, 0, 2, 1",
        "2,1.851851851851852e+17,1.4142135623730951,19.5,3, 3, 0, 0, 3, 3, 0, 0, 6, 0, 3, 0",
        "3,3.703703703703704e+17,1.7320508075688772,29.25,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0",
        "4,6.296296296296297e+17,2.0,39.0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0",
        "5,9.629629629629631e+17,2.23606797749979,48.75,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0",
        "6,1.3703703703703706e+18,2.449489742783178,58.5,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0",
        "7,1.851851851851852e+18,2.6457513110645907,68.25,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0",
        "8,2.4074074074074076e+18,2.8284271247461903,78.0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0",
        "9,3.037037037037037e+18,3.0,87.75,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0"]

        for idx, line in enumerate(data):
            if idx ==1:
                self.assertEqual(line.strip("\n").split(",")[:-1], expected[idx].split(",")[:-1]) #don't check githash as this will change
            else:
                self.assertEqual(line.strip("\n"), expected[idx])

        with open("DSBclusterSize_test.csv", "r") as f:
            dataDSB = f.readlines()

        expectedDSB = ["Filename,EnergyMeV,LET,numEvtIntersectingVolume,chromatinVolume,numBP,sugarPosFilename,SimGitHash,ClusteringGitHash",
                "testExtras.root,43.875,100.0,10,2.6999999999999997e-20,147,sugarPosTest.bin,abc123,M",
                "EventNum,DoseGy,PathLength_nm,Direct_1_SBperDSBcluster,Direct_2_SBperDSBcluster,Direct_3_SBperDSBcluster,Direct_4_SBperDSBcluster,Direct_5_SBperDSBcluster,Direct_6_SBperDSBcluster,Direct_7_SBperDSBcluster,Direct_8_SBperDSBcluster,Direct_9_SBperDSBcluster,Direct_10_SBperDSBcluster,Indirect_1_SBperDSBcluster,Indirect_2_SBperDSBcluster,Indirect_3_SBperDSBcluster,Indirect_4_SBperDSBcluster,Indirect_5_SBperDSBcluster,Indirect_6_SBperDSBcluster,Indirect_7_SBperDSBcluster,Indirect_8_SBperDSBcluster,Indirect_9_SBperDSBcluster,Indirect_10_SBperDSBcluster,Hybrid_1_SBperDSBcluster,Hybrid_2_SBperDSBcluster,Hybrid_3_SBperDSBcluster,Hybrid_4_SBperDSBcluster,Hybrid_5_SBperDSBcluster,Hybrid_6_SBperDSBcluster,Hybrid_7_SBperDSBcluster,Hybrid_8_SBperDSBcluster,Hybrid_9_SBperDSBcluster,Hybrid_10_SBperDSBcluster,Mixed_1_SBperDSBcluster,Mixed_2_SBperDSBcluster,Mixed_3_SBperDSBcluster,Mixed_4_SBperDSBcluster,Mixed_5_SBperDSBcluster,Mixed_6_SBperDSBcluster,Mixed_7_SBperDSBcluster,Mixed_8_SBperDSBcluster,Mixed_9_SBperDSBcluster,Mixed_10_SBperDSBcluster,Total_1_SBperDSBcluster,Total_2_SBperDSBcluster,Total_3_SBperDSBcluster,Total_4_SBperDSBcluster,Total_5_SBperDSBcluster,Total_6_SBperDSBcluster,Total_7_SBperDSBcluster,Total_8_SBperDSBcluster,Total_9_SBperDSBcluster,Total_10_SBperDSBcluster,Total_DSBdistances",
                "0,3.703703703703704e+16,0.0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, []",
                "1,7.407407407407408e+16,1.0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, []",
                "2,1.851851851851852e+17,1.4142135623730951,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, []",
                "3,3.703703703703704e+17,1.7320508075688772,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, []",
                "4,6.296296296296297e+17,2.0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, []",
                "5,9.629629629629631e+17,2.23606797749979,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, []",
                "6,1.3703703703703706e+18,2.449489742783178,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, []",
                "7,1.851851851851852e+18,2.6457513110645907,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, []",
                "8,2.4074074074074076e+18,2.8284271247461903,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, []",
                "9,3.037037037037037e+18,3.0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, []"]


        for idx, line in enumerate(dataDSB):
            if idx ==1:
                self.assertEqual(line.strip("\n").split(",")[:-1], expectedDSB[idx].split(",")[:-1]) #don't check githash as this will change
            else:
                self.assertEqual(line.strip("\n"), expectedDSB[idx])

    def test_runClustering_Photon(self):
        """Test full runClustering function with example root files for a photon simulation
        """
        # create EventEdep TTree
        input_tree = TTree("EventEdep","EventEdep")
        eventNum = array('i', [0])
        input_tree.Branch('EventNum', eventNum, 'EventNum/I')

        Edep_J = array('d', [0])
        input_tree.Branch('Edep_J', Edep_J, 'Edep_J/D')

        PathLengthChromatin = array('d', [0])
        input_tree.Branch('PathLengthChromatin', PathLengthChromatin, 'PathLengthChromatin/D')

        evts = [10, 21, 33]
        for i in range(3):
            eventNum[0]=evts[i]
            Edep_J[0] = float(i**2)+1
            PathLengthChromatin[0] = (i**0.5)

            input_tree.Fill()

        # Create info TTree

        finfo = TTree("Info","Info")

        NumIntersecting = array('i', [0])
        finfo.Branch('NumIntersecting', NumIntersecting, 'NumIntersecting/I')

        ChromatinVolume_m3 = array('d', [0])
        finfo.Branch('ChromatinVolume_m3', ChromatinVolume_m3, 'ChromatinVolume_m3/D')

        NumBasepairs = array('i', [0])
        finfo.Branch('NumBasepairs', NumBasepairs, 'NumBasepairs/I')

        GitHash = bytearray("abc123", "utf-8")
        finfo.Branch('GitHash', GitHash, 'GitHash/C')

        NumIntersecting[0] = 1849
        ChromatinVolume_m3[0] = 300e-9**3
        NumBasepairs[0] = 147

        finfo.Fill()

        fDirect = TTree("Direct","Direct")

        EventNum = array('i', [0])
        fDirect.Branch('EventNum', EventNum, 'EventNum/I')

        eDep_eV = array('d', [0])
        fDirect.Branch('eDep_eV', eDep_eV, 'eDep_eV/D')

        x = array('d', [0])
        fDirect.Branch('x', x, 'x/D')

        y = array('d', [0])
        fDirect.Branch('y', y, 'y/D')
    
        z = array('d', [0])
        fDirect.Branch('z', z, 'z/D')

        edep = [50]*16
        evt = [10, 10, 10, 10, 10, 10, 21, 21, 21, 21, 21, 21, 21, 33, 33, 33]

        xs = [10, 10, 10, 10, 10, 10, 20, 20, 20, 10, 20, 20, 10, 20, 20, 10]
        ys = [4, 3.6, -2, -8.8, -2.2, -10, 2.6, 7.8, 1, 6.2, 1.4, -8.4, 5.6, -5, -1.4, 3.4]
        zs = [10]*16

        for i in range(16):
            EventNum[0] = evt[i]
            eDep_eV[0] = edep[i]
            x[0] = xs[i]
            y[0] = ys[i]
            z[0] = zs[i]

            fDirect.Fill()

        outHistFile = TFile.Open("testPhoton.root" ,"RECREATE")
        dir = outHistFile.mkdir("ntuple")
        dir.cd()
        fDirect.Write()
        finfo.Write()
        input_tree.Write()
        outHistFile.Close()

        runClustering(filename = "test.root", 
                      outputFilename = "testPhoton.csv",
                      fEMinDamage = 5, 
                      fEMaxDamage = 37.5, 
                      probIndirect = 1, 
                      sugarPosFilename = "sugarPosTest.bin", 
                      simulationType = "photon", 
                      filenamePhoton = "testPhoton.root", 
                      continuous = True, 
                      part1_CopyNum = False, 
                      part1_particleSource = False)

        with open("testPhoton.csv", "r") as f:
            data = f.readlines()
        
        #  Expected output files
        expected = ["Filename,EnergyMeV,LET,numEvtIntersectingVolume,chromatinVolume,numBP,sugarPosFilename,GitHash,ClusteringGitHash",
        "test.root,N/A,N/A,1849,2.6999999999999997e-20,147,sugarPosTest.bin,abc123,M",
        "EventNum,DoseGy,PathLength_nm,MeanKEMeV,TotalSBdirect,SSBdirect,cSSBdirect,DSBdirect,TotalSBindirect,SSBindirect,cSSBindirect,DSBindirect,TotalSBtotal,SSBtotal,cSSBtotal,DSBtotal",
        "10,7.407407407407408e+16,N/A,N/A,6, 0, 3, 0, 6, 0, 3, 0, 12, 0, 3, 0",
        "21,1.4814814814814816e+17,N/A,N/A,7, 1, 1, 1, 7, 1, 1, 1, 14, 0, 2, 1",
        "33,3.703703703703704e+17,N/A,N/A,3, 3, 0, 0, 3, 3, 0, 0, 6, 0, 3, 0",
        "34,3.703703703703704e+17,N/A,N/A,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0",
        "35,6.296296296296297e+17,N/A,N/A,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0",
        "36,9.629629629629631e+17,N/A,N/A,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0",
        "37,1.3703703703703706e+18,N/A,N/A,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0",
        "38,1.851851851851852e+18,N/A,N/A,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0",
        "39,2.4074074074074076e+18,N/A,N/A,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0",
        "40,3.037037037037037e+18,N/A,N/A,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0"]

        for idx, line in enumerate(data):
            if idx ==1:
                self.assertEqual(line.strip("\n").split(",")[:-1], expected[idx].split(",")[:-1]) #don't check githash as this will change
            else:
                self.assertEqual(line.strip("\n"), expected[idx])

        with open("DSBclusterSize_testPhoton.csv", "r") as f:
            dataDSB = f.readlines()

        expectedDSB = ["Filename,EnergyMeV,LET,numEvtIntersectingVolume,chromatinVolume,numBP,sugarPosFilename,SimGitHash,ClusteringGitHash",
                "test.root,N/A,N/A,1849,2.6999999999999997e-20,147,sugarPosTest.bin,abc123,M",
                "EventNum,DoseGy,PathLength_nm,Direct_1_SBperDSBcluster,Direct_2_SBperDSBcluster,Direct_3_SBperDSBcluster,Direct_4_SBperDSBcluster,Direct_5_SBperDSBcluster,Direct_6_SBperDSBcluster,Direct_7_SBperDSBcluster,Direct_8_SBperDSBcluster,Direct_9_SBperDSBcluster,Direct_10_SBperDSBcluster,Indirect_1_SBperDSBcluster,Indirect_2_SBperDSBcluster,Indirect_3_SBperDSBcluster,Indirect_4_SBperDSBcluster,Indirect_5_SBperDSBcluster,Indirect_6_SBperDSBcluster,Indirect_7_SBperDSBcluster,Indirect_8_SBperDSBcluster,Indirect_9_SBperDSBcluster,Indirect_10_SBperDSBcluster,Hybrid_1_SBperDSBcluster,Hybrid_2_SBperDSBcluster,Hybrid_3_SBperDSBcluster,Hybrid_4_SBperDSBcluster,Hybrid_5_SBperDSBcluster,Hybrid_6_SBperDSBcluster,Hybrid_7_SBperDSBcluster,Hybrid_8_SBperDSBcluster,Hybrid_9_SBperDSBcluster,Hybrid_10_SBperDSBcluster,Mixed_1_SBperDSBcluster,Mixed_2_SBperDSBcluster,Mixed_3_SBperDSBcluster,Mixed_4_SBperDSBcluster,Mixed_5_SBperDSBcluster,Mixed_6_SBperDSBcluster,Mixed_7_SBperDSBcluster,Mixed_8_SBperDSBcluster,Mixed_9_SBperDSBcluster,Mixed_10_SBperDSBcluster,Total_1_SBperDSBcluster,Total_2_SBperDSBcluster,Total_3_SBperDSBcluster,Total_4_SBperDSBcluster,Total_5_SBperDSBcluster,Total_6_SBperDSBcluster,Total_7_SBperDSBcluster,Total_8_SBperDSBcluster,Total_9_SBperDSBcluster,Total_10_SBperDSBcluster,Total_DSBdistances",
                "10,7.407407407407408e+16,N/A,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, []",
                "21,1.4814814814814816e+17,N/A,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, []",
                "33,3.703703703703704e+17,N/A,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, []",
                "34,3.703703703703704e+17,N/A,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, []",
                "35,6.296296296296297e+17,N/A,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, []",
                "36,9.629629629629631e+17,N/A,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, []",
                "37,1.3703703703703706e+18,N/A,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, []",
                "38,1.851851851851852e+18,N/A,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, []",
                "39,2.4074074074074076e+18,N/A,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, []",
                "40,3.037037037037037e+18,N/A,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, []"]


        for idx, line in enumerate(dataDSB):
            if idx ==1:
                self.assertEqual(line.strip("\n").split(",")[:-1], expectedDSB[idx].split(",")[:-1]) #don't check githash as this will change
            else:
                self.assertEqual(line.strip("\n"), expectedDSB[idx])

    def test_runClustering_decayAllParticles(self):
        """Test full runClustering function with example root files for a decay simulation, including all particles in the clustering
        """
        runClustering(filename = "test.root", 
                      outputFilename = "testDecay.csv",
                      fEMinDamage = 5, 
                      fEMaxDamage = 37.5, 
                      probIndirect = 1, 
                      sugarPosFilename = "sugarPosTest.bin", 
                      simulationType = "decay", 
                      filenamePhoton = False, 
                      continuous = True, 
                      part1_CopyNum = 2, 
                      part1_particleSource = False)
        
        with open("testDecay.csv", "r") as f:
            data = f.readlines()
        
        expected = ["Filename,EnergyMeV,LET,numEvtIntersectingVolume,chromatinVolume,numBP,sugarPosFilename,GitHash,ClusteringGitHash",
        "test.root,N/A,N/A,3,2.6999999999999997e-20,147,sugarPosTest.bin,abc123,M",
        "EventNum,DoseGy,PathLength_nm,MeanKEMeV,TotalSBdirect,SSBdirect,cSSBdirect,DSBdirect,TotalSBindirect,SSBindirect,cSSBindirect,DSBindirect,TotalSBtotal,SSBtotal,cSSBtotal,DSBtotal",
        "10,3.703703703703704e+16,N/A,N/A,1, 1, 0, 0, 1, 1, 0, 0, 2, 0, 1, 0",
        "21,7.407407407407408e+16,N/A,N/A,2, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0",
        "33,1.851851851851852e+17,N/A,N/A,1, 1, 0, 0, 1, 1, 0, 0, 2, 0, 1, 0"]

        for idx, line in enumerate(data):
            if idx ==1:
                self.assertEqual(line.strip("\n").split(",")[:-1], expected[idx].split(",")[:-1]) #don't check githash as this will change
            else:
                self.assertEqual(line.strip("\n"), expected[idx])

        with open("DSBclusterSize_testDecay.csv", "r") as f:
            dataDSB = f.readlines()

        expectedDSB = ["Filename,EnergyMeV,LET,numEvtIntersectingVolume,chromatinVolume,numBP,sugarPosFilename,SimGitHash,ClusteringGitHash",
                "test.root,N/A,N/A,3,2.6999999999999997e-20,147,sugarPosTest.bin,abc123,M",
                "EventNum,DoseGy,PathLength_nm,Direct_1_SBperDSBcluster,Direct_2_SBperDSBcluster,Direct_3_SBperDSBcluster,Direct_4_SBperDSBcluster,Direct_5_SBperDSBcluster,Direct_6_SBperDSBcluster,Direct_7_SBperDSBcluster,Direct_8_SBperDSBcluster,Direct_9_SBperDSBcluster,Direct_10_SBperDSBcluster,Indirect_1_SBperDSBcluster,Indirect_2_SBperDSBcluster,Indirect_3_SBperDSBcluster,Indirect_4_SBperDSBcluster,Indirect_5_SBperDSBcluster,Indirect_6_SBperDSBcluster,Indirect_7_SBperDSBcluster,Indirect_8_SBperDSBcluster,Indirect_9_SBperDSBcluster,Indirect_10_SBperDSBcluster,Hybrid_1_SBperDSBcluster,Hybrid_2_SBperDSBcluster,Hybrid_3_SBperDSBcluster,Hybrid_4_SBperDSBcluster,Hybrid_5_SBperDSBcluster,Hybrid_6_SBperDSBcluster,Hybrid_7_SBperDSBcluster,Hybrid_8_SBperDSBcluster,Hybrid_9_SBperDSBcluster,Hybrid_10_SBperDSBcluster,Mixed_1_SBperDSBcluster,Mixed_2_SBperDSBcluster,Mixed_3_SBperDSBcluster,Mixed_4_SBperDSBcluster,Mixed_5_SBperDSBcluster,Mixed_6_SBperDSBcluster,Mixed_7_SBperDSBcluster,Mixed_8_SBperDSBcluster,Mixed_9_SBperDSBcluster,Mixed_10_SBperDSBcluster,Total_1_SBperDSBcluster,Total_2_SBperDSBcluster,Total_3_SBperDSBcluster,Total_4_SBperDSBcluster,Total_5_SBperDSBcluster,Total_6_SBperDSBcluster,Total_7_SBperDSBcluster,Total_8_SBperDSBcluster,Total_9_SBperDSBcluster,Total_10_SBperDSBcluster,Total_DSBdistances",
                "10,3.703703703703704e+16,N/A,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, []",
                "21,7.407407407407408e+16,N/A,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, []",
                "33,1.851851851851852e+17,N/A,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, []",
                    ]
                
        for idx, line in enumerate(dataDSB):
            if idx ==1:
                self.assertEqual(line.strip("\n").split(",")[:-1], expectedDSB[idx].split(",")[:-1]) #don't check githash as this will change
            else:
                self.assertEqual(line.strip("\n"), expectedDSB[idx])

    def test_runClustering_decayOneParticle(self):
        """Test full runClustering function with example root files for a decay simulation, including one particle in the clustering
        """
        runClustering(filename = "test.root", 
                      outputFilename = "testDecayAlpha.csv",
                      fEMinDamage = 5, 
                      fEMaxDamage = 37.5, 
                      probIndirect = 1, 
                      sugarPosFilename = "sugarPosTest.bin", 
                      simulationType = "decay", 
                      filenamePhoton = False, 
                      continuous = True, 
                      part1_CopyNum = 2, 
                      part1_particleSource = ["alphaRa224", "alphaRn220", "alphaPo216", "alphaBi212", "alphaPo212"]
                      )
        
        with open("testDecayAlpha.csv", "r") as f:
            data = f.readlines()
        
        expected = ["Filename,EnergyMeV,LET,numEvtIntersectingVolume,chromatinVolume,numBP,sugarPosFilename,GitHash,ClusteringGitHash",
        "test.root,N/A,N/A,2,2.6999999999999997e-20,147,sugarPosTest.bin,abc123,M",
        "EventNum,DoseGy,PathLength_nm,MeanKEMeV,TotalSBdirect,SSBdirect,cSSBdirect,DSBdirect,TotalSBindirect,SSBindirect,cSSBindirect,DSBindirect,TotalSBtotal,SSBtotal,cSSBtotal,DSBtotal",
        "10,4.1481481481481485e+17,N/A,N/A,1, 1, 0, 0, 1, 1, 0, 0, 2, 0, 1, 0",
        "21,7.037037037037038e+17,N/A,N/A,1, 1, 0, 0, 1, 1, 0, 0, 2, 0, 1, 0"]

        for idx, line in enumerate(data):
            if idx ==1:
                self.assertEqual(line.strip("\n").split(",")[:-1], expected[idx].split(",")[:-1]) #don't check githash as this will change
            else:
                self.assertEqual(line.strip("\n"), expected[idx])

        with open("DSBclusterSize_testDecayAlpha.csv", "r") as f:
            dataDSB = f.readlines()

        expectedDSB = ["Filename,EnergyMeV,LET,numEvtIntersectingVolume,chromatinVolume,numBP,sugarPosFilename,SimGitHash,ClusteringGitHash",
                "test.root,N/A,N/A,2,2.6999999999999997e-20,147,sugarPosTest.bin,abc123,M",
                "EventNum,DoseGy,PathLength_nm,Direct_1_SBperDSBcluster,Direct_2_SBperDSBcluster,Direct_3_SBperDSBcluster,Direct_4_SBperDSBcluster,Direct_5_SBperDSBcluster,Direct_6_SBperDSBcluster,Direct_7_SBperDSBcluster,Direct_8_SBperDSBcluster,Direct_9_SBperDSBcluster,Direct_10_SBperDSBcluster,Indirect_1_SBperDSBcluster,Indirect_2_SBperDSBcluster,Indirect_3_SBperDSBcluster,Indirect_4_SBperDSBcluster,Indirect_5_SBperDSBcluster,Indirect_6_SBperDSBcluster,Indirect_7_SBperDSBcluster,Indirect_8_SBperDSBcluster,Indirect_9_SBperDSBcluster,Indirect_10_SBperDSBcluster,Hybrid_1_SBperDSBcluster,Hybrid_2_SBperDSBcluster,Hybrid_3_SBperDSBcluster,Hybrid_4_SBperDSBcluster,Hybrid_5_SBperDSBcluster,Hybrid_6_SBperDSBcluster,Hybrid_7_SBperDSBcluster,Hybrid_8_SBperDSBcluster,Hybrid_9_SBperDSBcluster,Hybrid_10_SBperDSBcluster,Mixed_1_SBperDSBcluster,Mixed_2_SBperDSBcluster,Mixed_3_SBperDSBcluster,Mixed_4_SBperDSBcluster,Mixed_5_SBperDSBcluster,Mixed_6_SBperDSBcluster,Mixed_7_SBperDSBcluster,Mixed_8_SBperDSBcluster,Mixed_9_SBperDSBcluster,Mixed_10_SBperDSBcluster,Total_1_SBperDSBcluster,Total_2_SBperDSBcluster,Total_3_SBperDSBcluster,Total_4_SBperDSBcluster,Total_5_SBperDSBcluster,Total_6_SBperDSBcluster,Total_7_SBperDSBcluster,Total_8_SBperDSBcluster,Total_9_SBperDSBcluster,Total_10_SBperDSBcluster,Total_DSBdistances",
                "10,4.1481481481481485e+17,N/A,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, []",
                "21,7.037037037037038e+17,N/A,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, []"]
                
        for idx, line in enumerate(dataDSB):
            if idx ==1:
                self.assertEqual(line.strip("\n").split(",")[:-1], expectedDSB[idx].split(",")[:-1]) #don't check githash as this will change
            else:
                self.assertEqual(line.strip("\n"), expectedDSB[idx])
        
if __name__ == '__main__':
    unittest.main()
