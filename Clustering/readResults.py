import csv
import numpy as np
from typing import Tuple


def readIN(filename:str) -> list:   
    """Read in output files from clustering

    Args:
        filename (str): filename of dataset to read


    Returns:
        list: data
    """
    with open(filename, "r") as f:
        data = csv.reader(f)
        
        data= [a for a in data]
    return  data

def combine(files: list, infoHeaderIndex: int = 0 , infoStartIndex: int = 1, dataHeaderIndex: int = 2, dataStartIndex: int = 3) -> list:
    """Combine multiple datasets from different simulations.

    Args:
        files (list): datasets to combine
        infoHeaderIndex (int, optional): starting index for info header. Defaults to 0.
        infoStartIndex (int, optional):  starting index for info values. Defaults to 1.
        dataHeaderIndex (int, optional): starting index for data header. Defaults to 2.
        dataStartIndex (int, optional): starting index for dataset. Defaults to 3.

    Returns:
        list: data
    """
    data = []

    numEvent = float(files[0][infoStartIndex][files[0][infoHeaderIndex].index("numEvtIntersectingVolume")])
    meanEnergy = files[0][infoStartIndex][files[0][infoHeaderIndex].index("EnergyMeV")].replace('p',".")
    if meanEnergy != "N/A":
        meanEnergy = float(meanEnergy)
    meanLET = files[0][infoStartIndex][files[0][infoHeaderIndex].index("LET")]
    if meanLET != "N/A":
        meanLET = float(meanLET)


    for i in range(1,len(files)): #compare all subsequent files to first


        # check volume and no sugars same
        assert(files[0][infoStartIndex][files[0][infoHeaderIndex].index("numBP") ] == files[i][infoStartIndex][files[i][infoHeaderIndex].index("numBP") ])
        assert(files[0][infoStartIndex][files[0][infoHeaderIndex].index("chromatinVolume") ] == files[i][infoStartIndex][files[i][infoHeaderIndex].index("chromatinVolume") ])

        # calculate new total EventNo
        numEvent += float(files[i][infoStartIndex][files[i][infoHeaderIndex].index("numEvtIntersectingVolume")])
        # Calculate mean LET and energy
        if meanEnergy != "N/A":
            meanEnergy += float(files[i][infoStartIndex][files[i][infoHeaderIndex].index("EnergyMeV")].replace('p',"."))
        if meanLET != "N/A":
            meanLET += float(files[i][infoStartIndex][files[i][infoHeaderIndex].index("LET")])

    if meanEnergy != "N/A":
        meanEnergy = meanEnergy/len(files)
    if meanLET != "N/A":
        meanLET = meanLET/len(files)

    data.append(files[0][infoHeaderIndex])
    data.append([
        "combinedFile",
        meanEnergy,
        meanLET,
        numEvent,
        files[0][infoStartIndex][files[0][infoHeaderIndex].index("chromatinVolume") ] ,
        files[0][infoStartIndex][files[0][infoHeaderIndex].index("numBP") ] ,
        files[0][infoStartIndex][files[0][infoHeaderIndex].index("sugarPosFilename") ] 
    ])
    data.append(files[0][dataHeaderIndex])

    for i in files:
        data.extend(i[dataStartIndex:])

    return data

def calcBreaksperPathLength(SBtype: str, data: list, infoHeaderIndex: int = 0, infoStartIndex: int = 1, dataHeaderIndex: int = 2, dataStartIndex: int =3) -> float:
    """Calculate the number of strand breaks of a given type per path length

    Args:
        SBtype (str): TotalSBdirect, SSBdirect, cSSBdirect, DSBdirect, TotalSBindirect, SSBindirect, cSSBindirect, DSBindirect, TotalSBtotal, SSBtotal, cSSBtotal, DSBtotal
        data (list): from readIN
        infoHeaderIndex (int, optional): starting index for info header. Defaults to 0.
        infoStartIndex (int, optional):  starting index for info values. Defaults to 1.
        dataHeaderIndex (int, optional): starting index for data header. Defaults to 2.
        dataStartIndex (int, optional): starting index for dataset. Defaults to 3.

    Returns:
        float: mean
    """

    totalSB = [float(a[data[dataHeaderIndex].index(SBtype)]) for a in data[dataStartIndex:]]
    totalPathLength = [float(a[data[dataHeaderIndex].index("PathLength_nm")]) for a in data[dataStartIndex:]]

    mean = np.sum(totalSB)/ np.sum(totalPathLength)

    return mean


def calcBreaksperDose(SBtype: str, data: list, infoHeaderIndex: int = 0, infoStartIndex: int = 1, dataHeaderIndex: int = 2, dataStartIndex: int = 3) -> float:
    """Calculate the number of stand breaks of a given type per Gy per number of giga base pairs
    Args:
        SBtype (str): TotalSBdirect, SSBdirect, cSSBdirect, DSBdirect, TotalSBindirect, SSBindirect, cSSBindirect, DSBindirect, TotalSBtotal, SSBtotal, cSSBtotal, DSBtotal
        data (list): from readIN
        infoHeaderIndex (int, optional):  starting index for info header. Defaults to 0.
        infoStartIndex (int, optional):  starting index for info values. Defaults to 1.
        dataHeaderIndex (int, optional): starting index for data header. Defaults to 2.
        dataStartIndex (int, optional): starting index for dataset. Defaults to 3.

    Returns:
        float: mean
    """

    totalSB = [float(a[data[dataHeaderIndex].index(SBtype)]) for a in data[dataStartIndex:]]
    totalDose = np.sum([float(a[data[dataHeaderIndex].index("DoseGy")]) for a in data[dataStartIndex:]])
    numGbp = float(data[infoStartIndex][data[infoHeaderIndex].index("numBP")])*1e-9

    mean = np.sum(totalSB)/totalDose/numGbp

    return mean


def calcDSBClusterSizeperPathLength(SBtype: str,  clustersize: int, data: list, infoHeaderIndex: int = 0, infoStartIndex: int = 1, dataHeaderIndex: int = 2, dataStartIndex: int =3) -> float:
    """Calculate mean DSB cluster size per path length

    Args:
        SBtype (str): Direct, Indirect, Total, Hybrid, Mixed
        clustersize (int): [2,3,4,5,6,7,8,9,10]
        data (list): from readIN        
        infoHeaderIndex (int, optional): starting index for info header. Defaults to 0.
        infoStartIndex (int, optional):  starting index for info values. Defaults to 1.
        dataHeaderIndex (int, optional): starting index for data header. Defaults to 2.
        dataStartIndex (int, optional): starting index for dataset. Defaults to 3.

    Returns:
        float: mean
    """

    totalDSBClusterSize = [float(a[data[dataHeaderIndex].index("{}_{}_SBperDSBcluster".format(SBtype,clustersize))]) for a in data[dataStartIndex:]]
    totalPathLength = [float(a[data[dataHeaderIndex].index("PathLength_nm")]) for a in data[dataStartIndex:]]

    mean = np.sum(totalDSBClusterSize)/ np.sum(totalPathLength)

    return mean


def calcMeanDSBClusterSize(SBtype: str, data: list, dataHeaderIndex: int = 2, dataStartIndex: int =3) -> float:
    """ Calculate mean DSB cluster size for a given source type

    Args:
        SBtype (str): Direct, Indirect, Total, Hybrid, Mixed
        data (list):from readIN
        dataHeaderIndex (int, optional): starting index for data header. Defaults to 2.
        dataStartIndex (int, optional): starting index for dataset. Defaults to 3.

    Returns:
        float: mean
    """
    clustersize = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10] #DSB cluster has at least 2 SB, up to 10 recorded

    totalPerClusterSize = np.zeros([len(clustersize)])
    totalPerClusterSizeSQ = np.zeros([len(clustersize)])

    for i in range(len(clustersize)): 
        totalPerClusterSize[i] += np.sum([float(a[data[dataHeaderIndex].index("{}_{}_SBperDSBcluster".format(SBtype,i+1))]) for a in data[dataStartIndex:]])
        totalPerClusterSizeSQ[i] += np.sum([float(a[data[dataHeaderIndex].index("{}_{}_SBperDSBcluster".format(SBtype,i+1))])**2 for a in data[dataStartIndex:]])


    DSBclustersize = np.sum([totalPerClusterSize[a-1]*a for a in clustersize])/np.sum([totalPerClusterSize[a-1] for a in clustersize])

    return DSBclustersize

def getDensity(data: list, infoHeaderIndex: int = 0 , infoStartIndex: int = 1) -> float:
    """Get base pair density

    Args:
        data (list): from readIN
        infoHeaderIndex (int, optional): starting index for info header. Defaults to 0.
        infoStartIndex (int, optional):  starting index for info values. Defaults to 1.

    Returns:
        float: density
    """
    return float(data[infoStartIndex][data[infoHeaderIndex].index("numBP") ])/(1e27*float(data[infoStartIndex][data[infoHeaderIndex].index("chromatinVolume") ]))

def getLET(data: list, infoHeaderIndex: int = 0 , infoStartIndex: int = 1) -> float:
    """Get mean LET

    Args:
        data (list): from readIN
        infoHeaderIndex (int, optional): starting index for info header. Defaults to 0.
        infoStartIndex (int, optional):  starting index for info values. Defaults to 1.

    Returns:
        float: mean LET
    """
    return float(data[infoStartIndex][data[infoHeaderIndex].index("LET") ])

def getEnergy(data: list, infoHeaderIndex: int = 0 , infoStartIndex: int = 1) -> float:
    """Get mean kinetic energy

    Args:
        data (list): from readIN
        infoHeaderIndex (int, optional): starting index for info header. Defaults to 0.
        infoStartIndex (int, optional):  starting index for info values. Defaults to 1.

    Returns:
        float: mean kinetic energy
    """
    return float(data[infoStartIndex][data[infoHeaderIndex].index("EnergyMeV") ])

def getBP(data: list, infoHeaderIndex: int = 0 , infoStartIndex: int = 1) -> float:
    """Get number of base pairs

    Args:
        data (list): from readIN
        infoHeaderIndex (int, optional): _description_. Defaults to 0.
        infoStartIndex (int, optional): _description_. Defaults to 1.

    Returns:
        float: number of base pairs
    """
    return float(data[infoStartIndex][data[infoHeaderIndex].index("numBP") ])

def calculateDsbSSbRatio(SBtype: str, data: list, dataHeaderIndex: int = 2, dataStartIndex: int =3) -> float:
    """ Calculate DSB to SSB ratio for a given source type

    Args:
        SBtype (str): direct, indirect, total
        data (list): from readIN
        dataHeaderIndex (int, optional): starting index for data header. Defaults to 2.
        dataStartIndex (int, optional): starting index for dataset. Defaults to 3.

    Returns:
        float: mean
    """

    ssb = calcBreaksperDose(f"SSB{SBtype}", data)
    cssb = calcBreaksperDose(f"cSSB{SBtype}", data)
    dsb = calcBreaksperDose(f"DSB{SBtype}", data)

    return dsb/(ssb+cssb)

def getMean(data: list, meanType: str, dataHeaderIndex: int = 2, dataStartIndex: int =3) -> float:
    """ Calculate mean given type

    Args:
        data (list): from readIN
        meanType (str): PathLength_nm, DoseGy, MeanKEMeV
        dataHeaderIndex (int, optional): starting index for data header. Defaults to 2.
        dataStartIndex (int, optional): starting index for dataset. Defaults to 3.

    Returns:
        float: mean
    """
        
    p = [
            float(a[data[dataHeaderIndex].index(meanType)]) for a in data[dataStartIndex:]
        ]
    return np.mean(p)

def calcComplexDSBperDose(SBtype: str, data: list, infoHeaderIndex: int = 0, infoStartIndex: int = 1, dataHeaderIndex: int = 2, dataStartIndex: int = 3) -> float:
    """ complex DSB is any DSB cluster containing more than just a simple pair of breaks one on each strand. cDSB is DSB+ and DSB++

    Args:
        SBtype (str): Direct, Indirect, Total, Hybrid, Mixed
        data (list):from readIN        
        infoHeaderIndex (int, optional): starting index for info header. Defaults to 0.
        infoStartIndex (int, optional):  starting index for info values. Defaults to 1.
        dataHeaderIndex (int, optional): starting index for data header. Defaults to 2.
        dataStartIndex (int, optional): starting index for dataset. Defaults to 3.

    Returns:
        float: mean
    """
    clustersize = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10] #DSB cluster has at least 2 SB, up to 10 recorded

    totalPerClusterSize = np.zeros([len(clustersize)])
    totalPerClusterSizeSQ = np.zeros([len(clustersize)])

    for i in range(len(clustersize)): 
        totalPerClusterSize[i] += np.sum([float(a[data[dataHeaderIndex].index("{}_{}_SBperDSBcluster".format(SBtype,i+1))]) for a in data[dataStartIndex:]])
        totalPerClusterSizeSQ[i] += np.sum([float(a[data[dataHeaderIndex].index("{}_{}_SBperDSBcluster".format(SBtype,i+1))])**2 for a in data[dataStartIndex:]])


    cDSBtotal = np.sum([totalPerClusterSize[a-1] for a in clustersize if (a!=1) and (a!=2) ])
    totalDose = np.sum([float(a[data[dataHeaderIndex].index("DoseGy")]) for a in data[dataStartIndex:]])
    numGbp = float(data[infoStartIndex][data[infoHeaderIndex].index("numBP")])*1e-9

    return cDSBtotal/totalDose/numGbp

def calcComplexDSBperPathLength(SBtype: str, data: list, infoHeaderIndex: int = 0, infoStartIndex: int = 1, dataHeaderIndex: int = 2, dataStartIndex: int = 3) -> float:
    """ complex DSB is any DSB cluster containing more than just a simple pair of breaks one on each strand. cDSB is DSB+ and DSB++

    Args:
        SBtype (str): Direct, Indirect, Total, Hybrid, Mixed
        data (list):from readIN        
        infoHeaderIndex (int, optional): starting index for info header. Defaults to 0.
        infoStartIndex (int, optional):  starting index for info values. Defaults to 1.
        dataHeaderIndex (int, optional): starting index for data header. Defaults to 2.
        dataStartIndex (int, optional): starting index for dataset. Defaults to 3.

    Returns:
        float: mean
    """
    clustersize = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10] #DSB cluster has at least 2 SB, up to 10 recorded

    totalPerClusterSize = np.zeros([len(clustersize)])
    totalPerClusterSizeSQ = np.zeros([len(clustersize)])

    for i in range(len(clustersize)): 
        totalPerClusterSize[i] += np.sum([float(a[data[dataHeaderIndex].index("{}_{}_SBperDSBcluster".format(SBtype,i+1))]) for a in data[dataStartIndex:]])
        totalPerClusterSizeSQ[i] += np.sum([float(a[data[dataHeaderIndex].index("{}_{}_SBperDSBcluster".format(SBtype,i+1))])**2 for a in data[dataStartIndex:]])

    cDSBtotal = np.sum([totalPerClusterSize[a-1] for a in clustersize if (a!=1) and (a!=2) ])
    totalPathLength = [float(a[data[dataHeaderIndex].index("PathLength_nm")]) for a in data[dataStartIndex:]]

    return cDSBtotal/np.sum(totalPathLength)

def calcSimpleDSBperDose(SBtype: str, data: list, infoHeaderIndex: int = 0, infoStartIndex: int = 1, dataHeaderIndex: int = 2, dataStartIndex: int = 3) -> float:
    """ simple DSB is any DSB cluster containing only 2 single strand breaks one on each side

    Args:
        SBtype (str): Direct, Indirect, Total
        data (list):from readIN
        dataHeaderIndex (int, optional): _description_. Defaults to 2.
        dataStartIndex (int, optional): _description_. Defaults to 3.

    Returns:
        float, float: mean, SEM
    """
    clustersize = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10] #DSB cluster has at least 2 SB, up to 10 recorded

    totalPerClusterSize = np.zeros([len(clustersize)])
    totalPerClusterSizeSQ = np.zeros([len(clustersize)])

    for i in range(len(clustersize)): 
        totalPerClusterSize[i] += np.sum([float(a[data[dataHeaderIndex].index("{}_{}_SBperDSBcluster".format(SBtype,i+1))]) for a in data[dataStartIndex:]])
        totalPerClusterSizeSQ[i] += np.sum([float(a[data[dataHeaderIndex].index("{}_{}_SBperDSBcluster".format(SBtype,i+1))])**2 for a in data[dataStartIndex:]])


    sDSBtotal = totalPerClusterSize[1]
    totalDose = np.sum([float(a[data[dataHeaderIndex].index("DoseGy")]) for a in data[dataStartIndex:]])
    numGbp = float(data[infoStartIndex][data[infoHeaderIndex].index("numBP")])*1e-9

    return sDSBtotal/totalDose/numGbp

def calcSimpleDSBperPathLength(SBtype: str, data: list, infoHeaderIndex: int = 0, infoStartIndex: int = 1, dataHeaderIndex: int = 2, dataStartIndex: int = 3) -> float:
    """ simple DSB is any DSB cluster containing only 2 single strand breaks one on each side

    Args:
        SBtype (str): Direct, Indirect, Total, Hybrid, Mixed
        data (list):from readIN        
        infoHeaderIndex (int, optional): starting index for info header. Defaults to 0.
        infoStartIndex (int, optional):  starting index for info values. Defaults to 1.
        dataHeaderIndex (int, optional): starting index for data header. Defaults to 2.
        dataStartIndex (int, optional): starting index for dataset. Defaults to 3.

    Returns:
        float: mean
    """
    clustersize = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10] #DSB cluster has at least 2 SB, up to 10 recorded

    totalPerClusterSize = np.zeros([len(clustersize)])

    for i in range(len(clustersize)): 
        totalPerClusterSize[i] += np.sum([float(a[data[dataHeaderIndex].index("{}_{}_SBperDSBcluster".format(SBtype,i+1))]) for a in data[dataStartIndex:]])

    sDSBtotal = totalPerClusterSize[1]
    totalPathLength = [float(a[data[dataHeaderIndex].index("PathLength_nm")]) for a in data[dataStartIndex:]]


    return sDSBtotal/np.sum(totalPathLength)

def calcDSBClusterDistribution(SBtype: str, data: list, infoHeaderIndex: int = 0, infoStartIndex: int = 1, dataHeaderIndex: int = 2, dataStartIndex: int =3) -> float:
    """ Calculate the mean distribution of DSB, for a given source type

    Args:
        SBtype (str): Direct, Indirect, Total, Hybrid, Mixed
        data (list):from readIN
        infoHeaderIndex (int, optional): starting index for info header. Defaults to 0.
        infoStartIndex (int, optional):  starting index for info values. Defaults to 1.
        dataHeaderIndex (int, optional): starting index for data header. Defaults to 2.
        dataStartIndex (int, optional): starting index for dataset. Defaults to 3.

    Returns:
        float: mean
    """
    clustersize = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10] #DSB cluster has at least 2 SB, up to 10 recorded

    totalPerClusterSize = np.zeros([len(clustersize)])
    totalDose = np.sum([float(a[data[dataHeaderIndex].index("DoseGy")]) for a in data[dataStartIndex:]])
    numGbp = float(data[infoStartIndex][data[infoHeaderIndex].index("numBP")])*1e-9

    for i in range(len(clustersize)): 
        totalPerClusterSize[i] += np.sum([float(a[data[dataHeaderIndex].index("{}_{}_SBperDSBcluster".format(SBtype,i+1))]) for a in data[dataStartIndex:]])/totalDose/numGbp

    return totalPerClusterSize