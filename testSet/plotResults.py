
import sys
import inspect 
import os
import matplotlib.pyplot as plt

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir+"/Clustering") 

from readResults import *

def plotResults(path, marker, ax):
    dataset = []

    refResults = {
        'sim1': {'LET': 57.48697519796609, 'direct': [42.575016881338016, 3.494336670649517], 'indirect': [199.16998329688937, 10.818988414855568], 'total': [241.74500017822737, 13.358468244524062]}, 
        'sim2': {'LET': 62.746045173918745, 'direct': [40.701809915177336, 4.211427714856065], 'indirect': [190.21008027987983, 12.224832705220308], 'total': [230.9118901950572, 15.980519772760886]}, 
        'sim3': {'LET': 88.55998604058813, 'indirect': [171.08610437236626, 11.24374029786037], 'total': [212.9916388479665, 12.879660780787333], 'direct': [41.905534475600206, 3.5710209779560884]}
        }

    for i in range(3):
        dataset.append(readIN(currentdir+path+f"sim{i+1}.csv"))

    # plot strand breaks vs LET & energy per dose
    LET=[]
    en = []
    total = []
    direct=[]
    indirect =[]

    checkVals = True
    valsDifference =[]

    for i in range(3):
        LET.append(getLET(dataset[i]))
        en.append(getEnergy(dataset[i]))
        a = calcBreaksperDose("TotalSBtotal",(dataset[i]))
        total.append(a)

        if (abs(a-refResults[f"sim{i+1}"]["total"][0])> refResults[f"sim{i+1}"]["total"][1] ):
            checkVals = False
            diff = abs(a-refResults[f"sim{i+1}"]["total"][0])-refResults[f"sim{i+1}"]["total"][1]
            valsDifference.append(f"sim{i+1} total outside of expected by {diff:.2f}")

        a = calcBreaksperDose("TotalSBdirect",(dataset[i]))
        direct.append(a)

        if (abs(a-refResults[f"sim{i+1}"]["direct"][0])> refResults[f"sim{i+1}"]["direct"][1] ):
            checkVals = False
            print("direct",i, a)
            diff = abs(a-refResults[f"sim{i+1}"]["direct"][0])-refResults[f"sim{i+1}"]["direct"][1]
            valsDifference.append(f"sim{i+1} direct outside of expected by {diff:.2f}")

        a = calcBreaksperDose("TotalSBindirect",(dataset[i]))
        indirect.append(a)

        if (abs(a-refResults[f"sim{i+1}"]["indirect"][0])> refResults[f"sim{i+1}"]["indirect"][1] ):
            checkVals = False
            diff = abs(a-refResults[f"sim{i+1}"]["indirect"][0])-refResults[f"sim{i+1}"]["indirect"][1]
            valsDifference.append(f"sim{i+1} indirect outside of expected by {diff:.2f}")


    ax.plot(LET, total, color = 'k', marker=marker, linestyle= 'None', label=f'Total SB - result')
    ax.plot(LET, indirect, color = 'firebrick',marker=marker, linestyle= 'None', label=f'Indirect SB - result')
    ax.plot(LET, direct, color = 'dimgrey',marker=marker,  linestyle= 'None', label=f'Direct SB - result')

    ax.errorbar([refResults[a]['LET'] for a in refResults], [refResults[a]['total'][0] for a in refResults], yerr = [refResults[a]['total'][1] for a in refResults], color = 'k', marker=".", linestyle= 'None', label=f'Total SB - reference')
    
    ax.errorbar([refResults[a]['LET'] for a in refResults], [refResults[a]['indirect'][0] for a in refResults], yerr = [refResults[a]['indirect'][1] for a in refResults], color = 'firebrick', marker=".", linestyle= 'None', label=f'Indirect SB - reference')

    ax.errorbar([refResults[a]['LET'] for a in refResults], [refResults[a]['direct'][0] for a in refResults], yerr = [refResults[a]['direct'][1] for a in refResults], color = 'dimgrey', marker=".", linestyle= 'None', label=f'Direct SB - reference')

    ax.annotate('Error bars represent the maximum variation from the mean seen for 10 sets of 1000 particles',
                xy=(0.3, 0.9), xytext=(0, 10),
                xycoords=('axes fraction', 'figure fraction'),
                textcoords='offset points',
                size=10, va='bottom', wrap = True)

    if checkVals == False:
        ax.annotate('The following points are outside the previously seen variation:',
            xy=(0.3, 0.85), xytext=(0, 0),
            xycoords=('axes fraction', 'figure fraction'),
            textcoords='offset points',
            size=10, va='bottom', wrap = True)
        for ind, val in enumerate(valsDifference):
            ax.annotate(val,
            xy=(0.3, 0.85 - (ind+1)*0.05), xytext=(0, 0),
            xycoords=('axes fraction', 'figure fraction'),
            textcoords='offset points',
            size=10, va='bottom', wrap = True)
    else:
        ax.annotate('All points are within this range',
            xy=(0.3, 0.9), xytext=(0, 0),
            xycoords=('axes fraction', 'figure fraction'),
            textcoords='offset points',
            size=10, va='bottom', wrap = True)


fig, ax = plt.subplots()
plotResults("/results/", "x", ax)
ax.set_xlabel('LET (keV/$\mu$m)',fontsize=11)
ax.set_ylabel('Number of strand breaks ($Gy^{-1} Gbp^{-1}$)',fontsize=11)
ax.legend(loc = 'center')
plt.savefig('results.png')