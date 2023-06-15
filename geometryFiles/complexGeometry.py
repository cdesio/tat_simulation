'''
Create more complex DNA structure using fractalDNA, convert to input format for RBE simulation.
'''

from fractaldna.dna_models import dnachain as dna
from mayavi import mlab
import numpy as np
from fractaldna.structure_models import voxelisation as v
from fractaldna.structure_models import hilbert as h

def chromatinFibrePosition(solenoidType, voxelheight = 75, radius = 11, nhistones = 33, histone_angle = 50):
    voxelheight = voxelheight *10 # convert to A
    radius = radius *10 # convert to A
    Angstrum2mm = 1e-7 #A positions convert to Angstrum2mm, default unit of geant

    if solenoidType == "straight":
        chain = dna.Solenoid(
            voxelheight = voxelheight, #A
            radius = radius, #A
            nhistones = nhistones,
            histone_angle = histone_angle,
            twist = False,
            chain = 0, # not used - provides chain number for molecularDNA
        )
    elif solenoidType == "turntwist":
        chain = dna.TurnedSolenoid(
            voxelheight = voxelheight, #A
            radius = radius, #A
            nhistones = nhistones,
            histone_angle = histone_angle,
            twist = True,
            chain = 0,
        )
    if solenoidType == "turn":
        chain = dna.TurnedSolenoid(
            voxelheight = voxelheight, #A
            radius = radius, #A
            nhistones = nhistones,
            histone_angle = histone_angle,
            twist = False,
            chain = 0,
        )

    chain.translate([0, 0, -voxelheight / 2.0])
    # MayaVI plots are best for visualisation here
    # plot = chain.to_strand_plot()

    # mlab.show()
    # plot.scene.save_jpg("single_solenoid_strand_plot.jpg")

    # plot = chain.to_line_plot()
    # mlab.show()

    # plot.scene.save_jpg("single_solenoid_line_plot.jpg")
    sugar = chain.to_frame()
    hist = chain.histones_to_frame()

    positionSugars = np.zeros((max(sugar.bp_idx)+1, 12))
    positionHistones = np.zeros((max(hist.histone_idx)+1, 6))

    for i in range(max(sugar.bp_idx)+1):
        atThatIdx = sugar[sugar.bp_idx==i]

        base0 = atThatIdx[atThatIdx.strand_idx==0].set_index("#name").drop(index = "Sugar").drop(index = "Phosphate")
        base1 = atThatIdx[atThatIdx.strand_idx==1].set_index("#name").drop(index = "Sugar").drop(index = "Phosphate")
        sugar0 = atThatIdx.loc[(atThatIdx["#name"]=="Sugar") & (atThatIdx["strand_idx"] ==0)]
        sugar1 = atThatIdx.loc[(atThatIdx["#name"]=="Sugar") & (atThatIdx["strand_idx"] ==1)]

        positionSugars[i][0:3] = [sugar0["pos_x"].values[0]*Angstrum2mm,sugar0["pos_y"].values[0]*Angstrum2mm,sugar0["pos_z"].values[0]*Angstrum2mm]
        positionSugars[i][3:6] = [sugar1["pos_x"].values[0]*Angstrum2mm,sugar1["pos_y"].values[0]*Angstrum2mm,sugar1["pos_z"].values[0]*Angstrum2mm]
        positionSugars[i][6:9] = [base0["pos_x"].values[0]*Angstrum2mm,base0["pos_y"].values[0]*Angstrum2mm,base0["pos_z"].values[0]*Angstrum2mm]
        positionSugars[i][9:12] = [base1["pos_x"].values[0]*Angstrum2mm,base1["pos_y"].values[0]*Angstrum2mm,base1["pos_z"].values[0]*Angstrum2mm]

    for i in range(max(hist.histone_idx)+1):
        histone = hist.loc[(hist["histone_idx"]==i)]

        positionHistones[i][:] = [histone.pos_x.values[0]*Angstrum2mm, histone.pos_y.values[0]*Angstrum2mm, histone.pos_z.values[0]*Angstrum2mm, 
            histone.rot_x.values[0], histone.rot_y.values[0], histone.rot_z.values[0]]

    return positionSugars, positionHistones

def complexGeometry(voxelheight, nn=2):
    nanometer = 1e-6 # to match the geant4 convention of mm being the default unit

    initial_string = "X"
    # Iterate it as required (here, nn=2)
    # for nn = 8, this takes about two hours on my 16GB RAM Mac
    # nn = 2
    iterated_lstring = h.iterate_lstring(initial_string)
    for _ in range(nn - 1):
        iterated_lstring = h.iterate_lstring(iterated_lstring)

    vf = v.VoxelisedFractal.fromLString(iterated_lstring, pbar=True)
    vf.center_fractal()
    # fig = vf.to_plot()
    # fig.savefig('plot.png')

    # fig = vf.to_pretty_plot()
    # fig.savefig('fractalprettyplot.png')

    simplePos = vf.to_frame()

    # scale to nm, based on one step being one chromatin fibre length which is 75 nm
    positions = simplePos.copy()
    positions.POS_X = simplePos.POS_X*voxelheight*nanometer
    positions.POS_Y = simplePos.POS_Y*voxelheight*nanometer
    positions.POS_Z = simplePos.POS_Z*voxelheight*nanometer

    boxSize = max(positions["POS_X"])*2/nanometer + voxelheight #size of box returned in nm (always cube)

    return positions,  boxSize
            
def transformSugarPositions(positionSugars, row):
    #  add rotation
    new = positionSugars.copy()

    alpha = row.EUL_PSI
    beta = row.EUL_THETA
    gamma = row.EUL_PHI 

    rotMatrix = np.array([
    [np.cos(beta)*np.cos(gamma),np.cos(gamma)*np.sin(beta)*np.sin(alpha)-np.sin(gamma)*np.cos(alpha),np.cos(alpha)*np.sin(beta)*np.cos(gamma)+np.sin(alpha)*np.sin(gamma)],
    [np.sin(gamma)*np.cos(beta),np.sin(alpha)*np.sin(beta)*np.sin(gamma)+np.cos(alpha)*np.cos(gamma),np.sin(gamma)*np.sin(beta)*np.cos(alpha)-np.cos(gamma)*np.sin(alpha)],
    [-1*np.sin(beta), np.cos(beta)*np.sin(alpha), np.cos(beta)*np.cos(alpha)]]).T

    # angle transform 
    for i in range(new.shape[0]):
        new[i,0:3] = np.matmul([positionSugars[i,0],positionSugars[i,1],positionSugars[i,2]],rotMatrix)
        new[i,3:6] = np.matmul([positionSugars[i,3],positionSugars[i,4],positionSugars[i,5]],rotMatrix)
        new[i,6:9] = np.matmul([positionSugars[i,6],positionSugars[i,7],positionSugars[i,8]],rotMatrix)
        new[i,9:12] = np.matmul([positionSugars[i,9],positionSugars[i,10],positionSugars[i,11]],rotMatrix)

    # translation
    new[:,0] += row.POS_X # sugar x
    new[:,1] += row.POS_Y # sugar y
    new[:,2] += row.POS_Z # sugar z
    new[:,3] += row.POS_X # sugar x
    new[:,4] += row.POS_Y # sugar y
    new[:,5] += row.POS_Z # sugar z
    new[:,6] += row.POS_X # base x
    new[:,7] += row.POS_Y # base y
    new[:,8] += row.POS_Z # base z
    new[:,9] += row.POS_X # base x
    new[:,10] += row.POS_Y # base y
    new[:,11] += row.POS_Z # base z

    return new

def transformHistonePositions(positionHistone, row):
    #  add rotation for more complex geometry
    new = positionHistone.copy()

    alpha = row.EUL_PSI
    beta = row.EUL_THETA
    gamma = row.EUL_PHI 

    rotMatrix = np.array([
    [np.cos(beta)*np.cos(gamma),np.cos(gamma)*np.sin(beta)*np.sin(alpha)-np.sin(gamma)*np.cos(alpha),np.cos(alpha)*np.sin(beta)*np.cos(gamma)+np.sin(alpha)*np.sin(gamma)],
    [np.sin(gamma)*np.cos(beta),np.sin(alpha)*np.sin(beta)*np.sin(gamma)+np.cos(alpha)*np.cos(gamma),np.sin(gamma)*np.sin(beta)*np.cos(alpha)-np.cos(gamma)*np.sin(alpha)],
    [-1*np.sin(beta), np.cos(beta)*np.sin(alpha), np.cos(beta)*np.cos(alpha)]]).T

    # angle transform 
    for i in range(new.shape[0]):
        new[i,0:3] = np.matmul([new[i,0],new[i,1],new[i,2]],rotMatrix)

    new[:,0] += row.POS_X # histone x
    new[:,1] += row.POS_Y # histone y
    new[:,2] += row.POS_Z # histone z

    # y and z reversed in RBE
    new[:,3] = positionHistone[:,3]
    new[:,5] = positionHistone[:,4] 
    new[:,4] = positionHistone[:,5]

    return new

# DNA structure parameters

voxelheight = 75 #nm
radius = 11 #nm
nhistones = 33
histone_angle = 50
n = 2 #number of iterations for fractal

# solenoid geometries are saved to save time
try:
    positionSugarsStraight = np.load(f"posSugarsStraight_H{voxelheight}nm_R{radius}nm_{nhistones}hist_{histone_angle}deg.npy")
    positionHistonesStraight = np.load(f"posHistonesStraight_H{voxelheight}nm_R{radius}nm_{nhistones}hist_{histone_angle}deg.npy")
    positionSugarsTurn = np.load(f"posSugarsTurn_H{voxelheight}nm_R{radius}nm_{nhistones}hist_{histone_angle}deg.npy")
    positionHistonesTurn = np.load(f"posHistonesTurn_H{voxelheight}nm_R{radius}nm_{nhistones}hist_{histone_angle}deg.npy")
    positionSugarsTurnTwist = np.load(f"posSugarsTurnTwist_H{voxelheight}nm_R{radius}nm_{nhistones}hist_{histone_angle}deg.npy")
    positionHistonesTurnTwist = np.load(f"posHistonesTurnTwist_H{voxelheight}nm_R{radius}nm_{nhistones}hist_{histone_angle}deg.npy")
    print("Saved solenoid geometry loaded")
except:
    positionSugarsStraight, positionHistonesStraight = chromatinFibrePosition("straight", voxelheight, radius, nhistones, histone_angle)
    positionSugarsTurn, positionHistonesTurn = chromatinFibrePosition("turn", voxelheight, radius, nhistones, histone_angle)
    positionSugarsTurnTwist, positionHistonesTurnTwist = chromatinFibrePosition("turntwist", voxelheight, radius, nhistones, histone_angle)

    np.save(f"posSugarsStraight_H{voxelheight}nm_R{radius}nm_{nhistones}hist_{histone_angle}deg", positionSugarsStraight)
    np.save(f"posHistonesStraight_H{voxelheight}nm_R{radius}nm_{nhistones}hist_{histone_angle}deg", positionHistonesStraight)
    np.save(f"posSugarsTurn_H{voxelheight}nm_R{radius}nm_{nhistones}hist_{histone_angle}deg", positionSugarsTurn)
    np.save(f"posHistonesTurn_H{voxelheight}nm_R{radius}nm_{nhistones}hist_{histone_angle}deg", positionHistonesTurn)
    np.save(f"posSugarsTurnTwist_H{voxelheight}nm_R{radius}nm_{nhistones}hist_{histone_angle}deg", positionSugarsTurnTwist)
    np.save(f"posHistonesTurnTwist_H{voxelheight}nm_R{radius}nm_{nhistones}hist_{histone_angle}deg", positionHistonesTurnTwist)
    print("Solenoid geometry saved")


chromatinFibrePositions, width = complexGeometry(voxelheight, n)
numSugars = 0
with open(f"sugarPos_n{n}_{int(width)}nm_H{voxelheight}nm_R{radius}nm_{nhistones}hist_{histone_angle}deg.csv", "w") as f:
    for ind, row in enumerate(chromatinFibrePositions.itertuples(index=True, name='Pandas')):
        if row.KIND == "straight":
            newSugarPositions = transformSugarPositions(positionSugarsStraight, row)
        elif row.KIND == "turn":
            newSugarPositions = transformSugarPositions(positionSugarsTurn, row)
        elif row.KIND == "turntwist":
            newSugarPositions = transformSugarPositions(positionSugarsTurnTwist, row)
        else:
            raise ValueError()

        for i in range(len(newSugarPositions)):
            f.write('\t'.join([str(a) for a in newSugarPositions[i]]))
            numSugars+=1
            if ind == len(chromatinFibrePositions)-1:
                if i < newSugarPositions.shape[0]-1:
                    f.write("\n")
            else:
                f.write("\n")


print("density = ", numSugars/(width)**3)
print("box size = ", (width), "nm")

with open(f"histonePos_n{n}_{int(width)}nm_H{voxelheight}nm_R{radius}nm_{nhistones}hist_{histone_angle}deg.csv", "w") as f:
    for ind, row in enumerate(chromatinFibrePositions.itertuples(index=True, name='Pandas')):
        if row.KIND == "straight":
            newHistonePositions = transformHistonePositions(positionHistonesStraight, row)
        elif row.KIND == "turn":
            newHistonePositions = transformHistonePositions(positionHistonesTurn, row)
        elif row.KIND == "turntwist":
            newHistonePositions = transformHistonePositions(positionHistonesTurnTwist, row)
        else:
            raise ValueError()

        for i in range(len(newHistonePositions)):
            f.write('\t'.join([str(a) for a in newHistonePositions[i]]))
            if ind == len(chromatinFibrePositions)-1:
                if i < newHistonePositions.shape[0]-1:
                    f.write("\n")
            else:
                f.write("\n")