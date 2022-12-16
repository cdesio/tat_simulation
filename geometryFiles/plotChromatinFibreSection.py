"""
Plot deoxyribose positions

B-DNA

"""
import numpy as np
import matplotlib.pyplot as plt
from math import pi

from mpl_toolkits.mplot3d import Axes3D

nanometer = 1e-6  # to match the geant4 convention of mm being the default unit


def plotNucleosome(x_0, y_0, z_0, alpha, beta, gamma):

    posSugar1 = []
    posSugar2 = []
    posBase1 = []
    posBase2 = []

    deoxyRadius = 0.29 * nanometer  # Geant4 molecule definition
    baseRadius = 0.3 * nanometer  # Geant4 molecule definition

    # https://www.sciencedirect.com/science/article/pii/B9780080571737500109
    # DNA Structure and Function 1994, Pages 1-57 Table 1.3
    # right handed helix
    DoubleHelixRadius = 1.0 * nanometer

    DoubleHelixPitch = 3.4 * nanometer  # the height for one rotation
    bpPerTurnDoubleHelix = 10.5
    distBetweenBP = DoubleHelixPitch / bpPerTurnDoubleHelix
    anglebetweenBP = 2 * pi / (bpPerTurnDoubleHelix)

    baseHelixRadius = DoubleHelixRadius - deoxyRadius - baseRadius

    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4598263/
    # 1.75 left handed superhelical turns
    numBPperNucleosome = 147
    nucleosomeHeight = 5.5 * nanometer
    nucleosomeRadius = 5.5 * nanometer  # nucleosome diameter is 11nm
    nuclosomePitch = (nucleosomeHeight) / 1.75
    angleNucleosome = 2 * pi / (numBPperNucleosome / 1.75)

    rotMatrix = np.array(
        [
            [
                np.cos(alpha) * np.cos(beta),
                np.cos(alpha) * np.sin(beta) * np.sin(gamma)
                - np.sin(alpha) * np.cos(gamma),
                np.cos(alpha) * np.sin(beta) * np.cos(gamma)
                + np.sin(alpha) * np.sin(gamma),
            ],
            [
                np.sin(alpha) * np.cos(beta),
                np.sin(alpha) * np.sin(beta) * np.sin(gamma)
                + np.cos(alpha) * np.cos(gamma),
                np.sin(alpha) * np.sin(beta) * np.cos(gamma)
                - np.cos(alpha) * np.sin(gamma),
            ],
            [
                -1 * np.sin(beta),
                np.cos(beta) * np.sin(gamma),
                np.cos(beta) * np.cos(gamma),
            ],
        ]
    )

    for n in range(1, numBPperNucleosome + 1):

        t = angleNucleosome * n
        R = (
            nucleosomeRadius - DoubleHelixRadius
        )  # not - deoxyRadius because the deoxy volumes overlap due to the small radius, increased to avoid this
        a = DoubleHelixRadius
        u = anglebetweenBP * n
        h = nuclosomePitch / (2 * pi)

        # Deoxyribose positions
        y1 = (
            R * np.cos(t)
            + a * np.cos(t) * np.cos(u)
            - (h * a * np.sin(t) * np.sin(u)) / np.sqrt(R ** 2 + h ** 2)
        )
        x1 = (
            R * np.sin(t)
            + a * np.sin(t) * np.cos(u)
            + (h * a * np.cos(t) * np.sin(u)) / np.sqrt(R ** 2 + h ** 2)
        )
        z1 = (
            h * t
            + R * a * np.sin(u) / np.sqrt(R ** 2 + h ** 2)
            - nucleosomeHeight / 2.0
        )

        u += pi  # second helix is 180 degrees out of phase

        y2 = (
            R * np.cos(t)
            + a * np.cos(t) * np.cos(u)
            - (h * a * np.sin(t) * np.sin(u)) / np.sqrt(R ** 2 + h ** 2)
        )
        x2 = (
            R * np.sin(t)
            + a * np.sin(t) * np.cos(u)
            + (h * a * np.cos(t) * np.sin(u)) / np.sqrt(R ** 2 + h ** 2)
        )
        z2 = (
            h * t
            + R * a * np.sin(u) / np.sqrt(R ** 2 + h ** 2)
            - nucleosomeHeight / 2.0
        )

        posSugar1.append(np.matmul([x1, y1, z1], rotMatrix))
        posSugar2.append(np.matmul([x2, y2, z2], rotMatrix))

        # Base positions
        a = baseHelixRadius

        y1 = (
            R * np.cos(t)
            + a * np.cos(t) * np.cos(u)
            - (h * a * np.sin(t) * np.sin(u)) / np.sqrt(R ** 2 + h ** 2)
        )
        x1 = (
            R * np.sin(t)
            + a * np.sin(t) * np.cos(u)
            + (h * a * np.cos(t) * np.sin(u)) / np.sqrt(R ** 2 + h ** 2)
        )
        z1 = (
            h * t
            + R * a * np.sin(u) / np.sqrt(R ** 2 + h ** 2)
            - nucleosomeHeight / 2.0
        )

        u += pi  # second helix is 180 degrees out of phase

        y2 = (
            R * np.cos(t)
            + a * np.cos(t) * np.cos(u)
            - (h * a * np.sin(t) * np.sin(u)) / np.sqrt(R ** 2 + h ** 2)
        )
        x2 = (
            R * np.sin(t)
            + a * np.sin(t) * np.cos(u)
            + (h * a * np.cos(t) * np.sin(u)) / np.sqrt(R ** 2 + h ** 2)
        )
        z2 = (
            h * t
            + R * a * np.sin(u) / np.sqrt(R ** 2 + h ** 2)
            - nucleosomeHeight / 2.0
        )

        posBase1.append(np.matmul([x1, y1, z1], rotMatrix))
        posBase2.append(np.matmul([x2, y2, z2], rotMatrix))

    posSugar1_x = [np.real(a[0]) + x_0 for a in posSugar1]
    posSugar1_y = [np.real(a[1]) + y_0 for a in posSugar1]
    posSugar1_z = [np.real(a[2]) + z_0 for a in posSugar1]

    posSugar2_x = [np.real(a[0]) + x_0 for a in posSugar2]
    posSugar2_y = [np.real(a[1]) + y_0 for a in posSugar2]
    posSugar2_z = [np.real(a[2]) + z_0 for a in posSugar2]

    posBase1_x = [np.real(a[0]) + x_0 for a in posBase1]
    posBase1_y = [np.real(a[1]) + y_0 for a in posBase1]
    posBase1_z = [np.real(a[2]) + z_0 for a in posBase1]

    posBase2_x = [np.real(a[0]) + x_0 for a in posBase2]
    posBase2_y = [np.real(a[1]) + y_0 for a in posBase2]
    posBase2_z = [np.real(a[2]) + z_0 for a in posBase2]

    return (
        posSugar1_x,
        posSugar1_y,
        posSugar1_z,
        posSugar2_x,
        posSugar2_y,
        posSugar2_z,
        posBase1_x,
        posBase1_y,
        posBase1_z,
        posBase2_x,
        posBase2_y,
        posBase2_z,
    )


fig = plt.figure(figsize=(12, 9))
ax = fig.add_subplot(111, projection="3d")
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4031381/
loopRadius = 11 * nanometer
numLayers = 24  # 300 nm z
layerSpacing = 12 * nanometer
numRows = 4
numColumns = 4
deoxyRadius = 0.29 * nanometer  # Geant4 molecule definition
numBPperNucleosome = 147
volumeWidth = 300  # x,y width nm

chromatinSpacingX = (
    volumeWidth * nanometer - deoxyRadius * (numRows - 1)
) / numRows  # - deoxyRadius to avoid overlap with the bounding box
chromatinHeight = volumeWidth * nanometer
chromatinSpacingY = (
    volumeWidth * nanometer - deoxyRadius * (numColumns - 1)
) / numColumns  # - deoxyRadius to avoid overlap with the bounding box

histonePositions = []
c = 0


with open("sugarPos_{}x{}_{}nm.csv".format(numRows, numColumns, volumeWidth), "w") as f:
    totalSugar = (
        ((numRows * numColumns) - numColumns // 2) * numLayers * 6 * numBPperNucleosome
    )
    for y in range(numColumns):
        yOffset = (y - (numColumns / 2 - 0.5)) * chromatinSpacingY
        for x in range(numRows):
            if y % 2 == 0:
                xOffset = (x - (numRows / 2 - 0.5)) * chromatinSpacingX
            else:
                xOffset = (x - (numRows / 2 - 1)) * chromatinSpacingX
                if x == numRows - 1:
                    continue
            for z in range(numLayers):
                for i in range(6):
                    # place 6 nucleosomes per layer
                    point = [
                        loopRadius,
                        0,
                        (z - (numLayers / 2 - 0.5)) * layerSpacing,
                    ]  # centre in z about 0

                    theta = (
                        (i + 0.5 * (z % 2) + 0.5 * (x % 2) + 0.5 * (y % 2)) * 2 * pi / 6
                    )

                    rotMatrix = np.array(
                        [
                            [np.cos(theta), -1 * np.sin(theta), 0],
                            [np.sin(theta), np.cos(theta), 0],
                            [0, 0, 1],
                        ]
                    )

                    point = np.matmul(point, rotMatrix)

                    histonePositions.append(
                        [
                            point[0] + xOffset,
                            point[1] + yOffset,
                            point[2],
                            0,
                            2 * pi - theta,
                            pi / 2,
                        ]
                    )

                    (
                        posSugar1_x,
                        posSugar1_y,
                        posSugar1_z,
                        posSugar2_x,
                        posSugar2_y,
                        posSugar2_z,
                        posBase1_x,
                        posBase1_y,
                        posBase1_z,
                        posBase2_x,
                        posBase2_y,
                        posBase2_z,
                    ) = plotNucleosome(
                        point[0],
                        point[1],
                        point[2],
                        ((i + 0.5 * z) * 2 * pi / 6),
                        2 * pi - theta,
                        pi / 2,
                    )  # each nucleosome is rotated about its axis
                    posSugar1_x = [a + xOffset for a in posSugar1_x]
                    posSugar1_y = [a + yOffset for a in posSugar1_y]
                    posSugar2_x = [a + xOffset for a in posSugar2_x]
                    posSugar2_y = [a + yOffset for a in posSugar2_y]
                    posBase1_x = [a + xOffset for a in posBase1_x]
                    posBase1_y = [a + yOffset for a in posBase1_y]
                    posBase2_x = [a + xOffset for a in posBase2_x]
                    posBase2_y = [a + yOffset for a in posBase2_y]

                    ax.scatter(posSugar1_x, posSugar1_y, posSugar1_z, color="r")
                    ax.scatter(posSugar2_x, posSugar2_y, posSugar2_z, color="b")
                    ax.scatter(posBase1_x, posBase1_y, posBase1_z, color="g")
                    ax.scatter(posBase2_x, posBase2_y, posBase2_z, color="y")

                    for j in range(len(posSugar1_x)):

                        f.write(
                            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
                                posSugar1_x[j],
                                posSugar1_y[j],
                                posSugar1_z[j],
                                posSugar2_x[j],
                                posSugar2_y[j],
                                posSugar2_z[j],
                                posBase1_x[j],
                                posBase1_y[j],
                                posBase1_z[j],
                                posBase2_x[j],
                                posBase2_y[j],
                                posBase2_z[j],
                            )
                        )
                        if c != totalSugar - 1:
                            f.write("\n")
                            c += 1

plt.show()


with open(
    "histonePositions_{}x{}_{}nm.csv".format(numRows, numColumns, volumeWidth), "w"
) as f:
    for i in range(len(histonePositions)):
        f.write(
            "{}\t{}\t{}\t{}\t{}\t{}".format(
                histonePositions[i][0],
                histonePositions[i][1],
                histonePositions[i][2],
                histonePositions[i][3],
                histonePositions[i][4],
                histonePositions[i][5],
            )
        )
        if i != len(histonePositions) - 1:
            f.write("\n")
