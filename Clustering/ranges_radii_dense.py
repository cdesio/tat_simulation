import numpy as np


def get_ranges_radii(ndiv_R = 80, ndiv_Z = 40, spacing=1, start_R = 10.5):
    div_theta = [int(np.ceil((2*np.pi * (start_R+r*spacing)/0.5)))
                for r in range(ndiv_R)]
    boxes_per_r = [int(np.ceil((2*np.pi * (start_R+r*spacing)/0.5)))
                * ndiv_Z for r in range(ndiv_R)]
    boxes_per_r.insert(0, 0)

    ranges_radii = {}
    for ri, (i, ip1) in enumerate(zip(np.cumsum(boxes_per_r), boxes_per_r[1:])):

        ranges_radii[ri] = range(i, i+ip1)
        # print(range(i, i+ip1), ip1, i+ip1)
    return ranges_radii, boxes_per_r


