import numpy as np
from icosphere import icosphere
from plyfile import PlyData, PlyElement
import sys


if len(sys.argv) == 3:

    ico = icosphere(int(sys.argv[1]))

    vert_raw = ico[0]

    vert_list = []
    for vert in vert_raw:
        vert_list.append((vert[0],vert[1],vert[2]))

    vertex = np.array(vert_list,
                          dtype=[('x', 'f4'), ('y', 'f4'),
                                              ('z', 'f4')])


    face_raw = ico[1]

    face_list = []
    for face in face_raw:
        face_list.append((face,0))

    face = np.array(face_list,
                           dtype=[('vertex_indices', 'i4', (3,)),
                                  ('useless', 'u1')])



    v = PlyElement.describe(vertex, 'vertex')
    f = PlyElement.describe(face,'face')

    PlyData([v,f],text=False).write(sys.argv[2])

