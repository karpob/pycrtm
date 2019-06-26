import numpy as np
import h5py 
with open('bt4.bin') as f:
    dataTb = np.fromfile(f, dtype='<f8' )
with open('emissivity4.bin') as f:
    emissivitySaved = np.fromfile(f, dtype='<f8' )
h5 = h5py.File('case4.h5','r+')
h5.create_dataset('Tb', data = dataTb)
h5.create_dataset('emissivity', data = emissivitySaved)
h5.close()



