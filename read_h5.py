#%%
import h5py

f = h5py.File('test_hl.h5', 'r')
# %%
import numpy as np
w = np.array(['abc', 'def'], dtype='S')
# f.create_dataset('test', data=w)
# f.create_dataset('test2', data=np.array(['abc', 'ab'], dtype='S'))
# f.create_dataset('test3', data=np.array(['abc'], dtype='S'))
f.create_dataset('test4', data='abc')
f.close()
# %%
