import numpy as np

a = np.array([[1,2,3,4],[4,3,2,1]])
print a
print a.shape
b = np.lib.pad(a,((5,9),(5,9)), 'constant', constant_values=0)
print b
print b.shape
