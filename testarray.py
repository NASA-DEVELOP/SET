from __future__ import division
import numpy
from numpy import *
import itertools

l_OC = array([[4,7,7,2],[5,7,3,5],[8,7,4,6],[9,9,7,7]])*.000000001
Chi = array([[4,7,7,2],[5,7,3,5],[8,7,4,6],[9,9,7,7]])*.000000001
u0 = array([[4,7,7,2],[5,7,3,5],[8,7,4,6],[9,9,7,7]])*.000000001
theta = array([[4,7,7,2],[5,7,3,5],[8,7,4,6],[9,9,7,7]])*.000000001
PropSumArray = zeros_like(l_OC)
for i,c,u,l,t in itertools.izip(nditer(PropSumArray, op_flags=['readwrite']),nditer(Chi, op_flags=['readwrite']),nditer(u0, op_flags=['readwrite']), nditer(l_OC, op_flags=['readwrite']),nditer(theta, op_flags=['readwrite'])):
	i[...] = c + u+ l + t
print "*************************Propogation Array*******************************"

savetxt("TesttestPropArray.txt", PropSumArray, fmt= "%.6e", delimiter= ',', newline=';')