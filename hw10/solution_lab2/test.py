#!/usr/bin/env python3
import numpy as np
def fun(x):
	return 3*(x>3)*(x<6)
x = np.arange(30)
for i in x:
	print(fun(i))