import numpy as np
import math
import matplotlib.pyplot as plt

def func(x, y):
	return math.cos(math.pi * x) * math.cos(math.pi * y)

N = 500
h = 1.0 / N
c = []

for j in range(0, N+1):
	cx = []	
	for i in range(0, N+1):
		cx.append(func(h*i, h*j))
	c.append(cx)

plt.matshow(c)

plt.colorbar()

plt.show()

