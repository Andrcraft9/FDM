import numpy as np
import matplotlib.pyplot as plt

f = open("result", "r")
N = int(f.readline())

c = []

for j in range(1, N):
	cx = []	
	for i in range(1, N):
		cx.append(float(f.readline()))
	c.append(cx)

plt.matshow(c)

plt.colorbar()

plt.show()

