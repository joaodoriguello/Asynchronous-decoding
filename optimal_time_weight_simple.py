import numpy as np
import lmfit
import matplotlib.pyplot as plt
import math as math



async0 = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
#time_weight = [1.28, 1.24, 1.23, 1.23, 1.16, 1.15, 1.06, 1.06, 1.0, 0.98, 1.0]
time_weight = [1.28, 1.24, 1.20, 1.19, 1.16, 1.13, 1.10, 1.09, 1.01, 0.98, 1.0]


plt.rcParams['mathtext.fontset'] = 'stix'

fig, ax = plt.subplots()
ax.plot(async0, time_weight, 'ko', markersize=5)
ax.set_xlabel(r'Synchronicity s', size=13)
ax.set_ylabel(r'Optimal time weight $\mathregular{w_{time}^{BG}}$', size=13)
ax.set_xlim(xmin=-0.02,xmax=1.02)
#ax.set_ylim(ymin=0.5,ymax=1.5)
ax.locator_params(axis='x', nbins=20)
ax.locator_params(axis='y', nbins=10)
ax.tick_params(direction='inout', labelright=1, right=1, labeltop=1, top=1)
ax.tick_params(axis="x", labelsize=11)
ax.tick_params(axis="y", labelsize=11)
#legend = ax.legend(loc='upper left', shadow=True)
plt.grid(linestyle = '--', linewidth = 0.3)
plt.savefig('///home/joao/optimal_time_weight_BG.png', bbox_inches='tight')
plt.show()
plt.close(fig)


