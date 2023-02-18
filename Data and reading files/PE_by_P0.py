import numpy as np
import lmfit
import matplotlib.pyplot as plt
import math as math


def data_analysis(file_name):

	PE_by_P0 = [[] for i in range(6)]

	with open(file_name, 'r') as f:
		data = f.readlines()

		for i in range(len(data)):

			if 'L=4' in data[i]:
				i += 2
				while '!' not in data[i]:
					PE_by_P0[0] += [float((data[i].rstrip('\n')).split(" ")[2])]
					i += 1

			if 'L=6' in data[i]:
				i += 2
				while '!' not in data[i]:
					PE_by_P0[1] += [float((data[i].rstrip('\n')).split(" ")[2])]
					i += 1

			if 'L=8' in data[i]:
				i += 2
				while '!' not in data[i]:
					PE_by_P0[2] += [float((data[i].rstrip('\n')).split(" ")[2])]
					i += 1

			if 'L=10' in data[i]:
				i += 2
				while '!' not in data[i]:
					PE_by_P0[3] += [float((data[i].rstrip('\n')).split(" ")[2])]
					i += 1

			if 'L=12' in data[i]:
				i += 2
				while '!' not in data[i]:
					PE_by_P0[4] += [float((data[i].rstrip('\n')).split(" ")[2])]
					i += 1

			if 'L=14' in data[i]:
				i += 2
				while '!' not in data[i]:
					PE_by_P0[5] += [float((data[i].rstrip('\n')).split(" ")[2])]
					i += 1



	return PE_by_P0


#######################################################################################################################################
#######################################################################################################################################

'''
    The 'data_analysis' function inputs two flags, either 'yes' or 'no'.

    CHANGE THE FIRST INPUT, WHICH IS THE FILE NAME, TO READ IT!
'''

PE_by_P0_synchronous = data_analysis('ratio_exp--new--async=1.0--time=2')
PE_by_P0_continuous = data_analysis('ratio_exp--new--async=0.0--time=2')


plt.rcParams['mathtext.fontset'] = 'stix'

fig, ax = plt.subplots()
#ax.plot([1,2,3,4,5,6,7,8,9,10], PE_by_P0_synchronous[0], color='0.60', linestyle='', marker='+', markersize=7, label='L = 4')
#ax.plot([1,2,3,4,5,6,7,8,9,10], PE_by_P0_synchronous[1], color='0.48', linestyle='', marker='x', markersize=6, label='L = 6')
#ax.plot([1,2,3,4,5,6,7,8,9,10], PE_by_P0_synchronous[2], color='0.36', linestyle='', marker='^', markersize=6, label='L = 8')
#ax.plot([1,2,3,4,5,6,7,8,9,10], PE_by_P0_synchronous[3], color='0.24', linestyle='', marker='D', markersize=5, label='L = 10')
#ax.plot([1,2,3,4,5,6,7,8,9,10], PE_by_P0_synchronous[4], color='0.12', linestyle='', marker='*', markersize=7, label='L = 12')
ax.plot([1,2,3,4,5,6,7,8,9,10], PE_by_P0_synchronous[5], color='0.00', linestyle='', marker='o', markersize=6, label='s = 1.0')

#ax.plot([1,2,3,4,5,6,7,8,9,10], PE_by_P0_continuous[0], color='r', linestyle='', marker='+', markersize=7, label='L = 4')
#ax.plot([1,2,3,4,5,6,7,8,9,10], PE_by_P0_continuous[1], color='r', linestyle='', marker='x', markersize=6, label='L = 6')
#ax.plot([1,2,3,4,5,6,7,8,9,10], PE_by_P0_continuous[2], color='r', linestyle='', marker='^', markersize=6, label='L = 8')
#ax.plot([1,2,3,4,5,6,7,8,9,10], PE_by_P0_continuous[3], color='r', linestyle='', marker='D', markersize=5, label='L = 10')
#ax.plot([1,2,3,4,5,6,7,8,9,10], PE_by_P0_continuous[4], color='r', linestyle='', marker='*', markersize=7, label='L = 12')
ax.plot([1,2,3,4,5,6,7,8,9,10], PE_by_P0_continuous[5], color='r', linestyle='', marker='o', markersize=6, label='s = 0.0')
ax.set_xlabel(r'$E$', size=18)
ax.set_ylabel(r'$\langle P_E/P_0\rangle$', size=18)
ax.locator_params(axis='x', nbins=10)
ax.locator_params(axis='y', nbins=10)
ax.tick_params(direction='inout', labelright=1, right=1, labeltop=1, top=1)
ax.tick_params(axis="x", labelsize=11)
ax.tick_params(axis="y", labelsize=11)
ax.set_xlim(xmin=0.75, xmax=10.15)
#ax.set_ylim(ymax=1.025)
legend = ax.legend(loc='upper right', prop={'size': 11}, shadow=True)
#plt.savefig('///home/joao/PE_by_P03.png', bbox_inches='tight')
plt.show()
plt.close(fig)



