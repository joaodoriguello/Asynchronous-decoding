import numpy as np
import lmfit
import matplotlib.pyplot as plt
import math as math


def data_analysis(data):

	synchronicity = []
	P2_or_P1_by_P0 = [[] for i in range(6)]

	with open(data, 'r') as f:
		for line in f:


			if 'synchronicity' in line:
				synchronicity += [float((line.rstrip('\n')).split("=")[1])]	

			if '# 4' in line:
				P2_or_P1_by_P0[0] += [float((line.rstrip('\n')).split(" ")[3])]

			if '# 6' in line:
				P2_or_P1_by_P0[1] += [float((line.rstrip('\n')).split(" ")[3])]

			if '# 8' in line:
				P2_or_P1_by_P0[2] += [float((line.rstrip('\n')).split(" ")[3])]

			if '# 10' in line:
				P2_or_P1_by_P0[3] += [float((line.rstrip('\n')).split(" ")[3])]

			if '# 12' in line:
				P2_or_P1_by_P0[4] += [float((line.rstrip('\n')).split(" ")[3])]

			if '# 14' in line:
				P2_or_P1_by_P0[5] += [float((line.rstrip('\n')).split(" ")[3])]




	return P2_or_P1_by_P0, synchronicity


#######################################################################################################################################
#######################################################################################################################################

'''
    The 'data_analysis' function inputs two flags, either 'yes' or 'no'.

    CHANGE THE FIRST INPUT, WHICH IS THE FILE NAME, TO READ IT!
'''

P1_by_P0, synchronicity = data_analysis('ratio_second_order--new--time=2')

plt.rcParams['mathtext.fontset'] = 'stix'

fig, ax = plt.subplots()
ax.plot(synchronicity, P1_by_P0[0], color='0.60', linestyle='', marker='+', markersize=7, label='L = 4')
ax.plot(synchronicity, P1_by_P0[1], color='0.48', linestyle='', marker='x', markersize=6, label='L = 6')
ax.plot(synchronicity, P1_by_P0[2], color='0.36', linestyle='', marker='^', markersize=6, label='L = 8')
ax.plot(synchronicity, P1_by_P0[3], color='0.24', linestyle='', marker='D', markersize=5, label='L = 10')
ax.plot(synchronicity, P1_by_P0[4], color='0.12', linestyle='', marker='*', markersize=7, label='L = 12')
ax.plot(synchronicity, P1_by_P0[5], color='0.00', linestyle='', marker='o', markersize=6, label='L = 14')
ax.set_xlabel(r'Synchronicity s', size=12)
ax.set_ylabel(r'$\langle P_1/P_0\rangle$', size=18)
ax.locator_params(axis='x', nbins=20)
ax.locator_params(axis='y', nbins=10)
ax.tick_params(direction='inout', labelright=1, right=1, labeltop=1, top=1)
ax.tick_params(axis="x", labelsize=11)
ax.tick_params(axis="y", labelsize=11)
ax.set_xlim(xmin=-0.025, xmax=1.025)
#ax.set_ylim(ymin=0.4, ymax=1.01)
legend = ax.legend(loc='lower left', prop={'size': 11}, shadow=True)
#plt.savefig('///home/joao/P1_by_P02.png', bbox_inches='tight')
plt.show()
plt.close(fig)



P2_by_P0, synchronicity = data_analysis('ratio_third_order--new--time=2')

plt.rcParams['mathtext.fontset'] = 'stix'

fig, ax = plt.subplots()
ax.plot(synchronicity, P2_by_P0[0], color='0.60', linestyle='', marker='+', markersize=7, label='L = 4')
ax.plot(synchronicity, P2_by_P0[1], color='0.48', linestyle='', marker='x', markersize=6, label='L = 6')
ax.plot(synchronicity, P2_by_P0[2], color='0.36', linestyle='', marker='^', markersize=6, label='L = 8')
ax.plot(synchronicity, P2_by_P0[3], color='0.24', linestyle='', marker='D', markersize=5, label='L = 10')
ax.plot(synchronicity, P2_by_P0[4], color='0.12', linestyle='', marker='*', markersize=7, label='L = 12')
ax.plot(synchronicity, P2_by_P0[5], color='0.00', linestyle='', marker='o', markersize=6, label='L = 14')
ax.set_xlabel(r'Synchronicity s', size=12)
ax.set_ylabel(r'$\langle P_2/P_0\rangle$', size=18)
ax.locator_params(axis='x', nbins=20)
ax.locator_params(axis='y', nbins=10)
ax.tick_params(direction='inout', labelright=1, right=1, labeltop=1, top=1)
ax.tick_params(axis="x", labelsize=11)
ax.tick_params(axis="y", labelsize=11)
ax.set_xlim(xmin=-0.025, xmax=1.025)
#ax.set_ylim(ymin=0.4, ymax=1.01)
legend = ax.legend(loc='lower left', prop={'size': 11}, shadow=True)
#plt.savefig('///home/joao/P2_or_P1_by_P02.png', bbox_inches='tight')
plt.show()
plt.close(fig)


