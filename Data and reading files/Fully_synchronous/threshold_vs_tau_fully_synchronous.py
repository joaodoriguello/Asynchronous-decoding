import numpy as np
import lmfit
import matplotlib.pyplot as plt
import math as math
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)


def fit_function(p, L, a=0.7, b=-10, c=-30., pth=0.015, v=1.0):
#    return a + b*np.power(p - pth, v) + c*np.power(p - pth, 2.*v)
    return a + b*(p - pth)*np.power(L, 1./v) + c*(p - pth)**2*np.power(L, 2./v)


#######################################################################################################################################

def display_fitting(L, physical_error, probability, result, synchronicity, tau):


	shape1 = ['ro', 'bo', 'go', 'yo']
	shape2 = ['r-', 'b-', 'g-', 'y-']

	count = 1
	for j in range(len(L)-1):
		if L[j] != L[j+1]:
			count += 1


	fig, ax = plt.subplots()
	for j in range(count):
		ax.plot(physical_error[int(len(physical_error)*j/count):int(len(physical_error)*(j+1)/count)], probability[int(len(physical_error)*j/count):int(len(physical_error)*(j+1)/count)], shape1[j])
		ax.plot(physical_error[int(len(physical_error)*j/count):int(len(physical_error)*(j+1)/count)], result.best_fit[int(len(physical_error)*j/count):int(len(physical_error)*(j+1)/count)], shape2[j], label='L = %d' % L[int(j*len(L)/count)])


	legend = ax.legend(loc='upper right', prop={'size': 12}, shadow=True)
	ax.set_xlabel('Physical error probability (%)', size=12)
	ax.set_ylabel('Decoding success probability (%)', size=12)
	ax.tick_params(axis="x", labelsize=11)
	ax.tick_params(axis="y", labelsize=11)
    	plt.title('async error = %s; tau = %s; threshold = %.4g %%' % (synchronicity, tau, 100*result.best_values['pth']))
#	plt.savefig('///home/joao/fig5-8b.png', bbox_inches='tight')
	plt.show()
	plt.close(fig)


#######################################################################################################################################

def data_analysis(data, flag = 'n'):

	synchronicity = []
	tau = []
	L = []
	physical_error = []
	probability = []
	L_aux = []
	physical_error_aux = []
	probability_aux = []


	with open(data, 'r') as f:
		for line in f:

			if 'tau' in line:
				tau += [float((line.rstrip('\n')).split("=")[1])]

			if 'synchronicity' in line:
				synchronicity = float((line.rstrip('\n')).split("=")[1])	

			if '#' in line:
				L_aux += [float((line.rstrip('\n')).split(" ")[1])]
				physical_error_aux += [100*float((line.rstrip('\n')).split(" ")[2])]
				probability_aux += [100*float((line.rstrip('\n')).split(" ")[3])]

			if '!' in line:
				L += [L_aux]
				physical_error += [physical_error_aux]
				probability += [probability_aux]
				L_aux = []
				physical_error_aux = []
				probability_aux = []



	threshold = []
	x_final = []
	mod = lmfit.Model(fit_function, independent_vars=['p', 'L'])
 
	for i in range(len(tau)):
		result = mod.fit(probability[i], p=physical_error[i], L=L[i])

		if flag == 'y':
			print(lmfit.fit_report(result.params, show_correl=False))
			display_fitting(L[i], physical_error[i], probability[i], result, synchronicity, time_weight, tau[i])


		if result.params['pth'].value > 0 and result.params['pth'].stderr/result.params['pth'].value < 0.07:
			x_final += [tau[i]]

		if synchronicity != 0:
			''' Rescaling of the threshold according to the asynchronism '''
			threshold += [0.5*(1 - math.pow(1 - 2*result.best_values['pth'], 1/synchronicity))]

		else:
			threshold += [result.best_values['pth']]
	    

	return threshold, x_final


#######################################################################################################################################
#######################################################################################################################################

'''
    The 'data_analysis' function inputs two flags, either 'yes' or 'no'. The flag decides if the threshold plots and the fitting results will be shown.

    CHANGE THE FIRST INPUT, WHICH IS THE FILE NAME, TO READ IT!
'''

threshold1, tau1 = data_analysis('3D_first_order--(10,12,14)', 'n')
threshold2, tau2 = data_analysis('3D_P0_plus_P1', 'y')

plt.rcParams['mathtext.fontset'] = 'stix'

fig, ax = plt.subplots()
ax.plot(tau1, np.array(threshold1), 'ko', markersize=5, label=r'$\Omega_0$')
ax.plot(tau2, np.array(threshold2), 'rx', markersize=6, label=r'$P_0+P_1$')
#ax.plot(tau3, 100.*np.array(threshold3), 'go', label='second order')
ax.set_xlabel(r'Degeneracy parameter $\tau$', size=13)
ax.set_ylabel('Threshold in physical error p (%)', size=13)
ax.locator_params(axis='x', nbins=20)
ax.locator_params(axis='y', nbins=20)
#ax.yaxis.set_major_locator(MultipleLocator(0.01))
#ax.yaxis.set_minor_locator(MultipleLocator(0.002))
ax.tick_params(direction='inout', labelright=1, right=1, labeltop=1, top=1)
ax.tick_params(axis="x", labelsize=11)
ax.tick_params(axis="y", labelsize=11)
ax.set_xlim(xmin=-0.03,xmax=2.03)
legend = ax.legend(loc='best', prop={'size': 13}, shadow=True)
plt.grid(linestyle = '--', linewidth = 0.3)
#plt.savefig('///home/joao/3D_tau.png', bbox_inches='tight')
plt.show()
plt.close(fig)



