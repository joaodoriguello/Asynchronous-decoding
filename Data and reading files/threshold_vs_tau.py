import numpy as np
import lmfit
import matplotlib.pyplot as plt
import math as math


def fit_function(p, L, a=0.7, b=-10, c=-30., pth=0.015, v=1.0):
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
		ax.plot(physical_error[len(physical_error)*j/count:len(physical_error)*(j+1)/count], probability[len(physical_error)*j/count:len(physical_error)*(j+1)/count], shape1[j])
		ax.plot(physical_error[len(physical_error)*j/count:len(physical_error)*(j+1)/count], result.best_fit[len(physical_error)*j/count:len(physical_error)*(j+1)/count], shape2[j], label='L = %s' % L[j*len(L)/count])


	legend = ax.legend(loc='upper right', shadow=True)
	ax.set_xlabel('physical error probability (%)')
	ax.set_ylabel('decoding success rate')
	plt.title('synchronicity = %s; tau = %s; threshold = %.4g %%' % (synchronicity, tau, 100*result.best_values['pth']))
	plt.show()
	plt.close(fig)


#######################################################################################################################################

def data_analysis(data, flag = 'no'):

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
				physical_error_aux += [float((line.rstrip('\n')).split(" ")[2])]
				probability_aux += [float((line.rstrip('\n')).split(" ")[3])]

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

		if flag == 'yes':
 			print(lmfit.fit_report(result.params, show_correl=False))
			display_fitting(L[i], physical_error[i], probability[i], result, synchronicity, tau[i])


		if result.params['pth'].value > 0 and result.params['pth'].stderr/result.params['pth'].value < 0.07:
			x_final += [tau[i]]
			threshold += [result.best_values['pth']]
	    

	return threshold, x_final


#######################################################################################################################################
#######################################################################################################################################

'''
    The 'data_analysis' function inputs two flags, either 'yes' or 'no'. The flag decides if the threshold plots and the fitting results will be shown.

    CHANGE THE FIRST INPUT, WHICH IS THE FILE NAME, TO READ IT!
'''

threshold1, tau1 = data_analysis('CG_decoder/Tau_function/CG_decoder_synchronicity=0.0', 'n')
threshold2, tau2 = data_analysis('CG_decoder/Tau_function/CG_decoder_first_order_synchronicity=0.0', 'n')
threshold3, tau3 = data_analysis('CG_decoder/Tau_function/CG_decoder_first_and_second_order_synchronicity=0.0', 'n')

plt.rcParams['mathtext.fontset'] = 'stix'

fig, ax = plt.subplots()
ax.plot(tau1, 100.*np.array(threshold1), 'm*', markersize=8, label='Contracted')
ax.plot(tau2, 100.*np.array(threshold2), 'c+', markersize=8, label=r'Contracted+$\Omega_0$')
ax.plot(tau3, 100.*np.array(threshold3), 'yx', markersize=8, label=r'Contracted+$\Omega_0$+$\Omega_1$')
#ax.plot([0.58,1.42], [1.688,1.688], 'g--', markersize=12, label='WC')
ax.set_xlabel(r'Degeneracy parameter $\tau$', size=13)
ax.set_ylabel('Threshold in physical error p (%)', size=13)
ax.locator_params(axis='x', nbins=10)
ax.locator_params(axis='y', nbins=5)
ax.tick_params(direction='inout', labelright=1, right=1, labeltop=1, top=1)
ax.tick_params(axis="x", labelsize=11)
ax.tick_params(axis="y", labelsize=11)
ax.set_xlim(xmin=0.58,xmax=1.42)
ax.set_ylim(ymax=1.701)
legend = ax.legend(loc='lower center', prop={'size': 13}, shadow=True)
plt.grid(linestyle = '--', linewidth = 0.3)
#plt.savefig('///home/joao/apbg_tau_dependence.png', bbox_inches='tight')
plt.show()
plt.close(fig)



