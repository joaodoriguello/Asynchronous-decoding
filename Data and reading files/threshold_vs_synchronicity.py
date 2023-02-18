import numpy as np
import lmfit
import matplotlib.pyplot as plt
import math as math


def fit_function(p, L, a=0.7, b=-10, c=-30., pth=0.01, v=1.0):
    return a + b*(p - pth)*np.power(L, 1./v) + c*(p - pth)**2*np.power(L, 2./v)


#######################################################################################################################################

def display_fitting(L, physical_error, probability, result, synchronicity):

	shape1 = ['ro', 'bo', 'go', 'yo']
	shape2 = ['r-', 'b-', 'g-', 'y-']

	count = 1
	for j in range(len(L)-1):
		if L[j] != L[j+1]:
			count += 1


	fig, ax = plt.subplots()
	for j in range(count):
		ax.plot(physical_error[int(len(physical_error)*j/count):int(len(physical_error)*(j+1)/count)], probability[int(len(physical_error)*j/count):int(len(physical_error)*(j+1)/count)], shape1[j])
		ax.plot(physical_error[int(len(physical_error)*j/count):int(len(physical_error)*(j+1)/count)], result.best_fit[int(len(physical_error)*j/count):int(len(physical_error)*(j+1)/count)], shape2[j], label='L = %s' % L[int(j*len(L)/count)])


	legend = ax.legend(loc='upper right', shadow=True)
	ax.set_xlabel('physical error probability (%)')
	ax.set_ylabel('decoding success rate')
	plt.title('synchronicity = %s; threshold = %.4g %%' % (synchronicity, 100*result.best_values['pth']))
	plt.show()
	plt.close(fig)


#######################################################################################################################################

def data_analysis(data, flag = 'n'):

	synchronicity = []
	L = []
	physical_error = []
	probability = []
	L_aux = []
	physical_error_aux = []
	probability_aux = []


	with open(data, 'r') as f:
		for line in f:

			if 'synchronicity' in line:
				synchronicity += [float((line.rstrip('\n')).split("=")[1])]
	

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

	for i in range(len(synchronicity)):
		result = mod.fit(probability[i], p=physical_error[i], L=L[i])

		if flag == 'y':
			print(lmfit.fit_report(result.params, show_correl=False))
			display_fitting(L[i], physical_error[i], probability[i], result, synchronicity[i])


		if result.params['pth'].value > 0 and result.params['pth'].stderr/result.params['pth'].value < 0.07:
			x_final += [synchronicity[i]]
			threshold += [result.best_values['pth']]
    

	return threshold, x_final



#######################################################################################################################################
#######################################################################################################################################

'''
    The 'data_analysis' function inputs two flags, either 'yes' or 'no'. The first flag decides if the threshold plots will be shown. The second flag decides if the fitting results will be shown.

    CHANGE THE FIRST INPUT, WHICH IS THE FILE NAME, TO READ IT!
'''


threshold0, synchronicity0 = data_analysis('BG_decoder/Synchronicity_function/BG_decoder', 'n')
threshold1, synchronicity1 = data_analysis('BG_decoder/Synchronicity_function/BG_decoder_degeneracy_tau=1', 'y')
threshold5, synchronicity5 = data_analysis('BG_decoder/Synchronicity_function/BG_decoder_degeneracy_tau=1 (old)', 'y')
threshold2, synchronicity2 = data_analysis('AP_decoder/Synchronicity_function/AP_decoder_optimal_time_weight', 'n')
threshold3, synchronicity3 = data_analysis('CG_decoder/Synchronicity_function/CG_decoder_first_and_second_order_degeneracy', 'n')
threshold4, synchronicity4 = data_analysis('CG_decoder/Synchronicity_function/CG_decoder', 'n')



plt.rcParams['mathtext.fontset'] = 'stix'
#plt.rcParams['font.family'] = 'STIXGeneral'

fig, ax = plt.subplots()
ax.plot(synchronicity0, 100.*np.array(threshold0), 'bo', markersize=5, label=r'BG')
ax.plot(synchronicity1, 100.*np.array(threshold1), 'k^', markersize=6, label=r'BG+$\Omega_0$')
ax.plot(synchronicity5, 100.*np.array(threshold5), 'm^', markersize=6, label=r'BG+$\Omega_0$')
#ax.plot(synchronicity2, 100.*np.array(threshold2), 'r+', markersize=7, label=r'AP')
#ax.plot(synchronicity3, 100.*np.array(threshold3), 'm*', markersize=5, label=r'Contracted+$\Omega_0$+$\Omega_1$')
ax.plot(synchronicity4, 100.*np.array(threshold4), 'gx', markersize=7, label=r'CG')
ax.set_xlabel('Synchronicity s', size=13)
ax.set_ylabel('Threshold in physical error p (%)', size=13)
ax.set_xlim(xmin=-0.01, xmax=1.01)
ax.set_ylim(ymin=1.1, ymax=3.0)
ax.locator_params(axis='x', nbins=20)
ax.locator_params(axis='y', nbins=10)
ax.tick_params(direction='inout', labelright=1, right=1, labeltop=1, top=1)
ax.tick_params(axis="x", labelsize=11)
ax.tick_params(axis="y", labelsize=11)
plt.grid(linestyle = '--', linewidth = 0.3)
legend = ax.legend(loc='upper left', prop={'size': 12}, shadow=True)
plt.savefig('///home/joao/decoders_comparison3.png', bbox_inches='tight')
plt.show()
plt.close(fig)

