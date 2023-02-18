import numpy as np
import lmfit
import matplotlib.pyplot as plt
import math as math


def fit_function(p, L, a=0.7, b=-10, c=-30., pth=0.015, v=1.0):
    return a + b*(p - pth)*np.power(L, 1./v) + c*(p - pth)**2*np.power(L, 2./v)


#######################################################################################################################################

def display_fitting(L, physical_error, probability, result, synchronicity, time_weight):

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
	plt.title('synchronicity = %s; time weight = %.5s; threshold = %.4g %%' % (synchronicity, time_weight, 100*result.best_values['pth']))
	plt.show()
	plt.close(fig)


#######################################################################################################################################

def data_analysis(data, flag = 'n'):

	time_weight = []
	L = []
	physical_error = []
	probability = []
	L_aux = []
	physical_error_aux = []
	probability_aux = []


	with open(data, 'r') as f:
		for line in f:

			if 'synchronicity' in line:
				synchronicity = float((line.rstrip('\n')).split("=")[1])

			if 'time_weight' in line:
				time_weight += [float((line.rstrip('\n')).split("=")[1])]	

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

	for i in range(len(time_weight)):
		result = mod.fit(probability[i], p=physical_error[i], L=L[i])

		if flag == 'y':
			print(lmfit.fit_report(result.params, show_correl=False))
			display_fitting(L[i], physical_error[i], probability[i], result, synchronicity, time_weight[i])


		if result.params['pth'].value > 0 and result.params['pth'].stderr/result.params['pth'].value < 0.07:
			x_final += [time_weight[i]]
			threshold += [result.best_values['pth']]
	    

	return threshold, x_final


#######################################################################################################################################
#######################################################################################################################################

'''
    The 'data_analysis' function inputs two flags, either 'yes' or 'no'. The first flag decides if the threshold plots will be shown. The second flag decides if the fitting results will be shown.

    CHANGE THE FIRST INPUT, WHICH IS THE FILE NAME, TO READ IT!
'''
threshold1, time_weight1 = data_analysis('AP_decoder/Time_weight_function/AP_decoder_synchronicity=0.0', 'n')

plt.rcParams['mathtext.fontset'] = 'stix'

fig, ax = plt.subplots()
ax.plot(time_weight1, 100.*np.array(threshold1), 'ko', markersize=5)
ax.set_xlabel(r'Time weight $\mathregular{w_{time}^{AP}}$', size=13)
ax.set_ylabel('Threshold in physical error p (%)', size=13)
ax.set_xlim(xmin=0.08,xmax=1.02)
ax.set_ylim(ymin=0.5,ymax=1.35)
ax.locator_params(axis='x', nbins=10)
ax.locator_params(axis='y', nbins=10)
ax.tick_params(direction='inout', labelright=1, right=1, labeltop=1, top=1)
ax.tick_params(axis="x", labelsize=11)
ax.tick_params(axis="y", labelsize=11)
#legend = ax.legend(loc='upper left', shadow=True)
plt.grid(linestyle = '--', linewidth = 0.3)
plt.savefig('///home/joao/average_decoder2.png', bbox_inches='tight')
plt.show()
plt.close(fig)


