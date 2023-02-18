import numpy as np


with open("simple--async=0.5", 'r') as f: #INPUT
	input_data = f.readlines()


output_data = []
synchronicity = 1.0

for i in range(len(input_data)):

	if 'synchronicity' in input_data[i]:
		synchronicity = float((input_data[i].rstrip('\n')).split("=")[1])
		output_data += [input_data[i]]
		
	elif '#' in input_data[i] and synchronicity != 0.0:  
		probability = float((input_data[i].rstrip('\n')).split(" ")[3])
		simulation_error = float((input_data[i].rstrip('\n')).split(" ")[2])
		physical_error = 0.5*(1 - pow(1 - 2*simulation_error, 1/synchronicity))
		output_data += [(input_data[i].rstrip('\n')).split(" ")[0] + ' ' + (input_data[i].rstrip('\n')).split(" ")[1] + ' ' + f'{physical_error:.6f}' + ' ' + f'{probability:.6f}' + '\n']
	
	elif '#' in input_data[i] and synchronicity == 0.0:  
		probability = float((input_data[i].rstrip('\n')).split(" ")[3])
		output_data += [(input_data[i].rstrip('\n')).split(" ")[0] + ' ' + (input_data[i].rstrip('\n')).split(" ")[1] + ' ' + (input_data[i].rstrip('\n')).split(" ")[2] + ' ' + f'{probability:.6f}' + '\n']
	else:
		output_data += [input_data[i]]




with open("simple--async=0.5--new", 'w') as f:
	f.writelines(output_data)

