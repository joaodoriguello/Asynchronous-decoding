import numpy as np


with open("3D_tau=0.0", 'r') as f: #INPUT
    data1 = f.readlines()

with open("3D_t=0--(10,12,14,16)", 'r') as f: #MODIFIER
    data2 = f.readlines()


for i in range(len(data1)):
    if 'precision' in data1[i]:
	precision1 = int((data1[i].rstrip('\n')).split("=")[1])
	j = 5
	precision2 = int((data2[j].rstrip('\n')).split("=")[1])
	data2[j] = 'precision=' + str(precision1 + precision2) + '\n'
	j += 1
	i += 1

	while '!' not in data1[i]:
	    data2[j] = (data1[i].rstrip('\n')).split(" ")[0] + ' ' + (data1[i].rstrip('\n')).split(" ")[1] + ' ' + (data1[i].rstrip('\n')).split(" ")[2] + ' ' + str( (precision1*float((data1[i].rstrip('\n')).split(" ")[3]) + precision2*float((data2[j].rstrip('\n')).split(" ")[3]))/(precision1 + precision2) ) + '\n'
	    i += 1
	    j += 1



with open("3D_t=0--(10,12,14,16)", 'w') as f:
    f.writelines(data2)

