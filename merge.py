import numpy as np


with open("simple--async=0.9", 'r') as f: #INPUT
    data1 = f.readlines()

with open("simple--async=0.9 (new)", 'r') as f: #MODIFIER
    data2 = f.readlines()


for i in range(len(data1)):
#	if 'async_error' in data1[i] and data1[i] in data2:
#	if 'tau' in data1[i] and data1[i] in data2:
	if 'time_weight' in data1[i] and data1[i] in data2:
#	if 'weight' in data1[i] and data1[i] in data2:
		precision1 = int((data1[i-1].rstrip('\n')).split("=")[1])
		j = data2.index(data1[i])
		precision2 = int((data2[j-1].rstrip('\n')).split("=")[1])
		data2[j-1] = 'precision=' + str(precision1 + precision2) + '\n'
		j += 1
		i += 1

		while '!' not in data1[i]:
	    		data2[j] = (data1[i].rstrip('\n')).split(" ")[0] + ' ' + (data1[i].rstrip('\n')).split(" ")[1] + ' ' + (data1[i].rstrip('\n')).split(" ")[2] + ' ' + str( (precision1*float((data1[i].rstrip('\n')).split(" ")[3]) + precision2*float((data2[j].rstrip('\n')).split(" ")[3]))/(precision1 + precision2) ) + '\n'
	    		i += 1
	    		j += 1



with open("simple--async=0.9 (new)", 'w') as f:
    f.writelines(data2)

