#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "omp.h"
#include "AP_decoder_continuous.h"

int main (void) {
	
	srand(time(NULL));

    	int precision = 4000;
    	int L[3] = {10, 12, 14};
	int time_factor = 2;
	double time_weight = 0.5;
	double p_min = 0.0123, p_max = 0.014;


	double result, physical_error[16];
	for (int i = 0; i < 16; i++) 
		physical_error[i] = p_min + i*(p_max - p_min)/15;

	
	for (int j = 0; j < 16; j++) {
    		for (int i = 0; i < 3; i++) {
			result = AP_decoder_continuous(L[i], physical_error[j], physical_error[j], L[i]*time_factor, precision, time_weight);
			printf("%d %lf %lf\n", L[i], physical_error[j], result);
		}
	}

    	return 0;
}

