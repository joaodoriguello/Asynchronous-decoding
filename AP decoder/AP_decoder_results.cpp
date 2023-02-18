#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "omp.h"
#include "AP_decoder.h"

int main (void) {

	srand(time(NULL));

    	int precision = 4000;
    	int L[3] = {10, 12, 14};
    	double synchronicity = 0.1;
	int time_factor = 2;
	double time_weight = 0.5;
	double p_min = 0.0095, p_max = 0.013;	


	double result, physical_error[16];
	for (int i = 0; i < 16; i++) 
		physical_error[i] = p_min + i*(p_max - p_min)/15;


	for (int j = 0; j < 16; j++) {
		for (int i = 0; i < 3; i++) {
			result = AP_decoder(L[i], 0.5*(1 - pow(1 - 2*physical_error[j], synchronicity)), physical_error[j], synchronicity, int(round(time_factor/synchronicity))*L[i], precision, time_weight);
			printf("%d %lf %lf\n", L[i], physical_error[j], result);
		}
	}

    	return 0;
}
