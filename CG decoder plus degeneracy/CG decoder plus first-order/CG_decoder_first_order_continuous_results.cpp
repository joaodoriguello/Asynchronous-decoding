#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "omp.h"
#include "CG_decoder_first_order_continuous.h"

int main (void) {
	
	srand(time(NULL));

    	int precision = 1000;
    	int L[3] = {10, 12, 14};
	double tau = 1.0;
	int time_factor = 3;
	double p_min = 0.0161, p_max = 0.018;


	double result, physical_error[16];
	for (int i = 0; i < 16; i++) 
		physical_error[i] = p_min + i*(p_max - p_min)/15;

	
	for (int j = 0; j < 16; j++) {
    		for (int i = 0; i < 3; i++) {
			result = CG_decoder_first_order_continuous(L[i], physical_error[j], physical_error[j], L[i]*time_factor, precision, tau);
			printf("%d %f %f\n", L[i], physical_error[j], result);
		}
	}

    	return 0;
}
