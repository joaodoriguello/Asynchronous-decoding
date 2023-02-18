#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "omp.h"
#include "CG_decoder_degeneracy_continuous.h"

int main (void) {
	
	srand(time(NULL));

    	int precision = 500;
    	int L[3] = {10, 12, 14};
	int time_factor = 2;
	double p_min = 0.016, p_max = 0.0177;
	double tau = 1.0;


	double result, physical_error[16];
	for (int i = 0; i < 16; i++) 
		physical_error[i] = p_min + i*(p_max - p_min)/15;


	for (int j = 0; j < 16; j++) {
    		for (int i = 0; i < 3; i++) {
			result = CG_decoder_degeneracy_continuous(L[i], physical_error[j], physical_error[j], L[i]*time_factor, precision, tau);
			printf("%d %f %f\n", L[i], physical_error[j], result);
		}
	}

    	return 0;
}
