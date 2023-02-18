#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "omp.h"
#include "P2_and_P1_by_P0.h"

int main (void) {

	srand(time(NULL));

    	int precision[5] = {250, 20, 5, 2, 1};
    	int L[5] = {6, 8, 10, 12, 14};
    	double synchronicity = 0.1;
	int time_factor = 2;
	double physical_error = 0.0017848;

	double P1_by_P0[5], P2_by_P0[5];
	long int average[5];
		
    	for (int i = 0; i < 5; i++) {
		std::tie(average[i], P1_by_P0[i], P2_by_P0[i]) = P2_and_P1_by_P0(L[i], physical_error, 0.5*(1 - pow(1 - 2*physical_error, 1/synchronicity)), synchronicity, int(round(time_factor/synchronicity))*L[i], precision[i]);
		printf("%d %ld %lf %lf\n", L[i], average[i], P1_by_P0[i]/average[i], P2_by_P0[i]/average[i]);	
	}

    	return 0;
}
