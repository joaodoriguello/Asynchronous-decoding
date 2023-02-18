#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "omp.h"
#include "P2_and_P1_by_P0_continuous.h"

int main (void) {

	srand(time(NULL));

    	int precision[5] = {250, 20, 5, 2, 1};
    	int L[5] = {6, 8, 10, 12, 14};
	int time_factor = 2;
	double physical_error = 0.016947;

	double P1_by_P0[5], P2_by_P0[5];
	long int average[5];
	
    	for (int i = 0; i < 5; i++) {
		std::tie(average[i], P1_by_P0[i], P2_by_P0[i]) = P2_and_P1_by_P0_continuous(L[i], physical_error, physical_error, L[i]*time_factor, precision[i]);
		printf("%d %ld %lf %lf\n", L[i], average[i], P1_by_P0[i]/average[i], P2_by_P0[i]/average[i]);
	}

    	return 0;
}

