#include <stdio.h>
#include "libCGP.h"

int main(void){

	struct parameters *params;
	struct chromosome *chromo;

	params = initialiseParameters();
	chromo = initialiseChromosome(params);

	printChromosome(params, chromo);
	
	mutateChromosome(params, chromo);
	
	printChromosome(params, chromo);
	/*
	setMu(params, 1);
	printf("\nmu: %d\n", getMu(params));
	*/
	
	return 1;
}
