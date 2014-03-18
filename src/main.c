#include <stdio.h>
#include "cgp.h"

int main(void){

	struct parameters *params;
	struct chromosome *chromo;
	
	float inputs[3] = {0,1,2};
	float outputs[2];
	
	params = initialiseParameters(3,6,2, 2);
	
	setFuctionSet(params, "add,add");
	printFuctionSet(params);
	
	setFuctionSet(params, "add");
	printFuctionSet(params);
	
	chromo = initialiseChromosome(params);
	
	
	printf("%f, %f, %f\n", inputs[0], inputs[1], inputs[2]   );
		
	printChromosome(params, chromo);
	/*
	mutateChromosome(params, chromo);
	printChromosome(params, chromo);
	*/
	
	/*
	setMu(params, 1);
	printf("\nmu: %d\n", getMu(params));
	*/
	
	
	executeChromosome(params, chromo, inputs, outputs);
	
	printf("%f %f\n", outputs[0], outputs[1] );
	
	
	
	freeChromosome(params, chromo);
	
	
	
	freeParameters(params);
	
	return 1;
}
