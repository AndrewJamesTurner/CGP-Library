#include <stdio.h>
#include <stdlib.h> 
#include "cgp.h"

#include "fullAdder.c"


int main(void){

	struct parameters *params;
	struct population *pop;
	struct data *dat;
		
	params = initialiseParameters(3,6,2, 2);
	setFuctionSet(params, "add");
	
	dat = initialiseDataFromFile("./example/fullAdder");
		
	pop = initialisePopulation(params);
	
	
	/*setFitnessFuction(params, fullAdder);*/
	
	
		
	/*printChromosome(params, chromo);*/
	/*
	mutateChromosome(params, chromo);
	printChromosome(params, chromo);
	*/
	
	/*
	setMu(params, 1);
	printf("\nmu: %d\n", getMu(params));
	*/
	
	
	/*executeChromosome(params, chromo, inputs, outputs);*/
	
	/*printf("%f %f\n", outputs[0], outputs[1] );*/
	
	
	
	/*freeChromosome(params, chromo);*/
	
	
	freePopulation(params, pop);
	freeData(dat);
	freeParameters(params);
	
	return 1;
}
