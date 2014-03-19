#include <stdio.h>
#include <stdlib.h> 
#include "cgp.h"




int main(void){

	struct parameters *params;
	struct population *pop;
	struct data *dat;
	float bestFit;
		
	params = initialiseParameters(3,10,2, 2);
	setFuctionSet(params, "and");
	
	dat = initialiseDataFromFile("./example/fullAdder");
		
	pop = initialisePopulation(params);
	
	
	
	bestFit =  evolvePopulation(params, pop, dat);
	
	printf("bestFit: %f\n", bestFit);
	
	/*setFitnessFuction(params, fullAdder, "fullAdder");*/
	
	
		
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
