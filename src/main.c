#include <stdio.h>
#include "cgp.h"

int main(void){

	struct parameters *params;
	struct chromosome *chromo;
	
	params = initialiseParameters(3,10,1, 2);
	
	setFuctionSet(params, "add,add");
	printFuctionSet(params);
	
	setFuctionSet(params, "add,sub ");
	printFuctionSet(params);
	
	chromo = initialiseChromosome(params);
	
	
	
	printChromosome(params, chromo);
	/*
	mutateChromosome(params, chromo);
	printChromosome(params, chromo);
	*/
	
	/*
	setMu(params, 1);
	printf("\nmu: %d\n", getMu(params));
	*/
	
	freeChromosome(params, chromo);
	freeParameters(params);
	
	return 1;
}
