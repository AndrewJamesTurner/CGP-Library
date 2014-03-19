#include <stdio.h>
#include <stdlib.h> 
#include "cgp.h"


int main(void){

	struct parameters *params;
	struct population *pop;
	struct data *trainingData;
	
	float bestFit;
	
	int numInputs = 3;
	int numNodes = 50;
	int numOutputs = 2;
	int nodeArity = 2;
	
	params = initialiseParameters(numInputs, numNodes, numOutputs, nodeArity);
	setFuctionSet(params, "and,or,nand,nor,xor");
	
	pop = initialisePopulation(params);
		
	trainingData = initialiseDataFromFile("./example/fullAdder");
		
	bestFit = evolvePopulation(params, pop, trainingData);
	
	printf("bestFit: %f\n", bestFit);
	
	
	freePopulation(params, pop);
	freeData(trainingData);
	freeParameters(params);
	
	
	
	
	/*setFitnessFuction(params, fullAdder, "fullAdder");*/
	
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
	
	
	
	
	return 1;
}
