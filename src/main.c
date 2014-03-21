#include <stdio.h>
#include <stdlib.h> 
#include "include/cgp.h" 


int main(void){

	struct parameters *params;
	struct population *pop;
	struct data *trainingData;
	struct chromosome *chromo;	
		
	int numInputs = 3;
	int numNodes = 10;
	int numOutputs = 2;
	int nodeArity = 2;
	
	float inputs[8][3] = {{0,0,0},{0,0,1},{0,1,0},{0,1,1},{1,0,0},{1,0,1},{1,1,0},{1,1,1}};
	float outputs[8][2] = {{0,0},{1,0},{1,0},{0,1},{1,0},{0,1},{0,1},{1,1}};	
		
	
	params = initialiseParameters(numInputs, numNodes, numOutputs, nodeArity);
			
	addNodeFunction(params, "and,or,nand,nor,xor");
	printFunctionSet(params);
	
	pop = initialisePopulation(params);
		

	/*trainingData = initialiseDataFromFile("./example/fullAdder");*/	
	trainingData = initialiseDataFromArrays(3,2,8, inputs[0], outputs[0]);
	
		
	evolvePopulation(params, pop, trainingData);
		
	
	chromo = getFittestChromosome(params, pop);
	printf("Best Fitness: %f found after %d generations. which used %d active nodes\n", getChromosomeFitness(chromo), getNumberOfGenerations(pop), getChromosomeNumActiveNodes(chromo));
	
	
	
	printChromosome(params,chromo);
	
	
	freePopulation(pop);
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
