#include <stdio.h>
#include <stdlib.h> 
#include <math.h>
#include "include/cgp.h" 

float symbolicEq1(float x){
	return powf(x,6.0) - 2*powf(x,4.0) + powf(x,2.0);
}


float symbolicRegression1(struct parameters *params, struct chromosome *chromo, struct data *dat){
	
	float i;
	float error = 0;
	float chromoInputs[1];	
	float chromoOutputs[1];	
						
	if(getNumInputs(params) != 1){
		printf("Error: The 'symbolicRegression1' fitness function requires one chromosome input; not %d\n", getNumInputs(params));
		exit(0);
	}				
		
	if(getNumOutputs(params) != 1){
		printf("Error: The 'symbolicRegression1' fitness function requires one chromosome output; not %d\n", getNumOutputs(params));
		exit(0);
	}		
					
	/* for each line in the truth table */				
	for(i=-5; i<=5; i=i+0.1){
		
		chromoInputs[0] = i;
		
		/* calculate the chromosome outputs for the set of inputs  */
		executeChromosome(chromo, chromoInputs, chromoOutputs);
		
		error = fabs(symbolicEq1(i) - chromoOutputs[0]);
	}				
					
	return error;
}



int main(void){

	struct parameters *params = NULL;
	struct results *rels = NULL;
	struct chromosome *chromo = NULL;
		
	int numInputs = 1;
	int numNodes = 20;
	int numOutputs = 1;
	int nodeArity = 2;
	
	int numGens = 5000;
	int numRuns = 10;
	
	params = initialiseParameters(numInputs, numNodes, numOutputs, nodeArity);
			
	/*setNumGenerations(params, numGens);		*/
			
	setTargetFitness(params, 0.1);		
			
	addNodeFunction(params, "add,sub,mul,div");
	
	setFitnessFunction(params, symbolicRegression1, "symBol1" );
	
	setMutationType(params, "probabilistic");
	setMutationRate(params, 0.3);
	
	printParameters(params);
	exit(0);
	
	rels = repeatCGP(params, NULL, numGens, numRuns);	
	
	/*	
	printf("\n");
	printChromosome(getChromosome(rels,0));
		
	saveChromosome(getChromosome(rels,0), "test.chromo");
	chromo = loadChromosome("test.chromo");
	printChromosome(chromo);		
		*/	
	free(chromo);		
	freeResults(rels);
	freeParameters(params);		
		
		
		
		
	/*printFunctionSet(params);*/	
		
	/*struct data *trainingData;*/
	/*struct chromosome *chromo = NULL;	*/	
		
		
		/*
	float inputs[8][3] = {{0,0,0},{0,0,1},{0,1,0},{0,1,1},{1,0,0},{1,0,1},{1,1,0},{1,1,1}};
	float outputs[8][2] = {{0,0},{1,0},{1,0},{0,1},{1,0},{0,1},{0,1},{1,1}};	
	*/	
		
			
	/*trainingData = initialiseDataFromFile("./example/fullAdder");*/	
	/*trainingData = initialiseDataFromArrays(3,2,8, inputs[0], outputs[0]);*/
	
	/*chromo = runCGP(params, NULL);
	printf("Chromo fitness: %f\n", getChromosomeFitness(chromo));
	*/				
			
			
			
	/*chromo = getFittestChromosome(params, pop);*/
	/*printf("Best Fitness: %f found after %d generations. which used %d active nodes\n", getChromosomeFitness(chromo), getNumberOfGenerations(pop), getChromosomeNumActiveNodes(chromo));*/
	
	
	
	/*printChromosome(params,chromo);*/
	
	/*freeChromosome(chromo);*/
	/*freePopulation(pop);*/
	/*freeData(trainingData);*/
	
	

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
