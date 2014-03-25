#include <stdio.h>


#include "../include/cgp.h" 


int main(void){

	struct parameters *params = NULL;
	struct data *trainingData = NULL;
	struct chromosome *chromo = NULL;
		
	int numInputs = 3;
	int numNodes = 50;
	int numOutputs = 2;
	int nodeArity = 2;
	
	int numGens = 10000;
	
	params = initialiseParameters(numInputs, numNodes, numOutputs, nodeArity);
			
	addNodeFunction(params, "and,or,not,xor, nor, nand");
		
	setUpdateFrequency(params, 250);	
		
	printParameters(params);
	
	trainingData = initialiseDataFromFile("./examples/fullAdder.data");
	
	
	chromo = runCGP(params, trainingData, numGens);	
	
	/*printChromosome(chromo);*/
	
	freeChromosome(chromo);		
	freeParameters(params);		
			
	return 1;
}
