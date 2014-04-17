#include <stdio.h>
#include "../include/cgp.h" 

int main(void){

	struct parameters *params = NULL;
	struct dataSet *trainingData = NULL;
	struct results *rels = NULL;

	int numInputs = 4;
	int numNodes = 50;
	int numOutputs = 1;
	int nodeArity = 2;

	int numGens = 500000;
	int numRuns = 20;
	
	
	params = initialiseParameters(numInputs, numNodes, numOutputs, nodeArity);

	addNodeFunction(params, "or,nor,and,nand");

	setMutationType(params, "point");

	setMutationRate(params, 0.03);

	/*printParameters(params);*/

	trainingData = initialiseDataSetFromFile("parity4bit.data");

	/*printDataSet(trainingData);*/
	
	/*runCGP(params, trainingData, numGens);*/

	rels = repeatCGP(params, trainingData, numGens, numRuns);

	freeDataSet(trainingData);
	freeResults(rels);
	freeParameters(params);

	return 0;
}
