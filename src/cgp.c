/*
	This file is part of CGP-Library
	Copyright (c) Andrew James Turner 2014 (andrew.turner@york.ac.uk)

    CGP-Library is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CGP-Library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with CGP-Library.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <float.h>

#include "../include/cgp.h"

#define FUNCTIONSETSIZE 50
#define FUNCTIONNAMELENGTH 11
#define FITNESSFUNCTIONNAMELENGTH 21
#define MUTATIONTYPENAMELENGTH 21
#define SELECTIONSCHEMENAMELENGTH 21
#define REPRODUCTIONSCHEMENAMELENGTH 21

/*
	Structure definitions
*/

struct parameters{
	int mu;
	int lambda;
	char evolutionaryStrategy;
	float mutationRate;
	float connectionWeightRange;
	int numInputs;
	int numNodes;
	int numOutputs;
	int arity;
	struct functionSet *funcSet;
	float targetFitness;
	int updateFrequency;
	void (*mutationType)(struct parameters *params, struct chromosome *chromo);
	char mutationTypeName[MUTATIONTYPENAMELENGTH];
	float (*fitnessFunction)(struct parameters *params, struct chromosome *chromo, struct dataSet *dat);
	char fitnessFunctionName[FITNESSFUNCTIONNAMELENGTH];
	void (*selectionScheme)(struct parameters *params, struct chromosome **parents, struct chromosome **candidateChromos, int numParents, int numCandidateChromos);
	char selectionSchemeName[SELECTIONSCHEMENAMELENGTH];
	void (*reproductionScheme)(struct parameters *params, struct chromosome **parents, struct chromosome **children, int numParents, int numChildren);
	char reproductionSchemeName[REPRODUCTIONSCHEMENAMELENGTH];
};

struct chromosome{
	int numInputs;
	int numOutputs;
	int numNodes;
	int numActiveNodes;
	int arity;
	struct node **nodes;
	int *outputNodes;
	int *activeNodes;
	float fitness;
	float *outputValues;
	struct functionSet *funcSet;
	float *nodeInputsHold;
	int generation;
};

struct node{
	int function;
	int *inputs;
	float *weights;
	int active;
	float output;
	int arity;
};

struct functionSet{
	int numFunctions;
	char functionNames[FUNCTIONSETSIZE][FUNCTIONNAMELENGTH];
	float (*functions[FUNCTIONSETSIZE])(const int numInputs, const float *inputs, const float *connectionWeights);
};

struct dataSet{
	int numSamples;
	int numInputs;
	int numOutputs;
	float **inputData;
	float **outputData;
};

struct results{
	int numRuns;
	struct chromosome **bestChromosomes;
};


/*
	Prototypes of functions used internally to CGP-Library
*/

/* chromosome functions */
static void setChromosomeActiveNodes(struct chromosome *chromo);
static void recursivelySetActiveNodes(struct chromosome *chromo, int nodeIndex);
static void sortChromosomeArray(struct chromosome **chromoArray, int numChromos);
static int cmpChromosome(const void *a, const void *b);

/* node functions */
static struct node *initialiseNode(int numInputs, int arity, int numFunctions, float connectionWeightRange, int nodePosition);
static void freeNode(struct node *n);
static void copyNode(struct node *nodeDest, struct node *nodeSrc);

/* getting gene value functions  */
static float getRandomConnectionWeight(float weightRange);
static int getRandomNodeInput(int numChromoInputs, int nodePosition);
static int getRandomFunction(int numFunctions);
static int getRandomChromosomeOutput(int numInputs, int numNodes);

/* function set functions */
static int addPresetFuctionToFunctionSet(struct parameters *params, char *functionName);
static void copyFuctionSet(struct functionSet *funcSetDest, struct functionSet *funcSetSrc);
static void printFunctionSet(struct parameters *params);

/* results functions */
struct results* initialiseResults(struct parameters *params, int numRuns);

/* mutation functions  */
static void probabilisticMutation(struct parameters *params, struct chromosome *chromo);
static void pointMutation(struct parameters *params, struct chromosome *chromo);
static void probabilisticMutationOnlyActive(struct parameters *params, 
struct chromosome *chromo);

/* selection scheme functions */
static void selectFittest(struct parameters *params, struct chromosome **parents, struct chromosome **candidateChromos, int numParents, int numCandidateChromos);

/* reproduction scheme functions */
static void mutateRandomParent(struct parameters *params, struct chromosome **parents, struct chromosome **children, int numParents, int numChildren);

/* fitness function */
static float supervisedLearning(struct parameters *params, struct chromosome *chromo, struct dataSet *data);

/* node functions defines in CGP-Library */
static float add(const int numInputs, const float *inputs, const float *connectionWeights);
static float sub(const int numInputs, const float *inputs, const float *connectionWeights);
static float mul(const int numInputs, const float *inputs, const float *connectionWeights);
static float divide(const int numInputs, const float *inputs, const float *connectionWeights);
static float and(const int numInputs, const float *inputs, const float *connectionWeights);
static float absolute(const int numInputs, const float *inputs, const float *connectionWeights);
static float squareRoot(const int numInputs, const float *inputs, const float *connectionWeights);
static float square(const int numInputs, const float *inputs, const float *connectionWeights);
static float cube(const int numInputs, const float *inputs, const float *connectionWeights);
static float exponential(const int numInputs, const float *inputs, const float *connectionWeights);
static float sine(const int numInputs, const float *inputs, const float *connectionWeights);
static float cosine(const int numInputs, const float *inputs, const float *connectionWeights);
static float tangent(const int numInputs, const float *inputs, const float *connectionWeights);
static float nand(const int numInputs, const float *inputs, const float *connectionWeights);
static float or(const int numInputs, const float *inputs, const float *connectionWeights);
static float nor(const int numInputs, const float *inputs, const float *connectionWeights);
static float xor(const int numInputs, const float *inputs, const float *connectionWeights);
static float xnor(const int numInputs, const float *inputs, const float *connectionWeights);
static float not(const int numInputs, const float *inputs, const float *connectionWeights);
static float sigmoid(const int numInputs, const float *inputs, const float *connectionWeights);
static float gaussian(const int numInputs, const float *inputs, const float *connectionWeights);
static float step(const int numInputs, const float *inputs, const float *connectionWeights);
static float softsign(const int numInputs, const float *inputs, const float *connectionWeights);
static float hyperbolicTangent(const int numInputs, const float *inputs, const float *connectionWeights);

/* other */
static float randFloat(void);
int randInt(int n);
static void sortIntArray(int *array, const int length);
static float sumWeigtedInputs(const int numInputs, const float *inputs, const float *connectionWeights);
static int cmpInt(const void * a, const void * b);


/*
	parameters function definitions
*/

/*
	Initialises a parameter struct with default values. These
	values can be individually changed via set functions.
*/
DLL_EXPORT struct parameters *initialiseParameters(const int numInputs, const int numNodes, const int numOutputs, const int arity){

	struct parameters *params;

	/* allocate memory for parameters */
	params = malloc(sizeof(struct parameters));

	/* Set default values */
	params->mu = 1;
	params->lambda = 4;
	params->evolutionaryStrategy = '+';
	params->mutationRate = 0.05;
	params->connectionWeightRange = 1;

	params->targetFitness = 0;

	params->updateFrequency = 1;

	setNumInputs(params, numInputs);
	setNumNodes(params, numNodes);
	setNumOutputs(params, numOutputs);
	setArity(params, arity);

	params->mutationType = probabilisticMutation;
	strncpy(params->mutationTypeName, "probabilistic", MUTATIONTYPENAMELENGTH);

	params->funcSet = malloc(sizeof(struct functionSet));
	params->funcSet->numFunctions = 0;

	params->fitnessFunction = supervisedLearning;
	strncpy(params->fitnessFunctionName, "supervisedLearning", FITNESSFUNCTIONNAMELENGTH);

	params->selectionScheme = selectFittest;
	strncpy(params->selectionSchemeName, "selectFittest", SELECTIONSCHEMENAMELENGTH);

	params->reproductionScheme = mutateRandomParent;
	strncpy(params->reproductionSchemeName, "mutateRandomParent", REPRODUCTIONSCHEMENAMELENGTH);

	/* Seed the random number generator */
	srand(time(NULL));

	return params;
}


/*
	Frees the memory associated with the given parameter structure
*/
DLL_EXPORT void freeParameters(struct parameters *params){

	/* attempt to prevent user double freeing */
	if(params == NULL){
		return;
	}

	free(params->funcSet);
	free(params);
}

/*
	prints the given parameters to the terminal
*/
DLL_EXPORT void printParameters(struct parameters *params){

	if(params == NULL){
		printf("Error: cannot print uninitialised parameters.\nTerminating CGP-Library.\n");
		exit(0);
	}

	printf("---------------------------------------------------\n");
	printf("                   Parameters                      \n");
	printf("---------------------------------------------------\n");
	printf("Evolutionary Strategy:\t\t(%d%c%d)-ES\n", params->mu, params->evolutionaryStrategy, params->lambda);
	printf("Inputs:\t\t\t\t%d\n", params->numInputs);
	printf("Nodes:\t\t\t\t%d\n", params->numNodes);
	printf("Outputs:\t\t\t%d\n", params->numOutputs);
	printf("Node Arity:\t\t\t%d\n", params->arity);
	printf("Connection weights range:\t+/- %f\n", params->connectionWeightRange);
	printf("Mutation Type:\t\t\t%s\n", params->mutationTypeName);
	printf("Mutation rate:\t\t\t%f\n", params->mutationRate);
	printf("Fitness Function:\t\t%s\n", params->fitnessFunctionName);
	printf("Target Fitness:\t\t\t%f\n", params->targetFitness);
	printf("Selection scheme:\t\t%s\n", params->selectionSchemeName);
	printf("Reproduction scheme:\t\t%s\n", params->reproductionSchemeName);
	printf("Update frequency:\t\t%d\n", params->updateFrequency);
	printFunctionSet(params);
	printf("---------------------------------------------------\n\n");
}


/*
	Sets the given function set to contain the per-set functions
	given in the char array. The function names must be comma separated
	and contain no spaces i.e. "and,or".
*/
DLL_EXPORT void addNodeFunction(struct parameters *params, char *functionNames){

	char *pch;
	char funcNames[FUNCTIONNAMELENGTH * FUNCTIONSETSIZE];

	/* make a local copy of the function names*/
	strncpy(funcNames, functionNames, FUNCTIONNAMELENGTH * FUNCTIONSETSIZE);

	/* get the first function name */
	pch = strtok(funcNames, ", ");

	/* while the function names char array contains function names */
	while (pch != NULL){

		/* add the named function to the function set */
		addPresetFuctionToFunctionSet(params, pch);

		/* get the next function name */
		pch = strtok(NULL, ", ");
	}

	/* if the function set is empty give warning */
	if(params->funcSet->numFunctions == 0){
		printf("Warning: No Functions added to function set.\n");
	}
}


/*
	Adds given node function to given function set with given name.
	Disallows exceeding the function set size.
*/
DLL_EXPORT void addNodeFunctionCustom(struct parameters *params, float (*function)(const int numInputs, const float *inputs, const float *weights), char *functionName){

	if(params->funcSet->numFunctions >= FUNCTIONSETSIZE){
		printf("Warning: functions set has reached maximum capacity (%d). Function '%s' not added.\n", FUNCTIONSETSIZE, functionName);
		return;
	}

	/* */
	params->funcSet->numFunctions++;

	/* */
	strncpy(params->funcSet->functionNames[params->funcSet->numFunctions-1], functionName, FUNCTIONNAMELENGTH);

	/* */
	params->funcSet->functions[params->funcSet->numFunctions-1] = function;
}


/*
	clears the given function set of functions
*/
DLL_EXPORT void clearFunctionSet(struct parameters *params){
	params->funcSet->numFunctions = 0;
}


/*
	sets num chromosome inputs in parameters
*/
DLL_EXPORT void setNumInputs(struct parameters *params, int numInputs){

	/* error checking */
	if(numInputs < 0){
		printf("Error: number of chromosome inputs cannot be negative; %d is invalid.\nTerminating CGP-Library.\n", numInputs);
		exit(0);
	}

	params->numInputs = numInputs;
}


/*
	sets num chromosome nodes in parameters
*/
DLL_EXPORT void setNumNodes(struct parameters *params, int numNodes){

	/* error checking */
	if(numNodes < 0){
		printf("Warning: number of chromosome nodes cannot be negative; %d is invalid.\nTerminating CGP-Library.\n", numNodes);
		exit(0);
	}

	params->numNodes = numNodes;
}

/*
	sets num chromosome outputs in parameters
*/
DLL_EXPORT void setNumOutputs(struct parameters *params, int numOutputs){

	/* error checking */
	if(numOutputs < 0){
		printf("Warning: number of chromosome outputs cannot be less than one; %d is invalid.\nTerminating CGP-Library.\n", numOutputs);
		exit(0);
	}

	params->numOutputs = numOutputs;
}


/*
	sets chromosome arity in parameters
*/
DLL_EXPORT void setArity(struct parameters *params, int arity){

	/* error checking */
	if(arity < 0){
		printf("Warning: node arity cannot be less than one; %d is invalid.\nTerminating CGP-Library.\n", arity);
		exit(0);
	}

	params->arity = arity;
}




/*
	Sets the mu value in given parameters to the new given value. If mu value
	is invalid a warning is displayed and the mu value is left unchanged.
*/
DLL_EXPORT void setMu(struct parameters *params, int mu){

	if(mu > 0){
		params->mu = mu;
	}
	else{
		printf("\nWarning: mu value '%d' is invalid. Mu value must have a value of one or greater. Mu value left unchanged as '%d'.\n", mu, params->mu);
	}
}


/*
	Sets the lambda value in given parameters to the new given value.
	If lambda value is invalid a warning is displayed and the lambda value
	is left unchanged.
*/
DLL_EXPORT void setLambda(struct parameters *params, int lambda){

	if(lambda > 0){
		params->lambda = lambda;
	}
	else{
		printf("\nWarning: lambda value '%d' is invalid. Lambda value must have a value of one or greater. Lambda value left unchanged as '%d'.\n", lambda, params->lambda);
	}
}


/*
	Sets the evolutionary strategy given in parameters to '+' or ','.
	If an invalid option is given a warning is displayed and the evolutionary
	strategy is left unchanged.
*/
DLL_EXPORT void setEvolutionaryStrategy(struct parameters *params, char evolutionaryStrategy){

	if(evolutionaryStrategy == '+' || evolutionaryStrategy == ','){
		params->evolutionaryStrategy = evolutionaryStrategy;
	}
	else{
		printf("\nWarning: the evolutionary strategy '%c' is invalid. The evolutionary strategy must be '+' or ','. The evolutionary strategy has been left unchanged as '%c'.\n", evolutionaryStrategy, params->evolutionaryStrategy);
	}
}


/*
	Sets the mutation rate given in parameters. IF an invalid mutation
	rate is given a warning is displayed and the mutation rate is left
	unchanged.
*/
DLL_EXPORT void setMutationRate(struct parameters *params, float mutationRate){

	if(mutationRate >= 0 && mutationRate <= 1){
		params->mutationRate = mutationRate;
	}
	else{
		printf("\nWarning: mutation rate '%f' is invalid. The mutation rate must be in the range [0,1]. The mutation rate has been left unchanged as '%f'.\n", mutationRate, params->mutationRate);
	}
}


/*
	Sets the connection weight range given in parameters.
*/
DLL_EXPORT void setConnectionWeightRange(struct parameters *params, float weightRange){

	params->connectionWeightRange = weightRange;
}


/*
	sets the fitness function to the fitnessFuction passed. If the fitnessFuction is NULL
	then the default supervisedLearning fitness function is used.
*/
DLL_EXPORT void setFitnessFunction(struct parameters *params, float (*fitnessFunction)(struct parameters *params, struct chromosome *chromo, struct dataSet *data), char *fitnessFunctionName){

	if(fitnessFunction == NULL){
		params->fitnessFunction = supervisedLearning;
		strncpy(params->fitnessFunctionName, "supervisedLearning", FITNESSFUNCTIONNAMELENGTH);
	}
	else{
		params->fitnessFunction = fitnessFunction;
		strncpy(params->fitnessFunctionName, fitnessFunctionName, FITNESSFUNCTIONNAMELENGTH);
	}
}



/*
	sets the selection scheme used to select the parents from the candidate chromosomes. If the selectionScheme is NULL
	then the default selectFittest selection scheme is used.
*/
DLL_EXPORT void setSelectionScheme(struct parameters *params, void (*selectionScheme)(struct parameters *params, struct chromosome **parents, struct chromosome **candidateChromos, int numParents, int numCandidateChromos), char *selectionSchemeName){

	if(selectionScheme == NULL){
		params->selectionScheme = selectFittest;
		strncpy(params->selectionSchemeName, "selectFittest", SELECTIONSCHEMENAMELENGTH);
	}
	else{
		params->selectionScheme = selectionScheme;
		strncpy(params->selectionSchemeName, selectionSchemeName, SELECTIONSCHEMENAMELENGTH);
	}
}


/*
	sets the reproduction scheme used to select the parents from the candidate chromosomes. If the reproductionScheme is NULL
	then the default mutateRandomParent selection scheme is used.
*/

DLL_EXPORT void setReproductionScheme(struct parameters *params, void (*reproductionScheme)(struct parameters *params, struct chromosome **parents, struct chromosome **children, int numParents, int numChildren), char *reproductionSchemeName){

	if(reproductionScheme == NULL){
		params->reproductionScheme = mutateRandomParent;
		strncpy(params->reproductionSchemeName, "mutateRandomParent", REPRODUCTIONSCHEMENAMELENGTH);
	}
	else{
		params->reproductionScheme = reproductionScheme;
		strncpy(params->reproductionSchemeName, reproductionSchemeName, REPRODUCTIONSCHEMENAMELENGTH);
	}
}

/*
	Sets the target fitness
*/
DLL_EXPORT void setTargetFitness(struct parameters *params, float targetFitness){
	params->targetFitness = targetFitness;
}


/*
	sets the mutation type in params
*/
DLL_EXPORT void setMutationType(struct parameters *params, char *mutationType){

	if(strncmp(mutationType, "probabilistic", MUTATIONTYPENAMELENGTH) == 0){

		params->mutationType = probabilisticMutation;
		strncpy(params->mutationTypeName, "probabilistic", MUTATIONTYPENAMELENGTH);
	}

	else if(strncmp(mutationType, "point", MUTATIONTYPENAMELENGTH) == 0){

		params->mutationType = pointMutation;
		strncpy(params->mutationTypeName, "point", MUTATIONTYPENAMELENGTH);
	}
	
	else if(strncmp(mutationType, "onlyActive", MUTATIONTYPENAMELENGTH) == 0){

		params->mutationType = probabilisticMutationOnlyActive;
		strncpy(params->mutationTypeName, "onlyActive", MUTATIONTYPENAMELENGTH);
	}
	
	else{
		printf("\nWarning: mutation type '%s' is invalid. The mutation type must be 'probabilistic' or 'point'. The mutation type has been left unchanged as '%s'.\n", mutationType, params->mutationTypeName);
	}
}


/*
	Sets the update frequency in generations 
*/
DLL_EXPORT void setUpdateFrequency(struct parameters *params, int updateFrequency){
	
	if(updateFrequency < 0){
		printf("Warning: update frequency of %d is invalid. Update frequency must be >= 0. Update frequency is left unchanged as %d.\n", updateFrequency, params->updateFrequency);
	}
	else{
		params->updateFrequency = updateFrequency;
	}
}


/*
	set the active nodes in the given chromosome
*/
static void setChromosomeActiveNodes(struct chromosome *chromo){

	int i;

	/* set the number of active nodes to zero */
	chromo->numActiveNodes = 0;

	/* reset the active nodes */
	for(i = 0; i < chromo->numNodes; i++){
		chromo->nodes[i]->active = 0;
	}

	/* start the recursive search for active nodes from the output nodes for the number of output nodes */
	for(i=0; i < chromo->numOutputs; i++){

		/* if the output connects to a chromosome input, skip */
		if(chromo->outputNodes[i] < chromo->numInputs){
			continue;
		}

		/* begin a recursive search for active nodes */
		recursivelySetActiveNodes(chromo, chromo->outputNodes[i]);
	}

	/* place active nodes in order */
	sortIntArray(chromo->activeNodes, chromo->numActiveNodes);
}

/*
	used by setActiveNodes to recursively search for active nodes
*/
static void recursivelySetActiveNodes(struct chromosome *chromo, int nodeIndex){

	int i;

	/* if the given node is an input, stop */
	if(nodeIndex < chromo->numInputs){
		return;
	}

	/* if the given node has already been flagged as active */
	if(chromo->nodes[nodeIndex - chromo->numInputs]->active == 1){
		return;
	}

	/* log the node as active */
	chromo->nodes[nodeIndex - chromo->numInputs]->active = 1;
	chromo->activeNodes[chromo->numActiveNodes] = nodeIndex - chromo->numInputs;
	chromo->numActiveNodes++;

	/* recursively log all the nodes to which the current nodes connect as active */
	for(i=0; i < chromo->arity; i++){
		recursivelySetActiveNodes(chromo, chromo->nodes[nodeIndex - chromo->numInputs]->inputs[i]);
	}
}

/*
	Switches the first chromosome with the last and then sorts the population.
*/
static void sortChromosomeArray(struct chromosome **chromoArray, int numChromos){

	qsort(chromoArray, numChromos, sizeof(struct chromosome *), cmpChromosome);
}

/*
	used by qsort in sortChromosomeArray
*/
static int cmpChromosome(const void *a, const void *b){
   
   const struct chromosome *chromoA = (* (struct chromosome **) a);
   const struct chromosome *chromoB = (* (struct chromosome **) b);
        
   if(chromoA->fitness < chromoB->fitness){
	   return -1;
   }    
   else if(chromoA->fitness == chromoB->fitness){
	   return 0;
   }
   else{
	   return 1;
   }  
}






















/*
	chromosome function definitions
*/
















/*
	save the given chromosome to a graphviz .dot file
	(www.graphviz.org/‎)
*/
DLL_EXPORT void saveChromosomeDot(struct chromosome *chromo, int weights, char *fileName){
	
	int i,j;
	FILE *fp;
	
	char colour[20];
	char weight[20];
	
	fp = fopen(fileName, "w");
	
	if(fp == NULL){
		return;
	}
	
	/* */
	fprintf(fp, "digraph NeuralNetwork {\n");
	
	/* landscape, square and centre */
	fprintf(fp, "rankdir=LR;\n");
	fprintf(fp, "size=\"4,3\";\n");	
	fprintf(fp, "center = true;\n");	
	
	/* for all the inputs */
	for(i=0; i < getNumChromosomeInputs(chromo); i++){
		
		fprintf(fp, "node%d [label=\"Input %d\", color=black, labelfontcolor=black, fontcolor=black];\n", i, i);
	}
	
	/* for all nodes */
	for(i=0; i<getNumChromosomeNodes(chromo); i++){
		
		if(chromo->nodes[i]->active == 1){
			strncpy(colour, "black", 20);
		}
		else{
			strncpy(colour, "lightgrey", 20);
		}
		
		fprintf(fp, "node%d [label=\"%s\", color=%s, labelfontcolor=%s, fontcolor=%s];\n", i+getNumChromosomeInputs(chromo), chromo->funcSet->functionNames[chromo->nodes[i]->function], colour, colour, colour);
		
		/* for each node input */
		for(j=0; j<getChromosomeNodeArity(chromo); j++){
			
			if(weights == 1){
				snprintf(weight, 20, "%f", chromo->nodes[i]->weights[j]);
			}
			else{
				snprintf(weight, 20, " ");
			}
			
			
			fprintf(fp, "node%d -> node%d [label=\"%s\", labelfontcolor=%s, fontcolor=%s, bold=true, color=%s];\n", chromo->nodes[i]->inputs[j], i+getNumChromosomeInputs(chromo), weight, colour, colour, colour);	
		}
	}
	
	for(i=0; i<getNumChromosomeOutputs(chromo); i++){
				
		fprintf(fp, "node%d [label=\"Output %d\", color=black, labelfontcolor=black, fontcolor=black];\n", i+getNumChromosomeInputs(chromo) + getNumChromosomeNodes(chromo), i);
		
		fprintf(fp, "node%d -> node%d [labelfontcolor=black, fontcolor=black, bold=true, color=black];\n", chromo->outputNodes[i], i + getNumChromosomeInputs(chromo) + getNumChromosomeNodes(chromo));	
	} 
	
	
	/* place inputs  on same line */
	fprintf(fp, "{ rank = source;");
	
	for(i=0; i < getNumChromosomeInputs(chromo); i++){
		fprintf(fp, " \"node%d\";", i);
	}
	fprintf(fp, " }\n");
	
		
	/* place outputs  on same line */
	fprintf(fp, "{ rank = max;");
	
	for(i = 0; i < getNumChromosomeOutputs(chromo); i++){
		fprintf(fp, "\"node%d\";", i + getNumChromosomeInputs(chromo) + getNumChromosomeNodes(chromo));
	}
	fprintf(fp, " }\n");
	
	
	
	/* last line of dot file */
	fprintf(fp, "}");
	
	fclose(fp);
}


/*
	Sets the random number seed
*/
DLL_EXPORT void setRandomNumberSeed(unsigned int seed){
	srand(seed);
}

/*
	Gets the number of chromosome inputs
*/
DLL_EXPORT int getNumChromosomeInputs(struct chromosome *chromo){
	return chromo->numInputs;
}

/*
	Gets the number of chromosome nodes
*/
DLL_EXPORT int getNumChromosomeNodes(struct chromosome *chromo){
	return chromo->numNodes;
}

DLL_EXPORT int getNumChromosomeOutputs(struct chromosome *chromo){
	return chromo->numOutputs;
}

DLL_EXPORT int getNumChromosomeActiveNodes(struct chromosome *chromo){
	return chromo->numActiveNodes;
}

DLL_EXPORT int getChromosomeNodeArity(struct chromosome *chromo){
	return chromo->arity;
}


DLL_EXPORT int getNumDataSetInputs(struct dataSet *data){
	return data->numInputs;
}

DLL_EXPORT int getNumDataSetOutputs(struct dataSet *data){
	return data->numOutputs;
}

DLL_EXPORT int getNumDataSetSamples(struct dataSet *data){
	return data->numSamples;
}

DLL_EXPORT float *getDataSetSampleInputs(struct dataSet *data, int sample){
	return data->inputData[sample];
}

DLL_EXPORT float getDataSetSampleInput(struct dataSet *data, int sample, int input){
	return data->inputData[sample][input];
}

DLL_EXPORT float *getDataSetSampleOutputs(struct dataSet *data, int sample){
	return data->outputData[sample];
}

DLL_EXPORT float getDataSetSampleOutput(struct dataSet *data, int sample, int output){
	return data->outputData[sample][output];
}


DLL_EXPORT float getChromosomeOutput(struct chromosome *chromo, int output){
	return chromo->outputValues[output];
}

/*
	Saves the given chromosome in a form which can be read in later
*/
DLL_EXPORT void saveChromosome(struct chromosome *chromo, char *fileName){

	int i,j;
	FILE *fp;

	/* create the chromosome file */
	fp = fopen(fileName, "w");

	/* ensure that the file was created correctly */
	if(fp == NULL){
		printf("Warning: cannot save chromosome to '%s'. Chromosome was not saved.\n", fileName);
		return;
	}

	/* save meta information */
	fprintf(fp, "numInputs,%d\n", chromo->numInputs);
	fprintf(fp, "numNodes,%d\n", chromo->numNodes);
	fprintf(fp, "numOutputs,%d\n", chromo->numOutputs);
	fprintf(fp, "arity,%d\n", chromo->arity);

	fprintf(fp, "fuctionSet");

	for(i=0; i<chromo->funcSet->numFunctions; i++){
		fprintf(fp, ",%s", chromo->funcSet->functionNames[i]);
	}
	fprintf(fp, "\n");

	/* save the chromosome structure */
	for(i=0; i<chromo->numNodes; i++){

		fprintf(fp, "%d\n", chromo->nodes[i]->function);

		for(j=0; j<chromo->arity; j++){
			fprintf(fp, "%d,%f\n", chromo->nodes[i]->inputs[j], chromo->nodes[i]->weights[j]);
		}
	}

	for(i=0; i<chromo->numOutputs; i++){
		fprintf(fp, "%d,", chromo->outputNodes[i]);
	}

	fclose(fp);
}


/*
	Returns a pointer to an initialised chromosome with values obeying the given parameters.
*/
DLL_EXPORT struct chromosome *initialiseChromosome(struct parameters *params){

	struct chromosome *chromo;
	int i;

	/* check that funcSet contains functions*/
	if(params->funcSet->numFunctions < 1){
		printf("Error: chromosome not initialised due to empty functionSet.\nTerminating CGP-Library.\n");
		exit(0);
	}

	/* allocate memory for chromosome */
	chromo = malloc(sizeof(struct chromosome));

	/* allocate memory for nodes */
	chromo->nodes = malloc(params->numNodes * sizeof(struct node));

	/* allocate memory for outputNodes matrix */
	chromo->outputNodes = malloc(params->numOutputs * sizeof(int));

	/* allocate memory for active nodes matrix */
	chromo->activeNodes = malloc(params->numNodes * sizeof(int));

	/* allocate memory for chromosome outputValues */
	chromo->outputValues = malloc(params->numOutputs * sizeof(float));

	/* Initialise each of the chromosomes nodes */
	for(i=0; i<params->numNodes; i++){
		chromo->nodes[i] = initialiseNode(params->numInputs, params->arity, params->funcSet->numFunctions, params->connectionWeightRange, i);
	}

	/* set each of the chromosomes outputs */
	for(i=0; i<params->numOutputs; i++){
		chromo->outputNodes[i] = getRandomChromosomeOutput(params->numInputs, params->numNodes);
	}

	/* Add all nodes to the active node matrix */
	for(i=0; i<params->numNodes; i++){
		chromo->activeNodes[i] = i;
	}

	/* set the number of inputs, nodes and outputs */
	chromo->numInputs = params->numInputs;
	chromo->numNodes = params->numNodes;
	chromo->numOutputs = params->numOutputs;
	chromo->arity = params->arity;

	/* set the number of active node to the number of nodes (all active) */
	chromo->numActiveNodes = params->numNodes;

	/* */
	chromo->fitness = -1;

	/*  */
	chromo->funcSet = malloc(sizeof(struct functionSet));

	/* */
	copyFuctionSet(chromo->funcSet, params->funcSet);

	/* set the active nodes in the newly generated chromosome */
	setChromosomeActiveNodes(chromo);


	chromo->nodeInputsHold = malloc(params->arity * sizeof(float));

	return chromo;
}

/*
	Reads in saved chromosomes 
*/
DLL_EXPORT struct chromosome* initialiseChromosomeFromFile(char *file){

	FILE *fp;
	struct chromosome *chromo;
	struct parameters *params;

	char *line, *record;
	char funcName[FUNCTIONNAMELENGTH];
	char buffer[1024];

	int i,j;
	int numInputs, numNodes, numOutputs, arity;


	/* create the chromosome file */
	fp = fopen(file, "r");

	/* ensure that the file was created correctly */
	if(fp == NULL){
		printf("Warning: cannot load chromosome: '%s'. Chromosome was not loaded.\n", file);
		return NULL;
	}

	/* get num inputs */
	line = fgets(buffer, sizeof(buffer), fp);
	if(line == NULL){/*error*/}
	record = strtok(line,",");
	record = strtok(NULL,",");
	numInputs = atoi(record);

	/* get num nodes */
	line = fgets(buffer, sizeof(buffer), fp);
	if(line == NULL){/*error*/}
	record = strtok(line,",");
	record = strtok(NULL,",");
	numNodes = atoi(record);

	/* get num outputs */
	line = fgets(buffer, sizeof(buffer), fp);
	if(line == NULL){/*error*/}
	record = strtok(line,",");
	record = strtok(NULL,",");
	numOutputs = atoi(record);

	/* get arity */
	line = fgets(buffer, sizeof(buffer), fp);
	if(line == NULL){/*error*/}
	record = strtok(line,",");
	record = strtok(NULL,",");
	arity = atoi(record);

	/* initialise parameters  */
	params = initialiseParameters(numInputs, numNodes, numOutputs, arity);

	/* get and set node functions */
	line = fgets(buffer, sizeof(buffer), fp);
	if(line == NULL){/*error*/}
	record = strtok(line,",\n");
	record = strtok(NULL,",\n");

	/* for each function name */
	while( record != NULL){

		strncpy(funcName, record, FUNCTIONNAMELENGTH);
		
		/* can only load functions defined within CGP-Library */
		if(addPresetFuctionToFunctionSet(params,funcName) == 0){
			printf("Error: cannot load chromosome which contains custom node functions.\n");
			printf("Terminating CGP-Library.\n");
			freeParameters(params);
			exit(0);
		}

		record = strtok(NULL,",\n");
	}

	/* initialise chromosome */
	chromo = initialiseChromosome(params);


	/* set the node parameters */
	for(i=0; i<numNodes; i++){

		/* get the function gene */
		line=fgets(buffer, sizeof(buffer), fp);
		record = strtok(line,",\n");
		chromo->nodes[i]->function = atoi(record);

		for(j=0; j<arity; j++){
			line=fgets(buffer, sizeof(buffer), fp);
			sscanf(line,"%d,%f", &chromo->nodes[i]->inputs[j], &chromo->nodes[i]->weights[j]);
		}
	}

	/* set the outputs */
	line=fgets(buffer, sizeof(buffer), fp);
	record = strtok(line,",\n");
	chromo->outputNodes[0] = atoi(record);

	for(i=1; i<numOutputs; i++){
		record = strtok(NULL,",\n");
		chromo->outputNodes[i] = atoi(record);
	}

	fclose(fp);

	freeParameters(params);

	return chromo;
}

/*
	Returns a pointer to an initialised chromosome with values obeying the given parameters.
*/
DLL_EXPORT struct chromosome *initialiseChromosomeFromChromosome(struct chromosome *chromo){

	struct chromosome *chromoNew;
	int i;

	/* check that funcSet contains functions*/
	if(chromo == NULL){
		printf("Error: cannot initialise chromosome from uninitialised chromosome.\nTerminating CGP-Library.\n");
		exit(0);
	}

	/* allocate memory for chromosome */
	chromoNew = malloc(sizeof(struct chromosome));

	/* allocate memory for nodes */
	chromoNew->nodes = malloc(chromo->numNodes * sizeof(struct node));

	/* allocate memory for outputNodes matrix */
	chromoNew->outputNodes = malloc(chromo->numOutputs * sizeof(int));

	/* allocate memory for active nodes matrix */
	chromoNew->activeNodes = malloc(chromo->numNodes * sizeof(int));

	/* allocate memory for chromosome outputValues */
	chromoNew->outputValues = malloc(chromo->numOutputs * sizeof(float));

	/* Initialise each of the chromosomes nodes */
	for(i=0; i<chromo->numNodes; i++){
		chromoNew->nodes[i] = initialiseNode(chromo->numInputs, chromo->arity, chromo->funcSet->numFunctions, 0, i);
		copyNode(chromoNew->nodes[i], chromo->nodes[i]);
	}

	/* set each of the chromosomes outputs */
	for(i=0; i<chromo->numOutputs; i++){
		chromoNew->outputNodes[i] = chromo->outputNodes[i];
	}

	/* set the number of inputs, nodes and outputs */
	chromoNew->numInputs = chromo->numInputs;
	chromoNew->numNodes = chromo->numNodes;
	chromoNew->numOutputs = chromo->numOutputs;
	chromoNew->arity = chromo->arity;

	
	/* */
	chromoNew->fitness = chromo->fitness;

	/* */
	chromoNew->generation = chromo->generation;

	/*  */
	chromoNew->funcSet = malloc(sizeof(struct functionSet));

	/* */
	copyFuctionSet(chromoNew->funcSet, chromo->funcSet);

	/* set the active nodes in the newly generated chromosome */
	setChromosomeActiveNodes(chromoNew);

	chromoNew->nodeInputsHold = malloc(chromo->arity * sizeof(float));

	return chromoNew;
}




/*
	initialises a results structure
*/
struct results* initialiseResults(struct parameters *params, int numRuns){

	int i;
	struct results *rels;

	rels = malloc(sizeof(struct results));
	rels->bestChromosomes = malloc(numRuns * sizeof(struct chromosome));


	rels->numRuns = numRuns;

	for(i=0; i<numRuns; i++){
		/* rels->bestChromosomes[i] = initialiseChromosome(params); */
	}

	return rels;
}

/*
	free a initialised results structure 
*/
DLL_EXPORT void freeResults(struct results *rels){

	int i;

	/* attempt to prevent user double freeing */
	if(rels == NULL){
		return;
	}

	for(i=0; i<rels->numRuns; i++){
		freeChromosome(rels->bestChromosomes[i]);
	}

	free(rels->bestChromosomes);
	free(rels);
}

/*
	saves results structure to file
*/
DLL_EXPORT void saveResults(struct results *rels, char *fileName){
	
	FILE *fp;
	int i;
	
	struct chromosome *chromoTemp;
	
	if(rels == NULL){
		printf("Warning: cannot save uninitialised results structure. Results not saved.\n");
		return;
	}
	
	fp = fopen(fileName, "w");
	
	if(fp == NULL){
		printf("Warning: cannot open '%s' and so cannot save results to that file. Results not saved.\n", fileName);
		return;
	}
	
	fprintf(fp, "Run,Fitness,Generations,Active Nodes\n");
		
	for(i=0; i<rels->numRuns; i++){
		
		chromoTemp = getChromosome(rels, i);
		
		fprintf(fp, "%d,%f,%d,%d\n", i, chromoTemp->fitness, chromoTemp->generation, chromoTemp->numActiveNodes);
			
		freeChromosome(chromoTemp);
	}
	
	fclose(fp);
}



DLL_EXPORT int getNumChromosomes(struct results *rels){
	return rels->numRuns;
}


/*
	returns the average number of chromosome active nodes from repeated
	run results specified in rels.
*/
DLL_EXPORT float getAverageActiveNodes(struct results *rels){

	int i;
	float avgActiveNodes = 0;
	struct chromosome *chromoTemp;


	for(i=0; i<getNumChromosomes(rels); i++){

		chromoTemp = rels->bestChromosomes[i];

		avgActiveNodes += getNumChromosomeActiveNodes(chromoTemp);
	}

	avgActiveNodes = avgActiveNodes / getNumChromosomes(rels);

	return avgActiveNodes;
}

/*
	returns the average chromosome fitness from repeated
	run results specified in rels.
*/
DLL_EXPORT float getAverageFitness(struct results *rels){

	int i;
	float avgFit = 0;
	struct chromosome *chromoTemp;


	for(i=0; i<getNumChromosomes(rels); i++){

		chromoTemp = rels->bestChromosomes[i];

		avgFit += getChromosomeFitness(chromoTemp);
	}

	avgFit = avgFit / getNumChromosomes(rels);

	return avgFit;
}

/*
	returns the number of generations used by each run  specified in rels.
*/
DLL_EXPORT float getAverageGenerations(struct results *rels){

	int i;
	float avgGens = 0;
	struct chromosome *chromoTemp;


	for(i=0; i<getNumChromosomes(rels); i++){

		chromoTemp = rels->bestChromosomes[i];

		avgGens += getChromosomeGenerations(chromoTemp);
	}

	avgGens = avgGens / getNumChromosomes(rels);

	return avgGens;
}

/*
	returns the number of generations required to find the given chromosome 
*/
DLL_EXPORT int getChromosomeGenerations(struct chromosome *chromo){
	return chromo->generation;
}

/*
	returns a pointer to a copy of the best chromosomes found on the given run in rels.
*/
DLL_EXPORT struct chromosome* getChromosome(struct results *rels, int run){

	struct chromosome *chromo;

	/* do some error checking */
	if(rels == NULL){
		printf("Error: cannot get best chromosome from uninitialised results.\nTerminating CGP-Library.\n");
		exit(0);
	}

	chromo = initialiseChromosomeFromChromosome(rels->bestChromosomes[run]);

	return chromo;
}


DLL_EXPORT struct results* repeatCGP(struct parameters *params, struct dataSet *data, int numGens, int numRuns){

	int i;
	struct results *rels;
	int updateFrequency = params->updateFrequency;

	/* set the update frequency so as to to so generational results */
	params->updateFrequency = 0;

	rels = initialiseResults(params, numRuns);

	printf("Run\tFitness\t\tGenerations\tActive Nodes\n");

	/* for each run */
	for(i=0; i<numRuns; i++){

		/* run cgp */
		rels->bestChromosomes[i] = runCGP(params, data, numGens);

		printf("%d\t%f\t%d\t\t%d\n", i, rels->bestChromosomes[i]->fitness, rels->bestChromosomes[i]->generation, rels->bestChromosomes[i]->numActiveNodes);
	}

	printf("----------------------------------------------------\n");
	printf("AVG\t%f\t%f\t%f\n", getAverageFitness(rels), getAverageGenerations(rels), getAverageActiveNodes(rels));
	printf("----------------------------------------------------\n\n");

	/* restore the original value for the update frequency */
	params->updateFrequency = updateFrequency;

	return rels;
}


DLL_EXPORT struct chromosome* runCGP(struct parameters *params, struct dataSet *data, int numGens){

	int i;
	int gen;

	/* */
	struct chromosome *chromo;
	struct chromosome *bestChromo;
	
	struct chromosome **parentChromos;
	struct chromosome **childrenChromos;

	/* storage for chromosomes used by selection scheme */
	struct chromosome **candidateChromos;  
	int numCandidateChromos;

	/* error checking */
	if(numGens < 0){
		printf("Error: %d generations is invalid. The number of generations must be >= 0.\n Terminating CGP-Library.\n", numGens);
		exit(0);
	}

	if(data != NULL && params->numInputs != data->numInputs){
		printf("Error: The number of inputs specified in the data (%d) does not match the number of inputs specified in the parameters (%d).\n", data->numInputs, params->numInputs);
		printf("Terminating CGP-Library.\n");
		exit(0);	
	}

	if(data != NULL && params->numOutputs != data->numOutputs){
		printf("Warning: The number of outputs specified in the data (%d) does not match the number of outputs specified in the parameters (%d).\n", data->numOutputs, params->numOutputs);
		printf("Terminating CGP-Library.\n");
		exit(0);
	}

	/* initialise parent chromosomes */
	parentChromos = malloc(params->mu * sizeof(struct chromosome *));

	for(i=0; i<params->mu; i++){
		parentChromos[i] = initialiseChromosome(params);
	}

	/* initialise children chromosomes */
	childrenChromos = malloc(params->lambda * sizeof(struct chromosome *));

	for(i=0; i<params->lambda; i++){
		childrenChromos[i] = initialiseChromosome(params);
	}

	/* initialise the chromosome to be returned */
	chromo = initialiseChromosome(params);

	/* determine the size of the Candidate Chromos based on the evolutionary Strategy */
	if(params->evolutionaryStrategy == '+'){
		numCandidateChromos = params->mu + params->lambda;
	}
	else if(params->evolutionaryStrategy == ','){
		numCandidateChromos = params->lambda;
	}
	else{
		printf("Error: the evolutionary strategy '%c' is not known.\nTerminating CGP-Library.\n", params->evolutionaryStrategy);
		exit(0);
	}

	/* initialise the candidateChromos */
	candidateChromos = malloc(numCandidateChromos * sizeof(struct chromosome *));

	for(i=0; i<numCandidateChromos; i++){
		candidateChromos[i] = initialiseChromosome(params);
	}

	/* set fitness of the parents */
	for(i=0; i<params->mu; i++){
		setChromosomeFitness(params, parentChromos[i], data);
	}
	
	/* show the user whats going on */
	if(params->updateFrequency != 0){
		printf("\n-- Starting CGP --\n\n");
		printf("Gen\tfitness\n");
	}

	/* for each generation */
	for(gen=0; gen <numGens; gen++){

		/* set fitness of the children of the population */
		for(i=0; i< params->lambda; i++){
			setChromosomeFitness(params, childrenChromos[i], data);
		}

		/* get best chromosome */
		bestChromo = parentChromos[0];

		for(i=1; i<params->mu; i++){

			if(parentChromos[i]->fitness <= bestChromo->fitness){
				bestChromo = parentChromos[i];
			}	
		}

		for(i=0; i<params->lambda; i++){

			if(childrenChromos[i]->fitness <= bestChromo->fitness){
				bestChromo = childrenChromos[i];
			}	
		}

		
		/* check termination conditions */
		if(getChromosomeFitness(bestChromo) <= params->targetFitness){

			if(params->updateFrequency != 0){
				printf("%d\t%f - Solution Found\n", gen, bestChromo->fitness);
			}

			break;
		}


		/* display progress to the user at the update frequency specified */
		if(params->updateFrequency != 0 && (gen % params->updateFrequency == 0 || gen >= numGens-1) ){
			printf("%d\t%f\n", gen, bestChromo->fitness);
		}
		

		/*
			Set the chromosomes which will be used by the selection scheme
			dependant upon the evolutionary strategy. i.e. '+' all are used
			by selection scheme, ',' only the children are.
			
			Note: the children are always placed before the parents when both are
			added to the candidate chromosomes. This ensures neutral genetic drift can 
			take place.
		*/
		for(i=0; i<numCandidateChromos; i++){

			if(i < params->lambda){
				copyChromosome(candidateChromos[i], childrenChromos[i] );
			}
			else{
				copyChromosome(candidateChromos[i], parentChromos[i - params->lambda] );
			}
		}

		/* select the parents from the candidateChromos */
		params->selectionScheme(params, parentChromos, candidateChromos, params->mu, numCandidateChromos);

		/* create the children from the parents */
		params->reproductionScheme(params, parentChromos, childrenChromos, params->mu, params->lambda);
	}

	/* deal with formatting for displaying progress */
    if(params->updateFrequency != 0){
        printf("\n");
    }

	/* get the fittest chromosome form the evolutionary run */
	bestChromo = parentChromos[0];

	for(i=1; i<params->mu; i++){

		if(parentChromos[i]->fitness <= bestChromo->fitness){
			bestChromo = parentChromos[i];
		}	
	}

	for(i=0; i<params->lambda; i++){

		if(childrenChromos[i]->fitness <= bestChromo->fitness){
			bestChromo = childrenChromos[i];
		}	
	}

	bestChromo->generation = gen;
	copyChromosome(chromo, bestChromo);

	

	for(i=0; i<params->mu; i++){
		freeChromosome(parentChromos[i]);
	}
	free(parentChromos);


	for(i=0; i<params->lambda; i++){
		freeChromosome(childrenChromos[i]);
	}
	free(childrenChromos);


	/* free the used chromosomes and population */
	for(i=0; i<numCandidateChromos; i++){
		freeChromosome(candidateChromos[i]);
	}
	free(candidateChromos);

	return chromo;
}


/* copies the contents of funcSetSrc to funcSetDest */
static void copyFuctionSet(struct functionSet *funcSetDest, struct functionSet *funcSetSrc){

	int i;

	funcSetDest->numFunctions = funcSetSrc->numFunctions;


	for(i=0; i<funcSetDest->numFunctions; i++){

		strncpy(funcSetDest->functionNames[i], funcSetSrc->functionNames[i], FUNCTIONNAMELENGTH);
		
		funcSetDest->functions[i] = funcSetSrc->functions[i];
	}
}





/*
	return the fitness of the given chromosome
*/
DLL_EXPORT float getChromosomeFitness(struct chromosome *chromo){
	return chromo->fitness;
}





/*
	Initialises data structure and assigns values of given file
*/
DLL_EXPORT struct dataSet *initialiseDataSetFromFile(char *file){

	int i;
	struct dataSet *data;
	FILE *fp;
	char *line, *record;
	char buffer[1024];
	int lineNum = -1;
	int col;

	/* attempt to open the given file */
	fp = fopen(file, "r");

	/* if the file cannot be found */
	if(fp == NULL){
		printf("Error: file '%s' cannot be found.\nTerminating CGP-Library.\n", file);
		exit(0);
	}

	/* initialise memory for data structure */
	data = malloc(sizeof(struct dataSet));

	/* for every line in the given file */
	while( (line=fgets(buffer, sizeof(buffer), fp)) != NULL){

		/* deal with the first line containing meta data */
		if(lineNum == -1){

			sscanf(line, "%d,%d,%d", &(data->numInputs), &(data->numOutputs), &(data->numSamples));

			data->inputData = malloc(data->numSamples * sizeof(float*));
			data->outputData = malloc(data->numSamples * sizeof(float*));

			for(i=0; i<data->numSamples; i++){
				data->inputData[i] = malloc(data->numInputs * sizeof(float));
				data->outputData[i] = malloc(data->numOutputs * sizeof(float));
			}
		}
		/* the other lines contain input output pairs */
		else{

			/* get the first value on the given line */
			record = strtok(line," ,\n");
			col = 0;

			/* until end of line */
			while(record != NULL){

				/* if its an input value */
				if(col < data->numInputs){
					data->inputData[lineNum][col] = atof(record);
				}

				/* if its an output value */
				else{

					data->outputData[lineNum][col - data->numInputs] = atof(record);
				}

				/* get the next value on the given line */
				record = strtok(NULL," ,\n");

				/* increment the current col index */
				col++;
			}
		}

		/* increment the current line index */
		lineNum++;
	}

	fclose(fp);

	return data;
}


/*
	Initialises data structure and assigns values of given arrays
	arrays must take the form
	inputs[numSamples][numInputs]
	outputs[numSamples][numOutputs]
*/
DLL_EXPORT struct dataSet *initialiseDataSetFromArrays(int numInputs, int numOutputs, int numSamples, float *inputs, float *outputs){

	int i,j;
	struct dataSet *data;

	/* initialise memory for data structure */
	data = malloc(sizeof(struct dataSet));

	data->numInputs = numInputs;
	data->numOutputs = numOutputs;
	data->numSamples = numSamples;

	data->inputData = malloc(data->numSamples * sizeof(float**));
	data->outputData = malloc(data->numSamples * sizeof(float**));

	for(i=0; i<data->numSamples; i++){

		data->inputData[i] = malloc(data->numInputs * sizeof(float));
		data->outputData[i] = malloc(data->numOutputs * sizeof(float));

		for(j=0; j<data->numInputs; j++){
			data->inputData[i][j] = inputs[(i*data->numInputs) + j];
		}

		for(j=0; j<data->numOutputs; j++){
			data->outputData[i][j] = outputs[(i*data->numOutputs) + j];
		}
	}

	return data;
}


/*
	frees given dataSet
*/
DLL_EXPORT void freeDataSet(struct dataSet *data){

	int i;

	/* attempt to prevent user double freeing */
	if(data == NULL){
		return;
	}

	for(i=0; i<data->numSamples; i++){
		free(data->inputData[i]);
		free(data->outputData[i]);
	}

	free(data->inputData);
	free(data->outputData);
	free(data);
}

/*
	saves dataset to file
*/
DLL_EXPORT void saveDataSet(struct dataSet *data, char *fileName){

	int i,j;
	FILE *fp;

	fp = fopen(fileName, "w");

	if(fp == NULL){
		printf("Warning: cannot save data set to %s. Data set was not saved.\n", fileName);
		return;
	}

	fprintf(fp, "%d,", data->numInputs);
	fprintf(fp, "%d,", data->numOutputs);
	fprintf(fp, "%d,", data->numSamples);
	fprintf(fp, "\n");


	for(i=0; i<data->numSamples; i++){

		for(j=0; j<data->numInputs; j++){
			fprintf(fp, "%f,", data->inputData[i][j]);
		}

		for(j=0; j<data->numOutputs; j++){
			fprintf(fp, "%f,", data->outputData[i][j]);
		}

		fprintf(fp, "\n");
	}

	fclose(fp);
}



/*
	prints the given data structure to the screen
*/
DLL_EXPORT void printDataSet(struct dataSet *data){

	int i,j;

	printf("DATA SET\n");
	printf("Inputs: %d, ", data->numInputs);
	printf("Outputs: %d, ", data->numOutputs);
	printf("Samples: %d\n", data->numSamples);

	for(i=0; i<data->numSamples; i++){

		for(j=0; j<data->numInputs; j++){
			printf("%f ", data->inputData[i][j]);
		}

		printf(" : ");

		for(j=0; j<data->numOutputs; j++){
			printf("%f ", data->outputData[i][j]);
		}

		printf("\n");
	}
}













/*
	returns mu value currently set in given parameters.
*/
DLL_EXPORT int getMu(struct parameters *params){
	return params->mu;
}










/*
*/
DLL_EXPORT int getNumInputs(struct parameters *params){
	return params->numInputs;
}

/*
*/
DLL_EXPORT int getNumOutputs(struct parameters *params){
	return params->numOutputs;
}




/*
	used as an interface to adding pre-set node functions
*/
static int addPresetFuctionToFunctionSet(struct parameters *params, char *functionName){

	int output = 1;

	if(strncmp(functionName, "add", FUNCTIONNAMELENGTH) == 0){
		addNodeFunctionCustom(params, add, "add");
	}
	else if(strncmp(functionName, "sub", FUNCTIONNAMELENGTH) == 0){
		addNodeFunctionCustom(params, sub, "sub");
	}
	else if(strncmp(functionName, "mul", FUNCTIONNAMELENGTH) == 0){
		addNodeFunctionCustom(params, mul, "mul");
	}
	else if(strncmp(functionName, "div", FUNCTIONNAMELENGTH) == 0){
		addNodeFunctionCustom(params, divide, "div");
	}
	else if(strncmp(functionName, "abs", FUNCTIONNAMELENGTH) == 0){
		addNodeFunctionCustom(params, absolute, "abs");
	}
	else if(strncmp(functionName, "sqrt", FUNCTIONNAMELENGTH) == 0){
		addNodeFunctionCustom(params, squareRoot, "sqrt");
	}
	else if(strncmp(functionName, "sq", FUNCTIONNAMELENGTH) == 0){
		addNodeFunctionCustom(params, square, "sq");
	}
	else if(strncmp(functionName, "cube", FUNCTIONNAMELENGTH) == 0){
		addNodeFunctionCustom(params, cube, "cube");
	}
	else if(strncmp(functionName, "exp", FUNCTIONNAMELENGTH) == 0){
		addNodeFunctionCustom(params, exponential, "exp");
	}
	else if(strncmp(functionName, "sin", FUNCTIONNAMELENGTH) == 0){
		addNodeFunctionCustom(params, sine, "sin");
	}
	else if(strncmp(functionName, "cos", FUNCTIONNAMELENGTH) == 0){
		addNodeFunctionCustom(params, cosine, "cos");
	}
	else if(strncmp(functionName, "tan", FUNCTIONNAMELENGTH) == 0){
		addNodeFunctionCustom(params, tangent, "tan");
	}


	else if(strncmp(functionName, "and", FUNCTIONNAMELENGTH) == 0){
		addNodeFunctionCustom(params, and, "and");
	}
	else if(strncmp(functionName, "nand", FUNCTIONNAMELENGTH) == 0){
		addNodeFunctionCustom(params, nand, "nand");
	}
	else if(strncmp(functionName, "or", FUNCTIONNAMELENGTH) == 0){
		addNodeFunctionCustom(params, or, "or");
	}
	else if(strncmp(functionName, "nor", FUNCTIONNAMELENGTH) == 0){
		addNodeFunctionCustom(params, nor, "nor");
	}
	else if(strncmp(functionName, "xor", FUNCTIONNAMELENGTH) == 0){
		addNodeFunctionCustom(params, xor, "xor");
	}
	else if(strncmp(functionName, "xnor", FUNCTIONNAMELENGTH) == 0){
		addNodeFunctionCustom(params, xnor, "xnor");
	}
	else if(strncmp(functionName, "not", FUNCTIONNAMELENGTH) == 0){
		addNodeFunctionCustom(params, not, "not");
	}

	else if(strncmp(functionName, "sig", FUNCTIONNAMELENGTH) == 0){
		addNodeFunctionCustom(params, sigmoid, "sig");
	}
	else if(strncmp(functionName, "gauss", FUNCTIONNAMELENGTH) == 0){
		addNodeFunctionCustom(params, gaussian, "gauss");
	}
	else if(strncmp(functionName, "step", FUNCTIONNAMELENGTH) == 0){
		addNodeFunctionCustom(params, step, "step");
	}
	else if(strncmp(functionName, "softsign", FUNCTIONNAMELENGTH) == 0){
		addNodeFunctionCustom(params, softsign, "soft");
	}
	else if(strncmp(functionName, "tanh", FUNCTIONNAMELENGTH) == 0){
		addNodeFunctionCustom(params, hyperbolicTangent, "tanh");
	}
	else{
		printf("Warning: function '%s' is not known and was not added.\n", functionName);
		output = 0;
	}

	return output;
}










/*
	Frees the memory associated with the given chromosome structure
*/
DLL_EXPORT void freeChromosome(struct chromosome *chromo){

	int i;

	/* attempt to prevent user double freeing */
	if(chromo == NULL){
		return;
	}

	for(i=0; i < chromo->numNodes; i++){
		freeNode(chromo->nodes[i]);
	}

	free(chromo->nodeInputsHold);
	free(chromo->funcSet);
	free(chromo->outputValues);
	free(chromo->nodes);
	free(chromo->outputNodes);
	free(chromo->activeNodes);
	free(chromo);
}


/*

*/
DLL_EXPORT void copyChromosome(struct chromosome *chromoDest, struct chromosome *chromoSrc){

	int i;

	/* error checking  */
	if(chromoDest->numInputs != chromoSrc->numInputs){
		printf("Error: cannot copy a chromosome to a chromosome of different dimensions. The number of chromosome inputs do not match.\n");
		printf("Terminating CGP-Library.\n");
		exit(0);
	}

	if(chromoDest->numNodes != chromoSrc->numNodes){
		printf("Error: cannot copy a chromosome to a chromosome of different dimensions. The number of chromosome nodes do not match.\n");
		printf("Terminating CGP-Library.\n");
		exit(0);
	}

	if(chromoDest->numOutputs != chromoSrc->numOutputs){
		printf("Error: cannot copy a chromosome to a chromosome of different dimensions. The number of chromosome outputs do not match.\n");
		printf("Terminating CGP-Library.\n");
		exit(0);
	}

	if(chromoDest->arity != chromoSrc->arity){
		printf("Error: cannot copy a chromosome to a chromosome of different dimensions. The arity of the chromosome nodes do not match.\n");
		printf("Terminating CGP-Library.\n");
		exit(0);
	}

	/* copy nodes */
	for(i=0; i<chromoSrc->numNodes; i++){
		copyNode(chromoDest->nodes[i],  chromoSrc->nodes[i]);
	}

	/* copy fuctionset */
	copyFuctionSet(chromoDest->funcSet, chromoSrc->funcSet);

	/* copy each of the chromosomes outputs */
	for(i=0; i<chromoSrc->numOutputs; i++){
		chromoDest->outputNodes[i] = chromoSrc->outputNodes[i];
	}

	/* copy the active node matrix */
	for(i=0; i<chromoSrc->numNodes; i++){
		chromoDest->activeNodes[i] = chromoSrc->activeNodes[i];
	}

	/* copy the number of active node */
	chromoDest->numActiveNodes = chromoSrc->numActiveNodes;

	/* copy the fitness */
	chromoDest->fitness = chromoSrc->fitness;

	/* copy generation */
	chromoDest->generation = chromoSrc->generation;
}


/*
*/
static void copyNode(struct node *nodeDest, struct node *nodeSrc){

	int i;

	/* copy the node's function */
	nodeDest->function = nodeSrc->function;

	/* copy active flag */
	nodeDest->active = nodeSrc->active;

	/* */
	nodeDest->arity = nodeSrc->arity;

	/* copy the nodes inputs and connection weights */
	for(i=0; i<nodeSrc->arity; i++){
		nodeDest->inputs[i] = nodeSrc->inputs[i];
		nodeDest->weights[i] = nodeSrc->weights[i];
	}
}

/*
	sets the fitness of the given chromosome
*/
DLL_EXPORT void setChromosomeFitness(struct parameters *params, struct chromosome *chromo, struct dataSet *data){

	float fitness;

	setChromosomeActiveNodes(chromo);

	fitness = params->fitnessFunction(params, chromo, data);

	chromo->fitness = fitness;
}


/*
	mutate Random parent reproduction method.
*/
static void mutateRandomParent(struct parameters *params, struct chromosome **parents, struct chromosome **children, int numParents, int numChildren){

	int i;

	/* for each child */
	for(i=0; i< numChildren; i++){

		/* set child as clone of random parent */
		copyChromosome(children[i], parents[randInt(numParents)]);

		/* mutate newly cloned child */
		mutateChromosome(params, children[i]);
	}
}


/*
	Selection scheme which selects the fittest members of the population
	to be the parents.
*/
static void selectFittest(struct parameters *params, struct chromosome **parents, struct chromosome **candidateChromos, int numParents, int numCandidateChromos ){

	int i;

	/*
	for(i=0; i<numCandidateChromos; i++){
		printf("%.0f, ", getChromosomeFitness(candidateChromos[i]));
	}

	printf("\n");
*/
	sortChromosomeArray(candidateChromos, numCandidateChromos);

/*

	for(i=0; i<numCandidateChromos; i++){
		printf("%.0f, ", getChromosomeFitness(candidateChromos[i]));
	}

	printf("\n");
	
	getchar();
*/
	for(i=0; i<numParents; i++){

		copyChromosome(parents[i], candidateChromos[i]);
	}
}


/*
	Executes the given chromosome
*/
DLL_EXPORT void executeChromosome(struct chromosome *chromo, const float *inputs){

	int i,j;
	int nodeInputLocation;
	int currentActiveNode;
	int currentActiveNodeFuction;

	const int numInputs = chromo->numInputs;
	const int numActiveNodes = chromo->numActiveNodes;
	const int numOutputs = chromo->numOutputs;
	const int arity = chromo->arity;

	/* error checking */
	if(chromo == NULL){
		printf("Error: cannot execute uninitialised chromosome.\n Terminating CGP-Library.\n");
		exit(0);
	}


	/* for all of the active nodes */
	for(i=0; i<numActiveNodes; i++){

		/* get the index of the current active node */
		currentActiveNode = chromo->activeNodes[i];

		/* for each of the active nodes inputs */
		for(j=0; j<arity; j++){

			/* gather the nodes inputs */
			nodeInputLocation = chromo->nodes[currentActiveNode]->inputs[j];

			if(nodeInputLocation < numInputs){
				chromo->nodeInputsHold[j] = inputs[nodeInputLocation];
			}
			else{
				chromo->nodeInputsHold[j] = chromo->nodes[nodeInputLocation - numInputs]->output;
			}
		}

		/* get the index of the active node under evaluation */
		currentActiveNodeFuction = chromo->nodes[currentActiveNode]->function;

		/* calculate the output of the active node under evaluation */
		chromo->nodes[currentActiveNode]->output = chromo->funcSet->functions[currentActiveNodeFuction](arity, chromo->nodeInputsHold, chromo->nodes[currentActiveNode]->weights);


		/* prevent float form going to inf and -inf */
		if(isinf(chromo->nodes[currentActiveNode]->output) != 0 ){

			if(chromo->nodes[currentActiveNode]->output > 0){
				chromo->nodes[currentActiveNode]->output = FLT_MAX;
			}
			else{
				chromo->nodes[currentActiveNode]->output = FLT_MIN;
			}
		}

		/* deal with floats becoming NAN */
		if(isnan(chromo->nodes[currentActiveNode]->output) != 0){
			chromo->nodes[currentActiveNode]->output = 0;
		}
	}

	/* Set the chromosome outputs */ 
	for(i=0; i<numOutputs; i++){

		if(chromo->outputNodes[i] < numInputs){
			chromo->outputValues[i] = inputs[chromo->outputNodes[i]];
		}
		else{
			chromo->outputValues[i] = chromo->nodes[chromo->outputNodes[i] - numInputs]->output;
		}
	}
	
	/*
	for(i=0; i<numInputs; i++){
		printf("%d, ", (int)inputs[i]);
	}
	
	printf("\n%d\n", (int)chromo->outputValues[0]);
	
	printChromosome(chromo);
	
	
	getchar();
	*/
}


/*
	returns a pointer to an initialised node. Initialised means that necessary
	memory has been allocated and values set.
*/
static struct node *initialiseNode(int numInputs, int arity, int numFunctions, float connectionWeightRange, int nodePosition){

	struct node *n;
	int i;

	/* allocate memory for node */
	n = malloc(sizeof(struct node));

	/* allocate memory for the node's inputs and connection weights */
	n->inputs = malloc(arity * sizeof(int));
	n->weights = malloc(arity * sizeof(float));

	/* set the node's function */
	n->function = getRandomFunction(numFunctions);

	/* set as active by default */
	n->active = 1;

	/* set the nodes inputs and connection weights */
	for(i=0; i<arity; i++){
		n->inputs[i] = getRandomNodeInput(numInputs,nodePosition);
		n->weights[i] = getRandomConnectionWeight(connectionWeightRange);
	}

	/* */
	n->output = 0;

	/* */
	n->arity = arity;

	return n;
}


/*

*/
static void freeNode(struct node *n){

	/* attempt to prevent user double freeing */
	if(n == NULL){
		return;
	}

	free(n->inputs);
	free(n->weights);
	free(n);
}

/*
	returns a random connection weight value
*/
static float getRandomConnectionWeight(float weightRange){
	return (randFloat() * 2 * weightRange) - weightRange;
}

/*
	returns a random function index
*/
static int getRandomFunction(int numFunctions){

	/* check that funcSet contains functions*/
	if(numFunctions <1){
		printf("Error: cannot assign the function gene a value as the Fuction Set is empty.\nTerminating CGP-Library.\n");
		exit(0);
	}

	return randInt(numFunctions);
}

/*
 returns a random input for the given node
*/
static int getRandomNodeInput(int numChromoInputs, int nodePosition){

	int input;

	input = randInt(numChromoInputs + nodePosition);

	return input;
}




/*
	returns a random chromosome output
*/
static int getRandomChromosomeOutput(int numInputs, int numNodes){

	int output;

	output = randInt(numInputs + numNodes);

	return output;
}




/*
	Prints the given chromosome to the screen
*/
DLL_EXPORT void printChromosome(struct chromosome *chromo){

	int i,j;

	/* set the active nodes in the given chromosome */
	setChromosomeActiveNodes(chromo);

	/* for all the chromo inputs*/
	for(i=0; i<chromo->numInputs; i++){
		printf("(%d):\tinput\n", i);
	}

	/* for all the hidden nodes */
	for(i = 0; i < chromo->numNodes; i++){

		/* print the node function */
		printf("(%d):\t%s\t", chromo->numInputs + i, chromo->funcSet->functionNames[chromo->nodes[i]->function]);

		/* for the arity of the node */
		for(j = 0; j < chromo->arity; j++){

			/* print the node input information */
			printf("%d,%+.1f\t", chromo->nodes[i]->inputs[j], chromo->nodes[i]->weights[j]);
		}

		/* Highlight active nodes */
		if(chromo->nodes[i]->active == 1){
			printf("*");
		}

		printf("\n");
	}

	/* for all of the outputs */
	printf("outputs: ");
	for(i = 0; i < chromo->numOutputs; i++){

		/* print the output node locations */
		printf("%d ", chromo->outputNodes[i]);
	}

	printf("\n\n");
}

/*
	Mutates the given chromosome using the mutation method described in parameters
*/
DLL_EXPORT void mutateChromosome(struct parameters *params, struct chromosome *chromo){

	params->mutationType(params, chromo);

	setChromosomeActiveNodes(chromo);
}


/*
	Conductions point mutation on the give chromosome. A predetermined
	number of chromosome genes are randomly selected and changed to
	a random valid allele. The number of mutations is the number of chromosome
	genes multiplied by the mutation rate. Each gene has equal probability
	of being selected.
*/
static void pointMutation(struct parameters *params, struct chromosome *chromo){

	int i;
	int numGenes;
	int numFunctionGenes, numInputGenes, numWeightGenes, numOutputGenes;
	int numGenesToMutate;
	int geneToMutate;
	int nodeIndex;
	int nodeInputIndex;

	/* get the number of each type of gene */
	numFunctionGenes = params->numNodes;
	numInputGenes = params->numNodes * params->arity;
	numWeightGenes = params->numNodes * params->arity;
	numOutputGenes = params->numOutputs;

	/* set the total number of chromosome genes */
	numGenes = numFunctionGenes + numInputGenes + numWeightGenes + numOutputGenes;

	/* calculate the number of genes to mutate */
	numGenesToMutate = (int)roundf(numGenes * params->mutationRate);

	/* for the number of genes to mutate */
	for(i=0; i<numGenesToMutate; i++){

		/* select a random gene */
		geneToMutate = randInt(numGenes);

		/* mutate function gene */
		if(geneToMutate < numFunctionGenes){

			nodeIndex = geneToMutate;

			chromo->nodes[nodeIndex]->function = getRandomFunction(chromo->funcSet->numFunctions);
		}

		/* mutate node input gene */
		else if(geneToMutate < numFunctionGenes + numInputGenes){

			nodeIndex = (int) ((geneToMutate - numFunctionGenes) / chromo->arity);
			nodeInputIndex = (geneToMutate - numFunctionGenes) % chromo->arity;

			chromo->nodes[nodeIndex]->inputs[nodeInputIndex] = getRandomNodeInput(chromo->numInputs, nodeIndex);
		}

		/* mutate connection weight */
		else if(geneToMutate < numFunctionGenes + numInputGenes + numWeightGenes){

			nodeIndex = (int) ((geneToMutate - numFunctionGenes - numInputGenes) / chromo->arity);
			nodeInputIndex = (geneToMutate - numFunctionGenes - numInputGenes) % chromo->arity;

			chromo->nodes[nodeIndex]->weights[nodeInputIndex] = getRandomConnectionWeight(params->connectionWeightRange);
		}

		/* mutate output gene */
		else{
			nodeIndex = geneToMutate - numFunctionGenes - numInputGenes - numWeightGenes;
			chromo->outputNodes[nodeIndex] = getRandomChromosomeOutput(chromo->numInputs, chromo->numNodes);
		}
	}
}


/*
	Conductions probabilistic mutation on the give chromosome. Each chromosome
	gene is changed to a random valid allele with a probability specified in
	parameters.
*/
static void probabilisticMutation(struct parameters *params, struct chromosome *chromo){

	int i,j;

	/* for every nodes in the chromosome */
	for(i=0; i<params->numNodes; i++){

		/* mutate the function gene */
		if(randFloat() <= params->mutationRate){
			chromo->nodes[i]->function = getRandomFunction(chromo->funcSet->numFunctions);
		}

		/* for every input to each chromosome */
		for(j=0; j<params->arity; j++){

			/* mutate the node input */
			if(randFloat() <= params->mutationRate){
				chromo->nodes[i]->inputs[j] = getRandomNodeInput(chromo->numInputs, i);
			}

			/* mutate the node connection weight */
			if(randFloat() <= params->mutationRate){
				chromo->nodes[i]->weights[j] = getRandomConnectionWeight(params->connectionWeightRange);
			}
		}
	}

	/* for every chromosome output */
	for(i=0; i<params->numOutputs; i++){

		/* mutate the chromosome output */
		if(randFloat() <= params->mutationRate){
			chromo->outputNodes[i] = getRandomChromosomeOutput(chromo->numInputs, chromo->numNodes);
		}
	}
}

/*
	Conductions probabilistic mutation on the active nodes in the give 
	chromosome. Each chromosome gene is changed to a random valid allele 
	with a probability specified in parameters.
*/
static void probabilisticMutationOnlyActive(struct parameters *params, struct chromosome *chromo){

	int i,j;

	/* for every nodes in the chromosome */
	for(i=0; i<params->numNodes; i++){

		/* skip inactive nodes */
		if(chromo->nodes[i]->active == 0){
			continue;
		}

		/* mutate the function gene */
		if(randFloat() <= params->mutationRate){
			chromo->nodes[i]->function = getRandomFunction(chromo->funcSet->numFunctions);
		}

		/* for every input to each chromosome */
		for(j=0; j<params->arity; j++){

			/* mutate the node input */
			if(randFloat() <= params->mutationRate){
				chromo->nodes[i]->inputs[j] = getRandomNodeInput(chromo->numInputs, i);
			}

			/* mutate the node connection weight */
			if(randFloat() <= params->mutationRate){
				chromo->nodes[i]->weights[j] = getRandomConnectionWeight(params->connectionWeightRange);
			}
		}
	}

	/* for every chromosome output */
	for(i=0; i<params->numOutputs; i++){

		/* mutate the chromosome output */
		if(randFloat() <= params->mutationRate){
			chromo->outputNodes[i] = getRandomChromosomeOutput(chromo->numInputs, chromo->numNodes);
		}
	}
}



/*
	Node function add. Returns the sum of the inputs.
*/
static float add(const int numInputs, const float *inputs, const float *connectionWeights){

	int i;
	float sum = 0;

	for(i=0; i<numInputs; i++){
		sum += inputs[i];
	}

	return sum;
}

/*
	Node function sub. Returns the first input minus all remaining inputs.
*/
static float sub(const int numInputs, const float *inputs, const float *connectionWeights){

	int i;
	float sum = inputs[0];

	for(i=1; i<numInputs; i++){
		sum -= inputs[i];
	}

	return sum;
}


/*
	Node function mul. Returns the multiplication of all the inputs.
*/
static float mul(const int numInputs, const float *inputs, const float *connectionWeights){

	int i;
	float multiplication = 1;

	for(i=0; i<numInputs; i++){
		multiplication *= inputs[i];
	}

	return multiplication;
}


/*
	Node function div. Returns the first input divided by the second input divided by the third input etc
*/
static float divide(const int numInputs, const float *inputs, const float *connectionWeights){

	int i;
	float divide = inputs[0];

	for(i=1; i<numInputs; i++){

		divide /= inputs[i];
	}

	return divide;
}

/*
	Node function abs. Returns the absolute of the first input
*/
static float absolute(const int numInputs, const float *inputs, const float *connectionWeights){

	return fabs(inputs[0]);
}

/*
	Node function sqrt.  Returns the square root of the first input
*/
static float squareRoot(const int numInputs, const float *inputs, const float *connectionWeights){

	return sqrtf(inputs[0]);
}

/*
	Node function squ.  Returns the square of the first input
*/
static float square(const int numInputs, const float *inputs, const float *connectionWeights){

	return powf(inputs[0],2);
}

/*
	Node function cub.  Returns the cube of the first input
*/
static float cube(const int numInputs, const float *inputs, const float *connectionWeights){

	return powf(inputs[0],3);
}

/*
	Node function exp.  Returns the exponential of the first input
*/
static float exponential(const int numInputs, const float *inputs, const float *connectionWeights){

	return expf(inputs[0]);
}


/*
	Node function sin.  Returns the sine of the first input
*/
static float sine(const int numInputs, const float *inputs, const float *connectionWeights){

	return sinf(inputs[0]);
}

/*
	Node function cos.  Returns the cosine of the first input
*/
static float cosine(const int numInputs, const float *inputs, const float *connectionWeights){

	return sinf(inputs[0]);
}

/*
	Node function tan.  Returns the tangent of the first input
*/
static float tangent(const int numInputs, const float *inputs, const float *connectionWeights){

	return tanf(inputs[0]);
}


/*
	Node function and. logical AND, returns '1' if all inputs are '1'
	else, '0'
*/
static float and(const int numInputs, const float *inputs, const float *connectionWeights){

	int i;
	
	for(i=0; i<numInputs; i++){

		if(inputs[i] == 0){
			return 0;
		}
	}

	return 1;
}


/*
	Node function and. logical NAND, returns '0' if all inputs are '1'
	else, '1'
*/
static float nand(const int numInputs, const float *inputs, const float *connectionWeights){

	int i;
	
	for(i=0; i<numInputs; i++){

		if(inputs[i] == 0){
			return 1;
		}
	}

	return 0;
}


/*
	Node function or. logical OR, returns '0' if all inputs are '0'
	else, '1'
*/
static float or(const int numInputs, const float *inputs, const float *connectionWeights){

	int i;
	
	for(i=0; i<numInputs; i++){

		if(inputs[i] == 1){
			return 1;
		}
	}

	return 0;
}


/*
	Node function nor. logical NOR, returns '1' if all inputs are '0'
	else, '0'
*/
static float nor(const int numInputs, const float *inputs, const float *connectionWeights){

	int i;
	
	for(i=0; i<numInputs; i++){

		if(inputs[i] == 1){
			return 0;
		}
	}

	return 1;
}


/*
	Node function xor. logical XOR, returns '1' iff one of the inputs is '1'
	else, '0'. AKA 'one hot'.
*/
static float xor(const int numInputs, const float *inputs, const float *connectionWeights){

	int i;
	int numOnes = 0;
	int out;

	for(i=0; i<numInputs; i++){

		if(inputs[i] == 1){
			numOnes++;
		}

		if(numOnes > 1){
			break;
		}
	}

	if(numOnes == 1){
		out = 1;
	}
	else{
		out = 0;
	}

	return out;
}

/*
	Node function xnor. logical XNOR, returns '0' iff one of the inputs is '1'
	else, '1'.
*/
static float xnor(const int numInputs, const float *inputs, const float *connectionWeights){

	int i;
	int numOnes = 0;
	int out;

	for(i=0; i<numInputs; i++){

		if(inputs[i] == 1){
			numOnes++;
		}

		if(numOnes > 1){
			break;
		}
	}

	if(numOnes == 1){
		out = 0;
	}
	else{
		out = 1;
	}

	return out;
}

/*
	Node function not. logical NOT, returns '1' if first input is '0', else '1'
*/
static float not(const int numInputs, const float *inputs, const float *connectionWeights){

	float out;

	if(inputs[0] == 0){
		out = 1;
	}
	else{
		out = 0;
	}

	return out;
}


/*
	Node function sigmoid. returns the sigmoid of the sum of weighted inputs.
	The specific sigmoid function used in the logistic function.
	range: [0,1]
*/
static float sigmoid(const int numInputs, const float *inputs, const float *connectionWeights){

	float weightedInputSum;
	float out;

	weightedInputSum = sumWeigtedInputs(numInputs, inputs, connectionWeights);

	out = 1 / (1 + expf(-weightedInputSum));

	return out;
}

/*
	Node function Gaussian. returns the Gaussian of the sum of weighted inputs.
	range: [0,1]
*/
static float gaussian(const int numInputs, const float *inputs, const float *connectionWeights){

	float weightedInputSum;
	float out;

	int centre = 0;
	int width = 1;

	weightedInputSum = sumWeigtedInputs(numInputs, inputs, connectionWeights);

	out = expf(-(powf(weightedInputSum - centre,2))/(2*powf(width,2)));

	return out;
}


/*
	Node function step. returns the step function of the sum of weighted inputs.
	range: [0,1]
*/
static float step(const int numInputs, const float *inputs, const float *connectionWeights){

	float weightedInputSum;
	float out;

	weightedInputSum = sumWeigtedInputs(numInputs, inputs, connectionWeights);

	if(weightedInputSum < 0){
		out = 0;
	}
	else{
		out = 1;
	}

	return out;
}


/*
	Node function step. returns the step function of the sum of weighted inputs.
	range: [-1,1]
*/
static float softsign(const int numInputs, const float *inputs, const float *connectionWeights){

	float weightedInputSum;
	float out;

	weightedInputSum = sumWeigtedInputs(numInputs, inputs, connectionWeights);

	out = weightedInputSum / (1 + fabs(weightedInputSum));

	return out;
}


/*
	Node function tanh. returns the tanh function of the sum of weighted inputs.
	range: [-1,1]
*/
static float hyperbolicTangent(const int numInputs, const float *inputs, const float *connectionWeights){

	float weightedInputSum;
	float out;

	weightedInputSum = sumWeigtedInputs(numInputs, inputs, connectionWeights);

	out = tanhf(weightedInputSum);

	return out;
}


/*
	Returns the sum of the weighted inputs.
*/
static float sumWeigtedInputs(const int numInputs, const float *inputs, const float *connectionWeights){

	int i;
	float weightedSum = 0;

	for(i=0; i<numInputs; i++){
		weightedSum += (inputs[i] * connectionWeights[i]);
	}

	return weightedSum;
}


/*
	removes the inactive nodes from the given chromosome
*/
DLL_EXPORT void removeInactiveNodes(struct chromosome *chromo){

	int i,j,k;

	int originalNumNodes = chromo->numNodes;

	/* set the active nodes */
	setChromosomeActiveNodes(chromo);

	/* for all nodes */
	for(i=0; i<chromo->numNodes-1; i++){

		/* if the node is inactive */
		if(chromo->nodes[i]->active == 0){

			/* set the node to be the next node */
			for(j=i; j<chromo->numNodes-1; j++){
				copyNode(chromo->nodes[j], chromo->nodes[j+1]);
			}

			/* */
			for(j=0; j<chromo->numNodes; j++){
				for(k=0; k<chromo->arity; k++){

					if(chromo->nodes[j]->inputs[k] >= i + chromo->numInputs){
						chromo->nodes[j]->inputs[k]--;
					}
				}
			}
			
			/* for the number of chromosome outputs */
			for(j=0; j<chromo->numOutputs; j++){

				if(chromo->outputNodes[j] >= i + chromo->numInputs){
					chromo->outputNodes[j]--;
				}
			}

			/* de-increment the number of nodes */
			chromo->numNodes--;

			/* made the newly assigned node be evaluated */
			i--;
		}
	}

	for(i=chromo->numNodes; i<originalNumNodes; i++){
		freeNode(chromo->nodes[i]);
	}

	if(chromo->nodes[chromo->numNodes-1]->active == 0){
		freeNode(chromo->nodes[chromo->numNodes-1]);
		chromo->numNodes--;
	}

	/* reallocate the memory associated with the chromosome */
	chromo->nodes = realloc(chromo->nodes, chromo->numNodes * sizeof(struct node));
	chromo->activeNodes = realloc(chromo->activeNodes, chromo->numNodes * sizeof(int));

	/* set the active nodes */
	setChromosomeActiveNodes(chromo);
}


/*
	The default fitness function used by CGP-Library.
	Simply assigns an error of the sum of the absolute differences between the target and actual outputs
*/
static float supervisedLearning(struct parameters *params, struct chromosome *chromo, struct dataSet *data){

	int i,j;
	float error = 0;

	/* error checking */
	if(getNumChromosomeInputs(chromo) != getNumDataSetInputs(data)){
		printf("Error: the number of chromosome inputs must match the number of inputs specified in the dataSet.\n");
		printf("Terminating CGP-Library.\n");
		exit(0);
	}

	if(getNumChromosomeOutputs(chromo) != getNumDataSetOutputs(data)){
		printf("Error: the number of chromosome outputs must match the number of outputs specified in the dataSet.\n");
		printf("Terminating CGP-Library.\n");
		exit(0);
	}

	/* for each sample in data */
	for(i=0 ; i<getNumDataSetSamples(data); i++){

		/* calculate the chromosome outputs for the set of inputs  */
		executeChromosome(chromo, getDataSetSampleInputs(data, i));

		/* for each chromosome output */
		for(j=0; j<getNumChromosomeOutputs(chromo); j++){

			error += fabs(getChromosomeOutput(chromo, j) - getDataSetSampleOutput(data, i, j));
		}
	}

	return error;
}


/*
	returns a random float between [0,1]
*/
static float randFloat(void){
	return (float)rand()/(float)RAND_MAX;
}

/*
	sort int array using qsort
*/
static void sortIntArray(int *array, const int length){

	qsort(array, length, sizeof(int), cmpInt);
}

/*
	used by qsort in sortIntArray 
*/
static int cmpInt(const void * a, const void * b){
   return ( *(int*)a - *(int*)b );
}



/*
	Prints the current functions in the function set to
	the terminal.
*/
static void printFunctionSet(struct parameters *params){

	int i;

	printf("Function Set:");

	for(i=0; i<params->funcSet->numFunctions; i++){
		printf(" %s", params->funcSet->functionNames[i]);
	}

	printf(" (%d)\n", params->funcSet->numFunctions);
}


/*
	random integer between zero and n without modulo bias.
	
	adapted from: http://zuttobenkyou.wordpress.com/2012/10/18/generating-random-numbers-without-modulo-bias/
*/
int randInt(int n){
	
	return rand() % n;
	
	/*
	int x;
	int randLimit;
	int randExcess;
	
	randExcess = (RAND_MAX % n) + 1;
	randLimit = RAND_MAX - randExcess;
	
	do{
		x = rand();
	}
	while (x > randLimit);
		
	return x % n;
	*/
}
