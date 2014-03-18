/*
	This file is part of CGP-Library, Andrew James Turner 2014.

    CGP-Library is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published 
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CGP-Library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with CGP-Library.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h> 
#include <string.h>
#include "cgp.h"

#define FUNCTIONSETSIZE 50
#define FUNCTIONNAMELENGTH 512

struct parameters{
	
	/* can use default values*/
	int mu;
	int lambda;
	char stratergy;
	float mutationRate;
	float connectionsWeightRange;
	unsigned int randSeed;
	int generations;
	
	/* need to be set by user */
	int numInputs;
	int numNodes;
	int numOutputs;
	int arity;
	
	float *nodeInputsHold;
	
	struct fuctionSet *funcSet;
	
	void (*mutationType)(struct parameters *params, struct chromosome *chromo);
	float (*fitnessFuction)(struct parameters *params, struct chromosome *chromo, struct data *dat);		
};

struct population{
	struct chromosome **chromosomes;
};

struct fuctionSet{
	int numFunctions;
	char functionNames[FUNCTIONSETSIZE][FUNCTIONNAMELENGTH];
	float (*functions[FUNCTIONSETSIZE])(const int numInputs, float *inputs, float *weights);	
};

struct chromosome{
	
	int numInputs;
	int numOutputs;
	
	struct node **nodes;
	int *outputNodes;
	int numActiveNodes;
	int *activeNodes;
	float fitness;
	
	float *outputValues;
};

struct node{
	int function;
	int *inputs;
	float *weights;	
	int active;
	float output;
};

struct data{
	int numSamples;
	int numInputs;
	int numOutputs;
	
	float **inputData;
	float **outputData;
};


/* Prototypes of functions used internally to CGP-Library */

/* node functions */
struct node *initialiseNode(struct parameters *params, int nodePosition);
void freeNode(struct parameters *params, struct node *n);

/* getting gene value functions  */
float getRandomConnectionWeight(struct parameters *params);
int getRandomNodeInput(struct parameters *params, int nodePosition);
int getRandomFunction(struct parameters *params);
int getRandomChromosomeOutput(struct parameters *params);

/* active node functions */
void setActiveNodes(struct parameters *params, struct chromosome *chromo);
void recursivelySetActiveNodes(struct parameters *params, struct chromosome *chromo, int nodeIndex);

/* mutation functions  */
void probabilisticMutation(struct parameters *params, struct chromosome *chromo);

/* function set functions */
void addFuctionToFunctionSet(struct fuctionSet *funcSet, char *functionName);
void freeFuctionSet(struct fuctionSet *fs);
float supervisedLearning(struct parameters *params, struct chromosome *chromo, struct data *dat);

/* node functions defines in CGP-Library */
float add(const int numInputs, float *inputs, float *weights);
float sub(const int numInputs, float *inputs, float *weights);

/* other */
float randFloat(void);
void bubbleSortInt(int *array, const int length);
void sortPopulation(struct parameters *params, struct population *pop);

/*
	Initialises data structure and assigns values of given file
*/
struct data *initialiseDataFromFile(char *file){
	
	int i;
	struct data *dat;
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
	dat = malloc(sizeof(struct data));
	
	/* for every line in the given file */
	while( (line=fgets(buffer, sizeof(buffer), fp)) != NULL){
	
		/* deal with the first line containing meta data */
		if(lineNum == -1){
						
			sscanf(line, "%d,%d,%d", &(dat->numInputs), &(dat->numOutputs), &(dat->numSamples));
						
			dat->inputData = malloc(dat->numSamples * sizeof(float**));
			dat->outputData = malloc(dat->numSamples * sizeof(float**));
		
			for(i=0; i<dat->numSamples; i++){
				dat->inputData[i] = malloc(dat->numInputs * sizeof(float));
				dat->outputData[i] = malloc(dat->numOutputs * sizeof(float));
			}			
		}
		/* the other lines contain input output pairs */
		else{
			
			/* get the first value on the given line */
			record = strtok(line,",");
			col = 0;
		
			/* until end of line */
			while(record != NULL){
							
				/* if its an input value */				
				if(col < dat->numInputs){
					dat->inputData[lineNum][col] = atof(record);
				}
				
				/* if its an output value */
				else{
					dat->outputData[lineNum][col - dat->numInputs] = atof(record);
				}
				
				/* get the next value on the given line */
				record = strtok(NULL,",");
		
				/* increment the current col index */
				col++;
			}
		}	
		
		/* increment the current line index */
		lineNum++;
	}

	fclose(fp);

	return dat;
}


/*

*/
void freeData(struct data *dat){
	
	int i;
	
	for(i=0; i<dat->numSamples; i++){
		free(dat->inputData[i]);
		free(dat->outputData[i]);
	}
	
	free(dat->inputData);
	free(dat->outputData);
	free(dat);
}


/*
	prints the given data structure to the screen
*/
void printData(struct data *dat){
	
	int i,j;
	
	printf("DATA SET\n");
	printf("Inputs: %d, ", dat->numInputs);
	printf("Outputs: %d, ", dat->numOutputs);
	printf("Samples: %d\n", dat->numSamples);
	
	for(i=0; i<dat->numSamples; i++){
		
		for(j=0; j<dat->numInputs; j++){
			printf("%f ", dat->inputData[i][j]);
		}
		
		printf(" : ");
		
		for(j=0; j<dat->numOutputs; j++){
			printf("%f ", dat->outputData[i][j]);
		}
		
		printf("\n");
	}
}


/*
	Initialises a parameter struct with default values. These 
	values can be individually changed via set functions.
*/
struct parameters *initialiseParameters(const int numInputs, const int numNodes, const int numOutputs, const int arity){
		
	struct parameters *params;
	
	/* allocate memory for parameters */
	params = malloc(sizeof(struct parameters));
		
	/* Set default values */	
	params->mu = 1;
	params->lambda = 4;
	params->stratergy = '+';
	params->mutationRate = 0.5;	
	params->randSeed = 123456789;
	params->connectionsWeightRange = 1;
	params->generations = 1000;
	
	
	params->arity = arity;
	params->numInputs = numInputs;
	params->numNodes = numNodes;
	params->numOutputs = numOutputs;
		
	params->mutationType = probabilisticMutation;
		
	params->funcSet = malloc(sizeof(struct fuctionSet));
	
	params->nodeInputsHold = malloc(params->arity * sizeof(float));
	
	params->fitnessFuction = supervisedLearning;
	
	
	
	/* Seed the random number generator */
	srand(params->randSeed);
	
	return params;
}

/*
	Frees the memory associated with the given parameter structure
*/
void freeParameters(struct parameters *params){
	
	free(params->nodeInputsHold);
	free(params->funcSet);
			
	free(params);
}



/*

*/
struct population *initialisePopulation(struct parameters *params){
	
	int i;
	
	struct population *pop;
	
	pop = malloc(sizeof(struct population));
	
	pop->chromosomes = malloc( (params->mu + params->lambda) * sizeof(struct chromosome) );
	
	for(i=0; i < params->mu + params->lambda; i++){
		pop->chromosomes[i] = initialiseChromosome(params);
	}
	
	return pop;
}

/*

*/
void freePopulation(struct parameters *params, struct population *pop){
	
	int i;
	
	for(i=0; i < params->mu + params->lambda; i++){
		freeChromosome(params, pop->chromosomes[i]);
	}
	
	free(pop->chromosomes);
	free(pop);
}

/*
	returns mu value currently set in given parameters.
*/
int getMu(struct parameters *params){
	return params->mu;
}

/*
	Sets the mu value in given parameters to the new given value. If mu value
	is invalid a warning is displayed and the mu value is left unchanged.
*/
void setMu(struct parameters *params, int mu){
	
	if(mu > 0){
		params->mu = mu;
	}
	else{
		printf("\nWarning: mu value '%d' is invalid. Mu value must have a value of one or greater. Mu value left unchanged as '%d'.\n", mu, params->mu);
	}
}


/*
*/
int getNumInputs(struct parameters *params){
	return params->numInputs;
}

/*
*/
int getNumOutputs(struct parameters *params){
	return params->numOutputs;
}

/*
*/
void setFitnessFuction(struct parameters *params, float (*fitnessFuction)(struct parameters *params, struct chromosome *chromo, struct data *dat)){
	params->fitnessFuction = fitnessFuction;
}

/*
	
*/
void setFuctionSet(struct parameters *params, char *functionNames){

	char *pch;
	char functions[FUNCTIONNAMELENGTH];
			
	/* */
	strcpy(functions, functionNames);
		
	/* */
	params->funcSet->numFunctions = 0;
	
	pch = strtok(functions, ",");
											
	while (pch != NULL){
		addFuctionToFunctionSet(params->funcSet, pch);
		pch = strtok(NULL, " -:,");
	}
		
		
	if(params->funcSet->numFunctions == 0){
		printf("Warning: No Functions added to function set.\n");
	}	
}

/*
	
*/
void addFuctionToFunctionSet(struct fuctionSet *funcSet, char *functionName){
	
	funcSet->numFunctions++;
		
	if(strcmp(functionName, "add") == 0){
		strcpy(funcSet->functionNames[funcSet->numFunctions-1], "add");
		funcSet->functions[funcSet->numFunctions-1] = add;
	}
	else if(strcmp(functionName, "sub") == 0){
		strcpy(funcSet->functionNames[funcSet->numFunctions-1], "sub");
		funcSet->functions[funcSet->numFunctions-1] = sub;
	}
	else{
		printf("Warning: function '%s' is not known and was not added.\n", functionName);
		funcSet->numFunctions--;
	}	
}

/*
	Prints the current functions in the function set to
	the terminal.  
*/
void printFuctionSet(struct parameters *params){

	int i;

	printf("Functions (%d):", params->funcSet->numFunctions);
	
	for(i=0; i<params->funcSet->numFunctions; i++){
		printf(" %s", params->funcSet->functionNames[i]);
	}
	
	printf("\n");
}

/*
	Returns a pointer to an initialised chromosome with values obeying the given parameters.
*/
struct chromosome *initialiseChromosome(struct parameters *params){
	
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
		chromo->nodes[i] = initialiseNode(params, i);
	}
		
	/* set each of the chromosomes outputs */
	for(i=0; i<params->numOutputs; i++){
		chromo->outputNodes[i] = getRandomChromosomeOutput(params);
	}
	
	/* Add all nodes to the active node matrix */
	for(i=0; i<params->numNodes; i++){
		chromo->activeNodes[i] = i;
	}
	
	/* set the number of inputs and outputs */
	chromo->numInputs = params->numInputs;
	chromo->numOutputs = params->numOutputs;
	
	/* set the number of active node to the number of nodes (all active) */
	chromo->numActiveNodes = params->numNodes;
	
	return chromo;
}

/*
	Frees the memory associated with the given chromosome structure
*/
void freeChromosome(struct parameters *params, struct chromosome *chromo){

	int i;
		
	for(i=0; i < params->numNodes; i++){
		
		freeNode(params, chromo->nodes[i]);
	}
	
	free(chromo->outputValues);
	free(chromo->nodes);
	free(chromo->outputNodes);
	free(chromo->activeNodes);	
	free(chromo);
}



/*
	Evolves the given population using the parameters specified in the given parameters. The data 
	structure can be used by the fitness function if required, otherwise set as NULL. If the data
	structure is set as NULL the fitness function will be reset to the default fitness function - 
	one suited to supervised learning. 
*/
float evolvePopulation(struct parameters *params, struct population *pop, struct data *dat){
	
	int i;
	int gen;
	
	if(dat != NULL){
		
	}
	
	/* set fitness of the mu members of the population */
	for(i=0; i<params->mu; i++){
		pop->chromosomes[i]->fitness = params->fitnessFuction(params, pop->chromosomes[i], dat);
	}
	
	/* for each generation */
	for(gen=0; gen<params->generations; gen++){
		
		/* set fitness of the lambda members of the population */
		for(i=params->mu; i<params->lambda; i++){
			pop->chromosomes[i]->fitness = params->fitnessFuction(params, pop->chromosomes[i], dat);
		}
		
		/* sort the population into fitness order */
		
		/* check termination conditions */
		
		/* generate the next generation */
		
	}
	
	
	return pop->chromosomes[0]->fitness;
}

/* */
void sortPopulation(struct parameters *params, struct population *pop){
	
	struct chromosome *chromoTmp;
	int i;
	int finished = 0;
		
	while(finished == 0){
		
		finished = 1;
		
		for(i=0; i < params->mu + params->lambda -1; i++){
			
			if(pop->chromosomes[i]->fitness > pop->chromosomes[i+1]->fitness){
				
				finished = 0;
				chromoTmp = pop->chromosomes[i];
				pop->chromosomes[i] = pop->chromosomes[i+1];
				pop->chromosomes[i+1] = chromoTmp;
			}
		}
	}
}


/*
	Executes the given chromosome with the given outputs and placed the outputs in 'outputs'.
*/
void executeChromosome(struct parameters *params, struct chromosome *chromo, float *inputs, float *outputs){
	
	int i,j;
	int nodeInputLocation;
	int currentActiveNode;
	int currentActiveNodeFuction;
		
	/* for all of the active nodes */
	for(i=0; i<chromo->numActiveNodes; i++){
		
		/* get the index of the current active node */
		currentActiveNode = chromo->activeNodes[i];
		
		/* for each of the active nodes inputs */
		for(j=0; j<params->arity; j++){
			
			/* gather the nodes inputs */
			nodeInputLocation = chromo->nodes[currentActiveNode]->inputs[j];
			
			if(nodeInputLocation < params->numInputs){
				params->nodeInputsHold[j] = inputs[nodeInputLocation];
			}
			else{
				params->nodeInputsHold[j] = chromo->nodes[nodeInputLocation - params->numInputs]->output;
			}
		}
		
		/* get the index of the active node under evaluation */
		currentActiveNodeFuction = chromo->nodes[currentActiveNode]->function;
		
		/* calculate the output of the active node under evaluation */
		chromo->nodes[currentActiveNode]->output = params->funcSet->functions[currentActiveNodeFuction](params->arity, params->nodeInputsHold, chromo->nodes[currentActiveNode]->weights);
	}
	
	/* Set the chromosome outputs */
	for(i=0; i<params->numOutputs; i++){
	
		if(chromo->outputNodes[i] < params->numInputs){
			outputs[i] = inputs[chromo->outputNodes[i]];
		}
		else{
			outputs[i] = chromo->nodes[chromo->outputNodes[i] - params->numInputs]->output;
		}
	}
}


/*
	returns a pointer to an initialised node. Initialised means that necessary
	memory has been allocated and values set.
*/
struct node *initialiseNode(struct parameters *params, int nodePosition){
		
	struct node *n;
	int i;
	
	/* allocate memory for node */
	n = malloc(sizeof(struct node));
	
	/* allocate memory for the node's inputs and connection weights */
	n->inputs = malloc(params->arity * sizeof(int));
	n->weights = malloc(params->arity * sizeof(float));	

	/* set the node's function */
	n->function = getRandomFunction(params);

	/* set as active by default */
	n->active = 1;

	/* set the nodes inputs and connection weights */
	for(i=0; i<params->arity; i++){
		n->inputs[i] = getRandomNodeInput(params,nodePosition);
		n->weights[i] = getRandomConnectionWeight(params);
	}

	return n;
}


/*

*/
void freeNode(struct parameters *params, struct node *n){
	
	free(n->inputs);
	free(n->weights);
	free(n);
}

/* 
	returns a random connection weight value
*/
float getRandomConnectionWeight(struct parameters *params){
	return (randFloat() * 2 * params->connectionsWeightRange) - params->connectionsWeightRange;
}

/*
	returns a random function index
*/
int getRandomFunction(struct parameters *params){
	
	/* check that funcSet contains functions*/
	if(params->funcSet->numFunctions <1){
		printf("Error: cannot assign the function gene a value as the Fuction Set is empty.\nTerminating CGP-Library.\n");
		exit(0);
	}
	
	return rand() % (params->funcSet->numFunctions);
}

/*
 returns a random input for the given node
*/
int getRandomNodeInput(struct parameters *params, int nodePosition){
	
	int input;
	
	input = rand() % (params->numInputs + nodePosition); 
	
	return input;
}
	
/* 
	set the active nodes in the given chromosome
*/
void setActiveNodes(struct parameters *params, struct chromosome *chromo){
	
	int i;	
	
	/* set the number of active nodes to zero */
	chromo->numActiveNodes = 0;
	
	/* reset the active nodes */
	for(i = 0; i < params->numNodes; i++){
		chromo->nodes[i]->active = 0;
	}
	
	/* start the recursive search for active nodes from the output nodes for the number of output nodes */
	for(i=0; i < params->numOutputs; i++){
			
		/* if the output connects to a chromosome input, skip */	
		if(chromo->outputNodes[i] < params->numInputs){
			continue; 
		}

		/* begin a recursive search for active nodes */
		recursivelySetActiveNodes(params, chromo, chromo->outputNodes[i]);
	}
	
	/* place active nodes in order */
	bubbleSortInt(chromo->activeNodes, chromo->numActiveNodes);
}	
	
/* 
	used by setActiveNodes to recursively search for active nodes
*/
void recursivelySetActiveNodes(struct parameters *params, struct chromosome *chromo, int nodeIndex){
 
	int i;	

	/* if the given node is an input, stop */
	if(nodeIndex < params->numInputs){
		return;
	}
	 
	/* if the given node has already been flagged as active */
	if(chromo->nodes[nodeIndex - params->numInputs]->active == 1){
		return;
	}
	
	/* log the node as active */
	chromo->nodes[nodeIndex - params->numInputs]->active = 1;
	chromo->activeNodes[chromo->numActiveNodes] = nodeIndex - params->numInputs;
	chromo->numActiveNodes++;			
					
	/* recursively log all the nodes to which the current nodes connect as active */
	for(i=0; i < params->arity; i++){
		recursivelySetActiveNodes(params, chromo, chromo->nodes[nodeIndex - params->numInputs]->inputs[i]);
	}
}
	
	
/*
	returns a random chromosome output
*/	
int getRandomChromosomeOutput(struct parameters *params){
	
	int output;
	
	output = rand() % (params->numInputs + params->numNodes);
	
	return output;
}
	
/*
	Prints the given chromosome to the screen
*/	 
void printChromosome(struct parameters *params, struct chromosome *chromo){

	int i,j;			
				
	/* set the active nodes in the given chromosome */
	setActiveNodes(params, chromo);
								
	/* for all the chromo inputs*/
	for(i=0; i<params->numInputs; i++){
		printf("(%d):\tinput\n", i);
	}
	
	/* for all the hidden nodes */
	for(i = 0; i < params->numNodes; i++){ 
	
		/* print the node function */
		printf("(%d):\t%d ", params->numInputs + i, chromo->nodes[i]->function);
		
		/* for the arity of the node */
		for(j = 0; j < params->arity; j++){

			/* print the node input information */
			printf("%d,%+.1f  ", chromo->nodes[i]->inputs[j], chromo->nodes[i]->weights[j]);
		}
		
		/* Highlight active nodes */
		if(chromo->nodes[i]->active == 1){
			printf("*");
		}
		
		printf("\n");
	}

	/* for all of the outputs */
	printf("outputs: ");
	for(i = 0; i < params->numOutputs; i++){
		
		/* print the output node locations */
		printf("%d ", chromo->outputNodes[i]);
	}
	
	printf("\n");
}	

/*
	Mutates the given chromosome using the mutation method described in parameters
*/
void mutateChromosome(struct parameters *params, struct chromosome *chromo){
	
	params->mutationType(params, chromo);
}

/*
	Conductions probabilistic mutation on the give chromosome. Each chromosome
	gene is changed to a random valid allele with a probability specified in
	parameters.
*/
void probabilisticMutation(struct parameters *params, struct chromosome *chromo){
	
	int i,j;
	
	/* for every nodes in the chromosome */
	for(i=0; i<params->numNodes; i++){
		
		/* mutate the function gene */
		if(randFloat() <= params->mutationRate){
			chromo->nodes[i]->function = getRandomFunction(params);
		}
		
		/* for every input to each chromosome */
		for(j=0; j<params->arity; j++){
			
			/* mutate the node input */
			if(randFloat() <= params->mutationRate){
				chromo->nodes[i]->inputs[j] = getRandomNodeInput(params, i);
			}
			
			/* mutate the node connection weight */
			if(randFloat() <= params->mutationRate){
				chromo->nodes[i]->weights[j] = getRandomConnectionWeight(params);
			}
		}
	}
	
	/* for every chromosome output */ 
	for(i=0; i<params->numOutputs; i++){
		
		/* mutate the chromosome output */
		if(randFloat() <= params->mutationRate){
			chromo->outputNodes[i] = getRandomChromosomeOutput(params);
		}
	}
}

/*
	Node function add. Returns the sum of the inputs. 
*/ 	
float add(const int numInputs, float *inputs, float *weights){
	
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
float sub(const int numInputs, float *inputs, float *weights){
	
	int i;
	float sum = inputs[0];
	
	for(i=1; i<numInputs; i++){
		sum -= inputs[i];
	}
	
	return sum;
}	
	
/*
*/
float supervisedLearning(struct parameters *params, struct chromosome *chromo, struct data *dat){

	int i,j;
	float error = 0;
	

	/* error checking */
	if(chromo->numInputs != dat->numInputs){
		printf("Error: the number of chromosome inputs specified in the chromosome must match the number of required inputs specified in the data.\n");
		printf("Terminating CGP-Library.\n");
		exit(0);
	}

	if(chromo->numOutputs != dat->numOutputs){
		printf("Error: the number of chromosome outputs specified in the chromosome must match the number of required outputs specified in the data.\n");
		printf("Terminating CGP-Library.\n");
		exit(0);
	}

	
	/* for each sample in data */
	for(i=0 ; i<dat->numSamples; i++){
	
		/* calculate the chromosome outputs for the set of inputs  */
		executeChromosome(params, chromo, dat->inputData[i], chromo->outputValues);
	
		/* for each chromosome output */
		for(j=0; j<chromo->numOutputs; j++){
			
			error += chromo->outputValues[j] - dat->outputData[i][j];
		}	
	}

	return error;
}
/* 
	returns a random float between [0,1]
*/
float randFloat(void){
	return (float)rand()/(float)RAND_MAX;
}

/*
	simple bad sort - replace with standard qsort...
*/
void bubbleSortInt(int *array, const int length){
	
	int i,tmp;
	int finished = 0;
		
	while(finished == 0){
		
		finished = 1;
		
		for(i=0; i < length-1; i++){
			
			if(array[i] > array[i+1]){
				
				finished = 0;
				tmp = array[i];
				array[i] = array[i+1];
				array[i+1] = tmp;
			}
		}
	}
}

