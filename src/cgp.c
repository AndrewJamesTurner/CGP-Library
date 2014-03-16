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
#include "libCGP.h"

struct parameters{
	
	int mu;
	int lambda;
	float mutationRate;
	float connectionsWeightRange;
	unsigned int randSeed;
	
	int numFuctions;
	int numInputs;
	int numNodes;
	int numOutputs;
	int arity;
	
	void (*mutationType)(struct parameters *params, struct chromosome *chromo);
};

struct chromosome{
	struct node **nodes;
	int *outputNodes;
	int numActiveNodes;
	int *activeNodes;
	int fitness;
};

struct node{
	int function;
	int *inputs;
	float *weights;	
	int active;
};


/* Prototypes of functions used internally to CGP-Library */
struct node *initialiseNode(struct parameters *params, int nodePosition);

float getRandomConnectionWeight(struct parameters *params);
int getRandomNodeInput(struct parameters *params, int nodePosition);
int getRandomFunction(struct parameters *params);
int getRandomChromosomeOutput(struct parameters *params);

void setActiveNodes(struct parameters *params, struct chromosome *chromo);
void recursivelySetActiveNodes(struct parameters *params, struct chromosome *chromo, int nodeIndex);

void probabilisticMutation(struct parameters *params, struct chromosome *chromo);

float randFloat(void);
void bubbleSortInt(int *array, const int length);

/*
	Initialises a parameter struct with default values. These
	values can be individually changed via set functions.
*/
struct parameters *initialiseParameters(void){
		
	struct parameters *params;
	
	/* allocate memory for parameters */
	params = malloc(sizeof(struct parameters));
		
	/* Set default values*/	
	params->mu = 1;
	params->lambda = 4;
	params->mutationRate = 0.5;	
	params->randSeed = 123456789;
	params->connectionsWeightRange = 1;
	params->arity = 2;
	params->numInputs = 2;
	params->numNodes = 10;
	params->numOutputs = 1;
	params->numFuctions = 4;
	
	params->mutationType = probabilisticMutation;
	
	
	/* Seed the random number generator */
	srand(params->randSeed);
	
	return params;
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
	Returns a pointer to an initialised chromosome with values obeying the given parameters.
*/
struct chromosome *initialiseChromosome(struct parameters *params){
	
	struct chromosome *chromo;
	int i;
	
	/* allocate memory for chromosome */
	chromo = malloc(sizeof(struct chromosome));
		
	/* allocate memory for nodes */
	chromo->nodes = malloc(params->numNodes * sizeof(struct node *));
	
	/* allocate memory for outputNodes matrix */
	chromo->outputNodes = malloc(params->numOutputs * sizeof(int));
	
	/* allocate memory for active nodes matrix */
	chromo->activeNodes = malloc(params->arity * sizeof(int));

	/* Initialise each of the chromosomes nodes */
	for(i=0; i<params->numNodes; i++){
		chromo->nodes[i] = initialiseNode(params, i);
	}
		
	/* set each of the chromosomes outputs */
	for(i=0; i<params->numOutputs; i++){
		chromo->outputNodes[i] = getRandomChromosomeOutput(params);
	}
	
	return chromo;
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

	/* set the nodes inputs and connection weights */
	for(i=0; i<params->arity; i++){
		n->inputs[i] = getRandomNodeInput(params,nodePosition);
		n->weights[i] = getRandomConnectionWeight(params);
	}

	return n;
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
	return rand() % (params->numFuctions);
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

