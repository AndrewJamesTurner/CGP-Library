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
	int numNodes;
	int arity;
	
	int numFuctions;
	int numInputs;
	
	float connectionsWeightRange;
};

struct chromosome{
	
	struct node *nodes;
	int *outputNodes;
	
	int numactiveNodes;
	int *activeNodes;
	
	int fitness;
	
};

struct node{
	
	int fuction;
	int *inputs;
	float *weights;	
};


/* 
	Prototypes of fuctions used interally to CGP-Library
*/
float getRandomConnectionWeight(struct parameters *params);
int getRandomNodeInput(struct parameters *params, int nodePosition);
int getRandomFuction(struct parameters *params);
int getRandomChromosomeOutput(struct parameters *params);
float randFloat(void);

/*
	Initialises a parameter structs with default values. These
	values can be indevidually chaged via set fuctions.
*/
struct parameters *initialiseParameters(void){
		
	struct parameters *params;
	
	/* allocate memory for parameters */
	params = malloc(sizeof(struct parameters));
		
	/* Set default values*/	
	params->mu = 1;
	params->lambda = 4;
	params->mutationRate = 0.03;	
	
	params->connectionsWeightRange = 1;
	
	return params;
}


int getMu(struct parameters *params){
	return params -> mu;
}


struct node *nodes;
	int *outputNodes;
	
	int numactiveNodes;
	int *activeNodes;
	
	int fitness;


/*
*/
struct chromosome *initialiseChromosome(struct parameters *params){
	
	struct chromosome *chromo;
	int i;
	
	/* allocate memory for chromosome */
	chromo = malloc(sizeof(struct chromosome));
	
	/* allocate memory for nodes and acitve nodes matrix */
	for(i=0; i<params->numNodes; i++){
		
	}

}


/*
	returns a pointer to an initilised node. Initilised means that nessassary
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

	/* set the node's fuction */
	n->fuction = getRandomFuction(params);

	/* set the nodes inputs and connnection weights */
	for(i=0; i<params->arity; i++){
		n->inputs[i] = getRandomNodeInput(params,i);
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
	returns a random fuction index
*/
int getRandomFuction(struct parameters *params){
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
	returns a random chromosome output
*/	
int getRandomChromosomeOutput(struct parameters *params){
	
	int output;
	
	output = rand() % (params->numInputs + params->numNodes);
	
	return output;
}
	
/* 
	returns a random flaot between [0,1]
*/
float randFloat(void){
	return (float)rand()/(float)RAND_MAX;
}




