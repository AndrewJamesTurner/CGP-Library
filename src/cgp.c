/*
	This file is part of CGP-Library
	Copyright (c) Andrew James Turner 2014 (andrew.turner@york.ac.uk)

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
#include <time.h>
#include <math.h>
#include <float.h>

#include "include/cgp.h" 

#define FUNCTIONSETSIZE 50
#define FUNCTIONNAMELENGTH 512
#define FITNESSFUNCTIONNAMELENGTH 512
#define MUTATIONTYPENAMELENGTH 512
#define SELECTIONSCHEMENAMELENGTH 512
#define REPRODUCTIONSCHEMENAMELENGTH 512

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
	float targetFitness;
	struct functionSet *funcSet;
	void (*mutationType)(struct parameters *params, struct chromosome *chromo);
	char mutationTypeName[MUTATIONTYPENAMELENGTH];
	float (*fitnessFunction)(struct parameters *params, struct chromosome *chromo, struct data *dat);	
	char fitnessFunctionName[FITNESSFUNCTIONNAMELENGTH];
	void (*selectionScheme)(struct parameters *params, struct chromosome **parents, struct chromosome **candidateChromos, int numCandidateChromos);
	char selectionSchemeName[SELECTIONSCHEMENAMELENGTH];
	void (*reproductionScheme)(struct parameters *params, struct population *pop);
	char reproductionSchemeName[REPRODUCTIONSCHEMENAMELENGTH];
	int updateFrequency;
};

struct population{
		
	int mu;
	int lambda;	
	struct chromosome **parents;
	struct chromosome **children;
	int trainedGenerations;
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

struct data{
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
static void copyChromosome(struct chromosome *chromoDest, struct chromosome *chromoSrc);

/* node functions */
static struct node *initialiseNode(struct parameters *params, int nodePosition);
static void freeNode(struct node *n);
static void copyNode(struct node *nodeDest, struct node *nodeSrc);

/* getting gene value functions  */
static float getRandomConnectionWeight(float weightRange);
static int getRandomNodeInput(int numChromoInputs, int nodePosition);
static int getRandomFunction(int numFunctions);
static int getRandomChromosomeOutput(int numInputs, int numNodes);

/* active node functions */
static void setActiveNodes(struct chromosome *chromo);
static void recursivelySetActiveNodes(struct chromosome *chromo, int nodeIndex);

/* function set functions */
static void addPresetFuctionToFunctionSet(struct parameters *params, char *functionName);

/* mutation functions  */
static void probabilisticMutation(struct parameters *params, struct chromosome *chromo);
static void pointMutation(struct parameters *params, struct chromosome *chromo);

/* selection scheme functions */
static void pickHighest(struct parameters *params, struct chromosome **parents, struct chromosome **candidateChromos, int numCandidateChromos);	

/* reproduction scheme functions */
static void mutateRandomParent(struct parameters *params, struct population *pop);

/* fitness function */
static float supervisedLearning(struct parameters *params, struct chromosome *chromo, struct data *dat);

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
static void bubbleSortInt(int *array, const int length);
static void sortChromosomeArray(struct chromosome **chromoArray, int numChromos);
static void copyFuctionSet(struct functionSet *funcSetDest, struct functionSet *funcSetSrc);

/* population functions */
struct population *initialisePopulation(struct parameters *params);
void freePopulation(struct population *pop);

struct chromosome *getFittestChromosome(struct population *pop);	
int getNumberOfGenerations(struct population *pop);
struct results* initialiseResults(struct parameters *params, int numRuns);
static float sumWeigtedInputs(const int numInputs, const float *inputs, const float *connectionWeights);


void saveChromosome(struct chromosome *chromo, char *file){
	
	int i,j;
	FILE *fp;
	
	/* create the chromosome file */
	fp = fopen(file, "w");
	
	/* ensure that the file was created correctly */
	if(fp == NULL){
		printf("Warning: cannot save chromosome to '%s'. Chromosome was not saved.\n", file);
		return;
	}
		
	fprintf(fp, "numInputs,%d\n", chromo->numInputs);
	fprintf(fp, "numNodes,%d\n", chromo->numNodes);
	fprintf(fp, "numOutputs,%d\n", chromo->numOutputs);
	fprintf(fp, "arity,%d\n", chromo->arity);
	
	fprintf(fp, "fuctionSet");
	
	for(i=0; i<chromo->funcSet->numFunctions; i++){
		fprintf(fp, ",%s", chromo->funcSet->functionNames[i]);	
	}
	fprintf(fp, "\n");
	
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

struct chromosome* loadChromosome(char *file){
	
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
	printf("numInputs: %d\n", numInputs);
	
	/* get num nodes */
	line = fgets(buffer, sizeof(buffer), fp);
	if(line == NULL){/*error*/}
	record = strtok(line,",");
	record = strtok(NULL,",");
	numNodes = atoi(record);
	printf("numNodes: %d\n", numNodes);
	
	/* get num outputs */
	line = fgets(buffer, sizeof(buffer), fp);
	if(line == NULL){/*error*/}
	record = strtok(line,",");
	record = strtok(NULL,",");
	numOutputs = atoi(record);
	printf("numOutputs: %d\n", numOutputs);
	
	/* get arity */
	line = fgets(buffer, sizeof(buffer), fp);
	if(line == NULL){/*error*/}
	record = strtok(line,",");
	record = strtok(NULL,",");
	arity = atoi(record);
	printf("arity: %d\n", arity);
	
	/* initialise parameters  */
	params = initialiseParameters(numInputs, numNodes, numOutputs, arity);
	
	/* get and set node functions */
	line = fgets(buffer, sizeof(buffer), fp);
	if(line == NULL){/*error*/}
	record = strtok(line,",\n");
	record = strtok(NULL,",\n");
	
	/* for each function name */
	while( record != NULL){
		
		strcpy(funcName, record);				
		addPresetFuctionToFunctionSet(params,funcName);
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
			
	free(params);
	
	return chromo;
}


void setTargetFitness(struct parameters *params, float targetFitness){

	/* error checking */
	
	
	params->targetFitness = targetFitness;
}

/*
*/
void setUpdateFrequency(struct parameters *params, int updateFrequency){
	params->updateFrequency = updateFrequency;
}


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

void freeResults(struct results *rels){
	
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
	returns the average number of chromosome active nodes from repeated 
	run results specified in rels.
*/
float getAverageActiveNodes(struct results *rels){
	
	int i;
	float avgActiveNodes = 0;
	struct chromosome *chromoTemp;
	
	
	for(i=0; i<rels->numRuns; i++){
		
		chromoTemp = rels->bestChromosomes[i];

		avgActiveNodes += getChromosomeActiveNodes(chromoTemp);
	}
	
	avgActiveNodes = avgActiveNodes / rels->numRuns;
	
	return avgActiveNodes;
}

/*
	returns the average chromosome fitness from repeated 
	run results specified in rels.
*/
float getAverageFitness(struct results *rels){
	
	int i;
	float avgFit = 0;
	struct chromosome *chromoTemp;
	
	
	for(i=0; i<rels->numRuns; i++){
		
		chromoTemp = rels->bestChromosomes[i];

		avgFit += getChromosomeFitness(chromoTemp);
	}
	
	avgFit = avgFit / rels->numRuns;
	
	return avgFit;
}

/*
	returns the number of generations used by each run  specified in rels.
*/
float getAverageGenerations(struct results *rels){
	
	int i;
	float avgGens = 0;
	struct chromosome *chromoTemp;
	
	
	for(i=0; i<rels->numRuns; i++){
		
		chromoTemp = rels->bestChromosomes[i];

		avgGens += getChromosomeGenerations(chromoTemp);
	}
	
	avgGens = avgGens / rels->numRuns;
	
	return avgGens;
}

int getChromosomeGenerations(struct chromosome *chromo){
	return chromo->generation;
}

/*
	returns a pointer to the best chromosomes found on the given run in rels.
*/
struct chromosome* getChromosome(struct results *rels, int run){
	
	/* do some error checking */
	
	return rels->bestChromosomes[run];
}


struct results* repeatCGP(struct parameters *params, struct data *dat, int numGens, int numRuns){
	
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
		rels->bestChromosomes[i] = runCGP(params, dat, numGens);
		
		printf("%d\t%f\t%d\t\t%d\n", i, rels->bestChromosomes[i]->fitness, rels->bestChromosomes[i]->generation, rels->bestChromosomes[i]->numActiveNodes);	
	}
	
	printf("----------------------------------------------------\n");
	printf("AVG\t%f\t%f\t%f\n", getAverageFitness(rels), getAverageGenerations(rels), getAverageActiveNodes(rels));	
	printf("----------------------------------------------------\n");
	
	/* restore the original value for the update frequency */	
	params->updateFrequency = updateFrequency;
	
	return rels;
}


struct chromosome* runCGP(struct parameters *params, struct data *dat, int numGens){
	
	int i;
	int gen;
	
	/* */
	struct chromosome *chromo;
	struct chromosome *bestChromo;
	
	/* */
	struct population *pop;
	
	/* storage for chromosomes used by selection scheme */
	struct chromosome **candidateChromos;
	int numCandidateChromos;
		
	/* */
	if(numGens < 0){
		printf("Error: %d generations is invalid. The number of generations must be >= 0. \n Terminating CGP-Library.\n", numGens);
	}	
		
	/* */
	chromo = initialiseChromosome(params);
		
	/* */
	pop = initialisePopulation(params);
		
	/* determine the size of the Candidate Chromos based on the evolutionary Strategy */	
	if(params->evolutionaryStrategy == '+'){
		numCandidateChromos = params->mu + params->lambda;
	}
	else{
		numCandidateChromos = params->lambda;
	}	
			
	/* initialise the candidateChromos */	
	candidateChromos = malloc(numCandidateChromos * sizeof(struct chromosome *));	
		
	for(i=0; i<numCandidateChromos; i++){
		candidateChromos[i] = initialiseChromosome(params);
	}
		
	/* if using '+' evolutionary strategy */
	if(params->evolutionaryStrategy == '+'){
		
		/* set fitness of the parents */
		for(i=0; i<params->mu; i++){
			setActiveNodes(pop->parents[i]);
			setChromosomeFitness(params, pop->parents[i], dat);
		}
	}
	
	if(params->updateFrequency != 0){
		
		printf("-- Starting CGP --\n\n");
		
		printf("Gen\tfit\n");
	}
	
	/* for each generation */
	for(gen=0; gen <numGens; gen++){
		
		/* set fitness of the children of the population */
		for(i=0; i< params->lambda; i++){
			setActiveNodes(pop->children[i]);
			setChromosomeFitness(params, pop->children[i], dat);
		}
					
		/* check termination conditions */
		bestChromo = getFittestChromosome(pop);
		if(bestChromo->fitness <= params->targetFitness){
			
			if(params->updateFrequency != 0){
				printf("%d\t%f - Solution Found\n", gen, bestChromo->fitness);
			}
			
			break;
		}			
				
			
		/* 
			Set the chromosomes which will be used by the selection scheme
			dependant upon the evolutionary strategy
		*/
		for(i=0; i<numCandidateChromos; i++){
			
			if(i < params->lambda){
				copyChromosome(candidateChromos[i], pop->children[i] );
			}
			else{
				copyChromosome(candidateChromos[i], pop->parents[i - params->lambda] );
			}
		}
						
		/* select the parents */		
		params->selectionScheme(params, pop->parents, candidateChromos, numCandidateChromos);
						
		
		
		/* */
		if(params->updateFrequency != 0 && (gen % params->updateFrequency == 0 || gen >= numGens-1) ){
			printf("%d\t%f\n", gen, pop->parents[0]->fitness);
		}
			
		/* create the children */
		params->reproductionScheme(params, pop);	
	}
	
		
	pop->trainedGenerations = gen;
	
	
	bestChromo = getFittestChromosome(pop);
	bestChromo->generation = gen;
	copyChromosome(chromo, bestChromo);
	
	
	
	for(i=0; i<numCandidateChromos; i++){
		freeChromosome(candidateChromos[i]);
	}
	free(candidateChromos);
	
	
	freePopulation(pop);
	
	
	return chromo;
}


/* */
static void copyFuctionSet(struct functionSet *funcSetDest, struct functionSet *funcSetSrc){
	
	int i;
	
	funcSetDest->numFunctions = funcSetSrc->numFunctions;
	

	for(i=0; i<funcSetDest->numFunctions; i++){
		
		strcpy(funcSetDest->functionNames[i], funcSetSrc->functionNames[i]);	
		funcSetDest->functions[i] = funcSetSrc->functions[i];
	}
}




/*
	Returns the number of generations the given population
	was ran for before terminating the search
*/
int getNumberOfGenerations(struct population *pop){
	return pop->trainedGenerations;
}




/*
	return the fitness of the given chromosome
*/
float getChromosomeFitness(struct chromosome *chromo){
	return chromo->fitness;
}

/*

*/
int getChromosomeActiveNodes(struct chromosome *chromo){
	return chromo->numActiveNodes;
}


/*
	Returns a pointer to the best chromosome in the given population
*/
struct chromosome *getFittestChromosome(struct population *pop){

	struct chromosome *bestChromo;
	float bestFitness;
	int i;
	
	/* set the best chromosome to be the first parent */
	bestChromo = pop->parents[0];
	bestFitness = pop->parents[0]->fitness;
	
	/* for all the parents except the first parent */
	for(i=1; i<pop->mu; i++){
		
		/* set this parent to be the best chromosome if it is fitter than the current best */
		if(pop->parents[i]->fitness < bestFitness){
			
			bestFitness = pop->parents[i]->fitness;
			bestChromo = pop->parents[i];
		}
	}
	
	/* for all the children */
	for(i=0; i<pop->lambda; i++){
		
		/* set this child to be the best chromosome if it is fitter than the current best */
		if(pop->children[i]->fitness < bestFitness){
			
			bestFitness = pop->children[i]->fitness;
			bestChromo = pop->children[i];
		}
	}
	
	return bestChromo;
}




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
struct data *initialiseDataFromArrays(int numInputs, int numOutputs, int numSamples, float *inputs, float *outputs){
	
	int i,j;
	struct data *dat;
	
	/* initialise memory for data structure */
	dat = malloc(sizeof(struct data));
	
	dat->numInputs = numInputs;
	dat->numOutputs = numOutputs;
	dat->numSamples = numSamples; 
	
	dat->inputData = malloc(dat->numSamples * sizeof(float**));
	dat->outputData = malloc(dat->numSamples * sizeof(float**));
		
	for(i=0; i<dat->numSamples; i++){
		
		dat->inputData[i] = malloc(dat->numInputs * sizeof(float));
		dat->outputData[i] = malloc(dat->numOutputs * sizeof(float));
	
		for(j=0; j<dat->numInputs; j++){
			dat->inputData[i][j] = inputs[(i*dat->numInputs) + j];
		}
	
		for(j=0; j<dat->numOutputs; j++){
			dat->outputData[i][j] = outputs[(i*dat->numOutputs) + j];
		}
	
	}			
	
	return dat;
}


/*

*/
void freeData(struct data *dat){
	
	int i;
	
	/* attempt to prevent user double freeing */	
	if(dat == NULL){
		return;
	}	
	
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
	params->evolutionaryStrategy = '+';
	params->mutationRate = 0.05;	
	params->connectionWeightRange = 1;
		
	params->targetFitness = 0;
	
	params->updateFrequency = 100;
	
	params->arity = arity;
	params->numInputs = numInputs;
	params->numNodes = numNodes;
	params->numOutputs = numOutputs;
		
	params->mutationType = probabilisticMutation;
	strcpy(params->mutationTypeName, "probabilistic"); 
		 
	params->funcSet = malloc(sizeof(struct functionSet));
	params->funcSet->numFunctions = 0;
	
	params->fitnessFunction = supervisedLearning;
	strcpy(params->fitnessFunctionName, "supervisedLearning");
	
	params->selectionScheme = pickHighest;
	strcpy(params->selectionSchemeName, "pickHighest");
		
	params->reproductionScheme = mutateRandomParent;
	strcpy(params->reproductionSchemeName, "mutateRandomParent"); 
	
	/* Seed the random number generator */
	srand(time(NULL));
	
	return params;
}




/*
	Frees the memory associated with the given parameter structure
*/
void freeParameters(struct parameters *params){
	
	/* attempt to prevent user double freeing */	
	if(params == NULL){
		return;
	}	
	
	free(params->funcSet);
			
	free(params);
}



/*
	Returns a pointer to an initialised population 
*/
struct population *initialisePopulation(struct parameters *params){
	
	int i;
				
	struct population *pop;
	
	pop = malloc(sizeof(struct population));
	
	pop->mu = params->mu;
	pop->lambda = params->lambda;
	
	pop->parents = malloc( pop->mu * sizeof(struct chromosome) );
	pop->children = malloc( pop->lambda * sizeof(struct chromosome) );
	
	for(i=0; i < pop->mu ; i++){
		pop->parents[i] = initialiseChromosome(params);
	}
	
	for(i=0; i < pop->lambda ; i++){
		pop->children[i] = initialiseChromosome(params);
	}
		
	pop->trainedGenerations = -1;
		
	return pop;
}

/*

*/
void freePopulation(struct population *pop){
	
	int i;
	
	/* attempt to prevent user double freeing */	
	if(pop == NULL){
		return;
	}	
	
	for(i=0; i < pop->mu; i++){
		freeChromosome(pop->parents[i]);
	}
	
	for(i=0; i < pop->lambda; i++){
		freeChromosome(pop->children[i]);
	}
	
	
	free(pop->parents);
	free(pop->children);
		
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
	Sets the lambda value in given parameters to the new given value. 
	If lambda value is invalid a warning is displayed and the lambda value 
	is left unchanged.
*/
void setLambda(struct parameters *params, int lambda){
	
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
void setEvolutionaryStrategy(struct parameters *params, char evolutionaryStrategy){
	
	if(evolutionaryStrategy == '+' || evolutionaryStrategy == ','){
		params->evolutionaryStrategy = evolutionaryStrategy;
	}
	else{
		printf("\nWarning: the evolutionary strategy '%c' is invalid. The evolutionary strategy must be '+' or ','. The evolutionary strategy has been left unchanged as '%c'.\n", evolutionaryStrategy, params->evolutionaryStrategy);
	}	
}



void setMutationType(struct parameters *params, char *mutationType){

	if(strcmp(mutationType, "probabilistic") == 0){
		
		params->mutationType = probabilisticMutation;
		strcpy(params->mutationTypeName, "probabilistic");
	}
	
	else if(strcmp(mutationType, "point") == 0){
		
		params->mutationType = pointMutation;
		strcpy(params->mutationTypeName, "point");
	}
	
	else{
		printf("\nWarning: mutation type '%s' is invalid. The mutation type must be 'probabilistic' or 'point'. The mutation type has been left unchanged as '%s'.\n", mutationType, params->mutationTypeName);
	}
}


/*
	Sets the mutation rate given in parameters. IF an invalid mutation
	rate is given a warning is displayed and the mutation rate is left
	unchanged.
*/
void setMutationRate(struct parameters *params, float mutationRate){

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
void setConnectionWeightRange(struct parameters *params, float weightRange){

	params->connectionWeightRange = weightRange;
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
	sets the fitness function to the fitnessFuction passed. If the fitnessFuction is NULL 
	then the default supervisedLearning fitness function is used. 
*/
void setFitnessFunction(struct parameters *params, float (*fitnessFunction)(struct parameters *params, struct chromosome *chromo, struct data *dat), char *fitnessFunctionName){
	
	if(fitnessFunction == NULL){
		params->fitnessFunction = supervisedLearning;
		strcpy(params->fitnessFunctionName, "supervisedLearning");
	}
	else{
		params->fitnessFunction = fitnessFunction;
		strcpy(params->fitnessFunctionName, fitnessFunctionName);
	}
}





/*
	Sets the given function set to contain the per-set functions 
	given in the char array. The function names must be comma separated 
	and contain no spaces i.e. "and,or".
*/
void addNodeFunction(struct parameters *params, char *functionNames){

	char *pch;
	char funcNames[FUNCTIONNAMELENGTH];
			
	/* make a local copy of the function names*/
	strcpy(funcNames, functionNames);
			
	/* get the first function name */
	pch = strtok(funcNames, ",");
			
	/* while the function names char array contains function names */										
	while (pch != NULL){
		
		/* add the named function to the function set */
		addPresetFuctionToFunctionSet(params, pch);
		
		/* get the next function name */
		pch = strtok(NULL, ",");
	}
		
	/* if the function set is empty give warning */	
	if(params->funcSet->numFunctions == 0){
		printf("Warning: No Functions added to function set.\n");
	}	
}

/*
	used as an interface to adding pre-set node functions
*/
static void addPresetFuctionToFunctionSet(struct parameters *params, char *functionName){
	
	if(strcmp(functionName, "add") == 0){
		addNodeFunctionCustom(params, add, "add");
	}
	else if(strcmp(functionName, "sub") == 0){
		addNodeFunctionCustom(params, sub, "sub");
	}
	else if(strcmp(functionName, "mul") == 0){
		addNodeFunctionCustom(params, mul, "mul");
	}
	else if(strcmp(functionName, "div") == 0){
		addNodeFunctionCustom(params, divide, "div");
	}
	else if(strcmp(functionName, "abs") == 0){
		addNodeFunctionCustom(params, absolute, "abs");
	}
	else if(strcmp(functionName, "sqrt") == 0){
		addNodeFunctionCustom(params, squareRoot, "sqrt");
	}
	else if(strcmp(functionName, "sq") == 0){
		addNodeFunctionCustom(params, square, "sq");
	}
	else if(strcmp(functionName, "cube") == 0){
		addNodeFunctionCustom(params, cube, "cube");
	}
	else if(strcmp(functionName, "exp") == 0){
		addNodeFunctionCustom(params, exponential, "exp");
	}
	else if(strcmp(functionName, "sin") == 0){
		addNodeFunctionCustom(params, sine, "sin");
	}
	else if(strcmp(functionName, "cos") == 0){
		addNodeFunctionCustom(params, cosine, "cos");
	}
	else if(strcmp(functionName, "tan") == 0){
		addNodeFunctionCustom(params, tangent, "tan");
	}
	
		
	else if(strcmp(functionName, "and") == 0){
		addNodeFunctionCustom(params, and, "and");
	}	
	else if(strcmp(functionName, "nand") == 0){
		addNodeFunctionCustom(params, nand, "nand");
	}
	else if(strcmp(functionName, "or") == 0){
		addNodeFunctionCustom(params, or, "or");
	}	
	else if(strcmp(functionName, "nor") == 0){
		addNodeFunctionCustom(params, nor, "nor");
	}	
	else if(strcmp(functionName, "xor") == 0){
		addNodeFunctionCustom(params, xor, "xor");
	}	
	else if(strcmp(functionName, "xnor") == 0){
		addNodeFunctionCustom(params, xnor, "xnor");
	}	
	else if(strcmp(functionName, "not") == 0){
		addNodeFunctionCustom(params, not, "not");
	}	
	
	else if(strcmp(functionName, "sig") == 0){
		addNodeFunctionCustom(params, sigmoid, "sig");
	}	
	else if(strcmp(functionName, "gauss") == 0){
		addNodeFunctionCustom(params, gaussian, "gauss");
	}	
	else if(strcmp(functionName, "step") == 0){
		addNodeFunctionCustom(params, step, "step");
	}	
	else if(strcmp(functionName, "softsign") == 0){
		addNodeFunctionCustom(params, softsign, "softsign");
	}	
	else if(strcmp(functionName, "tanh") == 0){
		addNodeFunctionCustom(params, hyperbolicTangent, "tanh");
	}	
	else{
		printf("Warning: function '%s' is not known and was not added.\n", functionName);
	}	
}

/*
*/
void clearFunctionSet(struct parameters *params){
	params->funcSet->numFunctions = 0;
}

/*
	Adds given node function to given function set with given name. 
	Disallows exceeding the function set size.
*/
void addNodeFunctionCustom(struct parameters *params, float (*function)(const int numInputs, const float *inputs, const float *weights), char *functionName){
	
	if(params->funcSet->numFunctions >= FUNCTIONSETSIZE){
		printf("Warning: functions set has reached maximum capacity (%d). Function '%s' not added.\n", FUNCTIONSETSIZE, functionName);
		return;
	}
	
	/* */
	params->funcSet->numFunctions++;
	
	/* */
	strcpy(params->funcSet->functionNames[params->funcSet->numFunctions-1], functionName);
	
	/* */
	params->funcSet->functions[params->funcSet->numFunctions-1] = function;
}


/*
	Prints the current functions in the function set to
	the terminal.  
*/
void printFunctionSet(struct parameters *params){

	int i;

	printf("Functions:");
	
	for(i=0; i<params->funcSet->numFunctions; i++){
		printf(" %s", params->funcSet->functionNames[i]);
	}
	
	printf(" (%d)\n", params->funcSet->numFunctions);
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
		
	/* */
	setActiveNodes(chromo);
	
	
	chromo->nodeInputsHold = malloc(params->arity * sizeof(float));
	
	return chromo;
}

/*
	Frees the memory associated with the given chromosome structure
*/
void freeChromosome(struct chromosome *chromo){

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
static void copyChromosome(struct chromosome *chromoDest, struct chromosome *chromoSrc){

	int i;
	
	/* copy num inputs, nodes, output and arity */
	chromoDest->numInputs = chromoSrc->numInputs;
	chromoDest->numNodes = chromoSrc->numNodes;
	chromoDest->numOutputs = chromoSrc->numOutputs;
	chromoDest->arity = chromoSrc->arity;
	
	
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
void setChromosomeFitness(struct parameters *params, struct chromosome *chromo, struct data *dat){
	
	float fitness;
	
	fitness = params->fitnessFunction(params, chromo, dat);
	
	chromo->fitness = fitness;
}







/*
	mutate Random parent reproduction method.
*/
static void mutateRandomParent(struct parameters *params, struct population *pop){
	
	int i;
	
	/* for each child */
	for(i=0; i< params->lambda; i++){
		
		/* set child as clone of random parent */
		copyChromosome(pop->children[i], pop->parents[rand() % params->mu]);
		
		/* mutate newly cloned child */
		params->mutationType(params, pop->children[i]);
	}
}



/*
	Selection scheme which selects the fittest members of the population
	to be the parents. For the '+' evolutionary strategy the parents are
	selected form the current parents and children. For the ',' 
	evolutionary strategy the parents are selected form only the children 
*/
static void pickHighest(struct parameters *params, struct chromosome **parents, struct chromosome **candidateChromos, int numCandidateChromos ){
	
	int i;
			
	sortChromosomeArray(candidateChromos, numCandidateChromos);	
			
	for(i=0; i<params->mu; i++){
		
		copyChromosome(parents[i], candidateChromos[i]);
	}	
}


/* 
	Switches the first chromosome with the last and then sorts the population.
*/
static void sortChromosomeArray(struct chromosome **chromoArray, int numChromos){
	
	struct chromosome *chromoTmp;
	int i;
	int finished = 0;
	
	/* 
		place first chromosome at the end of the population.
		has the effect of always choosing new blood allowing
		for neural genetic drift to take place.
	*/
	chromoTmp = chromoArray[0];
	chromoArray[0] = chromoArray[numChromos -1];
	chromoArray[numChromos -1] = chromoTmp;
	
	/* bubble sort population */	
	while(finished == 0){
		
		finished = 1;
		
		for(i=0; i < numChromos -1; i++){
			
			if(chromoArray[i]->fitness > chromoArray[i+1]->fitness){
				
				finished = 0;
				chromoTmp = chromoArray[i];
				chromoArray[i] = chromoArray[i+1];
				chromoArray[i+1] = chromoTmp;
			}
		}
	}
}


/*
	Executes the given chromosome with the given outputs and placed the outputs in 'outputs'.
*/
void executeChromosome(struct chromosome *chromo, float *inputs, float *outputs){
	
	int i,j;
	int nodeInputLocation;
	int currentActiveNode;
	int currentActiveNodeFuction;
		
	/* for all of the active nodes */
	for(i=0; i<chromo->numActiveNodes; i++){
		
		/* get the index of the current active node */
		currentActiveNode = chromo->activeNodes[i];
		
		/* for each of the active nodes inputs */
		for(j=0; j<chromo->arity; j++){
			
			/* gather the nodes inputs */
			nodeInputLocation = chromo->nodes[currentActiveNode]->inputs[j];
			
			if(nodeInputLocation < chromo->numInputs){
				chromo->nodeInputsHold[j] = inputs[nodeInputLocation];
			}
			else{
				chromo->nodeInputsHold[j] = chromo->nodes[nodeInputLocation - chromo->numInputs]->output;
			}
		}
		
		/* get the index of the active node under evaluation */
		currentActiveNodeFuction = chromo->nodes[currentActiveNode]->function;
		
		/* calculate the output of the active node under evaluation */
		chromo->nodes[currentActiveNode]->output = chromo->funcSet->functions[currentActiveNodeFuction](chromo->arity, chromo->nodeInputsHold, chromo->nodes[currentActiveNode]->weights);
	
		  
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
	for(i=0; i<chromo->numOutputs; i++){
	
		if(chromo->outputNodes[i] < chromo->numInputs){
			outputs[i] = inputs[chromo->outputNodes[i]];
		}
		else{
			outputs[i] = chromo->nodes[chromo->outputNodes[i] - chromo->numInputs]->output;
		}
	}
}


/*
	returns a pointer to an initialised node. Initialised means that necessary
	memory has been allocated and values set.
*/
static struct node *initialiseNode(struct parameters *params, int nodePosition){
		
	struct node *n;
	int i;
	
	/* allocate memory for node */
	n = malloc(sizeof(struct node));
	
	/* allocate memory for the node's inputs and connection weights */
	n->inputs = malloc(params->arity * sizeof(int));
	n->weights = malloc(params->arity * sizeof(float));	

	/* set the node's function */
	n->function = getRandomFunction(params->funcSet->numFunctions);

	/* set as active by default */
	n->active = 1;

	/* set the nodes inputs and connection weights */
	for(i=0; i<params->arity; i++){
		n->inputs[i] = getRandomNodeInput(params->numInputs,nodePosition);
		n->weights[i] = getRandomConnectionWeight(params->connectionWeightRange);
	}
	
	/* */
	n->output = 0;
	
	/* */
	n->arity = params->arity;
	
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
	
	return rand() % (numFunctions);
}

/*
 returns a random input for the given node
*/
static int getRandomNodeInput(int numChromoInputs, int nodePosition){
	
	int input;
	
	input = rand() % (numChromoInputs + nodePosition); 
	
	return input;
}
	
/* 
	set the active nodes in the given chromosome
*/
static void setActiveNodes(struct chromosome *chromo){
	
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
	bubbleSortInt(chromo->activeNodes, chromo->numActiveNodes);
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
	returns a random chromosome output
*/	
static int getRandomChromosomeOutput(int numInputs, int numNodes){
	
	int output;
	
	output = rand() % (numInputs + numNodes);
	
	return output;
}
	
	
void printParameters(struct parameters *params){

	printf("\t-- Parameters --\n");
	
	printf("evolutionaryStrategy: (%d%c%d)-ES\n", params->mu, params->evolutionaryStrategy, params->lambda);
	printFunctionSet(params);

	printf("Mutation Type: %s\n", params->mutationTypeName);
	printf("Mutation rate: %f\n", params->mutationRate);
	printf("Connection weights range: +/- %f\n", params->connectionWeightRange);
	printf("Fitness Function: %s\n", params->fitnessFunctionName);
	printf("Target Fitness: %f\n", params->targetFitness);
	printf("Selection scheme: %s\n", params->selectionSchemeName);
	printf("Reproduction scheme: %s\n", params->reproductionSchemeName);

	printf("Update frequency: %d\n", params->updateFrequency);

	printf("------------------------------------------\n\n");
/*
	
	int numInputs;
	int numNodes;
	int numOutputs;
	int arity;	
	*/
}
	
	
/*
	Prints the given chromosome to the screen
*/	 
void printChromosome(struct chromosome *chromo){

	int i,j;			
				
	/* set the active nodes in the given chromosome */
	setActiveNodes(chromo);
		
	printf("\n");	
								
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
	
	printf("\n");
}	

/*
	Mutates the given chromosome using the mutation method described in parameters
*/
void mutateChromosome(struct parameters *params, struct chromosome *chromo){
	
	params->mutationType(params, chromo);
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
		geneToMutate = rand() % numGenesToMutate;
		
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
	float out = 1;
	
	for(i=0; i<numInputs; i++){
		
		if(inputs[i] == 0){
			out = 0;
			break;
		}
	}
	
	return out;
}	


/*
	Node function and. logical NAND, returns '0' if all inputs are '1'
	else, '1'
*/ 	
static float nand(const int numInputs, const float *inputs, const float *connectionWeights){
		
	int i;
	float out = 0;
	
	for(i=0; i<numInputs; i++){
		
		if(inputs[i] == 0){
			out = 1;
			break;
		}
	}
	
	return out;
}	


/*
	Node function or. logical OR, returns '0' if all inputs are '0'
	else, '1'
*/ 	
static float or(const int numInputs, const float *inputs, const float *connectionWeights){
		
	int i;
	float out = 0;
	
	for(i=0; i<numInputs; i++){
		
		if(inputs[i] == 1){
			out = 1;
			break;
		}
	}
	
	return out;
}


/*
	Node function nor. logical NOR, returns '1' if all inputs are '0'
	else, '0'
*/ 	
static float nor(const int numInputs, const float *inputs, const float *connectionWeights){
		
	int i;
	float out = 1;
	
	for(i=0; i<numInputs; i++){
		
		if(inputs[i] == 1){
			out = 0;
			break;
		}
	}
	
	return out;
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
void removeInactiveNodes(struct chromosome *chromo){
	
	int i,j,k;

	int originalNumNodes = chromo->numNodes;
		
	/* set the active nodes */
	setActiveNodes(chromo);
	
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
			
			
			/* */
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
	setActiveNodes(chromo);
}


/*
*/
static float supervisedLearning(struct parameters *params, struct chromosome *chromo, struct data *dat){

	int i,j;
	float error = 0;
	float temp;
		
	/* error checking */
	if(chromo->numInputs != dat->numInputs){
		printf("Error: the number of chromosome inputs must match the number of inputs specified in the data.\n");
		printf("Terminating CGP-Library.\n");
		exit(0);
	}

	if(chromo->numOutputs != dat->numOutputs){
		printf("Error: the number of chromosome outputs must match the number of outputs specified in the data.\n");
		printf("Terminating CGP-Library.\n");
		exit(0);
	}

	/* for each sample in data */
	for(i=0 ; i<dat->numSamples; i++){
	
		/* calculate the chromosome outputs for the set of inputs  */
		executeChromosome(chromo, dat->inputData[i], chromo->outputValues);
	
		/* for each chromosome output */
		for(j=0; j<chromo->numOutputs; j++){
									
			temp = (chromo->outputValues[j] - dat->outputData[i][j]);
			
			if(temp >= 0){
				error += temp;
			}
			else{
				error -= temp;
			}
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
	simple bad sort - replace with standard qsort...
*/
static void bubbleSortInt(int *array, const int length){
	
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

