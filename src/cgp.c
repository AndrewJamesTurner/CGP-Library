/*
	This file is part of CGP-Library
	Copyright (c) Andrew James Turner 2014, 2015 (andrew.turner@york.ac.uk)

    CGP-Library is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CGP-Library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with CGP-Library. If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <float.h>


#include "cgp.h"

/*
	Hard limits on the size of the function set
	and the names of various functions.
	(could make the function set size dynamic)
*/
#define FUNCTIONSETSIZE 50
#define FUNCTIONNAMELENGTH 11
#define FITNESSFUNCTIONNAMELENGTH 21
#define MUTATIONTYPENAMELENGTH 21
#define SELECTIONSCHEMENAMELENGTH 21
#define REPRODUCTIONSCHEMENAMELENGTH 21

/*
	Structure definitions
*/

struct parameters {
	int mu;
	int lambda;
	char evolutionaryStrategy;
	double mutationRate;
	double recurrentConnectionProbability;
	double connectionWeightRange;
	int numInputs;
	int numNodes;
	int numOutputs;
	int arity;
	struct functionSet *funcSet;
	double targetFitness;
	int updateFrequency;
	int shortcutConnections;
	void (*mutationType)(struct parameters *params, struct chromosome *chromo);
	char mutationTypeName[MUTATIONTYPENAMELENGTH];
	double (*fitnessFunction)(struct parameters *params, struct chromosome *chromo, struct dataSet *dat);
	char fitnessFunctionName[FITNESSFUNCTIONNAMELENGTH];
	void (*selectionScheme)(struct parameters *params, struct chromosome **parents, struct chromosome **candidateChromos, int numParents, int numCandidateChromos);
	char selectionSchemeName[SELECTIONSCHEMENAMELENGTH];
	void (*reproductionScheme)(struct parameters *params, struct chromosome **parents, struct chromosome **children, int numParents, int numChildren);
	char reproductionSchemeName[REPRODUCTIONSCHEMENAMELENGTH];
	int numThreads;
};

struct chromosome {
	int numInputs;
	int numOutputs;
	int numNodes;
	int numActiveNodes;
	int arity;
	struct node **nodes;
	int *outputNodes;
	int *activeNodes;
	double fitness;
	double *outputValues;
	struct functionSet *funcSet;
	double *nodeInputsHold;
	int generation;
};

struct node {
	int function;
	int *inputs;
	double *weights;
	int active;
	double output;
	int maxArity;
	int actArity;
};

struct functionSet {
	int numFunctions;
	char functionNames[FUNCTIONSETSIZE][FUNCTIONNAMELENGTH];
	int maxNumInputs[FUNCTIONSETSIZE];
	double (*functions[FUNCTIONSETSIZE])(const int numInputs, const double *inputs, const double *connectionWeights);
};

struct dataSet {
	int numSamples;
	int numInputs;
	int numOutputs;
	double **inputData;
	double **outputData;
};

struct results {
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
static void getBestChromosome(struct chromosome **parents, struct chromosome **children, int numParents, int numChildren, struct chromosome *best);
static void saveChromosomeLatexRecursive(struct chromosome *chromo, int index, FILE *fp);

/* node functions */
static struct node *initialiseNode(int numInputs, int numNodes, int arity, int numFunctions, double connectionWeightRange, double recurrentConnectionProbability, int nodePosition);
static void freeNode(struct node *n);
static void copyNode(struct node *nodeDest, struct node *nodeSrc);

/* getting gene value functions  */
static double getRandomConnectionWeight(double weightRange);
static int getRandomNodeInput(int numChromoInputs, int numNodes, int nodePosition, double recurrentConnectionProbability);
static int getRandomFunction(int numFunctions);
static int getRandomChromosomeOutput(int numInputs, int numNodes, int shortcutConnections);

/* function set functions */
static int addPresetFunctionToFunctionSet(struct parameters *params, char const *functionName);
static void copyFunctionSet(struct functionSet *funcSetDest, struct functionSet *funcSetSrc);
static void printFunctionSet(struct parameters *params);

/* results functions */
struct results* initialiseResults(struct parameters *params, int numRuns);

/* mutation functions  */
static void probabilisticMutation(struct parameters *params, struct chromosome *chromo);
static void pointMutation(struct parameters *params, struct chromosome *chromo);
static void pointMutationANN(struct parameters *params, struct chromosome *chromo);
static void probabilisticMutationOnlyActive(struct parameters *params, struct chromosome *chromo);
static void singleMutation(struct parameters *params, struct chromosome *chromo);

/* selection scheme functions */
static void selectFittest(struct parameters *params, struct chromosome **parents, struct chromosome **candidateChromos, int numParents, int numCandidateChromos);

/* reproduction scheme functions */
static void mutateRandomParent(struct parameters *params, struct chromosome **parents, struct chromosome **children, int numParents, int numChildren);

/* fitness function */
static double supervisedLearning(struct parameters *params, struct chromosome *chromo, struct dataSet *data);

/* node functions defines in CGP-Library */
static double _add(const int numInputs, const double *inputs, const double *connectionWeights);
static double _sub(const int numInputs, const double *inputs, const double *connectionWeights);
static double _mul(const int numInputs, const double *inputs, const double *connectionWeights);
static double _divide(const int numInputs, const double *inputs, const double *connectionWeights);
static double _and(const int numInputs, const double *inputs, const double *connectionWeights);
static double _absolute(const int numInputs, const double *inputs, const double *connectionWeights);
static double _squareRoot(const int numInputs, const double *inputs, const double *connectionWeights);
static double _square(const int numInputs, const double *inputs, const double *connectionWeights);
static double _cube(const int numInputs, const double *inputs, const double *connectionWeights);
static double _power(const int numInputs, const double *inputs, const double *connectionWeights);
static double _exponential(const int numInputs, const double *inputs, const double *connectionWeights);
static double _sine(const int numInputs, const double *inputs, const double *connectionWeights);
static double _cosine(const int numInputs, const double *inputs, const double *connectionWeights);
static double _tangent(const int numInputs, const double *inputs, const double *connectionWeights);
static double _randFloat(const int numInputs, const double *inputs, const double *connectionWeights);
static double _constOne(const int numInputs, const double *inputs, const double *connectionWeights);
static double _constZero(const int numInputs, const double *inputs, const double *connectionWeights);
static double _constPI(const int numInputs, const double *inputs, const double *connectionWeights);
static double _nand(const int numInputs, const double *inputs, const double *connectionWeights);
static double _or(const int numInputs, const double *inputs, const double *connectionWeights);
static double _nor(const int numInputs, const double *inputs, const double *connectionWeights);
static double _xor(const int numInputs, const double *inputs, const double *connectionWeights);
static double _xnor(const int numInputs, const double *inputs, const double *connectionWeights);
static double _not(const int numInputs, const double *inputs, const double *connectionWeights);
static double _wire(const int numInputs, const double *inputs, const double *connectionWeights);
static double _sigmoid(const int numInputs, const double *inputs, const double *connectionWeights);
static double _gaussian(const int numInputs, const double *inputs, const double *connectionWeights);
static double _step(const int numInputs, const double *inputs, const double *connectionWeights);
static double _softsign(const int numInputs, const double *inputs, const double *connectionWeights);
static double _hyperbolicTangent(const int numInputs, const double *inputs, const double *connectionWeights);

/* other */
static double randDecimal(void);
static int randInt(int n);
static double sumWeigtedInputs(const int numInputs, const double *inputs, const double *connectionWeights);
static void sortIntArray(int *array, const int length);
static void sortDoubleArray(double *array, const int length);
static int cmpInt(const void * a, const void * b);
static int cmpDouble(const void * a, const void * b);
static double medianInt(const int *anArray, const int length);
static double medianDouble(const double *anArray, const int length);


/*
	parameters function definitions
*/

/*
	Initialises a parameter struct with default values. These
	values can be individually changed via set functions.
*/
DLL_EXPORT struct parameters *initialiseParameters(const int numInputs, const int numNodes, const int numOutputs, const int arity) {

	struct parameters *params;

	/* allocate memory for parameters */
	params = (struct parameters*)malloc(sizeof(struct parameters));

	/* Set default values */
	params->mu = 1;
	params->lambda = 4;
	params->evolutionaryStrategy = '+';
	params->mutationRate = 0.05;
	params->recurrentConnectionProbability = 0.0;
	params->connectionWeightRange = 1;
	params->shortcutConnections = 1;

	params->targetFitness = 0;

	params->updateFrequency = 1;

	setNumInputs(params, numInputs);
	setNumNodes(params, numNodes);
	setNumOutputs(params, numOutputs);
	setArity(params, arity);

	params->mutationType = probabilisticMutation;
	strncpy(params->mutationTypeName, "probabilistic", MUTATIONTYPENAMELENGTH);

	params->funcSet = (struct functionSet*)malloc(sizeof(struct functionSet));
	params->funcSet->numFunctions = 0;

	params->fitnessFunction = supervisedLearning;
	strncpy(params->fitnessFunctionName, "supervisedLearning", FITNESSFUNCTIONNAMELENGTH);

	params->selectionScheme = selectFittest;
	strncpy(params->selectionSchemeName, "selectFittest", SELECTIONSCHEMENAMELENGTH);

	params->reproductionScheme = mutateRandomParent;
	strncpy(params->reproductionSchemeName, "mutateRandomParent", REPRODUCTIONSCHEMENAMELENGTH);

	params->numThreads = 1;

	/* Seed the random number generator */
	srand(time(NULL));

	return params;
}


/*
	Frees the memory associated with the given parameter structure
*/
DLL_EXPORT void freeParameters(struct parameters *params) {

	/* attempt to prevent user double freeing */
	if (params == NULL) {
		printf("Warning: double freeing of parameters prevented.\n");
		return;
	}

	free(params->funcSet);
	free(params);
}

/*
	prints the given parameters to the terminal
*/
DLL_EXPORT void printParameters(struct parameters *params) {

	if (params == NULL) {
		printf("Error: cannot print uninitialised parameters.\nTerminating CGP-Library.\n");
		exit(0);
	}

	printf("-----------------------------------------------------------\n");
	printf("                       Parameters                          \n");
	printf("-----------------------------------------------------------\n");
	printf("Evolutionary Strategy:\t\t\t(%d%c%d)-ES\n", params->mu, params->evolutionaryStrategy, params->lambda);
	printf("Inputs:\t\t\t\t\t%d\n", params->numInputs);
	printf("Nodes:\t\t\t\t\t%d\n", params->numNodes);
	printf("Outputs:\t\t\t\t%d\n", params->numOutputs);
	printf("Node Arity:\t\t\t\t%d\n", params->arity);
	printf("Connection weights range:\t\t+/- %f\n", params->connectionWeightRange);
	printf("Mutation Type:\t\t\t\t%s\n", params->mutationTypeName);
	printf("Mutation rate:\t\t\t\t%f\n", params->mutationRate);
	printf("Recurrent Connection Probability:\t%f\n", params->recurrentConnectionProbability);
	printf("Shortcut Connections:\t\t\t%d\n", params->shortcutConnections);
	printf("Fitness Function:\t\t\t%s\n", params->fitnessFunctionName);
	printf("Target Fitness:\t\t\t\t%f\n", params->targetFitness);
	printf("Selection scheme:\t\t\t%s\n", params->selectionSchemeName);
	printf("Reproduction scheme:\t\t\t%s\n", params->reproductionSchemeName);
	printf("Update frequency:\t\t\t%d\n", params->updateFrequency);
	printf("Threads:\t\t\t%d\n", params->numThreads);
	printFunctionSet(params);
	printf("-----------------------------------------------------------\n\n");
}


/*
	Adds the give pre-defined functions to the given function set. The
	functions must be given in the char array. The function names must
	be comma separated and contain no spaces i.e. "and,or".
*/
DLL_EXPORT void addNodeFunction(struct parameters *params, char const *functionNames) {

	char *pch;
	char functionNamesAsArray[FUNCTIONNAMELENGTH * FUNCTIONSETSIZE];

	/* make a local copy of the function names*/
	strncpy(functionNamesAsArray, functionNames, FUNCTIONNAMELENGTH * FUNCTIONSETSIZE);

	/* get the first function name */
	pch = strtok(functionNamesAsArray, ", ");

	/* while the function names char array contains function names */
	while (pch != NULL) {

		/* add the named function to the function set */
		addPresetFunctionToFunctionSet(params, pch);

		/* get the next function name */
		pch = strtok(NULL, ", ");
	}

	/* if the function set is empty give warning */
	if (params->funcSet->numFunctions == 0) {
		printf("Warning: No Functions added to function set.\n");
	}
}


/*
	Adds given node function to given function set with given name.
	Disallows exceeding the function set size.
*/
DLL_EXPORT void addCustomNodeFunction(struct parameters *params, double (*function)(const int numInputs, const double *inputs, const double *weights), char const *functionName, int maxNumInputs) {

	if (params->funcSet->numFunctions >= FUNCTIONSETSIZE) {
		printf("Warning: functions set has reached maximum capacity (%d). Function '%s' not added.\n", FUNCTIONSETSIZE, functionName);
		return;
	}

	/* set the function name as the given function name */
	strncpy(params->funcSet->functionNames[params->funcSet->numFunctions], functionName, FUNCTIONNAMELENGTH);

	/* set the number of function inputs as the given number of function inputs */
	params->funcSet->maxNumInputs[params->funcSet->numFunctions] = maxNumInputs;

	/* add the given function to the function set */
	params->funcSet->functions[params->funcSet->numFunctions] = function;

	params->funcSet->numFunctions++;
}


/*
	used as an interface to adding pre-set node functions.
	returns one if successful, zero otherwise.
*/
static int addPresetFunctionToFunctionSet(struct parameters *params, char const *functionName) {

	int successfullyAdded = 1;

	/* Symbolic functions */

	if (strncmp(functionName, "add", FUNCTIONNAMELENGTH) == 0) {
		addCustomNodeFunction(params, _add, "add", -1);
	}
	else if (strncmp(functionName, "sub", FUNCTIONNAMELENGTH) == 0) {
		addCustomNodeFunction(params, _sub, "sub", -1);
	}
	else if (strncmp(functionName, "mul", FUNCTIONNAMELENGTH) == 0) {
		addCustomNodeFunction(params, _mul, "mul", -1);
	}
	else if (strncmp(functionName, "div", FUNCTIONNAMELENGTH) == 0) {
		addCustomNodeFunction(params, _divide, "div", -1);
	}
	else if (strncmp(functionName, "abs", FUNCTIONNAMELENGTH) == 0) {
		addCustomNodeFunction(params, _absolute, "abs", 1);
	}
	else if (strncmp(functionName, "sqrt", FUNCTIONNAMELENGTH) == 0) {
		addCustomNodeFunction(params, _squareRoot, "sqrt", 1);
	}
	else if (strncmp(functionName, "sq", FUNCTIONNAMELENGTH) == 0) {
		addCustomNodeFunction(params, _square, "sq", 1);
	}
	else if (strncmp(functionName, "cube", FUNCTIONNAMELENGTH) == 0) {
		addCustomNodeFunction(params, _cube, "cube", 1);
	}
	else if (strncmp(functionName, "pow", FUNCTIONNAMELENGTH) == 0) {
		addCustomNodeFunction(params, _power, "pow", 2);
	}
	else if (strncmp(functionName, "exp", FUNCTIONNAMELENGTH) == 0) {
		addCustomNodeFunction(params, _exponential, "exp", 1);
	}
	else if (strncmp(functionName, "sin", FUNCTIONNAMELENGTH) == 0) {
		addCustomNodeFunction(params, _sine, "sin", 1);
	}
	else if (strncmp(functionName, "cos", FUNCTIONNAMELENGTH) == 0) {
		addCustomNodeFunction(params, _cosine, "cos", 1);
	}
	else if (strncmp(functionName, "tan", FUNCTIONNAMELENGTH) == 0) {
		addCustomNodeFunction(params, _tangent, "tan", 1);
	}


	/* Boolean logic gates */

	else if (strncmp(functionName, "and", FUNCTIONNAMELENGTH) == 0) {
		addCustomNodeFunction(params, _and, "and", -1);
	}
	else if (strncmp(functionName, "nand", FUNCTIONNAMELENGTH) == 0) {
		addCustomNodeFunction(params, _nand, "nand", -1);
	}
	else if (strncmp(functionName, "or", FUNCTIONNAMELENGTH) == 0) {
		addCustomNodeFunction(params, _or, "or", -1);
	}
	else if (strncmp(functionName, "nor", FUNCTIONNAMELENGTH) == 0) {
		addCustomNodeFunction(params, _nor, "nor", -1);
	}
	else if (strncmp(functionName, "xor", FUNCTIONNAMELENGTH) == 0) {
		addCustomNodeFunction(params, _xor, "xor", -1);
	}
	else if (strncmp(functionName, "xnor", FUNCTIONNAMELENGTH) == 0) {
		addCustomNodeFunction(params, _xnor, "xnor", -1);
	}
	else if (strncmp(functionName, "not", FUNCTIONNAMELENGTH) == 0) {
		addCustomNodeFunction(params, _not, "not", 1);
	}

	/* Neuron functions */

	else if (strncmp(functionName, "sig", FUNCTIONNAMELENGTH) == 0) {
		addCustomNodeFunction(params, _sigmoid, "sig", -1);
	}
	else if (strncmp(functionName, "gauss", FUNCTIONNAMELENGTH) == 0) {
		addCustomNodeFunction(params, _gaussian, "gauss", -1);
	}
	else if (strncmp(functionName, "step", FUNCTIONNAMELENGTH) == 0) {
		addCustomNodeFunction(params, _step, "step", -1);
	}
	else if (strncmp(functionName, "softsign", FUNCTIONNAMELENGTH) == 0) {
		addCustomNodeFunction(params, _softsign, "soft", -1);
	}
	else if (strncmp(functionName, "tanh", FUNCTIONNAMELENGTH) == 0) {
		addCustomNodeFunction(params, _hyperbolicTangent, "tanh", -1);
	}

	/* other */

	else if (strncmp(functionName, "rand", FUNCTIONNAMELENGTH) == 0) {
		addCustomNodeFunction(params, _randFloat, "rand", 0);
	}
	else if (strncmp(functionName, "1", FUNCTIONNAMELENGTH) == 0) {
		addCustomNodeFunction(params, _constOne, "1", 0);
	}
	else if (strncmp(functionName, "0", FUNCTIONNAMELENGTH) == 0) {
		addCustomNodeFunction(params, _constZero, "0", 0);
	}
	else if (strncmp(functionName, "pi", FUNCTIONNAMELENGTH) == 0) {
		addCustomNodeFunction(params, _constPI, "pi", 0);
	}
	else if (strncmp(functionName, "wire", FUNCTIONNAMELENGTH) == 0) {
		addCustomNodeFunction(params, _wire, "wire", 1);
	}

	else {
		printf("Warning: function '%s' is not known and was not added.\n", functionName);
		successfullyAdded = 0;
	}

	return successfullyAdded;
}


/*
	clears the given function set of functions
*/
DLL_EXPORT void clearFunctionSet(struct parameters *params) {
	params->funcSet->numFunctions = 0;
}


/*
	sets num chromosome inputs in parameters
*/
DLL_EXPORT void setNumInputs(struct parameters *params, int numInputs) {

	/* error checking */
	if (numInputs <= 0) {
		printf("Error: number of chromosome inputs cannot be less than one; %d is invalid.\nTerminating CGP-Library.\n", numInputs);
		exit(0);
	}

	params->numInputs = numInputs;
}


/*
	sets num chromosome nodes in parameters
*/
DLL_EXPORT void setNumNodes(struct parameters *params, int numNodes) {

	/* error checking */
	if (numNodes < 0) {
		printf("Warning: number of chromosome nodes cannot be negative; %d is invalid.\nTerminating CGP-Library.\n", numNodes);
		exit(0);
	}

	params->numNodes = numNodes;
}


/*
	sets num chromosome outputs in parameters
*/
DLL_EXPORT void setNumOutputs(struct parameters *params, int numOutputs) {

	/* error checking */
	if (numOutputs < 0) {
		printf("Warning: number of chromosome outputs cannot be less than one; %d is invalid.\nTerminating CGP-Library.\n", numOutputs);
		exit(0);
	}

	params->numOutputs = numOutputs;
}


/*
	sets chromosome arity in parameters
*/
DLL_EXPORT void setArity(struct parameters *params, int arity) {

	/* error checking */
	if (arity < 0) {
		printf("Warning: node arity cannot be less than one; %d is invalid.\nTerminating CGP-Library.\n", arity);
		exit(0);
	}

	params->arity = arity;
}


/*
	Sets the mu value in given parameters to the new given value. If mu value
	is invalid a warning is displayed and the mu value is left unchanged.
*/
DLL_EXPORT void setMu(struct parameters *params, int mu) {

	if (mu > 0) {
		params->mu = mu;
	}
	else {
		printf("\nWarning: mu value '%d' is invalid. Mu value must have a value of one or greater. Mu value left unchanged as '%d'.\n", mu, params->mu);
	}
}


/*
	Sets the lambda value in given parameters to the new given value.
	If lambda value is invalid a warning is displayed and the lambda value
	is left unchanged.
*/
DLL_EXPORT void setLambda(struct parameters *params, int lambda) {

	if (lambda > 0) {
		params->lambda = lambda;
	}
	else {
		printf("\nWarning: lambda value '%d' is invalid. Lambda value must have a value of one or greater. Lambda value left unchanged as '%d'.\n", lambda, params->lambda);
	}
}


/*
	Sets the evolutionary strategy given in parameters to '+' or ','.
	If an invalid option is given a warning is displayed and the evolutionary
	strategy is left unchanged.
*/
DLL_EXPORT void setEvolutionaryStrategy(struct parameters *params, char evolutionaryStrategy) {

	if (evolutionaryStrategy == '+' || evolutionaryStrategy == ',') {
		params->evolutionaryStrategy = evolutionaryStrategy;
	}
	else {
		printf("\nWarning: the evolutionary strategy '%c' is invalid. The evolutionary strategy must be '+' or ','. The evolutionary strategy has been left unchanged as '%c'.\n", evolutionaryStrategy, params->evolutionaryStrategy);
	}
}


/*
	Sets the mutation rate given in parameters. If an invalid mutation
	rate is given a warning is displayed and the mutation rate is left
	unchanged.
*/
DLL_EXPORT void setMutationRate(struct parameters *params, double mutationRate) {

	if (mutationRate >= 0 && mutationRate <= 1) {
		params->mutationRate = mutationRate;
	}
	else {
		printf("\nWarning: mutation rate '%f' is invalid. The mutation rate must be in the range [0,1]. The mutation rate has been left unchanged as '%f'.\n", mutationRate, params->mutationRate);
	}
}


/*
	Sets the recurrent connection probability given in parameters. If an invalid
	value is given a warning is displayed and the value is left	unchanged.
*/
DLL_EXPORT void setRecurrentConnectionProbability(struct parameters *params, double recurrentConnectionProbability) {

	if (recurrentConnectionProbability >= 0 && recurrentConnectionProbability <= 1) {
		params->recurrentConnectionProbability = recurrentConnectionProbability;
	}
	else {
		printf("\nWarning: recurrent connection probability '%f' is invalid. The recurrent connection probability must be in the range [0,1]. The recurrent connection probability has been left unchanged as '%f'.\n", recurrentConnectionProbability, params->recurrentConnectionProbability);
	}
}


/*
	Sets the whether shortcut connections are used. If an invalid
	value is given a warning is displayed and the value is left	unchanged.
*/
DLL_EXPORT void setShortcutConnections(struct parameters *params, int shortcutConnections) {

	if (shortcutConnections == 0 || shortcutConnections == 1) {
		params->shortcutConnections = shortcutConnections;
	}
	else {
		printf("\nWarning: shortcut connection '%d' is invalid. The shortcut connections takes values 0 or 1. The shortcut connection has been left unchanged as '%d'.\n", shortcutConnections, params->shortcutConnections);
	}
}


/*
	Sets the connection weight range given in parameters.
*/
DLL_EXPORT void setConnectionWeightRange(struct parameters *params, double weightRange) {

	params->connectionWeightRange = weightRange;
}


/*
	sets the fitness function to the fitnessFunction passed. If the fitnessFunction is NULL
	then the default supervisedLearning fitness function is used.
*/
DLL_EXPORT void setCustomFitnessFunction(struct parameters *params, double (*fitnessFunction)(struct parameters *params, struct chromosome *chromo, struct dataSet *data), char const *fitnessFunctionName) {

	if (fitnessFunction == NULL) {
		params->fitnessFunction = supervisedLearning;
		strncpy(params->fitnessFunctionName, "supervisedLearning", FITNESSFUNCTIONNAMELENGTH);
	}
	else {
		params->fitnessFunction = fitnessFunction;
		strncpy(params->fitnessFunctionName, fitnessFunctionName, FITNESSFUNCTIONNAMELENGTH);
	}
}



/*
	sets the selection scheme used to select the parents from the candidate chromosomes. If the selectionScheme is NULL
	then the default selectFittest selection scheme is used.
*/
DLL_EXPORT void setCustomSelectionScheme(struct parameters *params, void (*selectionScheme)(struct parameters *params, struct chromosome **parents, struct chromosome **candidateChromos, int numParents, int numCandidateChromos), char const *selectionSchemeName) {

	if (selectionScheme == NULL) {
		params->selectionScheme = selectFittest;
		strncpy(params->selectionSchemeName, "selectFittest", SELECTIONSCHEMENAMELENGTH);
	}
	else {
		params->selectionScheme = selectionScheme;
		strncpy(params->selectionSchemeName, selectionSchemeName, SELECTIONSCHEMENAMELENGTH);
	}
}


/*
	sets the reproduction scheme used to select the parents from the candidate chromosomes. If the reproductionScheme is NULL
	then the default mutateRandomParent selection scheme is used.
*/

DLL_EXPORT void setCustomReproductionScheme(struct parameters *params, void (*reproductionScheme)(struct parameters *params, struct chromosome **parents, struct chromosome **children, int numParents, int numChildren), char const *reproductionSchemeName) {

	if (reproductionScheme == NULL) {
		params->reproductionScheme = mutateRandomParent;
		strncpy(params->reproductionSchemeName, "mutateRandomParent", REPRODUCTIONSCHEMENAMELENGTH);
	}
	else {
		params->reproductionScheme = reproductionScheme;
		strncpy(params->reproductionSchemeName, reproductionSchemeName, REPRODUCTIONSCHEMENAMELENGTH);
	}
}


/*
	Sets the target fitness
*/
DLL_EXPORT void setTargetFitness(struct parameters *params, double targetFitness) {

	params->targetFitness = targetFitness;
}


/*
	sets the mutation type in params
*/
DLL_EXPORT void setMutationType(struct parameters *params, char const *mutationType) {

	if (strncmp(mutationType, "probabilistic", MUTATIONTYPENAMELENGTH) == 0) {

		params->mutationType = probabilisticMutation;
		strncpy(params->mutationTypeName, "probabilistic", MUTATIONTYPENAMELENGTH);
	}

	else if (strncmp(mutationType, "point", MUTATIONTYPENAMELENGTH) == 0) {

		params->mutationType = pointMutation;
		strncpy(params->mutationTypeName, "point", MUTATIONTYPENAMELENGTH);
	}

	else if (strncmp(mutationType, "pointANN", MUTATIONTYPENAMELENGTH) == 0) {

		params->mutationType = pointMutationANN;
		strncpy(params->mutationTypeName, "pointANN", MUTATIONTYPENAMELENGTH);
	}

	else if (strncmp(mutationType, "onlyActive", MUTATIONTYPENAMELENGTH) == 0) {

		params->mutationType = probabilisticMutationOnlyActive;
		strncpy(params->mutationTypeName, "onlyActive", MUTATIONTYPENAMELENGTH);
	}

	else if (strncmp(mutationType, "single", MUTATIONTYPENAMELENGTH) == 0) {

		params->mutationType = singleMutation;
		strncpy(params->mutationTypeName, "single", MUTATIONTYPENAMELENGTH);
	}

	else {
		printf("\nWarning: mutation type '%s' is invalid. The mutation type must be 'probabilistic' or 'point'. The mutation type has been left unchanged as '%s'.\n", mutationType, params->mutationTypeName);
	}
}


/*
	Sets the update frequency in generations
*/
DLL_EXPORT void setUpdateFrequency(struct parameters *params, int updateFrequency) {

	if (updateFrequency < 0) {
		printf("Warning: update frequency of %d is invalid. Update frequency must be >= 0. Update frequency is left unchanged as %d.\n", updateFrequency, params->updateFrequency);
	}
	else {
		params->updateFrequency = updateFrequency;
	}
}




/*
	sets num chromosome inputs in parameters
*/
DLL_EXPORT void setNumThreads(struct parameters *params, int numThreads) {

	/* error checking */
	if (numThreads <= 0) {
		printf("Error: number threads cannot be less than one; %d is invalid. The number threads is left unchanged as %d.\n", numThreads, numThreads);
	}

	params->numThreads = numThreads;
}


/*
	chromosome function definitions
*/


/*
	Returns a pointer to an initialised chromosome with values obeying the given parameters.
*/
DLL_EXPORT struct chromosome *initialiseChromosome(struct parameters *params) {

	struct chromosome *chromo;
	int i;

	/* check that funcSet contains functions */
	if (params->funcSet->numFunctions < 1) {
		printf("Error: chromosome not initialised due to empty functionSet.\nTerminating CGP-Library.\n");
		exit(0);
	}

	/* allocate memory for chromosome */
	chromo = (struct chromosome*)malloc(sizeof(struct chromosome));

	/* allocate memory for nodes */
	chromo->nodes = (struct node**)malloc(params->numNodes * sizeof(struct node*));

	/* allocate memory for outputNodes matrix */
	chromo->outputNodes = (int*)malloc(params->numOutputs * sizeof(int));

	/* allocate memory for active nodes matrix */
	chromo->activeNodes = (int*)malloc(params->numNodes * sizeof(int));

	/* allocate memory for chromosome outputValues */
	chromo->outputValues = (double*)malloc(params->numOutputs * sizeof(double));

	/* Initialise each of the chromosomes nodes */
	for (i = 0; i < params->numNodes; i++) {
		chromo->nodes[i] = initialiseNode(params->numInputs, params->numNodes, params->arity, params->funcSet->numFunctions, params->connectionWeightRange, params->recurrentConnectionProbability, i);
	}

	/* set each of the chromosomes outputs */
	for (i = 0; i < params->numOutputs; i++) {
		chromo->outputNodes[i] = getRandomChromosomeOutput(params->numInputs, params->numNodes, params->shortcutConnections);
	}

	/* set the number of inputs, nodes and outputs */
	chromo->numInputs = params->numInputs;
	chromo->numNodes = params->numNodes;
	chromo->numOutputs = params->numOutputs;
	chromo->arity = params->arity;

	/* set the number of active node to the number of nodes (all active) */
	chromo->numActiveNodes = params->numNodes;

	/* set the fitness to initial value */
	chromo->fitness = -1;

	/* copy the function set from the parameters to the chromosome */
	chromo->funcSet = (struct functionSet*)malloc(sizeof(struct functionSet));
	copyFunctionSet(chromo->funcSet, params->funcSet);

	/* set the active nodes in the newly generated chromosome */
	setChromosomeActiveNodes(chromo);

	/* used interally when exicuting chromosome */
	chromo->nodeInputsHold = (double*)malloc(params->arity * sizeof(double));

	return chromo;
}


/*
	Reads in saved chromosomes
*/
DLL_EXPORT struct chromosome* initialiseChromosomeFromFile(char const *file) {

	int i, j;

	FILE *fp;
	struct chromosome *chromo;
	struct parameters *params;

	char *line, *record;
	char funcName[FUNCTIONNAMELENGTH];
	char buffer[1024];

	int numInputs, numNodes, numOutputs, arity;

	/* open the chromosome file */
	fp = fopen(file, "r");

	/* ensure that the file was opened correctly */
	if (fp == NULL) {
		printf("Warning: cannot open chromosome: '%s'. Chromosome was not open.\n", file);
		return NULL;
	}

	/* get num inputs */
	line = fgets(buffer, sizeof(buffer), fp);
	if (line == NULL) {/*error*/}
	record = strtok(line, ",");
	record = strtok(NULL, ",");
	numInputs = atoi(record);

	/* get num nodes */
	line = fgets(buffer, sizeof(buffer), fp);
	if (line == NULL) {/*error*/}
	record = strtok(line, ",");
	record = strtok(NULL, ",");
	numNodes = atoi(record);

	/* get num outputs */
	line = fgets(buffer, sizeof(buffer), fp);
	if (line == NULL) {/*error*/}
	record = strtok(line, ",");
	record = strtok(NULL, ",");
	numOutputs = atoi(record);

	/* get arity */
	line = fgets(buffer, sizeof(buffer), fp);
	if (line == NULL) {/*error*/}
	record = strtok(line, ",");
	record = strtok(NULL, ",");
	arity = atoi(record);

	/* initialise parameters  */
	params = initialiseParameters(numInputs, numNodes, numOutputs, arity);

	/* get and set node functions */
	line = fgets(buffer, sizeof(buffer), fp);
	if (line == NULL) {/*error*/}
	record = strtok(line, ",\n");
	record = strtok(NULL, ",\n");

	/* for each function name */
	while ( record != NULL) {

		strncpy(funcName, record, FUNCTIONNAMELENGTH);

		/* can only load functions defined within CGP-Library */
		if (addPresetFunctionToFunctionSet(params, funcName) == 0) {
			printf("Error: cannot load chromosome which contains custom node functions.\n");
			printf("Terminating CGP-Library.\n");
			freeParameters(params);
			exit(0);
		}

		record = strtok(NULL, ",\n");
	}

	/* initialise a chromosome beased on the parameters accociated with given chromosome */
	chromo = initialiseChromosome(params);

	/* set the node parameters */
	for (i = 0; i < numNodes; i++) {

		/* get the function gene */
		line = fgets(buffer, sizeof(buffer), fp);
		record = strtok(line, ",\n");
		chromo->nodes[i]->function = atoi(record);

		for (j = 0; j < arity; j++) {
			line = fgets(buffer, sizeof(buffer), fp);
			sscanf(line, "%d,%lf", &chromo->nodes[i]->inputs[j], &chromo->nodes[i]->weights[j]);
		}
	}

	/* set the outputs */
	line = fgets(buffer, sizeof(buffer), fp);
	record = strtok(line, ",\n");
	chromo->outputNodes[0] = atoi(record);

	for (i = 1; i < numOutputs; i++) {
		record = strtok(NULL, ",\n");
		chromo->outputNodes[i] = atoi(record);
	}

	fclose(fp);
	freeParameters(params);

	/* set the active nodes in the copied chromosome */
	setChromosomeActiveNodes(chromo);

	return chromo;
}


/*
	Returns a pointer to an initialised chromosome with values obeying the given parameters.
*/
DLL_EXPORT struct chromosome *initialiseChromosomeFromChromosome(struct chromosome *chromo) {

	struct chromosome *chromoNew;
	int i;

	/* check that funcSet contains functions*/
	if (chromo == NULL) {
		printf("Error: cannot initialise chromosome from uninitialised chromosome.\nTerminating CGP-Library.\n");
		exit(0);
	}

	/* allocate memory for chromosome */
	chromoNew = (struct chromosome*)malloc(sizeof(struct chromosome));

	/* allocate memory for nodes */
	chromoNew->nodes = (struct node**)malloc(chromo->numNodes * sizeof(struct node*));

	/* allocate memory for outputNodes matrix */
	chromoNew->outputNodes = (int*)malloc(chromo->numOutputs * sizeof(int));

	/* allocate memory for active nodes matrix */
	chromoNew->activeNodes = (int*)malloc(chromo->numNodes * sizeof(int));

	/* allocate memory for chromosome outputValues */
	chromoNew->outputValues = (double*)malloc(chromo->numOutputs * sizeof(double));

	/* Initialise each of the chromosomes nodes */
	for (i = 0; i < chromo->numNodes; i++) {
		chromoNew->nodes[i] = initialiseNode(chromo->numInputs, chromo->numNodes, chromo->arity, chromo->funcSet->numFunctions, 0, 0, i);
		copyNode(chromoNew->nodes[i], chromo->nodes[i]);
	}

	/* set each of the chromosomes outputs */
	for (i = 0; i < chromo->numOutputs; i++) {
		chromoNew->outputNodes[i] = chromo->outputNodes[i];
	}

	/* set the number of inputs, nodes and outputs */
	chromoNew->numInputs = chromo->numInputs;
	chromoNew->numNodes = chromo->numNodes;
	chromoNew->numOutputs = chromo->numOutputs;
	chromoNew->arity = chromo->arity;


	/* copy over the chromsosme fitness */
	chromoNew->fitness = chromo->fitness;

	/* copy over the number of gnerations to find a solution */
	chromoNew->generation = chromo->generation;

	/* copy over the functionset */
	chromoNew->funcSet = (struct functionSet*)malloc(sizeof(struct functionSet));
	copyFunctionSet(chromoNew->funcSet, chromo->funcSet);

	/* set the active nodes in the newly generated chromosome */
	setChromosomeActiveNodes(chromoNew);

	/* used internally by exicute chromosome */
	chromoNew->nodeInputsHold = (double*)malloc(chromo->arity * sizeof(double));

	return chromoNew;
}


/*
	Frees the memory associated with the given chromosome structure
*/
DLL_EXPORT void freeChromosome(struct chromosome *chromo) {

	int i;

	/* attempt to prevent user double freeing */
	if (chromo == NULL) {
		printf("Warning: double freeing of chromosome prevented.\n");
		return;
	}

	for (i = 0; i < chromo->numNodes; i++) {
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
	Prints the given chromosome to the screen
*/
DLL_EXPORT void printChromosome(struct chromosome *chromo, int weights) {

	int i, j;

	/* error checking */
	if (chromo == NULL) {
		printf("Error: chromosome has not been initialised and cannot be printed.\n");
		return;
	}

	/* set the active nodes in the given chromosome */
	setChromosomeActiveNodes(chromo);

	/* for all the chromo inputs*/
	for (i = 0; i < chromo->numInputs; i++) {
		printf("(%d):\tinput\n", i);
	}

	/* for all the hidden nodes */
	for (i = 0; i < chromo->numNodes; i++) {

		/* print the node function */
		printf("(%d):\t%s\t", chromo->numInputs + i, chromo->funcSet->functionNames[chromo->nodes[i]->function]);

		/* for the arity of the node */
		for (j = 0; j < getChromosomeNodeArity(chromo, i); j++) {

			/* print the node input information */
			if (weights == 1) {
				printf("%d,%+.1f\t", chromo->nodes[i]->inputs[j], chromo->nodes[i]->weights[j]);
			}
			else {
				printf("%d ", chromo->nodes[i]->inputs[j]);
			}
		}

		/* Highlight active nodes */
		if (chromo->nodes[i]->active == 1) {
			printf("*");
		}

		printf("\n");
	}

	/* for all of the outputs */
	printf("outputs: ");
	for (i = 0; i < chromo->numOutputs; i++) {

		/* print the output node locations */
		printf("%d ", chromo->outputNodes[i]);
	}

	printf("\n\n");
}



/*
	Executes the given chromosome
*/
DLL_EXPORT void executeChromosome(struct chromosome *chromo, const double *inputs) {

	int i, j;
	int nodeInputLocation;
	int currentActiveNode;
	int currentActiveNodeFunction;
	int nodeArity;

	const int numInputs = chromo->numInputs;
	const int numActiveNodes = chromo->numActiveNodes;
	const int numOutputs = chromo->numOutputs;

	/* error checking */
	if (chromo == NULL) {
		printf("Error: cannot execute uninitialised chromosome.\n Terminating CGP-Library.\n");
		exit(0);
	}

	/* for all of the active nodes */
	for (i = 0; i < numActiveNodes; i++) {

		/* get the index of the current active node */
		currentActiveNode = chromo->activeNodes[i];

		/* get the arity of the current node */
		nodeArity = chromo->nodes[currentActiveNode]->actArity;

		/* for each of the active nodes inputs */
		for (j = 0; j < nodeArity; j++) {

			/* gather the nodes input locations */
			nodeInputLocation = chromo->nodes[currentActiveNode]->inputs[j];

			if (nodeInputLocation < numInputs) {
				chromo->nodeInputsHold[j] = inputs[nodeInputLocation];
			}
			else {
				chromo->nodeInputsHold[j] = chromo->nodes[nodeInputLocation - numInputs]->output;
			}
		}

		/* get the functionality of the active node under evaluation */
		currentActiveNodeFunction = chromo->nodes[currentActiveNode]->function;

		/* calculate the output of the active node under evaluation */
		chromo->nodes[currentActiveNode]->output = chromo->funcSet->functions[currentActiveNodeFunction](nodeArity, chromo->nodeInputsHold, chromo->nodes[currentActiveNode]->weights);


		/* deal with doubles becoming NAN */
		if (isnan(chromo->nodes[currentActiveNode]->output) != 0) {
			chromo->nodes[currentActiveNode]->output = 0;
		}

		/* prevent double form going to inf and -inf */
		else if (isinf(chromo->nodes[currentActiveNode]->output) != 0 ) {

			if (chromo->nodes[currentActiveNode]->output > 0) {
				chromo->nodes[currentActiveNode]->output = DBL_MAX;
			}
			else {
				chromo->nodes[currentActiveNode]->output = DBL_MIN;
			}
		}
	}

	/* Set the chromosome outputs */
	for (i = 0; i < numOutputs; i++) {

		if (chromo->outputNodes[i] < numInputs) {
			chromo->outputValues[i] = inputs[chromo->outputNodes[i]];
		}
		else {
			chromo->outputValues[i] = chromo->nodes[chromo->outputNodes[i] - numInputs]->output;
		}
	}
}

/*
	used to access the chromosome outputs after executeChromosome
	has been called
*/
DLL_EXPORT double getChromosomeOutput(struct chromosome *chromo, int output) {

	if (output < 0 || output > chromo->numOutputs) {
		printf("Error: output less than or greater than the number of chromosome outputs. Called from getChromosomeOutput.\n");
		exit(0);
	}

	return chromo->outputValues[output];
}



/*
	used to access the chromosome node values after executeChromosome
	has been called
*/
DLL_EXPORT double getChromosomeNodeValue(struct chromosome *chromo, int node) {
	if (node < 0 || node > chromo->numNodes) {
		printf("Error: node less than or greater than the number of nodes  in chromosome. Called from getChromosomeNodeValue.\n");
		exit(0);
	}

	return chromo->nodes[node]->output;
}


/*
	returns whether the specified node is active in the given chromosome
*/
DLL_EXPORT int isNodeActive(struct chromosome *chromo, int node) {

	if (node < 0 || node > chromo->numNodes) {
		printf("Error: node less than or greater than the number of nodes  in chromosome. Called from isNodeActive.\n");
		exit(0);
	}

	return chromo->nodes[node]->active;
}


/*
	Saves the given chromosome in a form which can be read in later
*/
DLL_EXPORT void saveChromosome(struct chromosome *chromo, char const *fileName) {

	int i, j;
	FILE *fp;

	/* create the chromosome file */
	fp = fopen(fileName, "w");

	/* ensure that the file was created correctly */
	if (fp == NULL) {
		printf("Warning: cannot save chromosome to '%s'. Chromosome was not saved.\n", fileName);
		return;
	}

	/* save meta information */
	fprintf(fp, "numInputs,%d\n", chromo->numInputs);
	fprintf(fp, "numNodes,%d\n", chromo->numNodes);
	fprintf(fp, "numOutputs,%d\n", chromo->numOutputs);
	fprintf(fp, "arity,%d\n", chromo->arity);

	fprintf(fp, "functionSet");

	for (i = 0; i < chromo->funcSet->numFunctions; i++) {
		fprintf(fp, ",%s", chromo->funcSet->functionNames[i]);
	}
	fprintf(fp, "\n");

	/* save the chromosome structure */
	for (i = 0; i < chromo->numNodes; i++) {

		fprintf(fp, "%d\n", chromo->nodes[i]->function);

		for (j = 0; j < chromo->arity; j++) {
			fprintf(fp, "%d,%f\n", chromo->nodes[i]->inputs[j], chromo->nodes[i]->weights[j]);
		}
	}

	for (i = 0; i < chromo->numOutputs; i++) {
		fprintf(fp, "%d,", chromo->outputNodes[i]);
	}

	fclose(fp);
}


/*
	save the given chromosome to a graphviz .dot file
	(www.graphviz.org/‎)
*/
DLL_EXPORT void saveChromosomeDot(struct chromosome *chromo, int weights, char const *fileName) {

	int i, j;
	FILE *fp;

	char colour[20];
	char weight[20];

	fp = fopen(fileName, "w");

	if (fp == NULL) {
		return;
	}

	/* */
	fprintf(fp, "digraph NeuralNetwork {\n");

	/* landscape, square and centre */
	fprintf(fp, "rankdir=LR;\n");
	fprintf(fp, "size=\"4,3\";\n");
	fprintf(fp, "center = true;\n");

	/* for all the inputs */
	for (i = 0; i < getNumChromosomeInputs(chromo); i++) {

		fprintf(fp, "node%d [label=\"(%d) Input\", color=black, labelfontcolor=black, fontcolor=black];\n", i, i);
	}

	/* for all nodes */
	for (i = 0; i < getNumChromosomeNodes(chromo); i++) {

		if (chromo->nodes[i]->active == 1) {
			strncpy(colour, "black", 20);
		}
		else {
			strncpy(colour, "lightgrey", 20);
		}

		fprintf(fp, "node%d [label=\"(%d) %s\", color=%s, labelfontcolor=%s, fontcolor=%s];\n", i + getNumChromosomeInputs(chromo), i + getNumChromosomeInputs(chromo), chromo->funcSet->functionNames[chromo->nodes[i]->function], colour, colour, colour);

		/* for each node input */
		for (j = 0; j < getChromosomeNodeArity(chromo, i); j++) {

			if (weights == 1) {
				snprintf(weight, 20, "%.2f", chromo->nodes[i]->weights[j]);
			}
			else {
				snprintf(weight, 20, " (%d)", j);
			}


			fprintf(fp, "node%d -> node%d [label=\"%s\", labelfontcolor=%s, fontcolor=%s, bold=true, color=%s];\n", chromo->nodes[i]->inputs[j], i + getNumChromosomeInputs(chromo), weight, colour, colour, colour);
		}
	}

	for (i = 0; i < getNumChromosomeOutputs(chromo); i++) {

		fprintf(fp, "node%d [label=\"Output %d\", color=black, labelfontcolor=black, fontcolor=black];\n", i + getNumChromosomeInputs(chromo) + getNumChromosomeNodes(chromo), i);

		fprintf(fp, "node%d -> node%d [labelfontcolor=black, fontcolor=black, bold=true, color=black];\n", chromo->outputNodes[i], i + getNumChromosomeInputs(chromo) + getNumChromosomeNodes(chromo));
	}


	/* place inputs  on same line */
	fprintf(fp, "{ rank = source;");

	for (i = 0; i < getNumChromosomeInputs(chromo); i++) {
		fprintf(fp, " \"node%d\";", i);
	}
	fprintf(fp, " }\n");


	/* place outputs  on same line */
	fprintf(fp, "{ rank = max;");

	for (i = 0; i < getNumChromosomeOutputs(chromo); i++) {
		fprintf(fp, "\"node%d\";", i + getNumChromosomeInputs(chromo) + getNumChromosomeNodes(chromo));
	}
	fprintf(fp, " }\n");


	/* last line of dot file */
	fprintf(fp, "}");

	fclose(fp);
}


/*
	save the given chromosome as a latex equation

	Only compatible with feed-forward networks
	Only fully compatible with custom node functions
*/
DLL_EXPORT void saveChromosomeLatex(struct chromosome *chromo, int weights, char const *fileName) {

	int output;
	int i;
	FILE *fp;

	fp = fopen(fileName, "w");

	if (fp == NULL) {
		return;
	}

	/* document header */
	fprintf(fp, "\\documentclass{article}\n");
	fprintf(fp, "\\begin{document}\n");

	for (output = 0; output < chromo->numOutputs; output++) {

		fprintf(fp, "\\begin{equation}\n");

		/* function inputs */
		if (chromo->numInputs == 0) {
			fprintf(fp, "f()=");
		}
		else {

			fprintf(fp, "f_%d(x_0", output);

			for (i = 1; i < chromo->numInputs; i++) {

				fprintf(fp, ",x_%d", i);
			}

			fprintf(fp, ")=");
		}

		saveChromosomeLatexRecursive(chromo, chromo->outputNodes[output], fp);

		fprintf(fp, "\n\\end{equation}");
	}


	/* document footer */
	fprintf(fp, "\n\\end{document}");

	fclose(fp);
}

/*
	used by saveChromosomeLatex
*/
static void saveChromosomeLatexRecursive(struct chromosome *chromo, int index, FILE *fp) {

	int i;

	if (index < chromo->numInputs) {
		fprintf(fp, "x_%d", index);
		return;
	}

	/* add */
	if (strncmp(chromo->funcSet->functionNames[chromo->nodes[index - chromo->numInputs]->function], "add", FUNCTIONNAMELENGTH) == 0 ) {

		fprintf(fp, "\\left(");

		saveChromosomeLatexRecursive(chromo, chromo->nodes[index - chromo->numInputs]->inputs[0], fp);

		for (i = 1; i < getChromosomeNodeArity(chromo, index - chromo->numInputs); i++) {

			fprintf(fp, " + ");

			saveChromosomeLatexRecursive(chromo, chromo->nodes[index - chromo->numInputs]->inputs[i], fp);
		}

		fprintf(fp, "\\right)");
	}


	/* sub */
	else if (strncmp(chromo->funcSet->functionNames[chromo->nodes[index - chromo->numInputs]->function], "sub", FUNCTIONNAMELENGTH) == 0 ) {

		fprintf(fp, "\\left(");

		saveChromosomeLatexRecursive(chromo, chromo->nodes[index - chromo->numInputs]->inputs[0], fp);

		for (i = 1; i < getChromosomeNodeArity(chromo, index - chromo->numInputs); i++) {

			fprintf(fp, " - ");

			saveChromosomeLatexRecursive(chromo, chromo->nodes[index - chromo->numInputs]->inputs[i], fp);
		}

		fprintf(fp, "\\right)");
	}

	/* mul */
	else if (strncmp(chromo->funcSet->functionNames[chromo->nodes[index - chromo->numInputs]->function], "mul", FUNCTIONNAMELENGTH) == 0 ) {

		fprintf(fp, "\\left(");

		saveChromosomeLatexRecursive(chromo, chromo->nodes[index - chromo->numInputs]->inputs[0], fp);

		for (i = 1; i < getChromosomeNodeArity(chromo, index - chromo->numInputs); i++) {

			fprintf(fp, " \\times ");

			saveChromosomeLatexRecursive(chromo, chromo->nodes[index - chromo->numInputs]->inputs[i], fp);
		}

		fprintf(fp, "\\right)");
	}

	/* div (change to frac)*/
	else if (strncmp(chromo->funcSet->functionNames[chromo->nodes[index - chromo->numInputs]->function], "div", FUNCTIONNAMELENGTH) == 0 ) {

		if (getChromosomeNodeArity(chromo, index - chromo->numInputs) == 1) {
			saveChromosomeLatexRecursive(chromo, chromo->nodes[index - chromo->numInputs]->inputs[0], fp);
		}
		else {

			for (i = 0; i < getChromosomeNodeArity(chromo, index - chromo->numInputs); i++) {

				if (i + 1 < getChromosomeNodeArity(chromo, index - chromo->numInputs)) {
					fprintf(fp, "\\frac{");
					saveChromosomeLatexRecursive(chromo, chromo->nodes[index - chromo->numInputs]->inputs[i], fp);
					fprintf(fp, "}{");
				}
				else if (i + 1 == getChromosomeNodeArity(chromo, index - chromo->numInputs) && getChromosomeNodeArity(chromo, index - chromo->numInputs) > 2) {
					saveChromosomeLatexRecursive(chromo, chromo->nodes[index - chromo->numInputs]->inputs[i], fp);
					fprintf(fp, "}}");
				}
				else {
					saveChromosomeLatexRecursive(chromo, chromo->nodes[index - chromo->numInputs]->inputs[i], fp);
					fprintf(fp, "}");
				}
			}
		}
	}

	/* abs */
	else if (strncmp(chromo->funcSet->functionNames[chromo->nodes[index - chromo->numInputs]->function], "abs", FUNCTIONNAMELENGTH) == 0 ) {

		fprintf(fp, " \\left|");

		saveChromosomeLatexRecursive(chromo, chromo->nodes[index - chromo->numInputs]->inputs[0], fp);

		fprintf(fp, " \\right|");

	}

	/* sqrt */
	else if (strncmp(chromo->funcSet->functionNames[chromo->nodes[index - chromo->numInputs]->function], "sqrt", FUNCTIONNAMELENGTH) == 0 ) {

		fprintf(fp, " \\sqrt{");

		saveChromosomeLatexRecursive(chromo, chromo->nodes[index - chromo->numInputs]->inputs[0], fp);

		fprintf(fp, " }");

	}


	/* sq */
	else if (strncmp(chromo->funcSet->functionNames[chromo->nodes[index - chromo->numInputs]->function], "sq", FUNCTIONNAMELENGTH) == 0 ) {

		fprintf(fp, " (");

		saveChromosomeLatexRecursive(chromo, chromo->nodes[index - chromo->numInputs]->inputs[0], fp);

		fprintf(fp, " )^2");

	}

	/* cube */
	else if (strncmp(chromo->funcSet->functionNames[chromo->nodes[index - chromo->numInputs]->function], "cube", FUNCTIONNAMELENGTH) == 0 ) {

		fprintf(fp, " (");

		saveChromosomeLatexRecursive(chromo, chromo->nodes[index - chromo->numInputs]->inputs[0], fp);

		fprintf(fp, " )^3");

	}

	/* exp */
	else if (strncmp(chromo->funcSet->functionNames[chromo->nodes[index - chromo->numInputs]->function], "exp", FUNCTIONNAMELENGTH) == 0 ) {

		fprintf(fp, " e^{");

		saveChromosomeLatexRecursive(chromo, chromo->nodes[index - chromo->numInputs]->inputs[0], fp);

		fprintf(fp, " }");

	}

	/* sin */
	else if (strncmp(chromo->funcSet->functionNames[chromo->nodes[index - chromo->numInputs]->function], "sin", FUNCTIONNAMELENGTH) == 0 ) {

		fprintf(fp, "\\sin(");

		saveChromosomeLatexRecursive(chromo, chromo->nodes[index - chromo->numInputs]->inputs[0], fp);

		fprintf(fp, " )");

	}

	/* cos */
	else if (strncmp(chromo->funcSet->functionNames[chromo->nodes[index - chromo->numInputs]->function], "cos", FUNCTIONNAMELENGTH) == 0 ) {

		fprintf(fp, " \\cos(");

		saveChromosomeLatexRecursive(chromo, chromo->nodes[index - chromo->numInputs]->inputs[0], fp);

		fprintf(fp, " )");

	}

	/* tan */
	else if (strncmp(chromo->funcSet->functionNames[chromo->nodes[index - chromo->numInputs]->function], "tan", FUNCTIONNAMELENGTH) == 0 ) {

		fprintf(fp, " \\tan(");

		saveChromosomeLatexRecursive(chromo, chromo->nodes[index - chromo->numInputs]->inputs[0], fp);

		fprintf(fp, " )");

	}

	/* rand */
	else if (strncmp(chromo->funcSet->functionNames[chromo->nodes[index - chromo->numInputs]->function], "rand", FUNCTIONNAMELENGTH) == 0 ) {

		fprintf(fp, " rand()");
	}

	/* pi */
	else if (strncmp(chromo->funcSet->functionNames[chromo->nodes[index - chromo->numInputs]->function], "pi", FUNCTIONNAMELENGTH) == 0 ) {

		fprintf(fp, "\\pi");
	}

	/* 0 */
	else if (strncmp(chromo->funcSet->functionNames[chromo->nodes[index - chromo->numInputs]->function], "0", FUNCTIONNAMELENGTH) == 0 ) {

		fprintf(fp, " 0");
	}

	/* 1 */
	else if (strncmp(chromo->funcSet->functionNames[chromo->nodes[index - chromo->numInputs]->function], "1", FUNCTIONNAMELENGTH) == 0 ) {

		fprintf(fp, " 1");
	}

	/* other */
	else {

		fprintf(fp, "%s(", chromo->funcSet->functionNames[chromo->nodes[index - chromo->numInputs]->function]);

		for (i = 0; i < getChromosomeNodeArity(chromo, index - chromo->numInputs); i++) {

			saveChromosomeLatexRecursive(chromo, chromo->nodes[index - chromo->numInputs]->inputs[i], fp);

			if (i < getChromosomeNodeArity(chromo, index - chromo->numInputs) - 1)
				fprintf(fp, ", ");
		}

		fprintf(fp, ")");
	}

}

DLL_EXPORT int compareChromosomes(struct chromosome *chromoA, struct chromosome *chromoB) {

	int i, j;

	/* ensure that the chromosomes don't point to NULL */
	if (chromoA == NULL || chromoB == NULL) {
		return 0;
	}

	/* Check the high level parameters */
	if (chromoA->numInputs != chromoB->numInputs) {
		return 0;
	}

	if (chromoA->numNodes != chromoB->numNodes) {
		return 0;
	}

	if (chromoA->numOutputs != chromoB->numOutputs) {
		return 0;
	}

	if (chromoA->arity != chromoB->arity) {
		return 0;
	}

	/* for each node*/
	for (i = 0; i < chromoA->numNodes; i++) {

		/* Check the function genes */
		if (chromoA->nodes[i]->function != chromoB->nodes[i]->function) {
			return 0;
		}

		/* for each node input */
		for (j = 0; j < chromoA->arity; j++) {

			/* Check the node inputs */
			if (chromoA->nodes[i]->inputs[j] != chromoB->nodes[i]->inputs[j]) {
				return 0;
			}
		}
	}

	/* for all of the outputs */
	for (i = 0; i < chromoA->numOutputs; i++) {

		/* Check the outputs */
		if (chromoA->outputNodes[i] != chromoB->outputNodes[i] ) {
			return 0;
		}
	}

	return 1;
}


DLL_EXPORT int compareChromosomesANN(struct chromosome *chromoA, struct chromosome *chromoB) {

	int i, j;

	/* ensure that the chromosomes don't point to NULL */
	if (chromoA == NULL || chromoB == NULL) {
		return 0;
	}

	/* Check the high level parameters */
	if (chromoA->numInputs != chromoB->numInputs) {
		return 0;
	}

	if (chromoA->numNodes != chromoB->numNodes) {
		return 0;
	}

	if (chromoA->numOutputs != chromoB->numOutputs) {
		return 0;
	}

	if (chromoA->arity != chromoB->arity) {
		return 0;
	}

	/* for each node*/
	for (i = 0; i < chromoA->numNodes; i++) {

		/* Check the function genes */
		if (chromoA->nodes[i]->function != chromoB->nodes[i]->function) {
			return 0;
		}

		/* for each node input */
		for (j = 0; j < chromoA->arity; j++) {

			/* Check the node inputs */
			if (chromoA->nodes[i]->inputs[j] != chromoB->nodes[i]->inputs[j]) {
				return 0;
			}

			/* Check the connection weights inputs */
			if (chromoA->nodes[i]->weights[j] != chromoB->nodes[i]->weights[j]) {
				return 0;
			}
		}
	}

	/* for all of the outputs */
	for (i = 0; i < chromoA->numOutputs; i++) {

		/* Check the outputs */
		if (chromoA->outputNodes[i] != chromoB->outputNodes[i] ) {
			return 0;
		}
	}

	return 1;
}






DLL_EXPORT int compareChromosomesActiveNodes(struct chromosome *chromoA, struct chromosome *chromoB) {

	int i, j;

	/* ensure that the chromosomes don't point to NULL */
	if (chromoA == NULL || chromoB == NULL) {
		return 0;
	}

	/* Check the high level parameters */
	if (chromoA->numInputs != chromoB->numInputs) {
		return 0;
	}

	if (chromoA->numNodes != chromoB->numNodes) {
		return 0;
	}

	if (chromoA->numOutputs != chromoB->numOutputs) {
		return 0;
	}

	if (chromoA->arity != chromoB->arity) {
		return 0;
	}

	/* for each node*/
	for (i = 0; i < chromoA->numNodes; i++) {

		/* if the node is active in both chromosomes */
		if (chromoA->nodes[i]->active == 1 && chromoB->nodes[i]->active == 1) {

			/* Check the function genes */
			if (chromoA->nodes[i]->function != chromoB->nodes[i]->function) {
				return 0;
			}

			/* for each node input */
			for (j = 0; j < chromoA->arity; j++) {

				/* Check the node inputs */
				if (chromoA->nodes[i]->inputs[j] != chromoB->nodes[i]->inputs[j]) {
					return 0;
				}				
			}
		}
		/* if the node is active in one chromosome */
		else if (chromoA->nodes[i]->active != chromoB->nodes[i]->active ) {
			return 0;
		}

		/* The node is inactive in both chromosomes */
		else {
			/* do nothing */
		}
	}

	/* for all of the outputs */
	for (i = 0; i < chromoA->numOutputs; i++) {

		/* Check the outputs */
		if (chromoA->outputNodes[i] != chromoB->outputNodes[i] ) {
			return 0;
		}
	}

	return 1;
}


DLL_EXPORT int compareChromosomesActiveNodesANN(struct chromosome *chromoA, struct chromosome *chromoB) {

	int i, j;

	/* ensure that the chromosomes don't point to NULL */
	if (chromoA == NULL || chromoB == NULL) {
		return 0;
	}

	/* Check the high level parameters */
	if (chromoA->numInputs != chromoB->numInputs) {
		return 0;
	}

	if (chromoA->numNodes != chromoB->numNodes) {
		return 0;
	}

	if (chromoA->numOutputs != chromoB->numOutputs) {
		return 0;
	}

	if (chromoA->arity != chromoB->arity) {
		return 0;
	}

	/* for each node*/
	for (i = 0; i < chromoA->numNodes; i++) {

		/* if the node is active in both chromosomes */
		if (chromoA->nodes[i]->active == 1 && chromoB->nodes[i]->active == 1) {

			/* Check the function genes */
			if (chromoA->nodes[i]->function != chromoB->nodes[i]->function) {
				return 0;
			}

			/* for each node input */
			for (j = 0; j < chromoA->arity; j++) {

				/* Check the node inputs */
				if (chromoA->nodes[i]->inputs[j] != chromoB->nodes[i]->inputs[j]) {
					return 0;
				}
				
				/* Check the connection weights inputs */
				if (chromoA->nodes[i]->weights[j] != chromoB->nodes[i]->weights[j]) {
					return 0;
				}
			}
		}
		/* if the node is active in one chromosome */
		else if (chromoA->nodes[i]->active != chromoB->nodes[i]->active ) {
			return 0;
		}

		/* The node is inactive in both chromosomes */
		else {
			/* do nothing */
		}
	}

	/* for all of the outputs */
	for (i = 0; i < chromoA->numOutputs; i++) {

		/* Check the outputs */
		if (chromoA->outputNodes[i] != chromoB->outputNodes[i] ) {
			return 0;
		}
	}

	return 1;
}




/*
	Mutates the given chromosome using the mutation method described in parameters
*/
DLL_EXPORT void mutateChromosome(struct parameters *params, struct chromosome *chromo) {

	params->mutationType(params, chromo);

	setChromosomeActiveNodes(chromo);
}


/*
	removes the inactive nodes from the given chromosome
*/
DLL_EXPORT void removeInactiveNodes(struct chromosome *chromo) {

	int i, j, k;

	int originalNumNodes = chromo->numNodes;

	/* set the active nodes */
	setChromosomeActiveNodes(chromo);

	/* for all nodes */
	for (i = 0; i < chromo->numNodes - 1; i++) {

		/* if the node is inactive */
		if (chromo->nodes[i]->active == 0) {

			/* set the node to be the next node */
			for (j = i; j < chromo->numNodes - 1; j++) {
				copyNode(chromo->nodes[j], chromo->nodes[j + 1]);
			}

			/* */
			for (j = 0; j < chromo->numNodes; j++) {
				for (k = 0; k < chromo->arity; k++) {

					if (chromo->nodes[j]->inputs[k] >= i + chromo->numInputs) {
						chromo->nodes[j]->inputs[k]--;
					}
				}
			}

			/* for the number of chromosome outputs */
			for (j = 0; j < chromo->numOutputs; j++) {

				if (chromo->outputNodes[j] >= i + chromo->numInputs) {
					chromo->outputNodes[j]--;
				}
			}

			/* de-increment the number of nodes */
			chromo->numNodes--;

			/* made the newly assigned node be evaluated */
			i--;
		}
	}

	for (i = chromo->numNodes; i < originalNumNodes; i++) {
		freeNode(chromo->nodes[i]);
	}

	if (chromo->nodes[chromo->numNodes - 1]->active == 0) {
		freeNode(chromo->nodes[chromo->numNodes - 1]);
		chromo->numNodes--;
	}

	/* reallocate the memory associated with the chromosome */
	chromo->nodes = (struct node**)realloc(chromo->nodes, chromo->numNodes * sizeof(struct node*));
	chromo->activeNodes = (int*)realloc(chromo->activeNodes, chromo->numNodes * sizeof(int));

	/* set the active nodes */
	setChromosomeActiveNodes(chromo);
}


/*
	sets the fitness of the given chromosome
*/
DLL_EXPORT void setChromosomeFitness(struct parameters *params, struct chromosome *chromo, struct dataSet *data) {

	double fitness;

	setChromosomeActiveNodes(chromo);

	resetChromosome(chromo);

	fitness = params->fitnessFunction(params, chromo, data);

	chromo->fitness = fitness;
}


/*
	reset the output values of all chromosome nodes to zero
*/
DLL_EXPORT void resetChromosome(struct chromosome *chromo) {

	int i;

	for (i = 0; i < chromo->numNodes; i++) {
		chromo->nodes[i]->output = 0;
	}
}

/*
	copies the contents of one chromosome to another. Provided the number of inputs, nodes, outputs and node arity are the same.
*/
DLL_EXPORT void copyChromosome(struct chromosome *chromoDest, struct chromosome *chromoSrc) {

	int i;

	/* error checking  */
	if (chromoDest->numInputs != chromoSrc->numInputs) {
		printf("Error: cannot copy a chromosome to a chromosome of different dimensions. The number of chromosome inputs do not match.\n");
		printf("Terminating CGP-Library.\n");
		exit(0);
	}

	if (chromoDest->numNodes != chromoSrc->numNodes) {
		printf("Error: cannot copy a chromosome to a chromosome of different dimensions. The number of chromosome nodes do not match.\n");
		printf("Terminating CGP-Library.\n");
		exit(0);
	}

	if (chromoDest->numOutputs != chromoSrc->numOutputs) {
		printf("Error: cannot copy a chromosome to a chromosome of different dimensions. The number of chromosome outputs do not match.\n");
		printf("Terminating CGP-Library.\n");
		exit(0);
	}

	if (chromoDest->arity != chromoSrc->arity) {
		printf("Error: cannot copy a chromosome to a chromosome of different dimensions. The arity of the chromosome nodes do not match.\n");
		printf("Terminating CGP-Library.\n");
		exit(0);
	}

	/* copy nodes and which are active */
	for (i = 0; i < chromoSrc->numNodes; i++) {
		copyNode(chromoDest->nodes[i],  chromoSrc->nodes[i]);
		chromoDest->activeNodes[i] = chromoSrc->activeNodes[i];
	}

	/* copy functionset */
	copyFunctionSet(chromoDest->funcSet, chromoSrc->funcSet);

	/* copy each of the chromosomes outputs */
	for (i = 0; i < chromoSrc->numOutputs; i++) {
		chromoDest->outputNodes[i] = chromoSrc->outputNodes[i];
	}

	/* copy the number of active node */
	chromoDest->numActiveNodes = chromoSrc->numActiveNodes;

	/* copy the fitness */
	chromoDest->fitness = chromoSrc->fitness;

	/* copy generation */
	chromoDest->generation = chromoSrc->generation;
}

/*
	Gets the number of chromosome inputs
*/
DLL_EXPORT int getNumChromosomeInputs(struct chromosome *chromo) {
	return chromo->numInputs;
}

/*
	Gets the number of chromosome nodes
*/
DLL_EXPORT int getNumChromosomeNodes(struct chromosome *chromo) {
	return chromo->numNodes;
}

/*
	Gets the number of chromosome active nodes
*/
DLL_EXPORT int getNumChromosomeActiveNodes(struct chromosome *chromo) {
	return chromo->numActiveNodes;
}

/*
	Gets the number of chromosome outputs
*/
DLL_EXPORT int getNumChromosomeOutputs(struct chromosome *chromo) {
	return chromo->numOutputs;
}

/*
	Gets the chromosome node arity
*/
DLL_EXPORT int getChromosomeNodeArity(struct chromosome *chromo, int index) {

	int chromoArity = chromo->arity;
	int maxArity = chromo->funcSet->maxNumInputs[chromo->nodes[index]->function];

	if (maxArity == -1) {
		return chromoArity;
	}
	else if (maxArity < chromoArity) {
		return maxArity;
	}
	else {
		return chromoArity;
	}
}

/*
	Gets the chromosome fitness
*/
DLL_EXPORT double getChromosomeFitness(struct chromosome *chromo) {
	return chromo->fitness;
}

/*
	Gets the number of active connections in the given chromosome
*/
DLL_EXPORT int getNumChromosomeActiveConnections(struct chromosome *chromo) {

	int i;
	int complexity = 0;

	for (i = 0; i < chromo->numActiveNodes; i++) {
		complexity += chromo->nodes[chromo->activeNodes[i]]->actArity;
	}

	return complexity;
}

/*
	Gets the number of generations required to find the given chromosome
*/
DLL_EXPORT int getChromosomeGenerations(struct chromosome *chromo) {
	return chromo->generation;
}


/*
	set the active nodes in the given chromosome
*/
static void setChromosomeActiveNodes(struct chromosome *chromo) {

	int i;

	/* error checking */
	if (chromo == NULL) {
		printf("Error: chromosome has not been initialised and so the active nodes cannot be set.\n");
		return;
	}

	/* set the number of active nodes to zero */
	chromo->numActiveNodes = 0;

	/* reset the active nodes */
	for (i = 0; i < chromo->numNodes; i++) {
		chromo->nodes[i]->active = 0;
	}

	/* start the recursive search for active nodes from the output nodes for the number of output nodes */
	for (i = 0; i < chromo->numOutputs; i++) {

		/* if the output connects to a chromosome input, skip */
		if (chromo->outputNodes[i] < chromo->numInputs) {
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
static void recursivelySetActiveNodes(struct chromosome *chromo, int nodeIndex) {

	int i;

	/* if the given node is an input, stop */
	if (nodeIndex < chromo->numInputs) {
		return;
	}

	/* if the given node has already been flagged as active */
	if (chromo->nodes[nodeIndex - chromo->numInputs]->active == 1) {
		return;
	}

	/* log the node as active */
	chromo->nodes[nodeIndex - chromo->numInputs]->active = 1;
	chromo->activeNodes[chromo->numActiveNodes] = nodeIndex - chromo->numInputs;
	chromo->numActiveNodes++;

	/* set the nodes actual arity*/
	chromo->nodes[nodeIndex - chromo->numInputs]->actArity = getChromosomeNodeArity(chromo, nodeIndex - chromo->numInputs);

	/* recursively log all the nodes to which the current nodes connect as active */
	for (i = 0; i < chromo->nodes[nodeIndex - chromo->numInputs]->actArity; i++) {
		recursivelySetActiveNodes(chromo, chromo->nodes[nodeIndex - chromo->numInputs]->inputs[i]);
	}
}


/*
	Sorts the given array of chromosomes by fitness, lowest to highest
	uses insertion sort (quickish and stable)
*/
static void sortChromosomeArray(struct chromosome **chromoArray, int numChromos) {

	int i, j;
	struct chromosome *chromoTmp;

	for (i = 0; i < numChromos; i++) {
		for (j = i; j < numChromos; j++) {

			if (chromoArray[i]->fitness > chromoArray[j]->fitness) {
				chromoTmp = chromoArray[i];
				chromoArray[i] = chromoArray[j];
				chromoArray[j] = chromoTmp;
			}
		}
	}
}



/*
	Dataset functions
*/


/*
	Initialises data structure and assigns values for given arrays
	arrays must take the form
	inputs[numSamples][numInputs]
	outputs[numSamples][numOutputs]
*/
DLL_EXPORT struct dataSet *initialiseDataSetFromArrays(int numInputs, int numOutputs, int numSamples, double *inputs, double *outputs) {

	int i, j;
	struct dataSet *data;

	/* initialise memory for data structure */
	data = (struct dataSet*)malloc(sizeof(struct dataSet));

	data->numInputs = numInputs;
	data->numOutputs = numOutputs;
	data->numSamples = numSamples;

	data->inputData = (double**)malloc(data->numSamples * sizeof(double*));
	data->outputData = (double**)malloc(data->numSamples * sizeof(double*));

	for (i = 0; i < data->numSamples; i++) {

		data->inputData[i] = (double*)malloc(data->numInputs * sizeof(double));
		data->outputData[i] = (double*)malloc(data->numOutputs * sizeof(double));

		for (j = 0; j < data->numInputs; j++) {
			data->inputData[i][j] = inputs[(i * data->numInputs) + j];
		}

		for (j = 0; j < data->numOutputs; j++) {
			data->outputData[i][j] = outputs[(i * data->numOutputs) + j];
		}
	}

	return data;
}


/*
	Initialises data structure and assigns values of given file
*/
DLL_EXPORT struct dataSet *initialiseDataSetFromFile(char const *file) {

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
	if (fp == NULL) {
		printf("Error: file '%s' cannot be found.\nTerminating CGP-Library.\n", file);
		exit(0);
	}

	/* initialise memory for data structure */
	data = (struct dataSet*)malloc(sizeof(struct dataSet));

	/* for every line in the given file */
	while ( (line = fgets(buffer, sizeof(buffer), fp)) != NULL) {

		/* deal with the first line containing meta data */
		if (lineNum == -1) {

			sscanf(line, "%d,%d,%d", &(data->numInputs), &(data->numOutputs), &(data->numSamples));

			data->inputData = (double**)malloc(data->numSamples * sizeof(double*));
			data->outputData = (double**)malloc(data->numSamples * sizeof(double*));

			for (i = 0; i < data->numSamples; i++) {
				data->inputData[i] = (double*)malloc(data->numInputs * sizeof(double));
				data->outputData[i] = (double*)malloc(data->numOutputs * sizeof(double));
			}
		}
		/* the other lines contain input output pairs */
		else {

			/* get the first value on the given line */
			record = strtok(line, " ,\n");
			col = 0;

			/* until end of line */
			while (record != NULL) {

				/* if its an input value */
				if (col < data->numInputs) {
					data->inputData[lineNum][col] = atof(record);
				}

				/* if its an output value */
				else {

					data->outputData[lineNum][col - data->numInputs] = atof(record);
				}

				/* get the next value on the given line */
				record = strtok(NULL, " ,\n");

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
	frees given dataSet
*/
DLL_EXPORT void freeDataSet(struct dataSet *data) {

	int i;

	/* attempt to prevent user double freeing */
	if (data == NULL) {
		printf("Warning: double freeing of dataSet prevented.\n");
		return;
	}

	for (i = 0; i < data->numSamples; i++) {
		free(data->inputData[i]);
		free(data->outputData[i]);
	}

	free(data->inputData);
	free(data->outputData);
	free(data);
}


/*
	prints the given data structure to the screen
*/
DLL_EXPORT void printDataSet(struct dataSet *data) {

	int i, j;

	printf("DATA SET\n");
	printf("Inputs: %d, ", data->numInputs);
	printf("Outputs: %d, ", data->numOutputs);
	printf("Samples: %d\n", data->numSamples);

	for (i = 0; i < data->numSamples; i++) {

		for (j = 0; j < data->numInputs; j++) {
			printf("%f ", data->inputData[i][j]);
		}

		printf(" : ");

		for (j = 0; j < data->numOutputs; j++) {
			printf("%f ", data->outputData[i][j]);
		}

		printf("\n");
	}
}


/*
	saves dataset to file
*/
DLL_EXPORT void saveDataSet(struct dataSet *data, char const *fileName) {

	int i, j;
	FILE *fp;

	fp = fopen(fileName, "w");

	if (fp == NULL) {
		printf("Warning: cannot save data set to %s. Data set was not saved.\n", fileName);
		return;
	}

	fprintf(fp, "%d,", data->numInputs);
	fprintf(fp, "%d,", data->numOutputs);
	fprintf(fp, "%d,", data->numSamples);
	fprintf(fp, "\n");


	for (i = 0; i < data->numSamples; i++) {

		for (j = 0; j < data->numInputs; j++) {
			fprintf(fp, "%f,", data->inputData[i][j]);
		}

		for (j = 0; j < data->numOutputs; j++) {
			fprintf(fp, "%f,", data->outputData[i][j]);
		}

		fprintf(fp, "\n");
	}

	fclose(fp);
}


/*
	returns the number of inputs for each sample in the given dataSet
*/
DLL_EXPORT int getNumDataSetInputs(struct dataSet *data) {
	return data->numInputs;
}


/*
	returns the number of outputs for each sample in the given dataSet
*/
DLL_EXPORT int getNumDataSetOutputs(struct dataSet *data) {
	return data->numOutputs;
}


/*
	returns the number of samples in the given dataSet
*/
DLL_EXPORT int getNumDataSetSamples(struct dataSet *data) {
	return data->numSamples;
}


/*
	returns the inputs of the given sample of the given dataSet
*/
DLL_EXPORT double *getDataSetSampleInputs(struct dataSet *data, int sample) {
	return data->inputData[sample];
}


/*
	returns the given input of the given sample of the given dataSet
*/
DLL_EXPORT double getDataSetSampleInput(struct dataSet *data, int sample, int input) {
	return data->inputData[sample][input];
}


/*
	returns the outputs of the given sample of the given dataSet
*/
DLL_EXPORT double *getDataSetSampleOutputs(struct dataSet *data, int sample) {
	return data->outputData[sample];
}


/*
	returns the given output of the given sample of the given dataSet
*/
DLL_EXPORT double getDataSetSampleOutput(struct dataSet *data, int sample, int output) {
	return data->outputData[sample][output];
}



/*
	Results Functions
*/


/*
	initialises a results structure
*/
struct results* initialiseResults(struct parameters *params, int numRuns) {

	struct results *rels;

	rels = (struct results*)malloc(sizeof(struct results));
	rels->bestChromosomes = (struct chromosome**)malloc(numRuns * sizeof(struct chromosome*));

	rels->numRuns = numRuns;

	/*
		Initialised chromosomes are returns from runCGP and stored in a results structure.
		Therefore they should not be initialised here.
	*/

	return rels;
}


/*
	free a initialised results structure
*/
DLL_EXPORT void freeResults(struct results *rels) {

	int i;

	/* attempt to prevent user double freeing */
	if (rels == NULL) {
		printf("Warning: double freeing of results prevented.\n");
		return;
	}

	for (i = 0; i < rels->numRuns; i++) {
		freeChromosome(rels->bestChromosomes[i]);
	}

	free(rels->bestChromosomes);
	free(rels);
}


/*
	saves results structure to file
*/
DLL_EXPORT void saveResults(struct results *rels, char const *fileName) {

	FILE *fp;
	int i;

	struct chromosome *chromoTemp;

	if (rels == NULL) {
		printf("Warning: cannot save uninitialised results structure. Results not saved.\n");
		return;
	}

	fp = fopen(fileName, "w");

	if (fp == NULL) {
		printf("Warning: cannot open '%s' and so cannot save results to that file. Results not saved.\n", fileName);
		return;
	}

	fprintf(fp, "Run,Fitness,Generations,Active Nodes\n");

	for (i = 0; i < rels->numRuns; i++) {

		chromoTemp = getChromosome(rels, i);

		fprintf(fp, "%d,%f,%d,%d\n", i, chromoTemp->fitness, chromoTemp->generation, chromoTemp->numActiveNodes);

		freeChromosome(chromoTemp);
	}

	fclose(fp);
}


/*
	Gets the number of chromosomes in the results structure
*/
DLL_EXPORT int getNumChromosomes(struct results *rels) {
	return rels->numRuns;
}


/*
	returns the average number of chromosome active nodes from repeated
	run results specified in rels.
*/
DLL_EXPORT double getAverageActiveNodes(struct results *rels) {

	int i;
	double avgActiveNodes = 0;
	struct chromosome *chromoTemp;

	for (i = 0; i < getNumChromosomes(rels); i++) {

		chromoTemp = rels->bestChromosomes[i];

		avgActiveNodes += getNumChromosomeActiveNodes(chromoTemp);
	}

	avgActiveNodes = avgActiveNodes / getNumChromosomes(rels);

	return avgActiveNodes;
}


/*
	returns the median number of chromosome active nodes from repeated
	run results specified in rels.
*/
DLL_EXPORT double getMedianActiveNodes(struct results *rels) {

	int i;
	double medActiveNodes = 0;

	int *array = (int*)malloc(getNumChromosomes(rels) * sizeof(int));

	for (i = 0; i < getNumChromosomes(rels); i++) {
		array[i] = getNumChromosomeActiveNodes(rels->bestChromosomes[i]);
	}

	medActiveNodes = medianInt(array, getNumChromosomes(rels));

	free(array);

	return medActiveNodes;
}


static double medianInt(const int *anArray, const int length) {

	int i;
	int *copyArray = (int*)malloc(length * sizeof(int));
	double median;

	/* make a copy of the array */
	for (i = 0; i < length; i++) {
		copyArray[i] = anArray[i];
	}

	/* sort the copy array */
	sortIntArray(copyArray, length);

	/* if even */
	if (length % 2 == 0) {
		median = (copyArray[(length / 2)] + copyArray[(length / 2) - 1] ) / 2;
	}

	/* if odd */
	else {
		median = copyArray[ (length - 1) / 2];
	}

	free(copyArray);

	return median;
}

static double medianDouble(const double *anArray, const int length) {

	int i;
	double *copyArray = (double*)malloc(length * sizeof(double));
	double median;

	/* make a copy of the array */
	for (i = 0; i < length; i++) {
		copyArray[i] = anArray[i];
	}

	/* sort the copy array */
	sortDoubleArray(copyArray, length);

	/* if even */
	if (length % 2 == 0) {
		median = (copyArray[(length / 2)] + copyArray[(length / 2) - 1] ) / 2;
	}

	/* if odd */
	else {
		median = copyArray[ (length - 1) / 2];
	}

	free(copyArray);

	return median;
}



/*
	returns the average chromosome fitness from repeated
	run results specified in rels.
*/
DLL_EXPORT double getAverageFitness(struct results *rels) {

	int i;
	double avgFit = 0;
	struct chromosome *chromoTemp;


	for (i = 0; i < getNumChromosomes(rels); i++) {

		chromoTemp = rels->bestChromosomes[i];

		avgFit += getChromosomeFitness(chromoTemp);
	}

	avgFit = avgFit / getNumChromosomes(rels);

	return avgFit;
}


/*
	returns the median chromosome fitness from repeated
	run results specified in rels.
*/
DLL_EXPORT double getMedianFitness(struct results *rels) {

	int i;
	double med = 0;

	double *array = (double*)malloc(getNumChromosomes(rels) * sizeof(double));

	for (i = 0; i < getNumChromosomes(rels); i++) {
		array[i] = getChromosomeFitness(rels->bestChromosomes[i]);
	}

	med = medianDouble(array, getNumChromosomes(rels));

	free(array);

	return med;
}



/*
	returns the average number of generations used by each run  specified in rels.
*/
DLL_EXPORT double getAverageGenerations(struct results *rels) {

	int i;
	double avgGens = 0;
	struct chromosome *chromoTemp;

	for (i = 0; i < getNumChromosomes(rels); i++) {

		chromoTemp = rels->bestChromosomes[i];

		avgGens += getChromosomeGenerations(chromoTemp);
	}

	avgGens = avgGens / getNumChromosomes(rels);

	return avgGens;
}


/*
	returns the median number of generations used by each run  specified in rels.
*/
DLL_EXPORT double getMedianGenerations(struct results *rels) {

	int i;
	double med = 0;

	int *array = (int*)malloc(getNumChromosomes(rels) * sizeof(int));

	for (i = 0; i < getNumChromosomes(rels); i++) {
		array[i] = getChromosomeGenerations(rels->bestChromosomes[i]);
	}

	med = medianInt(array, getNumChromosomes(rels));

	free(array);

	return med;
}



/*
	returns a pointer to a copy of the best chromosomes found on the given run in rels.
*/
DLL_EXPORT struct chromosome* getChromosome(struct results *rels, int run) {

	struct chromosome *chromo;

	/* do some error checking */
	if (rels == NULL) {
		printf("Error: cannot get best chromosome from uninitialised results.\nTerminating CGP-Library.\n");
		exit(0);
	}

	chromo = initialiseChromosomeFromChromosome(rels->bestChromosomes[run]);

	return chromo;
}





/*
	CGP Functions
*/


/*
	Other Functions
*/








/*
	Mutation Methods
*/



/*
	Conductions point mutation on the give chromosome. A predetermined
	number of chromosome genes are randomly selected and changed to
	a random valid allele. The number of mutations is the number of chromosome
	genes multiplied by the mutation rate. Each gene has equal probability
	of being selected.

	DO NOT USE WITH ANN
*/
static void pointMutation(struct parameters *params, struct chromosome *chromo) {

	int i;
	int numGenes;
	int numFunctionGenes, numInputGenes, numOutputGenes;
	int numGenesToMutate;
	int geneToMutate;
	int nodeIndex;
	int nodeInputIndex;

	/* get the number of each type of gene */
	numFunctionGenes = params->numNodes;
	numInputGenes = params->numNodes * params->arity;
	numOutputGenes = params->numOutputs;

	/* set the total number of chromosome genes */
	numGenes = numFunctionGenes + numInputGenes + numOutputGenes;

	/* calculate the number of genes to mutate */
	numGenesToMutate = (int)roundf(numGenes * params->mutationRate);

	/* for the number of genes to mutate */
	for (i = 0; i < numGenesToMutate; i++) {

		/* select a random gene */
		geneToMutate = randInt(numGenes);

		/* mutate function gene */
		if (geneToMutate < numFunctionGenes) {

			nodeIndex = geneToMutate;

			chromo->nodes[nodeIndex]->function = getRandomFunction(chromo->funcSet->numFunctions);
		}

		/* mutate node input gene */
		else if (geneToMutate < numFunctionGenes + numInputGenes) {

			nodeIndex = (int) ((geneToMutate - numFunctionGenes) / chromo->arity);
			nodeInputIndex = (geneToMutate - numFunctionGenes) % chromo->arity;

			chromo->nodes[nodeIndex]->inputs[nodeInputIndex] = getRandomNodeInput(chromo->numInputs, chromo->numNodes, nodeIndex, params->recurrentConnectionProbability);
		}

		/* mutate output gene */
		else {
			nodeIndex = geneToMutate - numFunctionGenes - numInputGenes;
			chromo->outputNodes[nodeIndex] = getRandomChromosomeOutput(chromo->numInputs, chromo->numNodes, params->shortcutConnections);
		}
	}
}


/*
	Same as pointMutation but also mutates weight genes. The reason this is separated is
	that point mutation should always mutate the same number of genes. When weight genes are not
	used many mutations will not do anything and so the number of actual mutations varies.
	- needs explaining better...
*/
static void pointMutationANN(struct parameters *params, struct chromosome *chromo) {

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
	for (i = 0; i < numGenesToMutate; i++) {

		/* select a random gene */
		geneToMutate = randInt(numGenes);

		/* mutate function gene */
		if (geneToMutate < numFunctionGenes) {

			nodeIndex = geneToMutate;

			chromo->nodes[nodeIndex]->function = getRandomFunction(chromo->funcSet->numFunctions);
		}

		/* mutate node input gene */
		else if (geneToMutate < numFunctionGenes + numInputGenes) {

			nodeIndex = (int) ((geneToMutate - numFunctionGenes) / chromo->arity);
			nodeInputIndex = (geneToMutate - numFunctionGenes) % chromo->arity;

			chromo->nodes[nodeIndex]->inputs[nodeInputIndex] = getRandomNodeInput(chromo->numInputs, chromo->numNodes, nodeIndex, params->recurrentConnectionProbability);
		}

		/* mutate connection weight */
		else if (geneToMutate < numFunctionGenes + numInputGenes + numWeightGenes) {

			nodeIndex = (int) ((geneToMutate - numFunctionGenes - numInputGenes) / chromo->arity);
			nodeInputIndex = (geneToMutate - numFunctionGenes - numInputGenes) % chromo->arity;

			chromo->nodes[nodeIndex]->weights[nodeInputIndex] = getRandomConnectionWeight(params->connectionWeightRange);
		}

		/* mutate output gene */
		else {
			nodeIndex = geneToMutate - numFunctionGenes - numInputGenes - numWeightGenes;
			chromo->outputNodes[nodeIndex] = getRandomChromosomeOutput(chromo->numInputs, chromo->numNodes, params->shortcutConnections);
		}
	}
}



/*
	Conductions a single active mutation on the give chromosome.

	DO NOT USE WITH ANN
*/
static void singleMutation(struct parameters *params, struct chromosome *chromo) {

	int numFunctionGenes, numInputGenes, numOutputGenes;
	int numGenes;
	int geneToMutate;
	int nodeIndex;
	int nodeInputIndex;

	int mutatedActive = 0;
	int previousGeneValue;
	int newGeneValue;

	/* get the number of each type of gene */
	numFunctionGenes = params->numNodes;
	numInputGenes = params->numNodes * params->arity;
	numOutputGenes = params->numOutputs;

	/* set the total number of chromosome genes */
	numGenes = numFunctionGenes + numInputGenes + numOutputGenes;

	/* while active gene not mutated */
	while (mutatedActive == 0) {

		/* select a random gene */
		geneToMutate = randInt(numGenes);

		/* mutate function gene */
		if (geneToMutate < numFunctionGenes) {

			nodeIndex = geneToMutate;

			previousGeneValue = chromo->nodes[nodeIndex]->function;

			chromo->nodes[nodeIndex]->function = getRandomFunction(chromo->funcSet->numFunctions);

			newGeneValue = chromo->nodes[nodeIndex]->function;

			if ((previousGeneValue != newGeneValue) && (chromo->nodes[nodeIndex]->active == 1)) {
				mutatedActive = 1;
			}

		}

		/* mutate node input gene */
		else if (geneToMutate < numFunctionGenes + numInputGenes) {

			nodeIndex = (int) ((geneToMutate - numFunctionGenes) / chromo->arity);
			nodeInputIndex = (geneToMutate - numFunctionGenes) % chromo->arity;

			previousGeneValue = chromo->nodes[nodeIndex]->inputs[nodeInputIndex];

			chromo->nodes[nodeIndex]->inputs[nodeInputIndex] = getRandomNodeInput(chromo->numInputs, chromo->numNodes, nodeIndex, params->recurrentConnectionProbability);

			newGeneValue = chromo->nodes[nodeIndex]->inputs[nodeInputIndex];

			if ((previousGeneValue != newGeneValue) && (chromo->nodes[nodeIndex]->active == 1)) {
				mutatedActive = 1;
			}
		}

		/* mutate output gene */
		else {
			nodeIndex = geneToMutate - numFunctionGenes - numInputGenes;

			previousGeneValue = chromo->outputNodes[nodeIndex];

			chromo->outputNodes[nodeIndex] = getRandomChromosomeOutput(chromo->numInputs, chromo->numNodes, params->shortcutConnections);

			newGeneValue = chromo->outputNodes[nodeIndex];

			if (previousGeneValue != newGeneValue) {
				mutatedActive = 1;
			}
		}
	}
}


/*
	Conductions probabilistic mutation on the given chromosome. Each chromosome
	gene is changed to a random valid allele with a probability specified in
	parameters.
*/
static void probabilisticMutation(struct parameters *params, struct chromosome *chromo) {

	int i, j;

	/* for every nodes in the chromosome */
	for (i = 0; i < params->numNodes; i++) {

		/* mutate the function gene */
		if (randDecimal() <= params->mutationRate) {
			chromo->nodes[i]->function = getRandomFunction(chromo->funcSet->numFunctions);
		}

		/* for every input to each chromosome */
		for (j = 0; j < params->arity; j++) {

			/* mutate the node input */
			if (randDecimal() <= params->mutationRate) {
				chromo->nodes[i]->inputs[j] = getRandomNodeInput(chromo->numInputs, chromo->numNodes, i, params->recurrentConnectionProbability);
			}

			/* mutate the node connection weight */
			if (randDecimal() <= params->mutationRate) {
				chromo->nodes[i]->weights[j] = getRandomConnectionWeight(params->connectionWeightRange);
			}
		}
	}

	/* for every chromosome output */
	for (i = 0; i < params->numOutputs; i++) {

		/* mutate the chromosome output */
		if (randDecimal() <= params->mutationRate) {
			chromo->outputNodes[i] = getRandomChromosomeOutput(chromo->numInputs, chromo->numNodes, params->shortcutConnections);
		}
	}
}

/*
	Conductions probabilistic mutation on the active nodes in the given
	chromosome. Each chromosome gene is changed to a random valid allele
	with a probability specified in parameters.
*/
static void probabilisticMutationOnlyActive(struct parameters *params, struct chromosome *chromo) {

	int i, j;
	int activeNode;

	/* for every active node in the chromosome */
	for (i = 0; i < chromo->numActiveNodes; i++) {

		activeNode = chromo->activeNodes[i];

		/* mutate the function gene */
		if (randDecimal() <= params->mutationRate) {
			chromo->nodes[activeNode]->function = getRandomFunction(chromo->funcSet->numFunctions);
		}

		/* for every input to each chromosome */
		for (j = 0; j < params->arity; j++) {

			/* mutate the node input */
			if (randDecimal() <= params->mutationRate) {
				chromo->nodes[activeNode]->inputs[j] = getRandomNodeInput(chromo->numInputs, chromo->numNodes, activeNode, params->recurrentConnectionProbability);
			}

			/* mutate the node connection weight */
			if (randDecimal() <= params->mutationRate) {
				chromo->nodes[activeNode]->weights[j] = getRandomConnectionWeight(params->connectionWeightRange);
			}
		}
	}

	/* for every chromosome output */
	for (i = 0; i < params->numOutputs; i++) {

		/* mutate the chromosome output */
		if (randDecimal() <= params->mutationRate) {
			chromo->outputNodes[i] = getRandomChromosomeOutput(chromo->numInputs, chromo->numNodes, params->shortcutConnections);
		}
	}
}



/*
	Sets the random number seed
*/
DLL_EXPORT void setRandomNumberSeed(unsigned int seed) {
	srand(seed);
}


/*
	repetitively applies runCGP to obtain average behaviour
*/
DLL_EXPORT struct results* repeatCGP(struct parameters *params, struct dataSet *data, int numGens, int numRuns) {

	int i;
	struct results *rels;
	int updateFrequency = params->updateFrequency;

	/* set the update frequency so as to to so generational results */
	params->updateFrequency = 0;

	rels = initialiseResults(params, numRuns);

	printf("Run\tFitness\t\tGenerations\tActive Nodes\n");

	/* for each run */
	#pragma omp parallel for default(none), shared(numRuns,rels,params,data,numGens), schedule(dynamic), num_threads(params->numThreads)
	for (i = 0; i < numRuns; i++) {

		/* run cgp */
		rels->bestChromosomes[i] = runCGP(params, data, numGens);

		printf("%d\t%f\t%d\t\t%d\n", i, rels->bestChromosomes[i]->fitness, rels->bestChromosomes[i]->generation, rels->bestChromosomes[i]->numActiveNodes);
	}

	printf("----------------------------------------------------\n");
	printf("MEAN\t%f\t%f\t%f\n", getAverageFitness(rels), getAverageGenerations(rels), getAverageActiveNodes(rels));
	printf("MEDIAN\t%f\t%f\t%f\n", getMedianFitness(rels), getMedianGenerations(rels), getMedianActiveNodes(rels));
	printf("----------------------------------------------------\n\n");

	/* restore the original value for the update frequency */
	params->updateFrequency = updateFrequency;

	return rels;
}


DLL_EXPORT struct chromosome* runCGP(struct parameters *params, struct dataSet *data, int numGens) {

	int i;
	int gen;

	/* bestChromo found using runCGP */
	struct chromosome *bestChromo;

	/* arrays of the parents and children */
	struct chromosome **parentChromos;
	struct chromosome **childrenChromos;

	/* storage for chromosomes used by selection scheme */
	struct chromosome **candidateChromos;
	int numCandidateChromos;

	/* error checking */
	if (numGens < 0) {
		printf("Error: %d generations is invalid. The number of generations must be >= 0.\n Terminating CGP-Library.\n", numGens);
		exit(0);
	}

	if (data != NULL && params->numInputs != data->numInputs) {
		printf("Error: The number of inputs specified in the dataSet (%d) does not match the number of inputs specified in the parameters (%d).\n", data->numInputs, params->numInputs);
		printf("Terminating CGP-Library.\n");
		exit(0);
	}

	if (data != NULL && params->numOutputs != data->numOutputs) {
		printf("Error: The number of outputs specified in the dataSet (%d) does not match the number of outputs specified in the parameters (%d).\n", data->numOutputs, params->numOutputs);
		printf("Terminating CGP-Library.\n");
		exit(0);
	}

	/* initialise parent chromosomes */
	parentChromos = (struct chromosome**)malloc(params->mu * sizeof(struct chromosome*));

	for (i = 0; i < params->mu; i++) {
		parentChromos[i] = initialiseChromosome(params);
	}

	/* initialise children chromosomes */
	childrenChromos = (struct chromosome**)malloc(params->lambda * sizeof(struct chromosome*));

	for (i = 0; i < params->lambda; i++) {
		childrenChromos[i] = initialiseChromosome(params);
	}

	/* intilise best chromosome */
	bestChromo = initialiseChromosome(params);

	/* determine the size of the Candidate Chromos based on the evolutionary Strategy */
	if (params->evolutionaryStrategy == '+') {
		numCandidateChromos = params->mu + params->lambda;
	}
	else if (params->evolutionaryStrategy == ',') {
		numCandidateChromos = params->lambda;
	}
	else {
		printf("Error: the evolutionary strategy '%c' is not known.\nTerminating CGP-Library.\n", params->evolutionaryStrategy);
		exit(0);
	}

	/* initialise the candidateChromos */
	candidateChromos = (struct chromosome**)malloc(numCandidateChromos * sizeof(struct chromosome*));

	for (i = 0; i < numCandidateChromos; i++) {
		candidateChromos[i] = initialiseChromosome(params);
	}

	/* set fitness of the parents */
	for (i = 0; i < params->mu; i++) {
		setChromosomeFitness(params, parentChromos[i], data);
	}

	/* show the user whats going on */
	if (params->updateFrequency != 0) {
		printf("\n-- Starting CGP --\n\n");
		printf("Gen\tfitness\n");
	}

	/* for each generation */
	for (gen = 0; gen < numGens; gen++) {

		/* set fitness of the children of the population */
		#pragma omp parallel for default(none), shared(params, childrenChromos,data), schedule(dynamic), num_threads(params->numThreads)
		for (i = 0; i < params->lambda; i++) {
			setChromosomeFitness(params, childrenChromos[i], data);
		}

		/* get best chromosome */
		getBestChromosome(parentChromos, childrenChromos, params->mu, params->lambda, bestChromo);

		/* check termination conditions */
		if (getChromosomeFitness(bestChromo) <= params->targetFitness) {

			if (params->updateFrequency != 0) {
				printf("%d\t%f - Solution Found\n", gen, bestChromo->fitness);
			}

			break;
		}

		/* display progress to the user at the update frequency specified */
		if (params->updateFrequency != 0 && (gen % params->updateFrequency == 0 || gen >= numGens - 1) ) {
			printf("%d\t%f\n", gen, bestChromo->fitness);
		}

		/*
			Set the chromosomes which will be used by the selection scheme
			dependant upon the evolutionary strategy. i.e. '+' all are used
			by the selection scheme, ',' only the children are.
		*/
		if (params->evolutionaryStrategy == '+') {

			/*
				Note: the children are placed before the parents to
				ensure 'new blood' is always selected over old if the
				fitness are equal.
			*/

			for (i = 0; i < numCandidateChromos; i++) {

				if (i < params->lambda) {
					copyChromosome(candidateChromos[i], childrenChromos[i] );
				}
				else {
					copyChromosome(candidateChromos[i], parentChromos[i - params->lambda] );
				}
			}
		}
		else if (params->evolutionaryStrategy == ',') {

			for (i = 0; i < numCandidateChromos; i++) {
				copyChromosome(candidateChromos[i], childrenChromos[i] );
			}
		}

		/* select the parents from the candidateChromos */
		params->selectionScheme(params, parentChromos, candidateChromos, params->mu, numCandidateChromos);

		/* create the children from the parents */
		params->reproductionScheme(params, parentChromos, childrenChromos, params->mu, params->lambda);
	}

	/* deal with formatting for displaying progress */
	if (params->updateFrequency != 0) {
		printf("\n");
	}

	/* copy the best best chromosome */
	bestChromo->generation = gen;
	/*copyChromosome(chromo, bestChromo);*/

	/* free parent chromosomes */
	for (i = 0; i < params->mu; i++) {
		freeChromosome(parentChromos[i]);
	}
	free(parentChromos);

	/* free children chromosomes */
	for (i = 0; i < params->lambda; i++) {
		freeChromosome(childrenChromos[i]);
	}
	free(childrenChromos);

	/* free the used chromosomes and population */
	for (i = 0; i < numCandidateChromos; i++) {
		freeChromosome(candidateChromos[i]);
	}
	free(candidateChromos);

	return bestChromo;
}

/*
	returns a pointer to the fittest chromosome in the two arrays of chromosomes

	loops through parents and then the children in order for the children to always be selected over the parents
*/
static void getBestChromosome(struct chromosome **parents, struct chromosome **children, int numParents, int numChildren, struct chromosome *best) {

	int i;
	struct chromosome *bestChromoSoFar;

	bestChromoSoFar = parents[0];

	for (i = 1; i < numParents; i++) {

		if (parents[i]->fitness <= bestChromoSoFar->fitness) {
			bestChromoSoFar = parents[i];
		}
	}

	for (i = 0; i < numChildren; i++) {

		if (children[i]->fitness <= bestChromoSoFar->fitness) {
			bestChromoSoFar = children[i];
		}
	}

	copyChromosome(best, bestChromoSoFar);
}


/*
	copies the contents of funcSetSrc to funcSetDest
*/
static void copyFunctionSet(struct functionSet *funcSetDest, struct functionSet *funcSetSrc) {

	int i;

	funcSetDest->numFunctions = funcSetSrc->numFunctions;

	for (i = 0; i < funcSetDest->numFunctions; i++) {
		strncpy(funcSetDest->functionNames[i], funcSetSrc->functionNames[i], FUNCTIONNAMELENGTH);
		funcSetDest->functions[i] = funcSetSrc->functions[i];
		funcSetDest->maxNumInputs[i] = funcSetSrc->maxNumInputs[i];
	}
}


/*
	copys the contents from the src node into dest node.
*/
static void copyNode(struct node *nodeDest, struct node *nodeSrc) {

	int i;

	/* copy the node's function */
	nodeDest->function = nodeSrc->function;

	/* copy active flag */
	nodeDest->active = nodeSrc->active;

	/* copy the node arity */
	nodeDest->maxArity = nodeSrc->maxArity;
	nodeDest->actArity = nodeSrc->actArity;

	/* copy the nodes inputs and connection weights */
	for (i = 0; i < nodeSrc->maxArity; i++) {
		nodeDest->inputs[i] = nodeSrc->inputs[i];
		nodeDest->weights[i] = nodeSrc->weights[i];
	}
}


/*
	mutate Random parent reproduction method.
*/
static void mutateRandomParent(struct parameters *params, struct chromosome **parents, struct chromosome **children, int numParents, int numChildren) {

	int i;

	/* for each child */
	for (i = 0; i < numChildren; i++) {

		/* set child as clone of random parent */
		copyChromosome(children[i], parents[randInt(numParents)]);

		/* mutate newly cloned child */
		mutateChromosome(params, children[i]);
	}
}


/*
	Selection scheme which selects the fittest members of the population
	to be the parents.

	The candidateChromos contains the current children followed by the
 	current parents. This means that using a stable sort to order
	candidateChromos results in children being selected over parents if
	their fitnesses are equal. A desirable property in CGP to facilitate
	neutral genetic drift.
*/
static void selectFittest(struct parameters *params, struct chromosome **parents, struct chromosome **candidateChromos, int numParents, int numCandidateChromos ) {

	int i;

	sortChromosomeArray(candidateChromos, numCandidateChromos);

	for (i = 0; i < numParents; i++) {
		copyChromosome(parents[i], candidateChromos[i]);
	}
}



/*
	returns a pointer to an initialised node. Initialised means that necessary
	memory has been allocated and values set.
*/
static struct node *initialiseNode(int numInputs, int numNodes, int arity, int numFunctions, double connectionWeightRange, double recurrentConnectionProbability, int nodePosition) {

	struct node *n;
	int i;

	/* allocate memory for node */
	n = (struct node*)malloc(sizeof(struct node));

	/* allocate memory for the node's inputs and connection weights */
	n->inputs = (int*)malloc(arity * sizeof(int));
	n->weights = (double*)malloc(arity * sizeof(double));

	/* set the node's function */
	n->function = getRandomFunction(numFunctions);

	/* set as active by default */
	n->active = 1;

	/* set the nodes inputs and connection weights */
	for (i = 0; i < arity; i++) {
		n->inputs[i] = getRandomNodeInput(numInputs, numNodes, nodePosition, recurrentConnectionProbability);
		n->weights[i] = getRandomConnectionWeight(connectionWeightRange);
	}

	/* set the output of the node to zero*/
	n->output = 0;

	/* set the arity of the node */
	n->maxArity = arity;

	return n;
}


/*
	Free memory associated with given node
*/
static void freeNode(struct node *n) {

	/* attempt to prevent user double freeing */
	if (n == NULL) {
		printf("Warning: double freeing of node prevented.\n");
		return;
	}

	free(n->inputs);
	free(n->weights);
	free(n);
}

/*
	returns a random connection weight value
*/
static double getRandomConnectionWeight(double weightRange) {
	return (randDecimal() * 2 * weightRange) - weightRange;
}

/*
	returns a random function index
*/
static int getRandomFunction(int numFunctions) {

	/* check that funcSet contains functions */
	if (numFunctions < 1) {
		printf("Error: cannot assign the function gene a value as the Function Set is empty.\nTerminating CGP-Library.\n");
		exit(0);
	}

	return randInt(numFunctions);
}

/*
	returns a random input for the given node
*/
static int getRandomNodeInput(int numChromoInputs, int numNodes, int nodePosition, double recurrentConnectionProbability) {

	int input;

	/* pick any ahdead nodes or the node itself */
	if (randDecimal() < recurrentConnectionProbability) {
		input = randInt(numNodes - nodePosition) + nodePosition + 1;
	}
	/* pick any previous node including inputs */
	else {
		input = randInt(numChromoInputs + nodePosition);
	}

	return input;
}


/*
	returns a random chromosome output
*/
static int getRandomChromosomeOutput(int numInputs, int numNodes, int shortcutConnections) {

	int output;

	if ( shortcutConnections == 1)
		output = randInt(numInputs + numNodes);
	else
		output = randInt(numNodes) + numInputs;


	return output;
}


/*
	Node function add. Returns the sum of all the inputs.
*/
static double _add(const int numInputs, const double *inputs, const double *connectionWeights) {

	int i;
	double sum = inputs[0];

	for (i = 1; i < numInputs; i++) {
		sum += inputs[i];
	}

	return sum;
}

/*
	Node function sub. Returns the first input minus all remaining inputs.
*/
static double _sub(const int numInputs, const double *inputs, const double *connectionWeights) {

	int i;
	double sum = inputs[0];

	for (i = 1; i < numInputs; i++) {
		sum -= inputs[i];
	}

	return sum;
}


/*
	Node function mul. Returns the multiplication of all the inputs.
*/
static double _mul(const int numInputs, const double *inputs, const double *connectionWeights) {

	int i;
	double multiplication = inputs[0];

	for (i = 1; i < numInputs; i++) {
		multiplication *= inputs[i];
	}

	return multiplication;
}


/*
	Node function div. Returns the first input divided by the second input divided by the third input etc
*/
static double _divide(const int numInputs, const double *inputs, const double *connectionWeights) {

	int i;
	double divide = inputs[0];

	for (i = 1; i < numInputs; i++) {
		divide /= inputs[i];
	}

	return divide;
}


/*
	Node function abs. Returns the absolute of the first input
*/
static double _absolute(const int numInputs, const double *inputs, const double *connectionWeights) {

	return fabs(inputs[0]);
}


/*
	Node function sqrt.  Returns the square root of the first input
*/
static double _squareRoot(const int numInputs, const double *inputs, const double *connectionWeights) {

	return sqrt(inputs[0]);
}


/*
	Node function squ.  Returns the square of the first input
*/
static double _square(const int numInputs, const double *inputs, const double *connectionWeights) {

	return pow(inputs[0], 2);
}


/*
	Node function cub.  Returns the cube of the first input
*/
static double _cube(const int numInputs, const double *inputs, const double *connectionWeights) {

	return pow(inputs[0], 3);
}


/*
	Node function power.  Returns the first output to the power of the second
*/
static double _power(const int numInputs, const double *inputs, const double *connectionWeights) {

	return pow(inputs[0], inputs[1]);
}

/*
	Node function exp.  Returns the exponential of the first input
*/
static double _exponential(const int numInputs, const double *inputs, const double *connectionWeights) {

	return exp(inputs[0]);
}


/*
	Node function sin.  Returns the sine of the first input
*/
static double _sine(const int numInputs, const double *inputs, const double *connectionWeights) {

	return sin(inputs[0]);
}

/*
	Node function cos.  Returns the cosine of the first input
*/
static double _cosine(const int numInputs, const double *inputs, const double *connectionWeights) {

	return cos(inputs[0]);
}

/*
	Node function tan.  Returns the tangent of the first input
*/
static double _tangent(const int numInputs, const double *inputs, const double *connectionWeights) {

	return tan(inputs[0]);
}


/*
	Node function one.  Always returns 1
*/
static double _constOne(const int numInputs, const double *inputs, const double *connectionWeights) {
	return 1;
}

/*
	Node function one.  Always returns 0
*/
static double _constZero(const int numInputs, const double *inputs, const double *connectionWeights) {
	return 0;
}

/*
	Node function one.  Always returns PI
*/
static double _constPI(const int numInputs, const double *inputs, const double *connectionWeights) {
	return M_PI;
}


/*
	Node function rand.  Returns a random number between minus one and positive one
*/
static double _randFloat(const int numInputs, const double *inputs, const double *connectionWeights) {

	double random;

	random = (double)rand() / (double)RAND_MAX;

	random = (random * 2) - 1;

	return random;
}

/*
	Node function and. logical AND, returns '1' if all inputs are '1'
	else, '0'
*/
static double _and(const int numInputs, const double *inputs, const double *connectionWeights) {

	int i;

	for (i = 0; i < numInputs; i++) {

		if (inputs[i] == 0) {
			return 0;
		}
	}

	return 1;
}

/*
	Node function and. logical NAND, returns '0' if all inputs are '1'
	else, '1'
*/
static double _nand(const int numInputs, const double *inputs, const double *connectionWeights) {

	int i;

	for (i = 0; i < numInputs; i++) {

		if (inputs[i] == 0) {
			return 1;
		}
	}

	return 0;
}


/*
	Node function or. logical OR, returns '0' if all inputs are '0'
	else, '1'
*/
static double _or(const int numInputs, const double *inputs, const double *connectionWeights) {

	int i;

	for (i = 0; i < numInputs; i++) {

		if (inputs[i] == 1) {
			return 1;
		}
	}

	return 0;
}


/*
	Node function nor. logical NOR, returns '1' if all inputs are '0'
	else, '0'
*/
static double _nor(const int numInputs, const double *inputs, const double *connectionWeights) {

	int i;

	for (i = 0; i < numInputs; i++) {

		if (inputs[i] == 1) {
			return 0;
		}
	}

	return 1;
}


/*
	Node function xor. logical XOR, returns '1' iff one of the inputs is '1'
	else, '0'. AKA 'one hot'.
*/
static double _xor(const int numInputs, const double *inputs, const double *connectionWeights) {

	int i;
	int numOnes = 0;
	int out;

	for (i = 0; i < numInputs; i++) {

		if (inputs[i] == 1) {
			numOnes++;
		}

		if (numOnes > 1) {
			break;
		}
	}

	if (numOnes == 1) {
		out = 1;
	}
	else {
		out = 0;
	}

	return out;
}

/*
	Node function xnor. logical XNOR, returns '0' iff one of the inputs is '1'
	else, '1'.
*/
static double _xnor(const int numInputs, const double *inputs, const double *connectionWeights) {

	int i;
	int numOnes = 0;
	int out;

	for (i = 0; i < numInputs; i++) {

		if (inputs[i] == 1) {
			numOnes++;
		}

		if (numOnes > 1) {
			break;
		}
	}

	if (numOnes == 1) {
		out = 0;
	}
	else {
		out = 1;
	}

	return out;
}

/*
	Node function not. logical NOT, returns '1' if first input is '0', else '1'
*/
static double _not(const int numInputs, const double *inputs, const double *connectionWeights) {

	double out;

	if (inputs[0] == 0) {
		out = 1;
	}
	else {
		out = 0;
	}

	return out;
}


/*
	Node function wire. simply acts as a wire returning the first input
*/
static double _wire(const int numInputs, const double *inputs, const double *connectionWeights) {

	double out;

	out = inputs[0];

	return out;
}


/*
	Node function sigmoid. returns the sigmoid of the sum of weighted inputs.
	The specific sigmoid function used in the logistic function.
	range: [0,1]
*/
static double _sigmoid(const int numInputs, const double *inputs, const double *connectionWeights) {

	double weightedInputSum;
	double out;

	weightedInputSum = sumWeigtedInputs(numInputs, inputs, connectionWeights);

	out = 1 / (1 + exp(-weightedInputSum));

	return out;
}

/*
	Node function Gaussian. returns the Gaussian of the sum of weighted inputs.
	range: [0,1]
*/
static double _gaussian(const int numInputs, const double *inputs, const double *connectionWeights) {

	double weightedInputSum;
	double out;

	int centre = 0;
	int width = 1;

	weightedInputSum = sumWeigtedInputs(numInputs, inputs, connectionWeights);

	out = exp(-(pow(weightedInputSum - centre, 2)) / (2 * pow(width, 2)));

	return out;
}


/*
	Node function step. returns the step function of the sum of weighted inputs.
	range: [0,1]
*/
static double _step(const int numInputs, const double *inputs, const double *connectionWeights) {

	double weightedInputSum;
	double out;

	weightedInputSum = sumWeigtedInputs(numInputs, inputs, connectionWeights);

	if (weightedInputSum < 0) {
		out = 0;
	}
	else {
		out = 1;
	}

	return out;
}


/*
	Node function step. returns the step function of the sum of weighted inputs.
	range: [-1,1]
*/
static double _softsign(const int numInputs, const double *inputs, const double *connectionWeights) {

	double weightedInputSum;
	double out;

	weightedInputSum = sumWeigtedInputs(numInputs, inputs, connectionWeights);

	out = weightedInputSum / (1 + fabs(weightedInputSum));

	return out;
}


/*
	Node function tanh. returns the tanh function of the sum of weighted inputs.
	range: [-1,1]
*/
static double _hyperbolicTangent(const int numInputs, const double *inputs, const double *connectionWeights) {

	double weightedInputSum;
	double out;

	weightedInputSum = sumWeigtedInputs(numInputs, inputs, connectionWeights);

	out = tanh(weightedInputSum);

	return out;
}


/*
	Returns the sum of the weighted inputs.
*/
static double sumWeigtedInputs(const int numInputs, const double *inputs, const double *connectionWeights) {

	int i;
	double weightedSum = 0;

	for (i = 0; i < numInputs; i++) {
		weightedSum += (inputs[i] * connectionWeights[i]);
	}

	return weightedSum;
}





/*
	The default fitness function used by CGP-Library.
	Simply assigns an error of the sum of the absolute differences between the target and actual outputs for all outputs over all samples
*/
static double supervisedLearning(struct parameters *params, struct chromosome *chromo, struct dataSet *data) {

	int i, j;
	double error = 0;

	/* error checking */
	if (getNumChromosomeInputs(chromo) != getNumDataSetInputs(data)) {
		printf("Error: the number of chromosome inputs must match the number of inputs specified in the dataSet.\n");
		printf("Terminating CGP-Library.\n");
		exit(0);
	}

	if (getNumChromosomeOutputs(chromo) != getNumDataSetOutputs(data)) {
		printf("Error: the number of chromosome outputs must match the number of outputs specified in the dataSet.\n");
		printf("Terminating CGP-Library.\n");
		exit(0);
	}

	/* for each sample in data */
	for (i = 0 ; i < getNumDataSetSamples(data); i++) {

		/* calculate the chromosome outputs for the set of inputs  */
		executeChromosome(chromo, getDataSetSampleInputs(data, i));

		/* for each chromosome output */
		for (j = 0; j < getNumChromosomeOutputs(chromo); j++) {

			error += fabs(getChromosomeOutput(chromo, j) - getDataSetSampleOutput(data, i, j));
		}
	}

	return error;
}


/*
	returns a random decimal between [0,1]
*/
static double randDecimal(void) {
	return (double)rand() / (double)RAND_MAX;
}

/*
	sort int array using qsort
*/
static void sortIntArray(int *array, const int length) {

	qsort(array, length, sizeof(int), cmpInt);
}

/*
	used by qsort in sortIntArray
*/
static int cmpInt(const void * a, const void * b) {
	return ( *(int*)a - * (int*)b );
}


/*
	sort double array using qsort
*/
static void sortDoubleArray(double *array, const int length) {

	qsort(array, length, sizeof(double), cmpDouble);
}

/*
	used by qsort in sortDoubleArray
*/
static int cmpDouble(const void * a, const void * b) {

	if ( *(double*)a < * (double*)b) {
		return -1;
	}
	if ( *(double*)a == *(double*)b ) {
		return 0;
	}
	else {
		return 1;
	}
}



/*
	Prints the current functions in the function set to
	the terminal.
*/
static void printFunctionSet(struct parameters *params) {

	int i;

	printf("Function Set:");

	for (i = 0; i < params->funcSet->numFunctions; i++) {
		printf(" %s", params->funcSet->functionNames[i]);
	}

	printf(" (%d)\n", params->funcSet->numFunctions);
}


/*
	random integer between zero and n without modulo bias.
	adapted from: http://zuttobenkyou.wordpress.com/2012/10/18/generating-random-numbers-without-modulo-bias/
*/
static int randInt(int n) {

	int x;
	int randLimit;
	int randExcess;

	if (n == 0) {
		return 0;
	}

	randExcess = (RAND_MAX % n) + 1;
	randLimit = RAND_MAX - randExcess;

	do {
		x = rand();
	}
	while (x > randLimit);

	return x % n;
}
