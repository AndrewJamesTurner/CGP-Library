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
#include "../src/cgp.h"

int main(void){

	struct parameters *params = NULL;
	struct dataSet *trainingData = NULL;
	struct results *rels = NULL;
	struct chromosome *chromo = NULL;

	int numInputs = 1;
	int numNodes = 15;
	int numOutputs = 1;
	int nodeArity = 2;

	int numGens = 10000;
	int numRuns = 10;
	
	double targetFitness = 0.1;
	int updateFrequency = 500;

	double averageFitness;

	params = initialiseParameters(numInputs, numNodes, numOutputs, nodeArity);

	addNodeFunction(params, "add,sub,mul,div,sin");

	setTargetFitness(params, targetFitness);

	setUpdateFrequency(params, updateFrequency);

	trainingData = initialiseDataSetFromFile("./examples/symbolic.data");

	rels = repeatCGP(params, trainingData, numGens, numRuns);

	averageFitness = getAverageFitness(rels);

	printf("The average chromosome fitness is: %f\n", averageFitness);

	chromo = getChromosome(rels, 4);
	
	printf("The best chromosome found on run 4:\n");
	
	printChromosome(chromo, 0);

	saveResults(rels, "results.csv");

	freeDataSet(trainingData);
	freeChromosome(chromo);
	freeResults(rels);
	freeParameters(params);

	return 0;
}
