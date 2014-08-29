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
#include <math.h>
#include "../src/cgp.h"

#define POPULATIONSIZE 5
#define NUMINPUTS 1
#define NUMNODES 8
#define NUMOUTPUTS 1
#define ARITY 2

int main(void){

	struct parameters *params = NULL;
	struct chromosome *chromoA = NULL;
	struct chromosome *chromoB = NULL;
	struct chromosome *chromoC = NULL;
	struct dataSet *trainingData = NULL;

	double testInputs[NUMINPUTS];

	params = initialiseParameters(NUMINPUTS, NUMNODES, NUMOUTPUTS, ARITY);
	addNodeFunction(params, "add,sub,mul,sq,cube,sin");

	trainingData = initialiseDataSetFromFile("./dataSets/symbolic.data");

	chromoA = initialiseChromosome(params);
	chromoB = initialiseChromosome(params);

	setChromosomeFitness(params, chromoA, trainingData);

	mutateChromosome(params,chromoA);

	copyChromosome(chromoB, chromoA);

	removeInactiveNodes(chromoB);

	printf("chromoA with inactive nodes.\n");
	printChromosome(chromoA, 0);

	printf("chromoB without inactive nodes.\n");
	printChromosome(chromoB, 1);

	saveChromosome(chromoB, "chromoB.chromo");

	chromoC = initialiseChromosomeFromFile("chromoB.chromo");

	testInputs[0] = 3;

	executeChromosome(chromoC, testInputs);

	printf("Applied input: %f\n", testInputs[0]);
	printf("Generated output: %f\n", getChromosomeOutput(chromoC, 0));

	freeChromosome(chromoA);
	freeChromosome(chromoB);
	freeChromosome(chromoC);
	freeDataSet(trainingData);
	freeParameters(params);

	return 0;
}
