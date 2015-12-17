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
#include <math.h>
#include "../src/cgp.h"


double fitnessFunction(struct parameters *params, struct chromosome *chromo, struct dataSet *data) {

	int i;
	double inputs[1] = {0.5};
	int numExec = 10000;

	for (i = 0; i < numExec; i++) {
		executeChromosome(chromo, inputs);
	}

	return 10; /*(double)(rand() % 10000); squareError / (getNumDataSetSamples(data) * getNumDataSetOutputs(data));*/
}


int main(void) {

	struct parameters *params = NULL;

	int numInputs = 1;
	int numNodes = 100;
	int numOutputs = 1;
	int arity = 2;

	int gens = 100;
	int runs = 100;

	params = initialiseParameters(numInputs, numNodes, numOutputs, arity);
	setRandomNumberSeed(123456789);
	addNodeFunction(params, "add,sub,mul,div,sin");
	setMutationRate(params, 1.0);
	setCustomFitnessFunction(params, fitnessFunction, "FF");
	setNumThreads(params, 4);

	printParameters(params);

	repeatCGP(params, NULL, gens, runs);
	/*runCGP(params, NULL, gens);*/

	freeParameters(params);

	return 0;
}









