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
	struct chromosome *chromo = NULL;

	int numInputs = 1;
	int numNodes = 15;
	int numOutputs = 1;
	int nodeArity = 2;

	int numGens = 100000;
	int updateFrequency = 500;

	double recurrentConnectionProbability = 0.10;

	params = initialiseParameters(numInputs, numNodes, numOutputs, nodeArity);

	addNodeFunction(params, "add,sub,mul,div");

	setUpdateFrequency(params, updateFrequency);

	setRecurrentConnectionProbability(params, recurrentConnectionProbability);

	printParameters(params);

	trainingData = initialiseDataSetFromFile("./dataSets/fibonacci.data");

	chromo = runCGP(params, trainingData, numGens);

	printChromosome(chromo, 0);

	freeDataSet(trainingData);
	freeChromosome(chromo);
	freeParameters(params);

	return 0;
}
