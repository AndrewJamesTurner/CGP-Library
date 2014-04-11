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
#include <math.h>
#include "../include/cgp.h"

float sinWave(struct parameters *params, struct chromosome *chromo, struct dataSet *data){

	float i;
		
	float error = 0;
	float range = 6;
	float stepSize = 0.5;
	
	float inputs[1];

	for(i=0; i<range; i += stepSize){

		inputs[0] = i;

		executeChromosome(chromo, inputs);

		error += fabs(getChromosomeOutput(chromo, 0) - sin(i));
	}

	return error;
}


int main(void){

	struct parameters *params = NULL;
	struct chromosome *chromo = NULL;

	int numInputs = 1;
	int numNodes = 20;
	int numOutputs = 1;
	int nodeArity = 5;

	int numGens = 25000;
	float targetFitness = 0.5;
	int updateFrequency = 500;
	
	float weightRange = 5;

	params = initialiseParameters(numInputs, numNodes, numOutputs, nodeArity);

	setTargetFitness(params, targetFitness);

	setUpdateFrequency(params, updateFrequency);
	
	setConnectionWeightRange(params, weightRange);

	setFitnessFunction(params, sinWave, "sinWave");

	addNodeFunction(params, "tanh, softsign");

	chromo = runCGP(params, NULL, numGens);

	printChromosome(chromo);

	freeChromosome(chromo);
	freeParameters(params);

	return 0;
}
