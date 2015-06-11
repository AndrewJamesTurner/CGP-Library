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
#include "cgp.h" 


double meanSquareError(struct parameters *params, struct chromosome *chromo, struct dataSet *data){

	int i,j;
	double squareError = 0;

	if(getNumChromosomeInputs(chromo) !=getNumDataSetInputs(data)){
		printf("Error: the number of chromosome inputs must match the number of inputs specified in the dataSet.\n");
		printf("Terminating.\n");
		exit(0);
	}

	if(getNumChromosomeOutputs(chromo) != getNumDataSetOutputs(data)){
		printf("Error: the number of chromosome outputs must match the number of outputs specified in the dataSet.\n");
		printf("Terminating.\n");
		exit(0);
	}

	for(i=0; i<getNumDataSetSamples(data); i++){

		executeChromosome(chromo, getDataSetSampleInputs(data, i));

		for(j=0; j<getNumChromosomeOutputs(chromo); j++){

			squareError += pow(getDataSetSampleOutput(data,i,j) - getChromosomeOutput(chromo,j), 2);
		}
	}

	return squareError / (getNumDataSetSamples(data) * getNumDataSetOutputs(data));
}


int main(void){

	struct parameters *params = NULL;

	int numInputs = 1;
	int numNodes = 20;
	int numOutputs = 1;
	int arity = 2;

	params = initialiseParameters(numInputs, numNodes, numOutputs, arity);

	setCustomFitnessFunction(params, meanSquareError, "MSE");

	printParameters(params);

	freeParameters(params);

	return 0;
}









