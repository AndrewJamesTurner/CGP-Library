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
#include "../include/cgp.h"

float hypotenuse(const int numInputs, const float *inputs, const float *connectionWeights){

	int i;
	float sumOfSqrs = 0;
	float hypt;

	for(i=0; i<numInputs; i++){
		sumOfSqrs += powf(inputs[i], 2);
	}

	hypt = sqrtf(sumOfSqrs);

	return hypt;
}


int main(void){

	struct parameters *params = NULL;

	int numInputs = 2;
	int numNodes = 10;
	int numOutputs = 1;
	int arity = 3;

	params = initialiseParameters(numInputs, numNodes, numOutputs, arity);

	addNodeFunction(params, "add,sub");

	addNodeFunctionCustom(params, hypotenuse, "hypt");

	printParameters(params);

	freeParameters(params);

	return 0;
}









