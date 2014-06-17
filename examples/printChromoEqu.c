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
#include "../include/cgp.h"

int main(void){

	struct parameters *params = NULL;
	
	struct chromosome *chromo = NULL;

	int numInputs = 2;
	int numNodes = 5;
	int numOutputs = 1;
	int nodeArity = 2;

	

	params = initialiseParameters(numInputs, numNodes, numOutputs, nodeArity);

	addNodeFunction(params, "add,sin,exp");

	
	chromo = initialiseChromosome(params);
	
	printChromosome(chromo, 0);

	saveChromosomeLatex(chromo, 0, "Tmp.tex");
	
	
	freeChromosome(chromo);
	freeParameters(params);

	return 0;
}
