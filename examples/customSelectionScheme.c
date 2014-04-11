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
#include "../include/cgp.h"

void tournament(struct parameters *params, struct chromosome **parents, struct chromosome **candidateChromos, int numParents, int numCandidateChromos){

	int i;
	
	struct chromosome *candidateA;
	struct chromosome *candidateB;
	
	for(i=0; i<numParents; i++){
			
		candidateA = candidateChromos[rand() % numCandidateChromos];
		candidateB = candidateChromos[rand() % numCandidateChromos];
		
		if(getChromosomeFitness(candidateA) <= getChromosomeFitness(candidateB)){
			copyChromosome(parents[i], candidateA);
		}
		else{
			copyChromosome(parents[i], candidateB);
		}		
	}
}

int main(void){
	
	struct parameters *params = NULL;
	struct dataSet *trainingData = NULL;
	struct chromosome *chromo = NULL;
	
	int numInputs = 1;
	int numNodes = 50;
	int numOutputs = 1;
	int arity = 2;
	
	float targetFitness = 0.1;
	int updateFrequency = 1000;
	int numGens = 10000;
	
	params = initialiseParameters(numInputs, numNodes, numOutputs, arity);
	
	addNodeFunction(params, "add,sub,mul,div,sin");
	
	setTargetFitness(params, targetFitness);

    setUpdateFrequency(params, updateFrequency); 
	
	setSelectionScheme(params, tournament, "tournament");
	
	trainingData = initialiseDataSetFromFile("./examples/symbolic.data");
	
	chromo = runCGP(params, trainingData, numGens);
	
	freeChromosome(chromo);
	freeDataSet(trainingData);
	freeParameters(params);
	
	return 0;
}



















