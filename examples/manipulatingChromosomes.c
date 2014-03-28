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

#define POPULATIONSIZE 5

int main(void){
	
	int i;
	
	struct parameters *params = NULL;
	struct chromosome *population[POPULATIONSIZE];
	struct chromosome *fittestChromosome;
	struct dataSet *trainingData = NULL;
	
	int numInputs = 1;
	int numNodes = 30;
	int numOutputs = 1;
	int nodeArity = 2;
	
	int targetFitness = 1;
	
	int maxGens = 5000;
	int gen;
	
	float testInputs[1];
	float testOutputs[1];
	
	params = initialiseParameters(numInputs, numNodes, numOutputs, nodeArity);
	
	addNodeFunction(params, "add,sub,mul,sq,cube,sin");
	setTargetFitness(params, targetFitness);
	setMutationType(params, "point");
	setMutationRate(params, 0.15);
	
	trainingData = initialiseDataSetFromFile("./examples/symbolic.data");
	
	for(i=0; i<POPULATIONSIZE; i++){
		population[i] = initialiseChromosome(params);
	}
	
	fittestChromosome = initialiseChromosome(params);
	
	// for the number of allowed generations
	for(gen=0; gen<maxGens; gen++){
		
		// set the fitnesses of the population of chromosomes
		for(i=0; i<POPULATIONSIZE; i++){
			setChromosomeFitness(params, population[i], trainingData);
		}
		
		// copy over the last chromosome to fittestChromosome
		copyChromosome(fittestChromosome, population[POPULATIONSIZE - 1]);
		
		// for all chromosomes except the last
		for(i=0; i<POPULATIONSIZE-1; i++){
			
			// copy ith chromosome to fittestChromosome if fitter
			if(getChromosomeFitness(population[i]) < getChromosomeFitness(fittestChromosome)){
				copyChromosome(fittestChromosome, population[i]);
			}
		}
		
		// termination condition
		if(getChromosomeFitness(fittestChromosome) <= targetFitness){
			break;
		}
				
		// set the first member of the population to be the fittest chromosome
		copyChromosome(population[0], fittestChromosome);
		
		// set remaining member of the population to be mutations of the
		// fittest chromosome
		for(i=1; i<POPULATIONSIZE; i++){
			
			copyChromosome(population[i], population[0]);
			mutateChromosome(params, population[i]);
		}
	}
	
	
	printf("Fittest chromosome Found has fitness: %f\n", getChromosomeFitness(fittestChromosome));
	
	printf("\nFittest chromosome");
	printChromosome(fittestChromosome);
	
	removeInactiveNodes(fittestChromosome);
		
	saveChromosome(fittestChromosome, "fittestChromosome.chromo");
	
	freeChromosome(fittestChromosome);
	
	fittestChromosome = initialiseChromosomeFromFile("fittestChromosome.chromo");
	
	printf("\nFittest chromosome without inactive nodes loaded from file.");
	printChromosome(fittestChromosome);
	
	testInputs[0] = 3;
	executeChromosome(fittestChromosome, testInputs, testOutputs);
	
	printf("\n\n");
	printf("Chromo Input: %f\n", testInputs[0]);
	printf("Chromo Output: %f\n", testOutputs[0]);
	
	
	for(i=0; i<POPULATIONSIZE; i++){
		freeChromosome(population[i]);
	}
	
	freeChromosome(fittestChromosome);
	freeDataSet(trainingData);
	freeParameters(params);
	
	return 1;
}
