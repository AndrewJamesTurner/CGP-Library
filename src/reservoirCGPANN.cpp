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

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cgp.h" 
#include "../Eigen/Dense"
#include "../Eigen/LU" 

using namespace Eigen;
using namespace std;

double reservoirFitnesses(struct parameters *params, struct chromosome *chromo, struct dataSet *data){

	int i,j,index;; 

	int numActiveNodes = getNumChromosomeActiveNodes(chromo);
	int numSamples = getNumDataSetSamples(data);
	int numOutputs = getNumDataSetOutputs(data);

	MatrixXd states(numSamples,numActiveNodes);
	MatrixXd desiredOutputs(numSamples,numOutputs);
	MatrixXd Wout;
		
	// set the reservoir state
	for(i = 0; i<numSamples; i++){

		executeChromosome(chromo, getDataSetSampleInputs(data, i));

		index = 0;

		for(j=0; j<getNumChromosomeNodes(chromo); j++){
		
			if(isNodeActive(chromo, j)){

				states(i,index) = getChromosomeNodeValue(chromo, j);

				index++;
			}
		}
	}

	// set the desired outputs
	for(i = 0; i < numSamples; i++){
		for( j=0; j<numOutputs; j++){

			desiredOutputs(i,j) = getDataSetSampleOutput(data, i, j);
		}
	}


	// calculate the output weights
	Wout = states.colPivHouseholderQr().solve(desiredOutputs);
	//Wout = states.jacobiSvd(ComputeThinU | ComputeThinV).solve(desiredOutputs);


	// return the error
	return (states*Wout - desiredOutputs).norm() / desiredOutputs.norm();
}


void saveChromosomeBehaviour(struct chromosome *chromo, struct dataSet *data){

	int i,j,index;; 

	FILE *fp;

	int numActiveNodes = getNumChromosomeActiveNodes(chromo);
	int numSamples = getNumDataSetSamples(data);
	int numOutputs = getNumDataSetOutputs(data);
	int numInputs = getNumDataSetInputs(data);

	MatrixXd states(numSamples,numActiveNodes);
	MatrixXd desiredOutputs(numSamples,numOutputs);
	MatrixXd actualOutputs;
	MatrixXd Wout;
		
	// set the reservoir state
	for(i = 0; i<numSamples; i++){

		executeChromosome(chromo, getDataSetSampleInputs(data, i));

		index = 0;

		for(j=0; j<getNumChromosomeNodes(chromo); j++){
		
			if(isNodeActive(chromo, j)){

				states(i,index) = getChromosomeNodeValue(chromo, j);

				index++;
			}
		}
	}

	// set the desired outputs
	for(i = 0; i < numSamples; i++){
		for( j=0; j<numOutputs; j++){
			desiredOutputs(i,j) = getDataSetSampleOutput(data, i, j);
		}
	}

	// calculate the output weights
	//Wout = states.colPivHouseholderQr().solve(desiredOutputs);
	Wout = states.jacobiSvd(ComputeThinU | ComputeThinV).solve(desiredOutputs);

	// calculate the actual outputs
	actualOutputs = states * Wout;

	fp = fopen("tmp.csv", "w");
	fprintf(fp, "inputs,DesiredOutputs,ActualOutputs,\n");
	
	for(i = 0; i<numSamples; i++){
		
		for(j=0; j<numInputs; j++){
			fprintf(fp, "%f,", getDataSetSampleInput(data, i, j));
		}

		for(j=0; j<numOutputs; j++){
			fprintf(fp, "%f,", desiredOutputs(i,j));
		}

		for(j=0; j<numOutputs; j++){
			fprintf(fp, "%f,", actualOutputs(i,j));
		}


		fprintf(fp, "\n");
	}
	
//	fprintf(fp, "Run,Fitness,Generations,Active Nodes\n");
	fclose(fp);
	
	
}



int main(void){

	struct parameters* params = NULL;
	struct chromosome* chromo = NULL;
	struct dataSet* trainingData = NULL;

	int numInputs = 1;
	int numNodes = 100;
	int numOutputs = 1;
	int arity = 5;

	int numGens = 10;


	// set up parameters
	params = initialiseParameters(numInputs, numNodes, numOutputs, arity);
	addNodeFunction(params, "sig"); // softsign
	setConnectionWeightRange(params, 5);
	setCustomFitnessFunction(params, reservoirFitnesses, "reservoir");
	setRecurrentConnectionProbability(params, 0.5);
	setMutationRate(params, 0.05);
	setShortcutConnections(params, 0);
	printParameters(params);
	
	// set up traning data
	trainingData = initialiseDataSetFromFile("src/sin2saw.csv");

	// run CGP
	chromo = runCGP(params, trainingData, numGens);
	
	saveChromosomeBehaviour(chromo, trainingData);
		
	freeDataSet(trainingData);
	freeParameters(params);
	freeChromosome(chromo);

	return 0;
}









