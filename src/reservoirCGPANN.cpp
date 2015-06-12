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
//#include <stdio.h>
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

		for(j=0; j<getNumChromosomeActiveNodes(chromo); j++){
		
			if(isNodeActive(chromo, j)){

				states(i,j) = getChromosomeNodeValue(chromo, j);

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


	cout << numActiveNodes << "\n\n";

	// calculate the output weights
	Wout = states.jacobiSvd(ComputeThinU | ComputeThinV).solve(desiredOutputs);


	

	cout << Wout << endl;


	exit(0);
	return 0;
}


int main(void){

	struct parameters* params = NULL;
	struct chromosome* chromo = NULL;
	struct dataSet* trainingData = NULL;

	int numInputs = 1;
	int numNodes = 10;
	int numOutputs = 1;
	int arity = 2;

	double testInputs[1];
	testInputs[0] = 1;

	//getChromosomeNodeValue

	params = initialiseParameters(numInputs, numNodes, numOutputs, arity);
	addNodeFunction(params, "sig");
	chromo = initialiseChromosome(params);
	
	trainingData = initialiseDataSetFromFile("src/sin2saw.csv");

	setCustomFitnessFunction(params, reservoirFitnesses, "reservoir");
	setShortcutConnections(params, 0);

	printParameters(params);

	setChromosomeFitness(params, chromo, trainingData);
	
	freeDataSet(trainingData);
	freeParameters(params);
	freeChromosome(chromo);

	return 0;
}









