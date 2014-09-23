#!/usr/bin/env python

#    This file is part of CGP-Library
#    Copyright (c) Andrew James Turner 2014 (andrew.turner@york.ac.uk)
#
#    CGP-Library is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Lesser General Public License as published
#    by the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    CGP-Library is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public License
#    along with CGP-Library.  If not, see <http://www.gnu.org/licenses/>.
 

import cgp

numInputs = 1;
numNodes = 15;
numOutputs = 1;
nodeArity = 2;

numGens = 10000;
targetFitness = 0.1;
updateFrequency = 500;

params = cgp.initialiseParameters(numInputs, numNodes, numOutputs, nodeArity);

cgp.addNodeFunction(params, "add,sub,mul,div,sin");

cgp.setTargetFitness(params, targetFitness);

cgp.setUpdateFrequency(params, updateFrequency);

cgp.printParameters(params);

trainingData = cgp.initialiseDataSetFromFile("../dataSets/symbolic.data");

chromo = cgp.runCGP(params, trainingData, numGens);

cgp.printChromosome(chromo, 0);

cgp.freeDataSet(trainingData);
cgp.freeChromosome(chromo);
cgp.freeParameters(params);

