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

/*
	Title: API 
	
	Description of all the CGP-Library functions and structures.	
*/

#ifndef CGPLIB
#define CGPLIB
				
/*
	Title: Structures
	
	Description of all the structures used by CGP-Library
*/		
		
	/*
		variable: parameters
		
		Stores general parameters used by the CGP-Library.
		
		The values stored in parameters are set to default values when initialised and can be altered though setter functions. The parameters variables can not be accessed directly; getters and Setters must be used.
		
		Defaults:
		
			Default parameter values.
		
			> mu:						1
			> lambda:					4
			> evolutionary strategy:	+ (mu+lambda)-ES
			> mutation rate:			0.05
			> connection weight range:	1 +/- 
			> maximum generations:		10000
			> update frequency:			100
			> mutation type:			probabilistic
			> fitness function: 		supervisedLearning
			> selection scheme:			pickHighest
			> reproduction scheme:		mutateRandomParent
			
		See Also:
			<initialiseParameters>, <freeParameters>, <printParameters>
	*/	
	struct parameters;
	
	/*
		variable: chromosome
		
		Stores CGP chromosome instances used by the CGP-Library.
		
		See Also:
			<initialiseChromosome>, <initialiseChromosomeFromFile> <freeChromosome>, <printChromosome>, <executeChromosome>, <mutateChromosome>
		
	*/
	struct chromosome;
	
	/*
		variable: dataSet
		
		Stores data which can be used by fitness functions when calculating a chromosomes fitness. Typically contains input output pairs of data used when applying CGP to supervised learning.
		
		See Also:
			<initialiseDataSetFromFile>, <initialiseDataSetFromArrays>, <freeDataSet>, <printDataSet>
	*/
	struct dataSet;	
		
	/*
		variable: results
		
		
		
	*/
	struct results;
	
	
	
/*
	Title: Parameters Functions
	
	Description of all the functions related to CGP-Library parameters
*/		
	
	/*
		Function: initialiseParameters

		Initialises parameters used throughout the CGP-Library. The inputs describe the structure of the chromosomes created when using <initialiseChromosome>, <runCGP> or <repeatCGP>.

		Parameters:
			numInputs - the number of chromosome inputs required. 
			numNodes - the number of chromosome nodes required.
			numOutputs - the number of chromosome outputs required.
			arity - the arity of each chromosome node required.

		Returns:
			A pointer to an initialised parameters structure.

		Example:
			
			Initialising parameters
			> struct parameters *params;
			>
			> int numInputs = 3;
			> int numNodes = 10;
			> int numOutputs = 2;
			> int nodeArity = 2;		
			>
			> params = initialiseParameters(numInputs, numNodes, numOutputs, nodeArity); 

		See Also:
			<freeParameters>, <printParameters>
	*/
	struct parameters *initialiseParameters(const int numInputs, const int numNodes, const int numOutputs, const int arity);
	
	/*
		Function: freeParameters

		Frees parameters instance.

		Parameters:
			params - pointer to parameters structure.
			
		See Also:
			<initialiseParameters>
	*/
	void freeParameters(struct parameters *params);
	
	/*
		Function: printParameters

		Prints the given parameters to the screen in a human readable format. 

		Parameters:
			params - pointer to parameters structure.
			
	*/
	void printParameters(struct parameters *params);
	
			
	/*
		Function: addNodeFunction

		Adds pre-defined node function(s) to the set of functions made available to chromosome nodes.  
	
		If one function name is given that function is added to the function set. If multiple node function's names are given then each must be separated by a ','.
		
		If a node function name is given which is not recognised, a warning is given and that function is not added to the function set.

		Parameters:
			params - pointer to parameters structure
			functionNames - the name(s) of the function(s) to add to the function set
		
		Node Functions:
		
			mathematical operations 
			
			- add 		- 	summation over all inputs.
			- sub		-	subtracts all but the first input from the first input
			- mul		-	multiplies all of the inputs
			- div		-	divides the first input by the second and then the third etc
			- abs		-	the absolute of the first input
			- sqrt		- 	the square root of the first input
			- sq		-	the square of the first input
			- cube		- 	the cube of the first input
			- exp 		- 	the exponential of the first input
			- sin		- 	the sine of the first input
			- cos 		-	the cosine of the first input
			- tan		-	the tangent of the first input   
			
			logic gates
			
			- and		-	returns '1' if all inputs are '1', else '0'
			- nand		-	returns '0' if all inputs are '1', else, '1'
			- or		-	returns '0' if all inputs are '0', else, '1'
			- nor		-	returns '1' if all inputs are '0', else, '0'
			- xor		-	returns '1' if only one of the inputs is '1', else, '0'
			- xnor		-	returns '0' if only one of the inputs is '1', else, '1'
			- not		-	returns '1' if first input is '0', else '1'
			
			neuron transfer/activation functions
			
			- sig		- 	the logistic sigmoid of the weighted sum of inputs. Output range [0,1]
			- gauss		-	the Gaussian of the weighted sum of inputs. Output range [0,1]
			- step		-	the heaviside step function of the weighted sum of inputs. Output range [0,1]
			- softsign	-	the softsign of the weighted sum of inputs. Output range [-1,1]
			- tanh		-	the hyperbolic tangent of the weighted sum of inputs. Output range [-1,1]
			
		Example:
			
			Add the node functions logical AND OR NAND NOR and XOR to the function set.
			(begin code)
			addNodeFunction(params, "and,or,nand,nor,xor");
			(end)
			
		See Also:
			<clearFunctionSet>, <addNodeFunctionCustom>, <printFunctionSet>
	*/
	void addNodeFunction(struct parameters *params, char *functionNames);
	
	
		/*
		Function: addNodeFunctionCustom

		Adds a custom node function to the set of functions made available to chromosome nodes.  
	
		The custom fitness function must take the form
		> float (*nodeFunctionName)(const int numInputs, const float *inputs, const float *weights)
	
		where the user replaces 'nodeFunctionName' with their own function name.
	
		Parameters:
			params - pointer to parameters structure
			function - the custom node function
			functionName - the name of the added function
			
		Example:
			
			custom node function, add
			
			(begin code)
			float add(const int numInputs, const float *inputs, const float *connectionWeights){
			 
				int i;
				float sum = 0;
			 
				for(i=0; i<numInputs; i++){
					sum += inputs[i];
				}
			 
				return sum;
			}	
			(end)
			
			Adding the new custom node function
			
			(begin code)
			addNodeFunctionCustom(params, add, "add");
			(end)
			
		See Also:
			<clearFunctionSet>, <addNodeFunction>, <printFunctionSet>
	*/
	void addNodeFunctionCustom(struct parameters *params, float (*function)(const int numInputs, const float *inputs, const float *weights), char *functionName);
	
	
	/*
		Function: clearFunctionSet

		Resets the function set to contain no functions.
	
		Parameters:
			params - pointer to parameters structure
						
		See Also:
			<addNodeFunction>, <addNodeFunctionCustom>, <printFunctionSet>
	*/
	void clearFunctionSet(struct parameters *params);
	
	
	/*
		Function: printFunctionSet

		Prints the contents of the current function set in a human readable form.  
			
		Parameters:
			params - pointer to parameters structure
						
		See Also:
			<addNodeFunction>, <addNodeFunctionCustom>, <clearFunctionSet>
	*/
	void printFunctionSet(struct parameters *params);
	
	
	/*
		Function: setMu

		Sets the mu value in the given parameters.
		
		The given mu value is also parsed to ensure a valid mu value. 
		mu values <1 are invalid. If an invalid mu value is give a 
		warning is displayed and the mu value is left unchanged.
		
		Parameters:
			params - pointer to parameters structure.
			mu - The value of mu to be set.

		See Also:
			<getMu>
	*/
	void setMu(struct parameters *params, int mu);
	
	
	/*
		Function: setLambda

		Sets the lambda value in the given parameters.
		
		The given lambda value is also parsed to ensure a valid lambda value. 
		lambda values <1 are invalid. If an invalid lambda value is give a 
		warning is displayed and the lambda value is left unchanged.
		
		Parameters:
			params - pointer to parameters structure.
			lambda - The value of lambda to be set.

		See Also:
			<getLambda>
	*/
	void setLambda(struct parameters *params, int lambda);
	
	
	/*
		Function: setEvolutionaryStrategy

		Sets the evolutionary strategy in the given parameters.
		
		The given evolutionary strategy is also parsed to ensure a valid evolutionary strategy. 
		Evolutionary strategies other than '+' and ',' are invalid. If an invalid evolutionary strategy is give a 
		warning is displayed and the evolutionary strategy is left unchanged.
				
		Parameters:
			params - pointer to parameters structure.
			evolutionaryStrategy - The evolutionary strategy to be set.

		See Also:
			<getEvolutionaryStrategy>
	*/
	void setEvolutionaryStrategy(struct parameters *params, char evolutionaryStrategy);
	

	/*
		Function: setMutationRate

		Sets the mutation rate in the given parameters.
		
		The given mutation rate is also parsed to ensure a valid mutation rate. 
		mutation rate <0 or >1 are invalid. If an invalid mutation rate is give a 
		warning is displayed and the mutation rate is left unchanged.
		
		Parameters:
			params - pointer to parameters structure.
			mutationRate - The value of the mutation rate to be set.

		See Also:
			<getMutationRate>
	*/
	void setMutationRate(struct parameters *params, float mutationRate);
	
	
	/*
		Function: setConnectionWeightRange

		Sets the connection weight range in the given parameters. (only used by NeuroEvolution)
				
		Parameters:
			params - pointer to parameters structure.
			weightRange - The connection weight range to be set. (the range is +/- weightRange)

		See Also:
			<getConnectionWeightRange>
	*/
	void setConnectionWeightRange(struct parameters *params, float weightRange);
	
	/*
		Function: setFitnessFunction

		Set custom fitness function. 
		
		The custom fitness function must take the form:
		> float functionName(struct parameters *params, struct chromosome *chromo, struct data *dat)
		
		If the fitnessFunction parameter is set as NULL, the fitness function will be reset to the default supervised learning. 

		Parameters:
			params - pointer to parameters structure.
			fitnessFunction - the custom fitness function 
			fitnessFunctionName - name of custom fitness function

		Example:
			
			Defining a custom fitness function, full adder
			
			(begin code)
			float fullAdder(struct parameters *params, struct chromosome *chromo, struct data *dat){
			
			int i;
			float error = 0;
			float chromoOutputs[8];	
			
			// full adder truth table
			float inputs[8][3] = {{0,0,0},{0,0,1},{0,1,0},{0,1,1},{1,0,0},{1,0,1},{1,1,0},{1,1,1}};
			float outputs[8][2] = {{0,0},{1,0},{1,0},{0,1},{1,0},{0,1},{0,1},{1,1}};
			 
			 	//for each line in the truth table 				
			 	for(i=0; i<8; i++){
			 		
			 		// calculate the chromosome outputs for the set of inputs  
			 		executeChromosome(chromo, inputs[i], chromoOutputs);
			 		
			 		// If the first chromosome outputs differ from the correct outputs increment the error 
			 		if(outputs[i][0] != chromoOutputs[0]){
			 			error++;
			 		}
			 		
			 		// If the second chromosome outputs differ from the correct outputs increment the error 
			 		if(outputs[i][1] != chromoOutputs[1]){
			 			error++;
			 		}
			 	}				
			 					
			 	return error;
			}
			(end)

			Setting the new custom fitness function as the fitness function to be used
			(begin code)
			struct parameters *params = initialiseParameters();
			setFitnessFuction(params, fullAdder, "fullAdder");
			(end)			
	*/
	void setFitnessFunction(struct parameters *params, float (*fitnessFunction)(struct parameters *params, struct chromosome *chromo, struct dataSet *data), char *fitnessFunctionName);
	

	/*
		Function: setNumInputs
		
			Sets the number of chromosome inputs in the given parameters.
		
			The given number of chromosome inputs is also parsed to ensure a valid number of chromosome inputs. 
			A number of chromosome inputs <0 is invalid. If an invalid number of chromosome inputs is give a 
			warning is displayed but the number of chromosome inputs *is changed*.
		
		Parameters:
			params - pointer to parameters structure.
			numInputs - The number of chromosome inputs to be set.

		See Also:
			<getNumInputs>
	*/
	void setNumInputs(struct parameters *params, int numInputs);


	/*
		Function: setNumNodes
		
			Sets the number of chromosome nodes in the given parameters.
		
			The given number of chromosome nodes is also parsed to ensure a valid number of chromosome nodes. 
			A number of chromosome nodes <0 is invalid. If an invalid number of chromosome nodes is give a 
			warning is displayed but the number of chromosome nodes *is changed*.
		
		Parameters:
			params - pointer to parameters structure.
			nodes - The number of chromosome nodes to be set.

		See Also:
			<getNumNodes>
	*/
	void setNumNodes(struct parameters *params, int numNodes);

	/*
		Function: setNumOutputs
		
			Sets the number of chromosome outputs in the given parameters.
		
			The given number of chromosome outputs is also parsed to ensure a valid number of chromosome outputs. 
			A number of chromosome outputs <1 is invalid. If an invalid number of chromosome outputs is give a 
			warning is displayed but the number of chromosome outputs *is changed*.
		
		Parameters:
			params - pointer to parameters structure.
			numOutputs - The number of chromosome outputs to be set.

		See Also:
			<getNumOutputs>
	*/
	void setNumOutputs(struct parameters *params, int numOutputs);


	/*
		Function: setArity
		
			Sets the arity of the chromosome nodes in the given parameters.
		
			The given arity for each chromosome node is also parsed to ensure a valid chromosome node arity. 
			A chromosome node arity <1 is invalid. If an invalid chromosome node arity is give a 
			warning is displayed but the chromosome node arity  *is changed*.
		
		Parameters:
			params - pointer to parameters structure.
			arity - The chromosome node arity to be set.

		See Also:
			<getArity>
	*/
	void setArity(struct parameters *params, int arity);


	/*
		Function: setTargetFitness
		
			Sets the target fitness used when running CGP. 
			
			In all cases lower fitness values are used to represent fitter chromosomes.
				
		Parameters:
			params - pointer to parameters structure.
			targetFitness - The target fitness to be set.

		See Also:
			<getTargetFitness>
	*/
	void setTargetFitness(struct parameters *params, float targetFitness);
	
	
	/*
		Function: setTargetFitness
		
			Sets the target fitness used when running CGP. 
			
			In all cases lower fitness values are used to represent fitter chromosomes.
				
		Parameters:
			params - pointer to parameters structure.
			targetFitness - The target fitness to be set.

		See Also:
			<getTargetFitness>
	*/
	
	
	/*
		Function: setMutationType
		
			Sets the mutation methods used when mutating chromosomes. 
			
			There are two mutation type options: the default "probabilistic" and "point". 
			If a invalid mutation type is given a warning is displayed and the mutation type is left unchanged.
				
		Parameters:
			params - pointer to parameters structure.
			mutationType - char array specifying the mutation type.

		See Also:
			<getMutationType>
	*/
	void setMutationType(struct parameters *params, char *mutationType);
	
	
	/*
		Function: setUpdateFrequency
		
			Sets the frequency of the updates to the user when using runCGP. 
			
			The update frequency represents the number of generations which elapse between showing the user the current best fitness.
			
		Note:	
			A value of '0' is a special case which causes not updates to be shown.
				
		Parameters:
			params - pointer to parameters structure.
			updateFrequency - update frequency in generations.

		See Also:
			<getUpdateFrequency> <runCGP>
	*/
	void setUpdateFrequency(struct parameters *params, int updateFrequency);
	
	
	
	/*
		Function: getMu

		Gets the value of mu current set in the given parameters.

		Parameters:
			params - pointer to parameters structure.

		Returns:
			The value of mu in the given parameters.

		See Also:
			<setMu>
	*/
	int getMu(struct parameters *params);
	
	
		/*
		Function: getNumInputs

		Gets the number of chromosome inputs current set in the given parameters.

		Parameters:
			params - pointer to parameters structure.

		Returns:
			The number of chromosome inputs current set in the given parameters

		See Also:
			<setNumInputs>
	*/
	int getNumInputs(struct parameters *params);
	
	
	/*
		Function: getNumOutputs

		Gets the number of chromosome output current set in the given parameters.

		Parameters:
			params - pointer to parameters structure.

		Returns:
			The number of chromosome output current set in the given parameters

		See Also:
			<setNumOutputs>
	*/
	int getNumOutputs(struct parameters *params);
	
	
/*
	Title: Chromosome Functions
	
	Description of all functions and structures relating to chromosomes
*/	
	
	/*
		Function: initialiseChromosome
			Initialises a chromosome based on the given parameters.

		Parameters:
			params - pointer to parameters structure

		Returns:
			A pointer to an initialised chromosome structure.
			
		See Also:
			<freeChromosome>, <initialiseChromosomeFromFile>, <printChromosome>
	*/
	struct chromosome *initialiseChromosome(struct parameters *params);
	
	
	/*
		Function: initialiseChromosomeFromFile
			Initialises a chromosome from a given previously saved chromosome.
			
		Note:
			Only chromosomes which use node functions defined by the CGP-library can be loaded. Chromosomes which use custom node functions cannot be loaded.
			
		Parameters:
			file - char array giving the location of the chromosome to be loaded.
		
		Returns:
			A pointer to an initialised chromosome structure.
			
		Example:
		
			(begin code)
			struct Chromosome *chromo;
			char *chromoFile = "location of Chromosome";
			
			chromo = loadChromosome(chromoFile);
			(end)
			
		See Also:
			<freeChromosome>, <saveChromosome>
	*/
	struct chromosome* initialiseChromosomeFromFile(char *file);
	
	
	/*
		Function: freeChromosome

			Frees chromosome instance.

		Parameters:
			chromo - pointer to chromosome structure.
			
		See Also:
			<initialiseChromosome>
	*/
	void freeChromosome(struct chromosome *chromo);
	
	
	/*
		Function: printChromosome
			Prints the given chromosome to the screen in a human readable format. 
			
		Parameters:
			chromo - pointer to chromosome structure.
	*/
	void printChromosome(struct chromosome *chromo);
	
	
	/*
		Function: executeChromosome
			Executes the given chromosome producing outputs to for the given inputs.
			
			The dimensions of the inputs and outputs arrays must match the dimensions of the chromosome inputs and outputs respectively. 
			
		Parameters:
			chromo - pointer to an initialised chromosome structure.
			inputs - array of floats used as inputs to the chromosome 
			outputs - array of floats used to store the calculated chromosome outputs
		
		Example:
		
			for a chromosome with three inputs and two outputs.
		
			(begin code)
			struct parameters *params;
			struct chromosome *chromo;
			
			float inputs[3] = {1,2,3};
			float outputs[2];
			
			params = initialiseParameters(3,10,2,2);
			chromo = initialiseChromosome(params);
			
			executeChromosome(chromo, inputs, outputs);
			(end)
	*/
	void executeChromosome(struct chromosome *chromo, float *inputs, float *outputs);
	
	
	/*
		Function: saveChromosome
			Saves the given chromosome to a file which can used to initialise new chromosomes.
			
		Node:
			Only chromosome which use node functions defined by the CGP-library can be loaded. Chromosomes which use custom node functions cannot be loaded.

			
		Parameters:
			chromo - pointer to chromosome structure.
			file - char array giving the location of the chromosome to be saved.
						
		See Also:
			<initialiseChromosomeFromFile>
	*/
	void saveChromosome(struct chromosome *chromo, char *file);
	
	
	/*
		Function: mutateChromosome
			Mutate the given chromosome using the mutation method described in the given parameters.
						
		Parameters:
			params - pointer to parameters structure
			chromo - pointer to chromosome structure.
			
		Example: 
		
			(begin code)
			struct parameters *params;
			struct chromosome *chromo;
			
			params = initialiseParameters(3,10,2,2);
			chromo = initialiseChromosome(params);
			
			mutateChromosome(params, chromo);
			(end)
	*/
	void mutateChromosome(struct parameters *params, struct chromosome *chromo);
	
	
	/*
		Function: removeInactiveNodes
			Removes all of the inactive nodes from the given chromosome.
		
		Note:
			This operation causes the number of nodes in the chromosome to change.
			
		Parameters:
			chromo - pointer to chromosome structure.
	*/
	void removeInactiveNodes(struct chromosome *chromo);
	
	/*
		Function: setChromosomeFitness
			Sets the fitness of the chromosome using the fitness function given in the parameters
			
		Parameters:
			params - pointer to parameters structure
			chromo - pointer to chromosome structure.
			data - pointer to the data used by the fitness function (NULL if not used)
						
		See Also:
			<getChromosomeFitness>
	*/
	void setChromosomeFitness(struct parameters *params, struct chromosome *chromo, struct dataSet *data);
	
	
	/*
		Function: getNumChromosomeInputs
			Gets the number of chromosome inputs
			
		Parameters:
			chromo - pointer to an initialised chromosome structure.
			
		Returns:
			Number number of chromosome inputs
			
		See Also:
			<initialiseChromosome>
	*/
	int getNumChromosomeInputs(struct chromosome *chromo);



	/*
		Function: getNumChromosomeNodes
			Gets the number of chromosome nodes
			
		Parameters:
			chromo - pointer to an initialised chromosome structure.
			
		Returns:
			Number number of chromosome nodes
			
		See Also:
			<initialiseChromosome>
	*/
	int getNumChromosomeNodes(struct chromosome *chromo);
	
	/*
		Function: getNumChromosomeActiveNodes
			Gets the number of chromosome active nodes
			
		Parameters:
			chromo - pointer to an initialised chromosome structure.
			
		Returns:
			Number number of chromosome active nodes
			
		See Also:
			<initialiseChromosome>
	*/
	int getNumChromosomeActiveNodes(struct chromosome *chromo);
	
	
	/*
		Function: getNumChromosomeOutputs
			Gets the number of chromosome outputs
			
		Parameters:
			chromo - pointer to an initialised chromosome structure.
			
		Returns:
			Number number of chromosome outputs
			
		See Also:
			<initialiseChromosome>
	*/
	int getNumChromosomeOutputs(struct chromosome *chromo);


	/*
		Function: getChromosomeNodeArity
			Gets the arity of the chromosome nodes
			
		Parameters:
			chromo - pointer to an initialised chromosome structure.
			
		Returns:
			Arity of chromosome nodes 
			
		See Also:
			<initialiseChromosome>
	*/
	int getChromosomeNodeArity(struct chromosome *chromo);
	
	
	/*
		Function: getChromosomeFitness
			Gets the fitness of the given chromosome 
			
		Parameters:
			chromo - pointer to initialised chromosome structure.
			
		Returns:
			The fitness of the given chromosome
			
		See Also:
			<setChromosomeFitness> <getChromosomeGenerations> <getChromosomeNumActiveNodes>
	*/
	float getChromosomeFitness(struct chromosome *chromo);
	
	
	/*
		Function: getChromosomeGenerations
			Gets the number of generations for which the given chromosome has been trained.
			
			If the chromosome has not been trained then -1 is returned.
			
		Parameters:
			chromo - pointer to initialised chromosome structure.
			
		Returns:
			Number of generations for which the given chromosome has been trained
			
		See Also:
			<getChromosomeFitness> <getNumChromosomeActiveNodes>
	*/
	int getChromosomeGenerations(struct chromosome *chromo);
	
	
	/*
		Function: copyChromosome
	*/
	void copyChromosome(struct chromosome *chromoDest, struct chromosome *chromoSrc);
	
	
/*
	Title: DataSet Functions
	
	Description of all functions and structures relating to data sets
*/
	
	/*
		Function: initialiseDataSetFromArrays

		Initialises a data structures using the given input output pairs.
		
		The data structures can be used by the fitness functions.

		Parameters:
			numInputs - number of inputs per data sample
			numOutputs - number of outputs per data sample
			numSamples - number of data samples 
			inputs - pointer to the inputs to be stored in the data structure
			outputs - pointer to the outputs to be stored in the data structure

		Returns:
			A pointer to an initialised dataSet structure.
			
		Example:
		
			(begin code)
			struct data *trainingData;
			
			int numInputs = 3;
			int numOutputs = 2;
			int numSamples = 8;
			
			float inputs[8][3] = {{0,0,0},{0,0,1},{0,1,0},{0,1,1},{1,0,0},{1,0,1},{1,1,0},{1,1,1}};
			float outputs[8][2] = {{0,0},{1,0},{1,0},{0,1},{1,0},{0,1},{0,1},{1,1}};	
			
			trainingData = initialiseDataFromArrays(numInputs, numOutputs, numSamples, inputs[0], outputs[0]);
			(end)

		See Also:
			<freeDataSet>, <initialiseDataSetFromFile>, <printDataSet>
	*/
	struct dataSet *initialiseDataSetFromArrays(int numInputs, int numOutputs, int numSamples, float *inputs, float *outputs);
	
	
	/*
		Function: initialiseDataSetFromFile

		Initialises a data structures using the given file.
		
		The data structures can be used by the fitness functions.

		Parameters:
			file - the location of the file to be loaded into the data structure

		Returns:
			A pointer to an initialised dataSet structure.
			
		Example:
		
			The file containing the data must take the form
			
			> 2,3,4
			> 1,2,3,4,5
			> 6,7,8,9,10
			> 11,12,13,14,15
			> 16,17,18,19,20
			
			Where the header describes the number of inputs, number of outputs and the number of samples. The rest of the file should give the input and then the output data of each sample on a new line.
		
			The file can then be used to initialise a data structure. 
		
			(begin code)
			struct data *trainingData;
			
			trainingData = initialiseDataFromFile("file");
			(end)

		See Also:
			<freeDataSet>, <initialiseDataSetFromArrays>, <printDataSet>
	*/
	struct dataSet *initialiseDataSetFromFile(char *file);
	
	
	/*  
		Function: freeDataSet 

		Frees dataSet instance.

		Parameters:
			data - pointer to dataSet structure.
			
		See Also:
			<initialiseDataSetFromArrays>, <initialiseDataSetFromFile>
	*/
	void freeDataSet(struct dataSet *data);
	
	/*
		Function: printDataSet

		Prints the input output pairs held by a dataSet structure in a human readable form.  

		Parameters:
			data - pointer to dataSet structure.
			
		See Also:
			<initialiseDataSetFromArrays>, <initialiseDataSetFromFile>, <freeDataSet>
	*/
	void printDataSet(struct dataSet *data);
	
	
	/*
		Function: saveDataSet

		Saves the given dataSet to a file which can be read using <initialiseDataSetFromFile>.

		Parameters:
			data - pointer to dataSet structure.
			fileName - char array giving the location of the dataSet to be saved.
			
		See Also:
			<initialiseDataSetFromFile>, <freeDataSet>
	*/
	void saveDataSet(struct dataSet *data, char *fileName);
	
	
	
/*
	Title: Results Functions
*/	
	
	
	/*  
		Function: freeResults 

		Frees <results> instance.

		Parameters:
			rels - pointer to results structure.
			
		See Also:
			<getChromosome>, <getAverageFitness>, <getAverageActiveNodes>, <getAverageGenerations>
	*/
	void freeResults(struct results *rels); 
	
	
	struct chromosome* getChromosome(struct results *rels, int run);
	
	
	float getAverageFitness(struct results *rels); 
	float getAverageActiveNodes(struct results *rels);
	float getAverageGenerations(struct results *rels);
	
	
/*
	Title: CGP Functions
*/
	
	/*
		Function: runCGP
	
		Applies CGP to the given task.
	
		Returns the best chromosome found after applying CGP to a given task using the specified parameters. Depending upon the update frequency given in parameters (<setUpdateFrequency>) the progress made by CGP will be displayed in the terminal.
	
		Note:
			As runCGP returns an initialised chromosome this should later be free'd using <freeChromosome>. 
	
		Parameters:
			params - pointer to parameters structure.
			data - pointer to dataSet structure.
			gens - the number of allowed generations before terminating the search.
			
		Returns:	
			A pointer to an initialised chromosome.
		
		Example:
		
			(begin code)
			struct parameters *params;
			struct dataSet *data;
			struct chromosome *chromo;
			
			params = initialiseParameters(a,b,c,d);
			addNodeFunction(params, "aaa,bbb,ccc");
			
			data = initialiseDataSetFromFile("file");
			
			chromo = runCGP(params, data, 100);
			
			freeParameters(params);
			freeDataSet(data);
			freeChromosome(chromo);
			(end)
		
			
		See Also:
			<repeatCGP>, <initialiseParameters>, <initialiseDataSetFromFile>, <freeChromosome>
	*/	
	struct chromosome* runCGP(struct parameters *params, struct dataSet *data, int numGens);
	
	
	/*
		Function: repeatCGP
	
		Repeatedly applies CGP to the given task.
	
		Returns a <results> structure containing the results of applying CGP to the given task multiple times.
		
		Note:
			As repeatCGP returns an initialised results structure this should later be free'd using <freeResults>. 
	
		Parameters: 
			params - pointer to parameters structure.
			data - pointer to dataSet structure.
			numGens - the number of allowed generations before terminating the search.
			numRuns - the number of times CGP will be applied to the given task
			
		Returns:	
			A pointer to an initialised results structure.
			
		Example:
		
			(begin code)
			struct parameters *params;
			struct dataSet *data;
			struct results *rels;
			 
			params = initialiseParameters(a,b,c,d);
			addNodeFunction(params, "aaa,bbb,ccc");
			
			data = initialiseDataSetFromFile("file");
			
			rels = repeatCGP(params, data, 100, 50);
			
			freeParameters(params);
			freeDataSet(data);
			freeResults(rels);
			(end)	
			
			
		See Also:
			<runCGP>, <initialiseParameters>, <initialiseDataSetFromFile>, <freeResults> 
	*/	
	struct results* repeatCGP(struct parameters *params, struct dataSet *data, int numGens, int numRuns);
	

#endif 
