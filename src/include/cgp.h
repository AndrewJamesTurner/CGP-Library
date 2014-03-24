/*
	This file is part of CGP-Library, Andrew James Turner 2014.

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
	
	Description of all the CGP-Library functions and variables.	
	
*/


#ifndef CGPLIB
#define CGPLIB
		
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
			<initialiseParameters>, <freeParameters>
	*/	
	struct parameters;
	struct chromosome;
	struct fuctionSet;
	struct data;	
	struct population;
	struct results;
		
	
	/*
		Function: initialiseParameters

		Initialises parameters used throughout the CGP-Library. The inputs describe the structure of the chromosomes created when using <initialiseChromosome> or <initialisePopulation>.

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
			<freeParameters>
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
		Function: setMu

		Sets the mu value in the given parameters.
		
		The given mu value is also parsed to ensure a valid mu value. 
		mu values <1 are invalid. If a invalid mu value is give a 
		warning is displayed and the mu value is left unchanged.
		
		Parameters:
			params - pointer to parameters structure.
			mu - The value of mu to be set.

		See Also:
			<getMu>
	*/
	void setMu(struct parameters *params, int mu);
	
	
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
			> float fullAdder(struct parameters *params, struct chromosome *chromo, struct data *dat){
			>
			> int i;
			> float error = 0;
			> float chromoOutputs[8];	
			> float inputs[8][3] = {{0,0,0},{0,0,1},{0,1,0},{0,1,1},{1,0,0},{1,0,1},{1,1,0},{1,1,1}};
			> float outputs[8][2] = {{0,0},{1,0},{1,0},{0,1},{1,0},{0,1},{0,1},{1,1}};
			> 
			> 	//for each line in the truth table 				
			> 	for(i=0; i<8; i++){
			> 		
			> 		// calculate the chromosome outputs for the set of inputs  
			> 		executeChromosome(params, chromo, inputs[i], chromoOutputs);
			> 		
			> 		// If the first chromosome outputs differ from the correct outputs increment the error 
			> 		if(outputs[i][0] != chromoOutputs[0]){
			> 			error++;
			> 		}
			> 		
			> 		// If the second chromosome outputs differ from the correct outputs increment the error 
			> 		if(outputs[i][1] != chromoOutputs[1]){
			> 			error++;
			> 		}
			> 	}				
			> 					
			> 	return error;
			> }

			Setting the new custom fitness function as the fitness function to be used
			> struct parameters *params = initialiseParameters();
			> setFitnessFuction(params, fullAdder, "fullAdder");			
	*/
	void setFitnessFunction(struct parameters *params, float (*fitnessFunction)(struct parameters *params, struct chromosome *chromo, struct data *dat), char *fitnessFunctionName);
	
	
	
	
	/*
		Function: initialiseDataFromArrays

		Initialises a data structures using the given input output pairs.
		
		The data structures can be used by the fitness functions.

		Parameters:
			numInputs - number of inputs per data sample
			numOutputs - number of outputs per data sample
			numSamples - number of data samples 
			inputs - pointer to the inputs to be stored in the data structure
			outputs - pointer to the outputs to be stored in the data structure

		Returns:
			A pointer to an initialised data structure.
			
		Example:
		
			> struct data *trainingData;
			>
			> int numInputs = 3;
			> int numOutputs = 2;
			> int numSamples = 8;
			>
			> float inputs[8][3] = {{0,0,0},{0,0,1},{0,1,0},{0,1,1},{1,0,0},{1,0,1},{1,1,0},{1,1,1}};
			> float outputs[8][2] = {{0,0},{1,0},{1,0},{0,1},{1,0},{0,1},{0,1},{1,1}};	
			>
			> trainingData = initialiseDataFromArrays(numInputs, numOutputs, numSamples, inputs[0], outputs[0]);

		See Also:
			<freeData>, <initialiseDataFromFile>, <printData>
	*/
	struct data *initialiseDataFromArrays(int numInputs, int numOutputs, int numSamples, float *inputs, float *outputs);
	
	
	/*
		Function: initialiseDataFromFile

		Initialises a data structures using the given file.
		
		The data structures can be used by the fitness functions.

		Parameters:
			file - the location of the file to be loaded into the data structure

		Returns:
			A pointer to an initialised data structure.
			
		Example:
		
			The file containing the data must take the form
			
			> 2,3,4
			> 1,2,3,4,5
			> 6,7,8,9,10
			> 11,12,13,14,15
			> 16,17,18,19,20
			
			Where the header describes the number of inputs, number of outputs and the number of samples. The rest of the file should give the input and then the output data of each sample on a new line.
		
			The file can then be used to initialise a data structure. 
		
			> struct data *trainingData;
			>
			> trainingData = initialiseDataFromFile("file");

		See Also:
			<freeData>, <initialiseDataFromArrays>, <printData>
	*/
	struct data *initialiseDataFromFile(char *file);
	
	
	/*
		Function: freeData

		Frees data instance.

		Parameters:
			dat - pointer to data structure.
			
		See Also:
			<initialiseDataFromArrays>, <initialiseDataFromFile>
	*/
	void freeData(struct data *dat);
	
	/*
		Function: printData

		Prints the input output pairs held by a data structure in a human readable form.  

		Parameters:
			dat - pointer to data structure.
			
		See Also:
			<initialiseDataFromArrays>, <initialiseDataFromFile>, <freeData>
	*/
	void printData(struct data *dat);
	
	
	
	/*
		Function: addNodeFunction

		Adds pre-made node function(s) to the set of functions made available to chromosome nodes.  
	
		If one function name is given that function is added to the function set. If multiple node function's names are given then each must be separated by a ','.
		
		If a node function name is given which is not recognised, a warning is given and that function is not added to the function set.

		Parameters:
			params - pointer to parameters structure
			functionNames - the name(s) of the function(s) to add to the function set
		
		Node Functions:
		
			mathematical operations 
			
			- add 	- 	summation over all inputs.
			- sub	-	subtracts all but the first input from the first input
			
			logic gates
			
			- and	-	returns '1' if all inputs are '1', else '0'
			- nand	-	returns '0' if all inputs are '1', else, '1'
			- or	-	returns '0' if all inputs are '0', else, '1'
			- nor	-	returns '1' if all inputs are '0', else, '0'
			- xor	-	returns '1' if only one of the inputs is '1', else, '0'
			- xnor	-	returns '0' if only one of the inputs is '1', else, '1'
			- not	-	returns '1' if first input is '0', else '1'
			
		Example:
			
			Add the node functions logical AND OR NAND NOR and XOR to the function set.
			> addNodeFunction(params, "and,or,nand,nor,xor");
			
		See Also:
			<clearFunctionSet>, <addNodeFunctionCustom>, <printFunctionSet>
	*/
	void addNodeFunction(struct parameters *params, char *functionNames);
	
	
		/*
		Function: addNodeFunctionCustom

		Adds a custom node function to the set of functions made available to chromosome nodes.  
	
		The custom fitness function must take the form:
		> float (*nodeFunctionName)(const int numInputs, const float *inputs, const float *weights)
	
		Parameters:
			params - pointer to parameters structure
			function - the custom node function
			functionName - the name of the added function
			
		Example:
			
			custom node function, add
			> float add(const int numInputs, const float *inputs, const float *connectionWeights){
			> 
			> 	int i;
			> 	float sum = 0;
			> 
			> 	for(i=0; i<numInputs; i++){
			> 		sum += inputs[i];
			> 	}
			> 
			> 	return sum;
			> }	
			
			Adding the new custom node function
			> addNodeFunctionCustom(params, add, "add");
			
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
		Function: initialiseChromosome

		Initialises chromosome based on the parameters given in params.

		Parameters:
			params - pointer to parameters structure

		Returns:
			A pointer to an initialised chromosome structure.
			
		See Also:
			<freeChromosome>, <printChromosome>
	*/
	struct chromosome *initialiseChromosome(struct parameters *params);
	
	/*
		Function: freeChromosome

		Frees chromosome instance.

		Parameters:
			chromo - pointer to chromosome structure.
			
		See Also:
			<initialiseChromosome>
	*/
	void freeChromosome(struct chromosome *chromo);
	
	
	
	void executeChromosome(struct chromosome *chromo, float *inputs, float *outputs);
	void mutateChromosome(struct parameters *params, struct chromosome *chromo);
	void printChromosome(struct chromosome *chromo);
	float getChromosomeFitness(struct chromosome *chromo);
	
	void setChromosomeFitness(struct parameters *params, struct chromosome *chromo, struct data *dat);
	int getChromosomeActiveNodes(struct chromosome *chromo);
	

	/*
		returns an initilised chromosome which shoyulod be freed by the user.
	*/	
	struct chromosome* runCGP(struct parameters *params, struct data *dat);
	
	
	struct results* repeatCGP(struct parameters *params, struct data *dat, int numRuns);
	
	void freeResults(struct results *rels);
	
	float getAverageFitness(struct results *rels);
	float getAverageActiveNodes(struct results *rels);
	float getAverageGenerations(struct results *rels);
	
	struct chromosome* getChromosome(struct results *rels, int run);
	
	int getChromosomeGenerations(struct chromosome *chromo);
	
	void setNumGenerations(struct parameters *params, int numGens);
	void setTargetFitness(struct parameters *params, float targetFitness);
	
	void saveChromosome(struct chromosome *chromo, char *file);
	struct chromosome* loadChromosome(char *file);
	
	
#endif 
