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

#if defined(_WIN32) && defined(NO_DLL)
	#define DLL_EXPORT
#elif  defined(_WIN32) && defined(BUILD_DLL)
    #define DLL_EXPORT __declspec(dllexport)
#elif defined(_WIN32) && !defined(BUILD_DLL)
    #define DLL_EXPORT __declspec(dllimport)
#else
    #define DLL_EXPORT
#endif



/*
	Title: Structures

	Description of all the structures used by CGP-Library../t
*/

	/*
		variable: parameters

		Stores general parameters used by the CGP-Library.

		The values stored in parameters are set to default values when initialised and can be altered though setter functions. The parameters variables can not be accessed directly; getters and Setters must be used.

		Defaults:

			Default parameter values.

			> mu:						1
			> lambda:					4
			> evolutionary strategy:	+
			> mutation rate:			0.05
			> connection weight range:	1
			> update frequency:			1
			> mutation type:			probabilistic
			> fitness function: 		supervisedLearning
			> selection scheme:			pickHighest
			> reproduction scheme:		mutateRandomParent

			For information on probabilistic and other available mutation types see <setMutationType>. The supervisedLearning fitness function assigns the sum of the absolute difference between all actual and target outputs. This fitness function can be changed using <setFitnessFunction>. The pickHighest selection scheme selects the fittest members of the population to become the parents each generation. The mutateRandomParent reproduction scheme create the children by mutating randomly selected parents.

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

		Stores the the best chromosome found on each run when using <repeatCGP>

		See Also:

			<repeatCGP>, <freeResults>, <getChromosome>, <getResultsAverageFitness>, <getResultsAverageActiveNodes>, <getResultsAverageGenerations>

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

			(begin code)
			struct parameters *params;

			int numInputs = 3;
			int numNodes = 10;
			int numOutputs = 2;
			int nodeArity = 2;

			params = initialiseParameters(numInputs, numNodes, numOutputs, nodeArity);
			(end)

		See Also:
			<freeParameters>, <printParameters>
	*/
	DLL_EXPORT struct parameters *initialiseParameters(const int numInputs, const int numNodes, const int numOutputs, const int arity);

	/*
		Function: freeParameters

		Frees parameters instance.

		Parameters:
			params - pointer to parameters structure.

		See Also:
			<initialiseParameters>
	*/
	DLL_EXPORT void freeParameters(struct parameters *params);

	/*
		Function: printParameters

		Prints the given parameters to the screen in a human readable format.

		Parameters:
			params - pointer to parameters structure.

	*/
	DLL_EXPORT void printParameters(struct parameters *params);


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
	DLL_EXPORT void addNodeFunction(struct parameters *params, char *functionNames);


		/*
		Function: addNodeFunctionCustom

		Adds a custom node function to the set of functions made available to chromosome nodes.

		The custom fitness function must take the form

		(begin code)
		float nodeFunctionName(const int numInputs, const float *inputs, const float *weights);
		(end)

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
	DLL_EXPORT void addNodeFunctionCustom(struct parameters *params, float (*function)(const int numInputs, const float *inputs, const float *weights), char *functionName);


	/*
		Function: clearFunctionSet

		Resets the function set to contain no functions.

		Parameters:
			params - pointer to parameters structure

		See Also:
			<addNodeFunction>, <addNodeFunctionCustom>, <printFunctionSet>
	*/
	DLL_EXPORT void clearFunctionSet(struct parameters *params);


	/*
		Function: printFunctionSet

		Prints the contents of the current function set in a human readable form.

		Parameters:
			params - pointer to parameters structure

		See Also:
			<addNodeFunction>, <addNodeFunctionCustom>, <clearFunctionSet>
	*/
	DLL_EXPORT void printFunctionSet(struct parameters *params);


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
	DLL_EXPORT void setMu(struct parameters *params, int mu);


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
	DLL_EXPORT void setLambda(struct parameters *params, int lambda);


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
	DLL_EXPORT void setEvolutionaryStrategy(struct parameters *params, char evolutionaryStrategy);


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
	DLL_EXPORT void setMutationRate(struct parameters *params, float mutationRate);


	/*
		Function: setConnectionWeightRange

		Sets the connection weight range in the given parameters. (only used by NeuroEvolution)

		Parameters:
			params - pointer to parameters structure.
			weightRange - The connection weight range to be set. (the range is +/- weightRange)

		See Also:
			<getConnectionWeightRange>
	*/
	DLL_EXPORT void setConnectionWeightRange(struct parameters *params, float weightRange);

	/*
		Function: setFitnessFunction

		Set custom fitness function.

		The custom fitness function prototype must take the form

		(begin code)
		float functionName(struct parameters *params, struct chromosome *chromo, struct dataSet *data);
		(end)

		If the fitnessFunction parameter is set as NULL, the fitness function will be reset to the default supervised learning fitness function; See <parameters>

		Parameters:
			params - pointer to parameters structure.
			fitnessFunction - the custom fitness function
			fitnessFunctionName - name of custom fitness function

		Example:

			Defining a custom fitness function, full adder. Note that the <dataSet> does not have to be used.

			(begin code)
			float fullAdder(struct parameters *params, struct chromosome *chromo, struct data *dat){

			int i;
			float error = 0;

			// full adder truth table
			float inputs[8][3] = {{0,0,0},{0,0,1},{0,1,0},{0,1,1},{1,0,0},{1,0,1},{1,1,0},{1,1,1}};
			float outputs[8][2] = {{0,0},{1,0},{1,0},{0,1},{1,0},{0,1},{0,1},{1,1}};

			 	//for each line in the truth table
			 	for(i=0; i<8; i++){

			 		// calculate the chromosome outputs for each set of inputs
			 		executeChromosome(chromo, inputs[i]);

			 		// If the first chromosome outputs differ from the correct outputs increment the error
			 		if(outputs[i][0] != getChromosomeOutput(chromo, 0) ){
			 			error++;
			 		}

			 		// If the second chromosome outputs differ from the correct outputs increment the error
			 		if(outputs[i][1] != getChromosomeOutput(chromo, 1) ){
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

		See Also:
			<executeChromosome>, <getChromosomeOutput>

	*/
	DLL_EXPORT void setFitnessFunction(struct parameters *params, float (*fitnessFunction)(struct parameters *params, struct chromosome *chromo, struct dataSet *data), char *fitnessFunctionName);


	/*
	Function: setSelectionScheme
	*/
	DLL_EXPORT void setSelectionScheme(struct parameters *params, void (*selectionScheme)(struct parameters *params, struct chromosome **parents, struct chromosome **candidateChromos, int numCandidateChromos), char *selectionSchemeName);


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
	DLL_EXPORT void setNumInputs(struct parameters *params, int numInputs);


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
	DLL_EXPORT void setNumNodes(struct parameters *params, int numNodes);

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
	DLL_EXPORT void setNumOutputs(struct parameters *params, int numOutputs);


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
	DLL_EXPORT void setArity(struct parameters *params, int arity);


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
	DLL_EXPORT void setTargetFitness(struct parameters *params, float targetFitness);


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
	DLL_EXPORT void setMutationType(struct parameters *params, char *mutationType);


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
	DLL_EXPORT void setUpdateFrequency(struct parameters *params, int updateFrequency);


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
	DLL_EXPORT struct chromosome *initialiseChromosome(struct parameters *params);


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
	DLL_EXPORT struct chromosome* initialiseChromosomeFromFile(char *file);


	/*
		Function: freeChromosome

			Frees chromosome instance.

		Parameters:
			chromo - pointer to chromosome structure.

		See Also:
			<initialiseChromosome>, <initialiseChromosomeFromFile>
	*/
	DLL_EXPORT void freeChromosome(struct chromosome *chromo);


	/*
		Function: printChromosome
			Prints the given chromosome to the screen in a human readable format.

		Parameters:
			chromo - pointer to chromosome structure.
	*/
	DLL_EXPORT void printChromosome(struct chromosome *chromo);


	/*
		Function: executeChromosome
			Executes the given chromosome.

			Executes the given chromosome with the given inputs. The dimensions of the inputs arrays must match the dimensions of the chromosome inputs. The chromosome outputs are then accessed using <getChromosomeOutput>.

		Parameters:
			chromo - pointer to an initialised chromosome structure.
			inputs - array of floats used as inputs to the chromosome

		Example:

			for a chromosome with three inputs and one outputs.

			(begin code)
			struct parameters *params;
			struct chromosome *chromo;

			float inputs[3] = {1,2,3};

			params = initialiseParameters(3,10,1,2);
			chromo = initialiseChromosome(params);

			executeChromosome(chromo, inputs);

			printf("Output: %f\n", getChromosomeOutput(chromo, 0));

			(end)

		See Also:
				<getChromosomeOutput>

	*/
	DLL_EXPORT void executeChromosome(struct chromosome *chromo, float *inputs);



	/*
		Function: getChromosomeOutput
			gets the outputs of the given chromosome.

			After a given chromosome has been executed using <executeChromosome> the chromosome outputs are made available using <getChromosomeOutput>.

			Parameters:
				chromo - pointer to an initialised chromosome structure.
				output - The index of the output to be retrieved

		See Also:
				<executeChromosome>
	*/
	DLL_EXPORT float getChromosomeOutput(struct chromosome *chromo, int output);



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
	DLL_EXPORT void saveChromosome(struct chromosome *chromo, char *file);


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
	DLL_EXPORT void mutateChromosome(struct parameters *params, struct chromosome *chromo);


	/*
		Function: removeInactiveNodes
			Removes all of the inactive nodes from the given chromosome.

		Note:
			This operation causes the number of nodes in the chromosome to change.

		Parameters:
			chromo - pointer to chromosome structure.
	*/
	DLL_EXPORT void removeInactiveNodes(struct chromosome *chromo);

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
	DLL_EXPORT void setChromosomeFitness(struct parameters *params, struct chromosome *chromo, struct dataSet *data);


	/*
		Function: copyChromosome

		Copies the contents of one chromoosme into the other.

		In order for copy chromosome to opperate correctly the dimensions of the chromosomes must match i.e. the number of inputs, nodes, outputs and node arity must be the same.

		Parameters:
			chromoDest - pointer to an initialised chromosome to be copied too
			chromoSrc - pointer to an initialised chromosome to be copied from

		See Also:
			<initialiseChromosome>

	*/
	DLL_EXPORT void copyChromosome(struct chromosome *chromoDest, struct chromosome *chromoSrc);


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
	DLL_EXPORT int getNumChromosomeInputs(struct chromosome *chromo);



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
	DLL_EXPORT int getNumChromosomeNodes(struct chromosome *chromo);

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
	DLL_EXPORT int getNumChromosomeActiveNodes(struct chromosome *chromo);


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
	DLL_EXPORT int getNumChromosomeOutputs(struct chromosome *chromo);


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
	DLL_EXPORT int getChromosomeNodeArity(struct chromosome *chromo);


	/*
		Function: getChromosomeFitness
			Gets the fitness of the given chromosome

		Parameters:
			chromo - pointer to initialised chromosome structure.

		Returns:
			The fitness of the given chromosome

		See Also:
			<setChromosomeFitness> <getChromosomeGenerations> <getNumChromosomeActiveNodes>
	*/
	DLL_EXPORT float getChromosomeFitness(struct chromosome *chromo);


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
	DLL_EXPORT int getChromosomeGenerations(struct chromosome *chromo);








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
	DLL_EXPORT struct dataSet *initialiseDataSetFromArrays(int numInputs, int numOutputs, int numSamples, float *inputs, float *outputs);


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
	DLL_EXPORT struct dataSet *initialiseDataSetFromFile(char *file);


	/*
		Function: freeDataSet

		Frees dataSet instance.

		Parameters:
			data - pointer to dataSet structure.

		See Also:
			<initialiseDataSetFromArrays>, <initialiseDataSetFromFile>
	*/
	DLL_EXPORT void freeDataSet(struct dataSet *data);

	/*
		Function: printDataSet

		Prints the input output pairs held by a dataSet structure in a human readable form.

		Parameters:
			data - pointer to dataSet structure.

		See Also:
			<initialiseDataSetFromArrays>, <initialiseDataSetFromFile>, <freeDataSet>
	*/
	DLL_EXPORT void printDataSet(struct dataSet *data);


	/*
		Function: saveDataSet

		Saves the given dataSet to a file which can be read using <initialiseDataSetFromFile>.

		Parameters:
			data - pointer to dataSet structure.
			fileName - char array giving the location of the dataSet to be saved.

		See Also:
			<initialiseDataSetFromFile>, <freeDataSet>
	*/
	DLL_EXPORT void saveDataSet(struct dataSet *data, char *fileName);


	/*
		Function: getNumDataSetInputs
			Gets the number of dataSet inputs

		Parameters:
			data - pointer to an initialised dataSet structure.

		Returns:
			Number number of dataSet inputs

		See Also:
			<getNumDataSetOutputs>, <getNumDataSetSamples>
	*/
	DLL_EXPORT int getNumDataSetInputs(struct dataSet *data);



	/*
		Function: getNumDataSetOutputs
			Gets the number of dataSet outputs

		Parameters:
			data - pointer to an initialised dataSet structure.

		Returns:
			Number number of dataSet outputs

		See Also:
			<getNumDataSetInputs>, <getNumDataSetSamples>
	*/
	DLL_EXPORT int getNumDataSetOutputs(struct dataSet *data);


	/*
		Function: getNumDataSetSamples
			Gets the number of samples in the given dataSet

		Parameters:
			data - pointer to an initialised dataSet structure.

		Returns:
			Number number of dataSet samples

		See Also:
			<getNumDataSetInputs>, <getNumDataSetOutputs>
	*/
	DLL_EXPORT int getNumDataSetSamples(struct dataSet *data);


	/*
		Function: getDataSetSampleInputs
			Gets the dataSet inputs for the given sample index

		Parameters:
			data - pointer to an initialised dataSet structure.
			sample - index of the sample inputs

		Returns:
			Pointer to array containing the sample inputs

		See Also:
			<getDataSetSampleInput>, <getDataSetSampleOutputs>, <getDataSetSampleOutput>
	*/
	DLL_EXPORT float *getDataSetSampleInputs(struct dataSet *data, int sample);


	/*
		Function: getDataSetSampleInput
			Gets the dataSet input for the given sample index and input index

		Parameters:
			data - pointer to an initialised dataSet structure.
			sample - index of the sample inputs
			input - index of the input for the given sample

		Returns:
			The given input value for the given input for the given sample

		See Also:
			<getDataSetSampleInputs>, <getDataSetSampleOutput>, <getDataSetSampleOutputs>
	*/
	DLL_EXPORT float getDataSetSampleInput(struct dataSet *data, int sample, int input);


	/*
		Function: getDataSetSampleOutputs
			Gets the dataSet outputs for the given sample index

		Parameters:
			data - pointer to an initialised dataSet structure.
			sample - index of the sample outputs

		Returns:
			Pointer to array containing the sample outputs

		See Also:
			<getDataSetSampleOutput>, <getDataSetSampleinputs>, <getDataSetSampleinput>
	*/
	DLL_EXPORT float *getDataSetSampleOutputs(struct dataSet *data, int sample);

	/*
		Function: getDataSetSampleOutput
			Gets the dataSet output for the given sample index and output index

		Parameters:
			data - pointer to an initialised dataSet structure.
			sample - index of the sample inputs
			output - index of the output for the given sample

		Returns:
			The given output value for the given output for the given sample

		See Also:
			<getDataSetSampleOutputs>, <getDataSetSampleInputt>, <getDataSetSampleInputss>
	*/
	DLL_EXPORT float getDataSetSampleOutput(struct dataSet *data, int sample, int output);


/*
	Title: Results Functions
*/

	/*
		Function: freeResults

		Frees <results> instance.

		Parameters:
			rels - pointer to results structure.

		See Also:
			<getChromosome>, <getResultsAverageFitness>, <getResultsAverageActiveNodes>, <getResultsAverageGenerations>
	*/
	DLL_EXPORT void freeResults(struct results *rels);

	/*
		Function: getChromosome
			Gets the best chromosome found on the given run in results.

		Parameters:
			rels - pointer to an initialised results structure.
			run - the run for which the chromosome is got

		Returns:
			Pointer to an initialised chromosome structure

		See Also:
			<repeatCGP>, <getResultsAverageFitness>, <getResultsAverageActiveNodes>, <getResultsAverageGenerations>
	*/
	DLL_EXPORT struct chromosome* getChromosome(struct results *rels, int run);
	
	/*
		Function: getNumChromosomes
			Gets number of chromosomes stored in the given results structure.

		Parameters:
			rels - pointer to an initialised results structure.
			
		Returns:
			The number of chromosomes stored in the given results structure

		See Also:
			<repeatCGP>, <getResultsAverageFitness>, <getResultsAverageActiveNodes>, <getResultsAverageGenerations>
	*/
	DLL_EXPORT int getNumChromosomes(struct results *rels);

	/*
		Function: getResultsAverageFitness
			Gets the average fitness of the best chromosome found for each run in results.

		Parameters:
			rels - pointer to an initialised results structure.

		Returns:
			The average fitness of the best chromosome found for each run in results.

		See Also:
			<repeatCGP>, <getResultsAverageActiveNodes>, <getResultsAverageGenerations>
	*/
	DLL_EXPORT float getResultsAverageFitness(struct results *rels);


	/*
		Function: getResultsAverageActiveNodes
			Gets the average number of active nodes of the best chromosome found for each run in results.

		Parameters:
			rels - pointer to an initialised results structure.

		Returns:
			The average number of active nodes of the best chromosome found for each run in results.

		See Also:
			<repeatCGP>, <getResultsAverageFitness>, <getResultsAverageGenerations>

	*/
	DLL_EXPORT float getResultsAverageActiveNodes(struct results *rels);

	/*
		Function: getResultsAverageGenerations
			Gets the average number generations required to find the best best chromosome for each run in results.

		Parameters:
			rels - pointer to an initialised results structure.

		Returns:
			The average number generations required to find the best chromosome found for each run in results.

		See Also:
			<repeatCGP>, <getResultsAverageFitness>, <getResultsAverageActiveNodes>

	*/
	DLL_EXPORT float getResultsAverageGenerations(struct results *rels);


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
	DLL_EXPORT struct chromosome* runCGP(struct parameters *params, struct dataSet *data, int numGens);


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
	DLL_EXPORT struct results* repeatCGP(struct parameters *params, struct dataSet *data, int numGens, int numRuns);


/*
	Title: Other
*/

	/*
		Function: setRandomNumberSeed
			Sets the seed used by the random number generator in CGP-Library.

			By default the current time is used as the random number seed. When a random number seed is specified the CGP-Library will produce the same results if used in the same way.

		Note:
			<setRandomNumberSeed> *must* be called *after* <initialiseParameters> otherwise the time will be used as the random number seed.

		Parameters:
			seed - the random number seed to be used

		See Also:
			<initialiseParameters>
	*/
	DLL_EXPORT void setRandomNumberSeed(unsigned int seed);




#endif
