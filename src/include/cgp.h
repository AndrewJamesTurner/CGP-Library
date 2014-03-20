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
			
			> struct parameters *params;
			>
			> int numInputs = 3;
			> int numNodes = 10;
			> int numOutputs = 2;
			> int nodeArity = 2;		
			>
			> params = initialiseParameters(numInputs, numNodes, numOutputs, nodeArity);
			>
			> freeParameters(params);
			

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
	
	int getMu(struct parameters *params);
	void setMu(struct parameters *params, int mu);
	
	int getNumInputs(struct parameters *params);
	int getNumOutputs(struct parameters *params);
	
	void setFitnessFuction(struct parameters *params, float (*fitnessFuction)(struct parameters *params, struct chromosome *chromo, struct data *dat), char *fitnessFuctionName);
	
	/* data functions */
	struct data *initialiseDataFromArrays(int numInputs, int numOutputs, int numSamples, float *inputs, float *outputs);
	struct data *initialiseDataFromFile(char *file);
	void freeData(struct data *dat);
	void printData(struct data *dat);
	
	/* function set functions */
	void addNodeFuction(struct parameters *params, char *functionNames);
	void addNodeFuctionCustom(struct parameters *params, float (*function)(const int numInputs, const float *inputs, const float *weights), char *functionName);
	void clearFuctionSet(struct parameters *params);
	void printFuctionSet(struct parameters *params);
	
	/* chromosome functions  */
	struct chromosome *initialiseChromosome(struct parameters *params);
	void freeChromosome(struct chromosome *chromo);
	void executeChromosome(struct parameters *params, struct chromosome *chromo, float *inputs, float *outputs);
	void mutateChromosome(struct parameters *params, struct chromosome *chromo);
	void printChromosome(struct chromosome *chromo);
	float getChromosomeFitness(struct chromosome *chromo);
	
	/* population functions */
	struct population *initialisePopulation(struct parameters *params);
	void freePopulation(struct parameters *params, struct population *pop);
	void evolvePopulation(struct parameters *params, struct population *pop, struct data *dat);
	struct chromosome *getFittestChromosome(struct parameters *params, struct population *pop);	
	int getNumberOfGenerations(struct population *pop);
		
		
	
	
	/*
		variable: parameters
		Stores general parameters used by CGP-Library.
		
		The structures variables should not and can not be accessed directly. Getters and Setters must be used.
		
		See Also:
			<initialiseParameters> <freeParameters>
	*/
	
	struct parameters;
	struct population;
	struct chromosome;
	struct fuctionSet;
	struct data;
	
	
	
	
	
#endif 
