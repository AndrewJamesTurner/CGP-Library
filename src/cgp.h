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

#ifndef CGPLIB
#define CGPLIB

	/* Structure definitions given via opaque pointers */
	struct parameters;
	struct population;
	struct chromosome;
	struct node;
	struct fuctionSet;
	
	
	/* */
	struct parameters *initialiseParameters(const int numInputs, const int numNodes, const int numOutputs, const int arity);
	
	/*
		getters and setters for the parameters. Getters return the current values
		stored in parameters. Setter set the values in parameters to new values. If
		invalid values are passed to the setters and warning is given and the parameters
		value remains unchanged. 
	*/
	
	void freeParameters(struct parameters *params);
	
	int getMu(struct parameters *params);
	void setMu(struct parameters *params, int mu);
	
	int getNumInputs(struct parameters *params);
	int getNumOutputs(struct parameters *params);
	
	/*
	
	*/
	void executeChromosome(struct parameters *params, struct chromosome *chromo, float *inputs, float *outputs);
	
	
	/*
	
	*/
	void setFuctionSet(struct parameters *params, char *functionNames);
	void printFuctionSet(struct parameters *params);
	
	struct population *initialisePopulation(struct parameters *params);
	
	/* */
	struct chromosome *initialiseChromosome(struct parameters *params);
	void freeChromosome(struct parameters *params, struct chromosome * chromo);
	
	void printChromosome(struct parameters *params, struct chromosome *chromo);
	void mutateChromosome(struct parameters *params, struct chromosome *chromo);
	
	void setFitnessFuction(struct parameters *params, float (*fitnessFuction)(struct parameters *params, struct chromosome *chromo));
	
	
#endif 
