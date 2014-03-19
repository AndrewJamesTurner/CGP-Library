#include <stdio.h>
#include <stdlib.h> 
#include "cgp.h"


float fullAdder(struct parameters *params, struct chromosome *chromo, struct data *dat){
	
	int i;
	float error = 0;
	float chromoOutputs[8];	
	float inputs[8][3] = {{0,0,0},{0,0,1},{0,1,0},{0,1,1},{1,0,0},{1,0,1},{1,1,0},{1,1,1}};
	float outputs[8][2] = {{0,0},{1,0},{1,0},{0,1},{1,0},{0,1},{0,1},{1,1}};
		
					
	if(getNumInputs(params) != 3){
		printf("Error: The 'fullAdder' fitness function requires three chromosome inputs; not %d\n", getNumInputs(params));
		printf("Terminating CGP-Library.\n");
		exit(0);
	}				
		
	if(getNumOutputs(params) != 2){
		printf("Error: The 'fullAdder' fitness function requires two chromosome outputs; not %d\n", getNumOutputs(params));
		printf("Terminating CGP-Library.\n");
		exit(0);
	}		
					
	/* for each line in the truth table */				
	for(i=0; i<8; i++){
		
		/* calculate the chromosome outputs for the set of inputs  */
		executeChromosome(params, chromo, inputs[i], chromoOutputs);
		
		/* If the chromosome outputs differ from the correct outputs increment the error */
		if(outputs[i][0] != chromoOutputs[0]){
			error++;
		}
		
		if(outputs[i][1] != chromoOutputs[1]){
			error++;
		}
	}				
					
	return error;
}



int main(void){

	struct parameters *params;
	struct population *pop;
	struct data *dat;
	float bestFit;
		
	params = initialiseParameters(3,60,2, 2);
	setFuctionSet(params, "and,or,xor,nand,nor");
	printFuctionSet(params);
	
	setFitnessFuction(params, fullAdder, "fullAdder");
	
	dat = initialiseDataFromFile("./example/fullAdder");
		
	pop = initialisePopulation(params);
	
	
	
	bestFit =  evolvePopulation(params, pop, dat);
	
	
	
	printf("bestFit: %f\n", bestFit);
	
	/*setFitnessFuction(params, fullAdder, "fullAdder");*/
	
	
		
	/*printChromosome(params, chromo);*/
	/*
	mutateChromosome(params, chromo);
	printChromosome(params, chromo);
	*/
	
	/*
	setMu(params, 1);
	printf("\nmu: %d\n", getMu(params));
	*/
	
	
	/*executeChromosome(params, chromo, inputs, outputs);*/
	
	/*printf("%f %f\n", outputs[0], outputs[1] );*/
	
	
	
	/*freeChromosome(params, chromo);*/
	
	
	freePopulation(params, pop);
	freeData(dat);
	freeParameters(params);
	
	return 1;
}
