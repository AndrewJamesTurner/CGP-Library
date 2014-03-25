
float fullAdder(struct parameters *params, struct chromosome *chromo, struct data *dat){
	
	int i;
	float error = 0;
	float chromoOutputs[8];	
	float inputs[8][3] = {{0,0,0},{0,0,1},{0,1,0},{0,1,1},{1,0,0},{1,0,1},{1,1,0},{1,1,1}};
	float outputs[8][2] = {{0,0},{1,0},{1,0},{0,1},{1,0},{0,1},{0,1},{1,1}};
		
					
	if(getNumInputs(params) != 3){
		printf("Error: The 'fullAdder' fitness function requires three chromosome inputs; not %d\n", getNumInputs(params));
		exit(0);
	}				
		
	if(getNumOutputs(params) != 2){
		printf("Error: The 'fullAdder' fitness function requires two chromosome outputs; not %d\n", getNumOutputs(params));
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
