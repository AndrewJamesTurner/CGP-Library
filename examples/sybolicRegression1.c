double symbolicEq1(double x){
	return pow(x,6.0) - 2*powf(x,4.0) + powf(x,2.0);
}


double symbolicRegression1(struct parameters *params, struct chromosome *chromo, struct data *dat){
	
	double i;
	double error = 0;
	double chromoInputs[1];	
	double chromoOutputs[1];	
						
	if(getNumInputs(params) != 1){
		printf("Error: The 'symbolicRegression1' fitness function requires one chromosome input; not %d\n", getNumInputs(params));
		exit(0);
	}				
		
	if(getNumOutputs(params) != 1){
		printf("Error: The 'symbolicRegression1' fitness function requires one chromosome output; not %d\n", getNumOutputs(params));
		exit(0);
	}		
					
	/* for each line in the truth table */				
	for(i=-5; i<=5; i=i+0.2){
		
		chromoInputs[0] = i;
		
		/* calculate the chromosome outputs for the set of inputs  */
		executeChromosome(params, chromo, chromoInputs, chromoOutputs);
		
		error = fabs(symbolicEq1(i) - chromoOutputs[0]);
	}				
					
	return error;
}
