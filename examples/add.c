
float add(const int numInputs, const double *inputs, const double *weights){
	
	int i;
	double sum = 0;
	
	for(i=0; i<numInputs; i++){
		sum += inputs[i];
	}
	
	return sum;
}	
