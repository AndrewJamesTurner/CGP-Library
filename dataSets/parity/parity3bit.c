
#include <math.h>
#include <cgp.h>


#define PARITY 3
#define SAMPLES (int)pow(2,PARITY)

void main(void){
	
	int a,b,c;
	int sum;
	int sample = 0;
	
	
	struct dataSet *data;
	
	float inputs[SAMPLES][PARITY];
	float outputs[SAMPLES][1];
	
	
	for(a=0; a<2; a++){
		for(b=0; b<2; b++){
			for(c=0; c<2; c++){
			
				inputs[sample][0] = a;
				inputs[sample][1] = b;
				inputs[sample][2] = c;
			
				sum = a+b+c;
				
				outputs[sample][0] = sum % 2;
				
				sample++;
			
			}
		}
	}
	
		
	data = initialiseDataSetFromArrays(PARITY, 1, SAMPLES, inputs[0], outputs[0]);
	
	printDataSet(data);
	
	saveDataSet(data, "parity3bit.data");
	
	freeDataSet(data);
	
}
