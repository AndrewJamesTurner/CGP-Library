#include <stdio.h>
#include "libCGP.h"

int main(void){

	struct parameters *params;

	params = initialiseParameters();
	
	printf("%d\n", getMu(params));
	
	
	return 1;
}
