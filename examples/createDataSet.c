/*
	This file is part of CGP-Library
	Copyright (c) Andrew James Turner 2014 (andrew.turner@york.ac.uk)

    CGP-Library is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published 
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CGP-Library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with CGP-Library.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <math.h>
#include "../include/cgp.h"  

#define NUMINPUTS 1
#define NUMOUTPUTS 1
#define NUMSAMPLES 101

/*
	Returns x^6 - 2x^4 + x^2
*/
float symbolicEq1(float x){
	return powf(x,6) - 2*powf(x,4) + powf(x,2);
}

int main(void){
	
	int i;
		
	struct dataSet *data;
	
	float inputs[NUMSAMPLES][NUMINPUTS];
	float outputs[NUMSAMPLES][NUMOUTPUTS];
	
	float inputTemp;
	float outputTemp;
			
	for(i=0; i<NUMSAMPLES; i++){
		
		inputTemp = (i * (10.0/(NUMSAMPLES-1))) - 5.0;
		outputTemp = symbolicEq1(inputTemp);
		
		inputs[i][0] = inputTemp;
		outputs[i][0] = outputTemp;
	}
	
	data = initialiseDataSetFromArrays(NUMINPUTS, NUMOUTPUTS, NUMSAMPLES, inputs[0], outputs[0]);
	
	saveDataSet(data, "symbolic.data");
	
	freeDataSet(data);
	
	return 1;
}
