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

#include <stdio.h>
#include <stdlib.h> 
#include "libCGP.h"

struct parameters{
		int mu;
		int lambda;
	};

/*
	Initialises a parameter structs with default values. These
	values can be indevidually chaged via set fuctions.
*/
struct parameters *initialiseParameters(void){
		
	struct parameters *params;
	
	params = malloc(sizeof(struct parameters));
		
	params->mu = 1;
	params->lambda = 4;
	
	return params;
}


int getMu(struct parameters *params){
	return params -> mu;
}

