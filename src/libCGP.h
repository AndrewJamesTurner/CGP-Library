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
	
	
	
	/* */
	struct parameters *initialiseParameters(void);
	
	/*
		getters and setters for the parameters.
	*/
	int getMu(struct parameters *params);
	
#endif 
