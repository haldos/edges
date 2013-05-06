// Copyright (C) 2011-2012, Haldo Spontón <haldos@fing.edu.uy>
// Copyright (C) 2011-2012, Juan Cardelino <juanc@fing.edu.uy>

/*	
	This program is free software: you can use, modify and/or
	redistribute it under the terms of the GNU General Public
	License as published by the Free Software Foundation, either
	version 3 of the License, or (at your option) any later
	version. You should have received a copy of this license along
	this program. If not, see <http://www.gnu.org/licenses/>.
	
	This program is free software: you can use, modify and/or
	redistribute it under the terms of the simplified BSD
	License. You should have received a copy of this license along
	this program. If not, see
	<http://www.opensource.org/licenses/bsd-license.html>.
	
	This program is provided for research and education only: you can
	use and/or modify it for these purposes, but you are not allowed
	to redistribute this work or derivative works in source or
	executable form. A license must be obtained from the patent right
	holders for any other use.
*/

//		+------------------------------------------+
//		| Generate a gaussian kernel (header file) |
//		| Implemented by Haldo Spontón		       |
//		+------------------------------------------+


double *gaussian_kernel(int n, float sigma);
// input:
//			int n			-	size of kernel
//			float sigma		-	standard deviation of gaussian function
// output:
//			double *xxxx	-	n*n kernel array

void free_gaussian_kernel(double* kernel);
// input:
//			double* kernel	-	n*n kernel array, generated with function
//								gaussian_kernel.
