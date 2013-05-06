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
//		| Convolve image with kernel (header file) |
//		| Implemented by Haldo Spontón		       |
//		+------------------------------------------+

double *conv2d(double *input, int w, int h, double *kernel, int n, int padding_method);
// inputs:
//			double *input		-	pointer to input image
//			int w, int h		-	width and height of input image
//			double *kernel		-	pointer to kernel
//			int n				-	size of kernel
//			int padding_method	-	padding method: 0 for zero-padding, 1 for image boundary reflection
// output:
//			double *xxxx		-	convolved image (size: (w+n-1)*(h+n-1) )

double *get_neighborhood(double *im, int pos, int n, int* offsets);
// inputs:
//			double *im		-	input image
//			int pos			-	position index
//			int n			-	size of neighborhood
//			int* offsets	-	offsets of the n*n neighbors (generated with 
//								function get_neighbors_offset).
// output:
//			double *yyyy	-	neighborhood array (n*n length)

int *get_neighbors_offset(int w, int n);
// inputs:
//			int w			-	width of the image
//			int n			-	size of neighborhood
// output:
//			int *zzzz		-	offsets of the n*n neighbors.

void free_neighborhood(double* neighborhood);
// input:
//			double* neighborhood	-	neighborhood array (n*n length)
//										produced by function get_neighborhood.

void free_neighbors_offsets(int* offsets);
// input:
//			int* offsets	-	offsets of the n*n neighbors generated
//								with function get_neighbors_offset.
