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

//		+----------------------------------------------+
//		| Edge detection algorithms (Roberts, Prewitt, |
//		| Sobel, Haralick, Marr-Hildreth)			   |
//		| Implemented by Haldo Spontón		           |
//		+----------------------------------------------+

// Roberts edge detector
float *edges_roberts(double *im, int w, int h, float threshold, int padding_method);
// inputs:
//			double *input		-	pointer to input image
//			int w, int h		-	width and height of input image
//			float threshold		-	threshold of edge detection
//			int padding_method	-	padding method for convolution
// output:
//			float *				-	pointer to output image

// Prewitt edge detector
float *edges_prewitt(double *im, int w, int h, float threshold, int padding_method);
// inputs:
//			double *input		-	pointer to input image
//			int w, int h		-	width and height of input image
//			float threshold		-	threshold of edge detection
//			int padding_method	-	padding method for convolution
// output:
//			float *				-	pointer to output image

// Sobel edge detector
float *edges_sobel(double *im, int w, int h, float threshold, int padding_method);
// inputs:
//			double *input		-	pointer to input image
//			int w, int h		-	width and height of input image
//			float threshold		-	threshold of edge detection
//			int padding_method	-	padding method for convolution
// output:
//			float *				-	pointer to output image

// Haralick edge detector
float *edges_haralick(double *im, int w, int h, float rhozero, int padding_method);
// inputs:
//			double *input		-	pointer to input image
//			int w, int h		-	width and height of input image
//			float rhozero		-	
//			int padding_method	-	padding method for convolution
// output:
//			float *				-	pointer to output image

// Marr-Hildreth edge detector with Gaussian and Laplacian kernels
float *edges_mh(double *im, int w, int h, float sigma, int n, float tzc, int padding_method);
// inputs:
//			double *input		-	pointer to input image
//			int w, int h		-	width and height of input image
//			float sigma			-	gaussian standard deviation
//			int n				-	kernel size
//			float tzc			-	threshold in zero-crossing
//			int padding_method	-	padding method for convolution
// output:
//			float *				-	pointer to output image

// Marr-Hildreth edge detector with Laplacian of Gaussian kernel
float *edges_mh_log(double *im, int w, int h, float sigma, int n, float tzc, int padding_method);
// inputs:
//			double *input		-	pointer to input image
//			int w, int h		-	width and height of input image
//			float sigma			-	gaussian standard deviation
//			int n				-	kernel size
//			float tzc			-	threshold in zero-crossing
//			int padding_method	-	padding method for convolution
// output:
//			float *				-	pointer to output image
