// Copyright (C) 2011-2012, Haldo Spontón <haldos@fing.edu.uy>
// Copyright (C) 2011-2012, Juan Cardelino <juanc@fing.edu.uy>

//		+------------------------------------------+
//		| Convolve image with kernel (header file) |
//		| Implemented by Haldo Spontón		       |
//		+------------------------------------------+

double *conv2d(double *input, int w, int h, double *kernel, int n);
// inputs:
//			double *input	-	pointer to input image
//			int w, int h	-	width and height of input image
//			double *kernel	-	pointer to kernel
//			int n			-	size of kernel
// output:
//			double *xxxx	-	convolved image (size: (w+n-1)*(h+n-1) )

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
