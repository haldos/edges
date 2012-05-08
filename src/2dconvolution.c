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

/** @file
* \file 2dconvolution.c
* \author Haldo Spontón <haldos@fing.edu.uy> & Juan Cardelino <juanc@fing.edu.uy>
* \date May, 2012
* \see ``Review of edge detectors´´ IPOL publication.
* \brief This file implements the required functions to compute the convolution between an input array and a kernel array. This function also manages the allocation and liberation of memory used by neighborhoods and index offset arrays.
*/

//  Software Guide : BeginLatex
//  This file contains the necessary functions to compute the convolution between an 
//  input array and a kernel (usually the input array is larger than the kernel). \\
//
//	Includes:
//  Software Guide : EndLatex 

// Software Guide : BeginCodeSnippet
#include <stdlib.h> // malloc, calloc, free
#include <stdio.h> // fprintf
// Software Guide : EndCodeSnippet

/*!
 \fn void free_neighborhood(double* neighborhood)
 \brief This function frees the memory occupied by an array containing a neighborhood of a pixel.
 @param *neighborhood Pointer to the neighborhood array.
 \return None.
*/
//  Software Guide : BeginLatex
//	\vspace{0.5cm}
//	\Large{Function \texttt{free\_gaussian\_kernel}} \\
// 
//	\normalsize
//	This function frees the memory allocated in the function \texttt{get\_neighborhood}. It receives as parameter the pointer to the array to be freed. \\
//  Software Guide : EndLatex
//  Software Guide : BeginCodeSnippet
void free_neighborhood(double* neighborhood){
//  Software Guide : EndCodeSnippet
	
//  Software Guide : BeginLatex
//	Free memory: \\
//  Software Guide : EndLatex
//  Software Guide : BeginCodeSnippet
	free(neighborhood);
//  Software Guide : EndCodeSnippet

}

/*!
 \fn void free_neighbors_offsets(int* offsets)
 \brief This function frees the memory occupied by an array containing a neighbors offset index of a pixel.
 @param *offsets Pointer to the neghbors offset array.
 \return None.
*/
//  Software Guide : BeginLatex
//	\vspace{0.5cm}
//	\Large{Function \texttt{free\_gaussian\_kernel}} \\
// 
//	\normalsize
//	This function frees the memory allocated in the function \texttt{get\_neighbors\_offset}. It receives as parameter the pointer to the array to be freed. \\
//  Software Guide : EndLatex
//  Software Guide : BeginCodeSnippet
void free_neighbors_offsets(int* offsets){
//  Software Guide : EndCodeSnippet

//  Software Guide : BeginLatex
//	Free memory: \\
//  Software Guide : EndLatex
//  Software Guide : BeginCodeSnippet
	free(offsets);
//  Software Guide : EndCodeSnippet

}

/*!
 \fn int *get_neighbors_offset(int w, int n)
 \brief This function calculates the offset index for neighbors of a pixel.
 @param w Image width.
 @param n Size of the neighborhood (n x n).
 \return Integer array of size nxn containing the offtests of each neighbor pixel.
*/
//  Software Guide : BeginLatex
//	\vspace{0.5cm}
//	\Large{Function \texttt{get\_neighbors\_offset}} \\
// 
//	\normalsize
//	This function computes the linear array index offsets required to access the neighbor pixels in 
//	a neighborhood of size $n\times n$. This computation requires to know the width $w$ of the image. \\
//  Software Guide : EndLatex
//  Software Guide : BeginCodeSnippet
int *get_neighbors_offset(int w, int n){
//  Software Guide : EndCodeSnippet
	
//  Software Guide : BeginLatex
//	Memory allocation for the index offsets array: \\
//  Software Guide : EndLatex
//  Software Guide : BeginCodeSnippet
	int *aux = malloc(n*n*sizeof(int));
//  Software Guide : EndCodeSnippet
	int i, delta_fila, delta_col;
//  Software Guide : BeginLatex
//	Computation of the index offsets array: \\
//  Software Guide : EndLatex
//  Software Guide : BeginCodeSnippet
	int imax = n*n;
	int dif_fila_col = (n-1)/2;
	for (i=0;i<imax;i++){
		delta_fila = (int)(i/n)-dif_fila_col;
		delta_col = i-n*((int)(i/n))-dif_fila_col;
		aux[i] = delta_fila*w+delta_col;
	}
//  Software Guide : EndCodeSnippet

//  Software Guide : BeginLatex
//	Return: \\
//  Software Guide : EndLatex
//  Software Guide : BeginCodeSnippet
	return aux;
//  Software Guide : EndCodeSnippet
}

/*!
 \fn double *get_neighborhood(double *im, int pos, int n, int* offsets)
 \brief This function returns the neighborhood of a pixel, using the neighbors offset computed in \link get_neighbors_offset \endlink.
 @param *im Image array.
 @param pos Pixel position in the image array.
 @param n Size of the neighborhood (n x n).
 @param *offsets Neigbors index offset array.
 \return Double array of size nxn containing the neighborhood of the pixel.
*/
//  Software Guide : BeginLatex
//	\vspace{0.5cm}
//	\Large{Function \texttt{get\_neighborhood}} \\
// 
//	\normalsize
//	This function returns an array of size $n\times n$ containing the neighbors of the 
//  pixel at the position \texttt{pos} of the image \texttt{im}. For this calculation, 
//	the index offsets array calculated with the function \texttt{get\_neighbors\_offset} 
//	is needed.
//  Software Guide : EndLatex
//  Software Guide : BeginCodeSnippet
double *get_neighborhood(double *im, int pos, int n, int* offsets){
//  Software Guide : EndCodeSnippet
	
//  Software Guide : BeginLatex
//	Memory allocation for the neighborhood array: \\
//  Software Guide : EndLatex
//  Software Guide : BeginCodeSnippet
	double *aux = malloc(n*n*sizeof(double));
//  Software Guide : EndCodeSnippet
	if (aux == NULL){
		fprintf(stderr, "Out of memory...\n");
		exit(EXIT_FAILURE);
	}

//  Software Guide : BeginLatex
//	Assigning values ​​of the neighborhood: \\
//  Software Guide : EndLatex
//  Software Guide : BeginCodeSnippet
	int i;
	int imax = n*n;
	for (i=0;i<imax;i++){
		aux[i] = im[pos+offsets[i]];
	}
//  Software Guide : EndCodeSnippet
	
//  Software Guide : BeginLatex
//	Return: \\
//  Software Guide : EndLatex
//  Software Guide : BeginCodeSnippet
	return aux;
//  Software Guide : EndCodeSnippet
}

/*!
 \fn double *conv2d(double *input, int w, int h, double *kernel, int n)
 \brief This function calculates the convolution of input image with the input kernel.
 @param *input Input image array.
 @param w Input image width.
 @param h Input image heigth.
 @param *kernel Kernel array used in the convolution.
 @param n Size of the kernel (n x n).
 \return Double array of size (w+n-1)x(h+n-1), output of the convolution between input image and kernel, and using zero-padding.
*/
//  Software Guide : BeginLatex
//	\vspace{0.5cm}
//	\Large{Function \texttt{conv2d}} \\
// 
//	\normalsize
//	This function calculates the convolution of image \texttt{input} (size $w\times h$) with the kernel \texttt{kernel} (size $n\times n$).
//	Returns an image of size $(w+n-1)\times(h+n-1)$ (due to zero padding).
//  Software Guide : EndLatex
//  Software Guide : BeginCodeSnippet
double *conv2d(double *input, int w, int h, double *kernel, int n){
//  Software Guide : EndCodeSnippet
	
	// Zero-padding
//  Software Guide : BeginLatex
//	Zero padding: reserve memory for an array of size $(w+2(n-1))\times(h+2(n-1))$. \\
//  Software Guide : EndLatex
//  Software Guide : BeginCodeSnippet
	int wx = (w+2*(n-1));
	int hx = (h+2*(n-1));
	double *aux = calloc(wx*hx,sizeof(double));
//  Software Guide : EndCodeSnippet
	if (aux == NULL){
		fprintf(stderr, "Out of memory...\n");
		exit(EXIT_FAILURE);
	}

//  Software Guide : BeginLatex
//	Fill in the values ​​of the image \texttt{aux}, centering the original image \texttt{input} on it. \\
//  Software Guide : EndLatex
//  Software Guide : BeginCodeSnippet
	int i,j,fila,col;
	int imax = wx*hx;
	for(i=0;i<imax;i++){
		fila = (int)(i/wx);
		col = i-(wx*fila);
		if ( (fila>=n-1)&&(col>=n-1)&&(fila<h+n-1)&&(col<w+n-1) ) {
			aux[i] = input[(col-n+1)+(w*(fila-n+1))];
		}
	}
//  Software Guide : EndCodeSnippet

//  Software Guide : BeginLatex
//	Reserve memory for the output array of size $(w+n-1)\times(h+n-1)$. \\
//  Software Guide : EndLatex
//  Software Guide : BeginCodeSnippet
	double *out = malloc((w+n-1)*(h+n-1)*sizeof(double));
//  Software Guide : EndCodeSnippet
	if (out == NULL){
		fprintf(stderr, "Out of memory...\n");
		exit(EXIT_FAILURE);
	}

//  Software Guide : BeginLatex
//	Compute the convolution. Most of the operations are intended to calculate the relative 
//	positions between the images \texttt{aux} and \texttt{out}. \\
//  Software Guide : EndLatex
//  Software Guide : BeginCodeSnippet
	double acum = 0;
	int pos;
	// Convolution
	imax = (w+n-1)*(h+n-1);
	int jmax = n*n;
	int *offsets = get_neighbors_offset(wx, n);
	int dif_fila_col = (n-1)/2;
	for(i=0;i<imax;i++){
		fila = (int)(i/(w+n-1));
		col = i-((w+n-1)*fila);
		fila += dif_fila_col;
		col += dif_fila_col;
		pos = wx*fila + col;
		// compute convolution:
		acum = 0;
		for (j=0;j<jmax;j++){
			acum += aux[pos+offsets[j]]*kernel[j];
		}
		out[i] = acum;
	}
//  Software Guide : EndCodeSnippet
	
//  Software Guide : BeginLatex
//	Free and return: \\
//  Software Guide : EndLatex
//  Software Guide : BeginCodeSnippet
	free(aux);
	free_neighbors_offsets(offsets);
	return out;
//  Software Guide : EndCodeSnippet

}
