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
* \file gaussian_kernel.c
* \author Haldo Spontón <haldos@fing.edu.uy> & Juan Cardelino <juanc@fing.edu.uy>
* \date May, 2012
* \see ``Review of edge detectors´´ IPOL publication.
* \brief Implements the generation of both Gaussian and LoG (Laplacian of a Gaussian) kernels. This function also manages the allocation and liberation of memory used by kernels.
*/

//  Software Guide : BeginLatex
//	This file implements the necessary functions for generating 
//	Gaussian and LoG (Laplacian of a Gaussian) kernels. It is also 
//	done here the alloc and free of the memory needed. \\
//
//	Includes:
//  Software Guide : EndLatex 

// Software Guide : BeginCodeSnippet
	#include <math.h> // exp
	#include <stdlib.h> // malloc, calloc, free
	#include <stdio.h> // fprintf
// Software Guide : EndCodeSnippet

//  Software Guide : BeginLatex
//	\vspace{0.5cm}
//	\Large{Function \texttt{gaussian\_kernel}} \\
// 
//	\normalsize
//	This function generates a Gaussian kernel of size $n\times n$ and standard deviation $\sigma$. \\
//  Software Guide : EndLatex
//  Software Guide : BeginCodeSnippet
double *gaussian_kernel(int n, float sigma){
//  Software Guide : EndCodeSnippet
	
//  Software Guide : BeginLatex
//	Memory allocation for the kernel: \\
//  Software Guide : EndLatex
//  Software Guide : BeginCodeSnippet
	double *kernel = malloc(n*n*sizeof(double));
//  Software Guide : EndCodeSnippet
	if (kernel == NULL){
		fprintf(stderr, "Out of memory...\n");
		exit(EXIT_FAILURE);
	}
	
	// Generation of the kernel
//  Software Guide : BeginLatex
//	We generate a normalized Gaussian kernel using the expression $e^{-\frac{x^2+y^2}{2\sigma^2}}$. \\
//  Software Guide : EndLatex
//  Software Guide : BeginCodeSnippet
	int i;
	int fila, col, x, y;
	double suma = 0;
	int imax = n*n;
	for(i=0;i<imax;i++){
		fila = (int)(i/n);
		col = i-(n*fila);
		y = ((int)(n/2))-fila;
		x = col-((int)(n/2));
		kernel[i] = exp(-(x*x + y*y)/(2*sigma*sigma));
		suma += kernel[i];
	}
	for(i=0;i<n*n;i++){
		kernel[i] = kernel[i]/suma;
	}
//  Software Guide : EndCodeSnippet
	
//  Software Guide : BeginLatex
//	Return: \\
//  Software Guide : EndLatex
//  Software Guide : BeginCodeSnippet
	return kernel;
//  Software Guide : EndCodeSnippet
}

//  Software Guide : BeginLatex
//	\vspace{0.5cm}
//	\Large{Function \texttt{free\_gaussian\_kernel}} \\
// 
//	\normalsize
//	This function frees the memory allocated in the function \texttt{gaussian\_kernel}. It receives as parameter the pointer to the array to be freed. \\
//  Software Guide : EndLatex
//  Software Guide : BeginCodeSnippet
void free_gaussian_kernel(double* kernel){
//  Software Guide : EndCodeSnippet

//  Software Guide : BeginLatex
//	Free memory: \\
//  Software Guide : EndLatex
//  Software Guide : BeginCodeSnippet
	free(kernel);
//  Software Guide : EndCodeSnippet

}

//  Software Guide : BeginLatex
//	\vspace{0.5cm}
//	\Large{Function \texttt{LoG\_kernel}} \\
// 
//	\normalsize
//	This function generates a Laplacian of a Gaussian kernel (LoG kernel) of size $n\times n$ and standard deviation $\sigma$. \\
//  Software Guide : EndLatex
//  Software Guide : BeginCodeSnippet
double *LoG_kernel(int n, float sigma){
//  Software Guide : EndCodeSnippet
	
//  Software Guide : BeginLatex
//	Memory allocation for the kernel: \\
//  Software Guide : EndLatex
//  Software Guide : BeginCodeSnippet
	double *kernel = malloc(n*n*sizeof(double));
//  Software Guide : EndCodeSnippet
	if (kernel == NULL){
		fprintf(stderr, "Out of memory...\n");
		exit(EXIT_FAILURE);
	}
	
	// Generation of the kernel
//  Software Guide : BeginLatex
//	We generate a Laplacian of a Gaussian kernel using the expression $\frac{x^2 + y^2 - 2\sigma^2}{\sigma^4}e^{-\frac{x^2+y^2}{2\sigma^2}}$. \\
//  Software Guide : EndLatex
//  Software Guide : BeginCodeSnippet
	int i;
	int fila, col, x, y;
	double suma = 0;
	int imax = n*n;
	for(i=0;i<imax;i++){
		fila = (int)(i/n);
		col = i-(n*fila);
		y = ((int)(n/2))-fila;
		x = col-((int)(n/2));
		kernel[i] = ( (x*x + y*y - 2*sigma*sigma)/(sigma*sigma*sigma*sigma) )
					*exp(-(x*x + y*y)/(2*sigma*sigma));
	}
//  Software Guide : EndCodeSnippet

	//for(i=0;i<n*n;i++){
	//	kernel[i] = kernel[i]/suma;
	//	fprintf(stderr, "kernel[%d] = %f \n", i, kernel[i]);
	//}
	
//  Software Guide : BeginLatex
//	Return: \\
//  Software Guide : EndLatex
//  Software Guide : BeginCodeSnippet
	return kernel;
//  Software Guide : EndCodeSnippet
}
