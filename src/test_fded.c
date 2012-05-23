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
* \file test_fded.c
* \author Haldo Spontón <haldos@fing.edu.uy> & Juan Cardelino <juanc@fing.edu.uy>
* \date May, 2012
* \see ``Review of edge detectors´´ IPOL publication.
* \brief Implements some first derivative edge detection algorithms: Roberts, Prewitt and Sobel. In brief, this program opens an image, and convolves it with several operators. Then thresholded images are generated for each of the methods.
*/

//  Software Guide : BeginLatex
//  First derivative edge detectors (Roberts, Prewitt and Sobel), main C file.\\
//  
//  Parameters:
//	\begin{itemize}
//		\item \texttt{input\_image} - Input image.
//		\item \texttt{threshold} - Threshold for gradient image ($0\leq th \leq 1$).
//	\end{itemize}
//
//	\textit{Note: Output images are saved with the filenames ``\texttt{roberts.png}'', ``\texttt{prewitt.png}'' and ``\texttt{sobel.png}''.} \\
//
//	Includes:
//  Software Guide : EndLatex 

// Software Guide : BeginCodeSnippet
	#include "iio.c"
	#include "2dconvolution.c"
	#include <time.h>
// Software Guide : EndCodeSnippet

//  Software Guide : BeginLatex
//  Macros:
//  Software Guide : EndLatex
// Software Guide : BeginCodeSnippet
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define THRESHOLD(x, th) (((x) > (th)) ? (255) : (0))
// Software Guide : EndCodeSnippet

/*!
 \fn int main_fded(int argc, char *argv[])
 \brief Main function of the Roberts, Prewitt and Sobel edge detection algorithms.
 @param input_image Input image filename.
 @param threshold Gradient threshold used to generate binary edges images.
 \return None.
 \ingroup fded
	\note The real name of this function is main. Function name was temporarily changed to the proper functioning of doxygen.
*/

//  Software Guide : BeginLatex
//	\vspace{0.5cm}
//	\Large{Main function} \\
//  Software Guide : EndLatex
// Software Guide : BeginCodeSnippet
int main(int argc, char *argv[]) {
// Software Guide : EndCodeSnippet

	if (argc != 3) {
		printf("Usage: %s input_image threshold\n", argv[0]);
	} else {

		// Execution time:
		double start = (double)clock();

		// Load input image (using iio)
//  Software Guide : BeginLatex
//	Load input image (using \textit{iio}): \\
//  Software Guide : EndLatex
// Software Guide : BeginCodeSnippet
		int w, h, pixeldim;
		float *im_orig = iio_read_image_float_vec(argv[1], &w, &h, &pixeldim);
// Software Guide : EndCodeSnippet
		fprintf(stderr, "Input image loaded:\t %dx%d image with %d channel(s).\n", w, h, pixeldim);

		// Grayscale conversion (if necessary)
//  Software Guide : BeginLatex
//	Grayscale conversion (if necessary): explained in \ref{app:marr-hildreth}. \\ \\
//  Software Guide : EndLatex
		double *im = malloc(w*h*sizeof(double));
		if (im == NULL){
			fprintf(stderr, "Out of memory...\n");
			exit(EXIT_FAILURE);
		}
			// allocate memory for the grayscale image <im>, output of the grayscale conversion
			// and correct allocation check.
		int z;
			// <z> is just an integer used as array index.
		int zmax = w*h;		// number of elements of <im>
		if (pixeldim==3){	// if the image is color (RGB, three channels)...
			for(z=0;z<zmax;z++){		// for each pixel in the image <im>, calculate the gray 
										// value according to the expression: 
										// I = ( 6968*R + 23434*G + 2366*B ) / 32768.
				im[z] =  (double)(6968*im_orig[3*z] + 23434*im_orig[3*z + 1] + 2366*im_orig[3*z + 2])/32768;
			}
			fprintf(stderr, "images converted to grayscale\n");
		} else {		// the image was originally grayscale...
			for(z=0;z<zmax;z++){
				im[z] = (double)im_orig[z];		// only assign the value of im_orig to im, casting to double.
			}
			fprintf(stderr, "images are already in grayscale\n");
		}

		// Define operators:
//  Software Guide : BeginLatex
//	Define the Roberts, Prewitt and Sobel operators: \\
//	(We use $3\times 3$ operators)
//	\begin{itemize}
//		\item	Roberts:
//				$$
//				R_1 = \begin{bmatrix} -1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 0 \end{bmatrix}
//				$$
//				$$
//				R_2 = \begin{bmatrix} 0 & -1 & 0 \\ 1 & 0 & 0 \\ 0 & 0 & 0 \end{bmatrix}
//				$$
//		\item	Prewitt:
//				$$
//				P_1 = \begin{bmatrix} -1 & -1 & -1 \\ 0 & 0 & 0 \\ 1 & 1 & 1 \end{bmatrix}
//				$$
//				$$
//				P_2 = \begin{bmatrix} -1 & 0 & 1 \\ -1 & 0 & 1 \\ -1 & 0 & 1 \end{bmatrix}
//				$$
//		\item	Sobel:
//				$$
//				S_1 = \begin{bmatrix} -1 & -2 & -1 \\ 0 & 0 & 0 \\ 1 & 2 & 1 \end{bmatrix}
//				$$
//				$$
//				S_2 = \begin{bmatrix} -1 & 0 & 1 \\ -2 & 0 & 2 \\ -1 & 0 & 1 \end{bmatrix}
//				$$
//	\end{itemize}
//  Software Guide : EndLatex
// Software Guide : BeginCodeSnippet
		double roberts_1[9] = {-1, 0, 0, 0, 1, 0, 0, 0, 0};		// ROBERTS
		double roberts_2[9] = { 0,-1, 0, 1, 0, 0, 0, 0, 0};		// OPERATORS
		//---------------------------------------------------------------------------------
		double prewitt_1[9] = {-1,-1,-1, 0, 0, 0, 1, 1, 1};		// PREWITT
		double prewitt_2[9] = {-1, 0, 1,-1, 0, 1,-1, 0, 1};		// OPERATORS
		//---------------------------------------------------------------------------------
		double sobel_1[9] = {-1,-2,-1, 0, 0, 0, 1, 2, 1};		// SOBEL
		double sobel_2[9] = {-1, 0, 1,-2, 0, 2,-1, 0, 1};		// OPERATORS
// Software Guide : EndCodeSnippet

		// Convolve images:
//  Software Guide : BeginLatex
//	The input image is convolved with the defined operatos, using the \texttt{conv2d} function in \texttt{2dconvolution.c}:
//  Software Guide : EndLatex
// Software Guide : BeginCodeSnippet
		double *im_r1 = conv2d(im, w, h, roberts_1, 3);
		double *im_r2 = conv2d(im, w, h, roberts_2, 3);
		double *im_p1 = conv2d(im, w, h, prewitt_1, 3);
		double *im_p2 = conv2d(im, w, h, prewitt_2, 3);
		double *im_s1 = conv2d(im, w, h, sobel_1, 3);
		double *im_s2 = conv2d(im, w, h, sobel_2, 3);
// Software Guide : EndCodeSnippet

		// Allocate memory for final images:
//	Software Guide : BeginLatex
//	Allocate memory for final images:
//  Software Guide : EndLatex
// Software Guide : BeginCodeSnippet
		float *im_roberts = malloc(w*h*sizeof(float));
		float *im_prewitt = malloc(w*h*sizeof(float));
		float *im_sobel = malloc(w*h*sizeof(float));
// Software Guide : EndCodeSnippet

//	Software Guide : BeginLatex
//	For each method, two images are obtained (one for each operator). Then the first derivative magnitude image is constructed using $M=\sqrt{g_x^2+g_y^2}$. \\
//	Also the absolute maximum value of the constructed images is computed, for each method. \\
//	Software Guide : EndLatex
// Software Guide : BeginCodeSnippet
		int i,j, fila, col;
		double max_r = 0;
		double max_p = 0;
		double max_s = 0;
		int imax = w*h;
		for (i=0;i<imax;i++){
			fila = (int)(i/w);
			col = i - w*fila + 1;
			fila += 1;
			j = col + (w+2)*fila;
			// Max in each case
			im_roberts[i] = sqrt(im_r1[j]*im_r1[j] + im_r2[j]*im_r2[j]);
			im_prewitt[i] = sqrt(im_p1[j]*im_p1[j] + im_p2[j]*im_p2[j]);
			im_sobel[i] = sqrt(im_s1[j]*im_s1[j] + im_s2[j]*im_s2[j]);
			// Absolute max
			max_r = MAX(max_r,im_roberts[i]);
			max_p = MAX(max_p,im_prewitt[i]);
			max_s = MAX(max_s,im_sobel[i]);
		}
// Software Guide : EndCodeSnippet

		// Thresholding
//	Software Guide : BeginLatex
//	Thresholded images of each method is created, using the THRESHOLD macro: \\
//	Software Guide : EndLatex
// Software Guide : BeginCodeSnippet
		float th = atof(argv[2]);
		for (i=0;i<imax;i++){
			im_roberts[i] = THRESHOLD(im_roberts[i],th*max_r);
			im_prewitt[i] = THRESHOLD(im_prewitt[i],th*max_p);
			im_sobel[i] = THRESHOLD(im_sobel[i],th*max_s);
		}
// Software Guide : EndCodeSnippet

		// Save outout images
//  Software Guide : BeginLatex
//	Save output image (using \textit{iio}): \\
//  Software Guide : EndLatex
// Software Guide : BeginCodeSnippet
		iio_save_image_float_vec("roberts.png", im_roberts, w, h, 1);
		iio_save_image_float_vec("prewitt.png", im_prewitt, w, h, 1);
		iio_save_image_float_vec("sobel.png", im_sobel, w, h, 1);
// Software Guide : EndCodeSnippet
			fprintf(stderr, "Roberts image saved to roberts.png.\n");
			fprintf(stderr, "Prewitt image saved to prewitt.png.\n");
			fprintf(stderr, "Sobel image saved to sobel.png.\n");

		// Free memory:
		free(im_orig);
		free(im);
		free(im_r1);
		free(im_r2);
		free(im_p1);
		free(im_p2);
		free(im_s1);
		free(im_s2);
		free(im_roberts);
		free(im_prewitt);
		free(im_sobel);

		fprintf(stderr, "Edge detection algorithms based on first derivative computation done.\n");

		// Execution time:
		double finish = (double)clock();
		double exectime = (finish - start)/CLOCKS_PER_SEC;
		fprintf(stderr, "execution time: %1.3f s.\n", exectime);

		return 0;

	} // else (argc)

} // end of the program
