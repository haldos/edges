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
* \file test_mh.c
* \author Haldo Spontón <haldos@fing.edu.uy> & Juan Cardelino <juanc@fing.edu.uy>
* \date May, 2012
* \see ``Review of edge detectors´´ IPOL publication.
* \brief Implements the Marr-Hildreth edge detection algorithm, using a gaussian kernel. In summary, this function opens the image, converts it to grayscale (if necessary), then the convolves it first with a Gaussian kernel and then with an operator that approximates the calculation of the Laplacian, and finally, find zero crossings in the resulting image.
*/

//  Software Guide : BeginLatex
//  Marr-Hildreth edge detector, main C file.\\
//  
//  Parameters:
//	\begin{itemize}
//		\item \texttt{input\_image} - Input image.
//		\item \texttt{sigma} - Standard deviation $\sigma$ of the Gaussian kernel.
//		\item \texttt{n} - Size $n$ of the Gaussian kernel ($n\times n$).
//		\item \texttt{tzc} - Threshold in the Zero-Crossing calculation ($0\leq t_{zc}\leq 1$).
//		\item \texttt{output\_image} - Output image (edges).
//	\end{itemize}
//
//	Includes:
//  Software Guide : EndLatex 

// Software Guide : BeginCodeSnippet
	#include "iio.c"
	#include "gaussian_kernel.c"
	#include "2dconvolution.c"
	#include <time.h>
// Software Guide : EndCodeSnippet

// Values for saving debug images.
#define SAVE_KERNEL false
#define SAVE_SMOOTHED_IMAGE false
#define SAVE_LAPLACIAN_IMAGE false

/*!
 \fn int main_mhg(int argc, char *argv[])
 \brief Main function of the Marr-Hildreth edge detection algorithm, using Gaussian kernel.
 @param input_image Filename of the input image.
 @param sigma Standard deviation of the Gaussian kernel.
 @param n Size of the kernel (n x n).
 @param tzc Threshold for the zero-crossing algorithm.
 @param output_image Filename of the output image.
 \return None.
 \ingroup mhg 
	\note The real name of this function is main. Function name was temporarily changed to the proper functioning of doxygen.
*/

//  Software Guide : BeginLatex
//	\vspace{0.5cm}
//	\Large{Main function} \\
//  Software Guide : EndLatex
// Software Guide : BeginCodeSnippet
int main(int argc, char *argv[]) {
// Software Guide : EndCodeSnippet

	if (argc != 6) {
		printf("Usage: %s input_image sigma n tzc output_image\n", argv[0]);
	} else {

		// Variable to measure execution time:
		double start = (double)clock();
	
		// Parameters
		float sigma = atof(argv[2]);
			// <sigma> is the standard deviation of the gaussian function used to
			// create the kernel.
		int n = atoi(argv[3]);
			// <n> is the size of the kernel (n*n).
		float tzc = atof(argv[4]);
			// <tzc> is the threshold of the zero-crossing method.
	
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
//	Grayscale conversion (if necessary): \\
//	First we allocate memory for the grayscale image \texttt{im}, with
//	the corresponding correct allocation check. Then we check the number of channels of the image: if
//	\texttt{pixeldim}=3, we assume RGB image and conversion is needed, else, we assume single channel image and
//	no conversion is required.\\
//	The computation of the gray intensity from RGB levels is:
//	$$
//	G = \frac{6968\times R + 23434\times G + 2366\times B}{32768} \\
//	$$
//  Software Guide : EndLatex
// Software Guide : BeginCodeSnippet
		double *im = malloc(w*h*sizeof(double));
		if (im == NULL){
			fprintf(stderr, "Out of memory...\n");
			exit(EXIT_FAILURE);
		}
		int z;
		int zmax = w*h;	
		if (pixeldim==3){
			for(z=0;z<zmax;z++){
				im[z] =  (double)(6968*im_orig[3*z] + 23434*im_orig[3*z + 1] 
													+ 2366*im_orig[3*z + 2])/32768;
			}
			fprintf(stderr, "images converted to grayscale\n");
		} else {
			for(z=0;z<zmax;z++){
				im[z] = (double)im_orig[z];
			}
			fprintf(stderr, "images are already in grayscale\n");
		}
// Software Guide : EndCodeSnippet

		// Generate gaussian kernel
		// see gaussian_kernel.c
//  Software Guide : BeginLatex
//	Generate Gaussian kernel using the \texttt{gaussian\_kernel} function in \texttt{gaussian\_kernel.c}: \\
//  Software Guide : EndLatex
// Software Guide : BeginCodeSnippet
		double *kernel = gaussian_kernel(n,sigma);
// Software Guide : EndCodeSnippet

		// Debug: save kernel to image (this is not part of the algorithm itself)
		if (SAVE_KERNEL){
			float *kernel_float = malloc(n*n*sizeof(float));
			if (kernel_float == NULL){
				fprintf(stderr, "Out of memory...\n");
				exit(EXIT_FAILURE);
			}			
			int i;
			int imax = n*n;
			for (i=0;i<imax;i++){
				kernel_float[i] = 5000*(float)kernel[i];
			}
			iio_save_image_float_vec("kernel.png", kernel_float, n, n, 1);
			free(kernel_float);
			fprintf(stderr, "kernel saved to kernel.png\n");
		}
		// end of save kernel image
		
		// Smooth input image with gaussian kernel
		// see 2dconvolution.c
		// <im_smoothed> is calculated convolving the grayscale image <im> with
		// the gaussian kernel previously generated.
//  Software Guide : BeginLatex
//	Smooth input image with the Gaussian kernel previously generated, using the \texttt{conv2d} function in \texttt{2dconvolution.c}: \\
//  Software Guide : EndLatex
// Software Guide : BeginCodeSnippet
		double *im_smoothed = conv2d(im, w, h, kernel, n);
// Software Guide : EndCodeSnippet

		// Debug: save smoothed image (this is not part of the algorithm itself)
		if (SAVE_SMOOTHED_IMAGE){
			//float *smoothed = malloc((w+n-1)*(h+n-1)*sizeof(float));
			float *smoothed = malloc(w*h*sizeof(float));
			if (smoothed == NULL){
				fprintf(stderr, "Out of memory...\n");
				exit(EXIT_FAILURE);
			}
			int i,j, fila, col;
			int imax = w*h;
			int dif_fila_col = (n-1)/2;
			for (i=0;i<imax;i++){
				fila = (int)(i/w);
				col = i - w*fila + dif_fila_col;
				fila += dif_fila_col;
				j = col + (w+n-1)*fila;
				smoothed[i] = (float)im_smoothed[j];
			}
			iio_save_image_float_vec("smoothed.png", smoothed, w, h, 1);
			free(smoothed);
			fprintf(stderr, "smoothed image saved to smoothed.png\n");
		}
		// end of save smoothed image

		// Laplacian of the smoothed image
//  Software Guide : BeginLatex
//	Computation of the Laplacian of the smoothed image: \\
//	We use a $3\times 3$ approximation of the laplacian operator: \\
//	$$
//	\begin{bmatrix}
//		1 &  1 & 1 \\
//		1 & -8 & 1 \\
//		1 &  1 & 1 
//	\end{bmatrix}
//	$$
//	Using the \texttt{conv2d} function we obtain the \texttt{laplacian} image. \\
//  Software Guide : EndLatex
// Software Guide : BeginCodeSnippet
		double operator[9] = {1, 1, 1, 1, -8, 1, 1, 1, 1};
		double *laplacian = conv2d(im_smoothed, w+n-1, h+n-1, operator, 3);
// Software Guide : EndCodeSnippet
		
		// calculate max absolute value of laplacian:
		// required for thresholding in zero-crossing (next)
//  Software Guide : BeginLatex
//	Now we calculate the maximum absolute value of the \texttt{laplacian} image. This value is required
//	for thresholding in zero-crossing calculation. \\
//  Software Guide : EndLatex
// Software Guide : BeginCodeSnippet
		double max_l = 0;
		int p;
		int pmax = (w+n+1)*(h+n+1);
		for (p=0;p<pmax;p++){
			if (abs(laplacian[p])>max_l){
				max_l = abs(laplacian[p]);
			}
		}
// Software Guide : EndCodeSnippet

		// Debug: save laplacian image (this is not part of the algorithm itself)
		if (SAVE_LAPLACIAN_IMAGE){
			//float *lapl = malloc((w+n+1)*(h+n+1)*sizeof(float));
			float *lapl = malloc(w*h*sizeof(float));
			if (lapl == NULL){
				fprintf(stderr, "Out of memory...\n");
				exit(EXIT_FAILURE);
			}
			int i,j, fila, col;
			int imax = w*h;
			int dif_fila_col = (n+1)/2;
			for (i=0;i<imax;i++){
				fila = (int)(i/w);
				col = i - w*fila + dif_fila_col;
				fila += dif_fila_col;
				j = col + (w+n+1)*fila;
				lapl[i] = (float)laplacian[j];
			}
			iio_save_image_float_vec("laplacian.png", lapl, w, h, 1);
			free(lapl);
			fprintf(stderr, "laplacian image saved to laplacian.png\n");
		}
		// end of save laplacian image

		// Zero-crossing
//  Software Guide : BeginLatex
//	Zero-crossing: \\
//	We explore the image, looking in every pixel a change of sign between neighboring opposite pixels.
// 	In every pixel $p$ we use the funcion \texttt{get\_neighborhood} (in \texttt{2dconvolution.c}) to get the
//	9 pixels of its neighborhood:
//	$$
//	\begin{bmatrix}
//		p_{up,left}		& p_{up,middle}		& p_{up,right}		\\
//		p_{middle,left}	& p					& p_{middle,right}	\\
//		p_{down,left}	& p_{down,middle}	& p_{down,right}	
//	\end{bmatrix}
//	$$
//	Then the pixel $p$ is marked as edge pixel if it occurs that:
//	\begin{itemize}
//	\item $sign(p_{up,left}) \neq sign(p_{down,right})$, or 
//	\item $sign(p_{up,middle}) \neq sign(p_{down,middle})$, or 
//	\item $sign(p_{up,right}) \neq sign(p_{down,left})$, or
//	\item $sign(p_{middle,left}) \neq sign(p_{middle,right})$. \\
//	\end{itemize}
//  Software Guide : EndLatex
// Software Guide : BeginCodeSnippet
		float *zero_cross = calloc(w*h,sizeof(float));		
						// this image will only content values 0 and 255
						// but we use float for saving using iio.
		if (zero_cross == NULL){
			fprintf(stderr, "Out of memory...\n");
			exit(EXIT_FAILURE);
		}
		int ind_en_lapl, fila, col;
		int *offsets = get_neighbors_offset(w+n+1, 3);
		pmax = w*h;
		int dif_fila_col = (n+1)/2;
		for (p=0;p<pmax;p++){
			fila = ((int)(p/w));
			col = p-(w*fila) + dif_fila_col;
			fila += dif_fila_col;
			ind_en_lapl = col + (w+n+1)*fila;
			double *n3 = get_neighborhood(laplacian, ind_en_lapl, 3, offsets);
			if ((n3[3]*n3[5]<0)&&(abs(n3[3]-n3[5])>(tzc*max_l))) {
				// horizontal sign change
				zero_cross[p] = 255;
			} else if ((n3[1]*n3[7]<0)&&(abs(n3[1]-n3[7])>(tzc*max_l))) {
					// vertical sign change
					zero_cross[p] = 255;
				} else if ((n3[2]*n3[6]<0)&&(abs(n3[2]-n3[6])>(tzc*max_l))) {
						// +45deg sign change
						zero_cross[p] = 255;
					} else if ((n3[0]*n3[8]<0)&&(abs(n3[0]-n3[8])>(tzc*max_l))) {
							// -45deg sign change
							zero_cross[p] = 255;
						}
			free_neighborhood(n3);
		}
		free_neighbors_offsets(offsets);
// Software Guide : EndCodeSnippet

		// Save output image
//  Software Guide : BeginLatex
//	Save output image (using \textit{iio}): \\
//  Software Guide : EndLatex
// Software Guide : BeginCodeSnippet
		iio_save_image_float_vec(argv[5], zero_cross, w, h, 1);
// Software Guide : EndCodeSnippet
		fprintf(stderr, "Output Image saved in %s:\t %dx%d image with %d channel(s).\n", argv[5], w, h, pixeldim);
	
		// Free memory
		free(zero_cross);
		free(im_orig);
		free(im);
		free(im_smoothed);
		free(laplacian);
		free_gaussian_kernel(kernel);

		fprintf(stderr, "marr-hildreth edge detector computation done.\n");

		// Execution time:
		double finish = (double)clock();
		double exectime = (finish - start)/CLOCKS_PER_SEC;
		fprintf(stderr, "execution time: %1.3f s.\n", exectime);		

		return 0;
	
	} // else (argc)
}

//  Software Guide : BeginLatex
//	\textit{Note: the main function in \texttt{test\_mh\_log.c} is essentially the same. The only difference is that we generate a LoG kernel (instead a Gaussian kernel) using the
//	\texttt{LoG\_kernel} function in \texttt{gaussian\_kernel.c}. Therefore we don't need the laplacian operator, so we only make one convolution.} \\
//  Software Guide : EndLatex
