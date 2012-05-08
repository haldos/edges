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
* \file test_mh_log.c
* \author Haldo Spontón <haldos@fing.edu.uy> & Juan Cardelino <juanc@fing.edu.uy>
* \date May, 2012
* \see ``Review of edge detectors´´ IPOL publication.
* \brief Implements the Marr-Hildreth edge detection algorithm, using a LoG (Laplacian of a Gaussian) kernel. In summary, this function opens the image, converts it to grayscale (if necessary), then the convolves it with a LoG kernel, and finally, find zero crossings in the resulting image.
*/

#include "iio.c"
#include "gaussian_kernel.c"
#include "2dconvolution.c"
#include <time.h>

// Values for saving debug images.
#define SAVE_KERNEL false
#define SAVE_SMOOTHED_IMAGE false

/*!
 \fn int main_mhl(int argc, char *argv[])
 \brief Main function of the Marr-Hildreth edge detection algorithm, using LoG kernel.
 @param input_image Filename of the input image.
 @param sigma Standard deviation of the Laplacian of a Gaussian kernel.
 @param n Size of the kernel (n x n).
 @param tzc Threshold for the zero-crossing algorithm.
 @param output_image Filename of the output image.
 \return None.
 \ingroup mhl
	\note The real name of this function is main. Function name was temporarily changed to the proper functioning of doxygen.
*/

int main(int argc, char *argv[]) {
	if (argc != 6) {
		printf("Usage: %s input_image_1 sigma n tzc output\n", argv[0]);
	} else {

		// Execution time:
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
		int w, h, pixeldim;
		float *im_orig = iio_read_image_float_vec(argv[1], &w, &h, &pixeldim);
		fprintf(stderr, "Input image loaded:\t %dx%d image with %d channel(s).\n", w, h, pixeldim);

		// Grayscale conversion (if necessary)
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

		// Generate gaussian kernel
		// double *kernel = gaussian_kernel(n,sigma);
		// modificacion 2: cambio el kernel gaussiano por el LoG kernel.
		// con esto aplico el kernel gaussiano y el laplaciano en una sola convolucion.
		double *kernel = LoG_kernel(n,sigma);

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
			iio_save_image_float_vec("kernel_log.png", kernel_float, n, n, 1);
			free(kernel_float);
			fprintf(stderr, "kernel saved to kernel_log.png\n");
		}
		// end of save kernel image.
		
		// Smooth input image with gaussian kernel
		double *im_smoothed = conv2d(im, w, h, kernel, n);

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
			iio_save_image_float_vec("smoothed_log.png", smoothed, w, h, 1);
			free(smoothed);
			fprintf(stderr, "smoothed image saved to smoothed_log.png\n");
		}
		// end of save smoothed image.

		double max_l = 0;
		int p;
		int pmax = (w+n-1)*(h+n-1);
		for (p=0;p<pmax;p++){
			if (abs(im_smoothed[p])>max_l){
				max_l = abs(im_smoothed[p]);
			}
		}

		// Zero-crossing
		float *zero_cross = calloc(w*h,sizeof(float));
		if (zero_cross == NULL){
			fprintf(stderr, "Out of memory...\n");
			exit(EXIT_FAILURE);
		}
		int ind_en_lapl, fila, col;
		int *offsets = get_neighbors_offset(w+n-1, 3);
		pmax = w*h;
		int dif_fila_col = (n-1)/2;
		for (p=0;p<pmax;p++){
			fila = ((int)(p/w));
			col = p-(w*fila) + dif_fila_col;
			fila += dif_fila_col;
			ind_en_lapl = col + (w+n-1)*fila;
			double *n3 = get_neighborhood(im_smoothed, ind_en_lapl, 3, offsets);
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

		// Save output image
		iio_save_image_float_vec(argv[5], zero_cross, w, h, 1);
		fprintf(stderr, "Output Image saved in %s:\t %dx%d image with %d channel(s).\n", argv[5], w, h, pixeldim);
	
		// Free memory
		free(zero_cross);
		free(im_orig);
		free(im);
		free(im_smoothed);
		free_gaussian_kernel(kernel);

		fprintf(stderr, "marr-hildreth edge detector computation done.\n");

		// Execution time:
		double finish = (double)clock();
		double exectime = (finish - start)/CLOCKS_PER_SEC;
		fprintf(stderr, "execution time: %1.3f s.\n", exectime);		

		return 0;
	
	} // else (argc)
}
