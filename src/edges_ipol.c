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

/*
	This file tests all the edge detection algorithms implemented in
	edge_detectors.c. It also leads with input/output and with the
	parameters handling.
*/

#include "edge_detectors.c"		// for edge detection algorithms
#include "iio.c"				// for image input/output
#include <time.h>				// for execution time measurement
#include <string.h>				// for string comparison

/* Display help function */
void displayhelp(void) {
	printf("edges_ipol: Main program to run the edge detection algorithms\ndescribed in the IPOL article called \"A REVIEW OF CLASSIC EDGE\nDETECTION ALGORITHMS\". The functions called are from\nedge_detectors.c.\n\nAuthor: Haldo Spontón <haldos@fing.edu.uy>\n\nUsage: edges_ipol [options] input.png.\n\nOptions:\n\t-r : Roberts edge detection algorithm.\n\t\tAfter -r you must specify a threshold\n\t\t(float) between 0 and 1.\n\t-p : Prewitt edge detection algorithm.\n\t\tAfter -p you must specify a threshold\n\t\t(float) between 0 and 1.\n\t-s : Sobel edge detection algorithm.\n\t\tAfter -s you must specify a threshold\n\t\t(float) between 0 and 1.\n\t-h : Haralick edge detection algorithm.\n\t\tAfter -h you must specify a value for\n\t\tthe Haralick condition threshold (float).\n\t-m : Marr-Hildreth (Gaussian) edge detection algorithm.\n\t\tAfter -m you must specify a value for\n\t\tthe standard deviation of the Gaussian\n\t\tfunction (float), the size of the Gaussian\n\t\tkernel (integer) ant the threshold for the\n\t\tzero-crossing detection (float).\n\t-l : Marr-Hildreth (LoG) edge detection algorithm.\n\t\tAfter -l you must specify a value for\n\t\tthe standard deviation of the LoG (Laplacian\n\t\tof Gaussian) function (float), the size of\n\t\tthe Gaussian kernel (integer) ant the threshold\n\t\tfor the zero-crossing detection (float).\n\t-H : Display help.\n\nInput image:\n\tMust be a (m x n x 1) grayscale png-format image.\n\nOutput images:\n\tThe program will export one image for each edge\n\tdetection algorithm selected. Their filenames are:\n\t\tim_roberts.png\n\t\tim_prewitt.png\n\t\tim_sobel.png\n\t\tim_haralick.png\n\t\tim_mh.png\n\t\tim_mhl.png\n");
}

/* Main function */
int main(int argc, char *argv[]) {	

	/* Start */
	printf("A review of classic edge detection algorithms\n");
	printf("IPOL 2013 - Haldo Spontón & Juan Cardelino\n");
	
	/* Options and parameters handling */
	int n;						// loop counter
	int argc_sobel = 0;			// 
	int argc_prewitt = 0;		// 
	int argc_roberts = 0;		// indexes
	int argc_haralick = 0;		// in argv[]
	int argc_mh = 0;			//
	int argc_mhl = 0;			//
	int padding_method = 1;		// Reflection of image boundary. Hard-coded. Can be changed to 0 (zero-padding)
	
	for ( n = 1; n < argc; n++ ) {			// scan through args
		switch ( (int)argv[n][0] ) {		// check for option character
			case '-':	switch ( (int)argv[n][1] ) {
							case 'r':	argc_roberts = n;
										break;
							case 'p':	argc_prewitt = n;
										break;
							case 's':	argc_sobel = n;
										break;
							case 'h':	argc_haralick = n;
										break;
							case 'm':	argc_mh = n;
										break;
							case 'l':	argc_mhl = n;
										break;
							case 'H':	/* Display help! */
										displayhelp();
										exit(1);
										break;
							default:	printf("Error: Invalid option -> %s. Valid options are -r, -p, -s, -h, -m and -l (-H for help).\n", argv[n]);
										exit(1);
										break;
						}
						/* Error message if an option is the last input parameter
						   or if the next parameter after an option is another option */
						if ( n==argc-1 ) {
							printf("Error: Missing parameter(s) for %s option.\n",argv[n]);
							exit(1);
							break;
						} else if ( argv[n+1][0]=='-' ) {
							printf("Error: Missing parameter(s) for %s option.\n",argv[n]);
							exit(1);
							break;
						}
						/* In the case of -m and -l, three parameters are needed */
						if ( (argv[n][1]=='m') || (argv[n][1]=='l') ) { 
							if ( n+3>argc-1 ) {
								printf("Error: Missing parameter(s) for %s option.\n",argv[n]);
								exit(1);
								break;
							} else if ( (argv[n+2][0]=='-') || (argv[n+3][0]=='-') ) {
								printf("Error: Missing parameter(s) for %s option.\n",argv[n]);
								exit(1);
								break;
							}
						}
			default: 	break;
		}
	}
	
	/* Check for total number of parameters */
	int nparam = 2;
	if (argc_roberts!=0) nparam+=2;
	if (argc_prewitt!=0) nparam+=2;
	if (argc_sobel!=0) nparam+=2;
	if (argc_haralick!=0) nparam+=2;
	if (argc_mh!=0) nparam+=4;
	if (argc_mhl!=0) nparam+=4;
	if (nparam!=argc) {
		printf("Error: Wrong number of arguments (%d instead of %d).\nUsage: %s [options] input.png\nFor help type %s -H.\n",argc-1,nparam-1,argv[0],argv[0]);
		exit(1);
	}
	
	/* Read parameters once checked */
	float	th_roberts, 	// threshold in Roberts algorithm
			th_prewitt, 	// threshold in Prewitt algorithm
			th_sobel, 		// threshold in Sobel algorithm
			rhozero,		// threshold for the Haralick edge condition
			sigma_m, 		// standard deviation of Gaussian kernel in Marr-Hildreth algorithm
			tzc_m, 			// threshold for zero-crossing in Marr-Hildreth algorithm
			sigma_l, 		// standard deviation of LoG kernel in Marr-Hildreth-LoG algorithm
			tzc_l;			// threshold for zero-crossing in Marr-Hildreth-LoG algorithm
	int		n_m, 			// kernel size in Marr-Hildreth algorithm
			n_l;			// kernel size in Marr-Hildreth-LoG algorithm
	if (argc_roberts!=0) th_roberts = atof(argv[argc_roberts+1]);
	if (argc_roberts!=0) th_prewitt = atof(argv[argc_prewitt+1]);
	if (argc_roberts!=0) th_sobel = atof(argv[argc_sobel+1]);
	if (argc_haralick!=0) rhozero = atof(argv[argc_haralick+1]);
	if (argc_mh!=0) {
		sigma_m = atof(argv[argc_mh+1]);
		n_m = atoi(argv[argc_mh+2]);
		tzc_m = atof(argv[argc_mh+3]);
	}
	if (argc_mhl!=0) {
		sigma_l = atof(argv[argc_mhl+1]);
		n_l = atoi(argv[argc_mhl+2]);
		tzc_l = atof(argv[argc_mhl+3]);
	}
	
	/* Display parameters */
	printf("\nPARAMETERS:\n");
	if (argc_roberts!=0) printf("\n\t--> Roberts selected.\n\tThreshold = %.2f.\n",th_roberts); 
		else printf("\n\t--> Roberts not selected.\n");
	if (argc_prewitt!=0) printf("\n\t--> Prewitt selected.\n\tThreshold = %.2f.\n",th_prewitt); 
		else printf("\n\t--> Prewitt not selected.\n");
	if (argc_sobel!=0) printf("\n\t--> Sobel selected.\n\tThreshold = %.2f.\n",th_sobel); 
		else printf("\n\t--> Sobel not selected.\n");
	if (argc_haralick!=0) printf("\n\t--> Haralick selected.\n\tRhozero = %.2f.\n",rhozero); 
		else printf("\n\t--> Haralick not selected.\n");
	if (argc_mh!=0) printf("\n\t--> Marr-Hildreth (Gaussian) selected.\n\tSigma = %.2f.\n\tN = %d.\n\tTZC = %.2f.\n",sigma_m,n_m,tzc_m); 
		else printf("\n\t--> Marr-Hildreth (Gaussian) not selected.\n");
	if (argc_mhl!=0) printf("\n\t--> Marr-Hildreth (LoG) selected.\n\tSigma = %.2f.\n\tN = %d.\n\tTZC = %.2f.\n",sigma_l,n_l,tzc_l); 
		else printf("\n\t--> Marr-Hildreth (LoG) not selected.\n");
	
	/* Load input image */
	int w, h, pixeldim;
	float *im = iio_read_image_float_vec(argv[argc-1], &w, &h, &pixeldim);
	double *imdouble = malloc(w*h*sizeof(double));
	printf("\nINPUT IMAGE:\n\n\t%s\n\tDimensions: %d x %d pixeles. Channels: %d.\n",argv[argc-1],w,h,pixeldim);
	for ( int i=0; i<w*h; i++ ) {
		imdouble[i] = (double)im[i];
	}
	
	/* Run the selected edge detection algorithms */
	printf("\nRUNNING SELECTED EDGE DETECTION ALGORITHMS...\n");
	
	/* Roberts edge detection algorithm */
	if (argc_roberts!=0) {
		printf("\n\tRunning Roberts edge detection algorithm...\n");
		double start = (double)clock();
		/* Processing here */
			float *im_roberts = edges_roberts(imdouble, w, h, th_roberts, padding_method);
			iio_save_image_float_vec("im_roberts.png", im_roberts, w, h, 1);
		/* End of processing */
		double finish = (double)clock();
		double exectime = (finish - start)/CLOCKS_PER_SEC;
		printf("\tDone! Execution time: %1.3f s.\n", exectime);
		printf("\tOutput image saved as im_roberts.png.\n");
	}
	
	/* Prewitt edge detection algorithm */
	if (argc_prewitt!=0) {
		printf("\n\tRunning Prewitt edge detection algorithm...\n");
		double start = (double)clock();
		/* Processing here */
			float *im_prewitt = edges_prewitt(imdouble, w, h, th_prewitt, padding_method);
			iio_save_image_float_vec("im_prewitt.png", im_prewitt, w, h, 1);
		/* End of processing */
		double finish = (double)clock();
		double exectime = (finish - start)/CLOCKS_PER_SEC;
		printf("\tDone! Execution time: %1.3f s.\n", exectime);
		printf("\tOutput image saved as im_prewitt.png.\n");
	}
	
	/* Sobel edge detection algorithm */
	if (argc_sobel!=0) {
		printf("\n\tRunning Sobel edge detection algorithm...\n");
		double start = (double)clock();
		/* Processing here */
			float *im_sobel = edges_sobel(imdouble, w, h, th_sobel, padding_method);
			iio_save_image_float_vec("im_sobel.png", im_sobel, w, h, 1);
		/* End of processing */
		double finish = (double)clock();
		double exectime = (finish - start)/CLOCKS_PER_SEC;
		printf("\tDone! Execution time: %1.3f s.\n", exectime);
		printf("\tOutput image saved as im_sobel.png.\n");
	}
	
	/* Haralick edge detection algorithm */
	if (argc_haralick!=0) {
		printf("\n\tRunning Haralick edge detection algorithm...\n");
		double start = (double)clock();
		/* Processing here */
			float *im_haralick = edges_haralick(imdouble, w, h, rhozero, padding_method);
			iio_save_image_float_vec("im_haralick.png", im_haralick, w, h, 1);
		/* End of processing */
		double finish = (double)clock();
		double exectime = (finish - start)/CLOCKS_PER_SEC;
		printf("\tDone! Execution time: %1.3f s.\n", exectime);
		printf("\tOutput image saved as im_haralick.png.\n");
	}
	
	/* Marr-Hildreth (Gaussian) edge detection algorithm */
	if (argc_mh!=0) {
		printf("\n\tRunning Marr-Hildreth (Gaussian) edge detection algorithm...\n");
		double start = (double)clock();
		/* Processing here */
			float *im_mh = edges_mh(imdouble, w, h, sigma_m, n_m, tzc_m, padding_method);
			iio_save_image_float_vec("im_mh.png", im_mh, w, h, 1);
		/* End of processing */
		double finish = (double)clock();
		double exectime = (finish - start)/CLOCKS_PER_SEC;
		printf("\tDone! Execution time: %1.3f s.\n", exectime);
		printf("\tOutput image saved as im_mh.png.\n");
	}
	
	/* Marr-Hildreth (LoG) edge detection algorithm */
	if (argc_mhl!=0) {
		printf("\n\tRunning Marr-Hildreth (LoG) edge detection algorithm...\n");
		double start = (double)clock();
		/* Processing here */
			float *im_mhl = edges_mh_log(imdouble, w, h, sigma_l, n_l, tzc_l, padding_method);
			iio_save_image_float_vec("im_mhl.png", im_mhl, w, h, 1);
		/* End of processing */
		double finish = (double)clock();
		double exectime = (finish - start)/CLOCKS_PER_SEC;
		printf("\tDone! Execution time: %1.3f s.\n", exectime);
		printf("\tOutput image saved as im_mhl.png.\n");
	}
	
	/* Memory free */
	free(im);
	free(imdouble);
}
