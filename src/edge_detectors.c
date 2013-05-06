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
	This file implements the different edge detection algorithms
	presented in the article "A review of classic edge detectors".
	All functions implemented here receive only a pointer to the 
	input image (grayscale) and the different parameters of each 
	method. All them also return a pointer to the output image.
*/

#include "2dconvolution.c"		// Include for 2D convolutions
#include "gaussian_kernel.c"	// for Gaussian and LoG kernels (Marr-Hildreth)
#include <math.h>				// for mathematical operations

// Usefull macros
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define THRESHOLD(x, th) (((x) > (th)) ? (255) : (0))

// Roberts edge detector
float *edges_roberts(double *im, int w, int h, float threshold, int padding_method){
	
	/* Define operators */
	double roberts_1[9] = {-1, 0, 0, 0, 1, 0, 0, 0, 0};		// ROBERTS
	double roberts_2[9] = { 0,-1, 0, 1, 0, 0, 0, 0, 0};		// OPERATORS
	for ( int z=0; z<9; z++ ) {								// NORMALIZATION
		roberts_1[z] /= 2;
		roberts_2[z] /= 2;
	}
	/* 3x3 operators are used  (although it is unnecessary) in order to simplify the 
	convolution application, and also to unify the three detectors of the first family. */
	
	/* Convolution with operators */
	double *im_r1 = conv2d(im, w, h, roberts_1, 3, padding_method);
	double *im_r2 = conv2d(im, w, h, roberts_2, 3, padding_method);

	/* Allocate memory for output image */
	float *im_roberts = malloc(w*h*sizeof(float));
	if (im_roberts == NULL){
		fprintf(stderr, "Out of memory...\n");
		exit(EXIT_FAILURE);
	}
	
	/* Two images are obtained (one for each operator). Then the gradient magnitude 
	image is constructed using $M=\sqrt{g_x^2+g_y^2}$. Also the absolute maximum 
	value of the constructed images is computed */
	int i,j, fila, col;
	double max_r = 0;
	int imax = w*h;
	for (i=0;i<imax;i++){
		fila = (int)(i/w);
		col = i - w*fila + 1;
		fila += 1;
		j = col + (w+2)*fila;
		im_roberts[i] = sqrt(im_r1[j]*im_r1[j] + im_r2[j]*im_r2[j]);
		max_r = MAX(max_r,im_roberts[i]);
	}
	
	/* Thresholding */
	for (i=0;i<imax;i++){
		im_roberts[i] = THRESHOLD(im_roberts[i],threshold*max_r);
	}
	
	/* Free memory */
	free(im_r1);
	free(im_r2);
	
	/* Output image */
	return im_roberts;
}

// Prewitt edge detector
float *edges_prewitt(double *im, int w, int h, float threshold, int padding_method){

	/* Define operators */
	double prewitt_1[9] = {-1,-1,-1, 0, 0, 0, 1, 1, 1};		// PREWITT
	double prewitt_2[9] = {-1, 0, 1,-1, 0, 1,-1, 0, 1};		// OPERATORS
	for ( int z=0; z<9; z++ ) {								// NORMALIZATION
		prewitt_1[z] /= 6;
		prewitt_2[z] /= 6;
	}
	
	/* Convolution with operators */
	double *im_p1 = conv2d(im, w, h, prewitt_1, 3, padding_method);
	double *im_p2 = conv2d(im, w, h, prewitt_2, 3, padding_method);

	/* Allocate memory for output image */
	float *im_prewitt = malloc(w*h*sizeof(float));
	if (im_prewitt == NULL){
		fprintf(stderr, "Out of memory...\n");
		exit(EXIT_FAILURE);
	}
	
	/* Two images are obtained (one for each operator). Then the gradient magnitude 
	image is constructed using $M=\sqrt{g_x^2+g_y^2}$. Also the absolute maximum 
	value of the constructed images is computed */
	int i,j, fila, col;
	double max_p = 0;
	int imax = w*h;
	for ( i=0; i<imax; i++ ) {
		fila = (int)(i/w);
		col = i - w*fila + 1;
		fila += 1;
		j = col + (w+2)*fila;
		im_prewitt[i] = sqrt(im_p1[j]*im_p1[j] + im_p2[j]*im_p2[j]);
		max_p = MAX(max_p,im_prewitt[i]);
	}
	
	/* Thresholding */
	for ( i=0; i<imax; i++ ) {
		im_prewitt[i] = THRESHOLD(im_prewitt[i],threshold*max_p);
	}
	
	/* Free memory */
	free(im_p1);
	free(im_p2);
	
	/* Output image */
	return im_prewitt;
}

// Sobel edge detector
float *edges_sobel(double *im, int w, int h, float threshold, int padding_method){

	/* Define operators */
	double sobel_1[9] = {-1,-2,-1, 0, 0, 0, 1, 2, 1};		// SOBEL
	double sobel_2[9] = {-1, 0, 1,-2, 0, 2,-1, 0, 1};		// OPERATORS
	for ( int z=0; z<9; z++ ) {								// NORMALIZATION
		sobel_1[z] /= 8;
		sobel_2[z] /= 8;
	}
	
	/* Convolution with operators */
	double *im_s1 = conv2d(im, w, h, sobel_1, 3, padding_method);
	double *im_s2 = conv2d(im, w, h, sobel_2, 3, padding_method);

	/* Allocate memory for output image */
	float *im_sobel = malloc(w*h*sizeof(float));
	if (im_sobel == NULL){
		fprintf(stderr, "Out of memory...\n");
		exit(EXIT_FAILURE);
	}
	
	/* Two images are obtained (one for each operator). Then the gradient magnitude 
	image is constructed using $M=\sqrt{g_x^2+g_y^2}$. Also the absolute maximum 
	value of the constructed images is computed */
	int i,j, fila, col;
	double max_s = 0;
	int imax = w*h;
	for ( i=0; i<imax; i++ ) {
		fila = (int)(i/w);
		col = i - w*fila + 1;
		fila += 1;
		j = col + (w+2)*fila;
		im_sobel[i] = sqrt(im_s1[j]*im_s1[j] + im_s2[j]*im_s2[j]);
		max_s = MAX(max_s,im_sobel[i]);
	}
	
	/* Thresholding */
	for ( i=0; i<imax; i++ ) {
		im_sobel[i] = THRESHOLD(im_sobel[i],threshold*max_s);
	}
	
	/* Free memory */
	free(im_s1);
	free(im_s2);
	
	/* Output image */
	return im_sobel;
}

// Haralick edge detector
float *edges_haralick(double *im, int w, int h, float rhozero, int padding_method){

	int i,j;	// loop counters
	
	/* Haralick's masks for computing k1 to k10
	New masks calculated by 2-d fitting using LS, with the function
	f(x,y) = k1 + k2*x + k3*y + k4*x² + k5*xy + k6*y² + k7*x³ + k8*x²y + k9*xy² + k10*y³. */
	double masks[10][25] = { {       425,   275,  225,  275,  425,   
									 275,   125,   75,  125,  275,   
                                     225,    75,   25,   75,  225,   
									 275,   125,   75,  125,  275,   
									 425,   275,  225,  275,  425},
								 { -2260,  -620,    0,  620, 2260, 
								   -1660,  -320,    0,  320, 1660, 
								   -1460,  -220,    0,  220, 1460,
                                   -1660,  -320,    0,  320, 1660, 
                                   -2260,  -620,    0,  620, 2260},
								 {  2260,  1660, 1460, 1660, 2260,   
								     620,   320,  220,  320,  620,
								       0,     0,    0,    0,    0,  
								    -620,  -320, -220, -320, -620, 
								   -2260, -1660,-1460,-1660,-2260},
								 {  1130,   620,  450,  620, 1130,   
								     830,   320,  150,  320,  830,   
								     730,   220,   50,  220,  730,   										 
								     830,   320,  150,  320,  830,  
								    1130,   620,  450,  620, 1130},
								 {  -400,  -200,    0,  200,  400,  
								    -200,  -100,    0,  100,  200,     
								       0,     0,    0,    0,    0,   
								     200,   100,    0, -100, -200,   
								     400,   200,    0, -200, -400},
								 {  1130,   830,  730,  830, 1130,   
								     620,   320,  220,  320,  620,   
								     450,   150,   50,  150,  450,   
								     620,   320,  220,  320,  620,  
								    1130,   830,  730,  830, 1130},
								 { -8260, -2180,    0, 2180, 8260, 
								   -6220, -1160,    0, 1160, 6220, 
								   -5540,  -820,    0,  820, 5540, 
								   -6220, -1160,    0, 1160, 6220, 
								   -8260, -2180,    0, 2180, 8260},
								 {  5640,  3600, 2920, 3600, 5640,  
								    1800,   780,  440,  780, 1800,     
								       0,     0,    0,    0,    0, 
								   -1800,  -780, -440, -780,-1800, 
								   -5640, -3600,-2920,-3600,-5640},
								 { -5640, -1800,    0, 1800, 5640, 
								   -3600,  -780,    0,  780, 3600, 
								   -2920,  -440,    0,  440, 2920, 
								   -3600,  -780,    0,  780, 3600, 
								   -5640, -1800,    0, 1800, 5640},
								 {  8260,  6220, 5540, 6220, 8260,  
								    2180,  1160,  820, 1160, 2180,     
								       0,     0,    0,    0,    0, 
								   -2180, -1160, -820,-1160,-2180, 
								   -8260, -6220, 5540,-6220,-8260   } };
	
	/* Initialise edges image */
	float *edges = calloc(w*h,sizeof(float));
	if (edges == NULL){
		fprintf(stderr, "Out of memory...\n");
		exit(EXIT_FAILURE);
	}
	
	/* Padding (zero-padding or reflection of image boundary) */
	int wx = (w+8);
	int hx = (h+8);
	double *aux = calloc(wx*hx,sizeof(double));
	if (aux == NULL){
		fprintf(stderr, "Out of memory...\n");
		exit(EXIT_FAILURE);
	}
	int fila,col;
	int imax = wx*hx;
	/* Zero-padding */
	if (padding_method == 0) {
		for(i=0;i<imax;i++){
			fila = (int)(i/wx);
			col = i-(wx*fila);	
			if ( (fila>=4)&&(col>=4)&&(fila<h+4)&&(col<w+4) ) {
				aux[i] = im[(col-4)+(w*(fila-4))];
			}
		}
	}
	/* Reflection of image boundary */
	if (padding_method == 1) {
		int fila_refl, col_refl;
		for(i=0;i<imax;i++){
			fila = (int)(i/wx);
			col = i-(wx*fila);
			if (fila<4) {
				fila_refl = 7 - fila;
				if (col<4) { 				//zone1
					col_refl = 7 - col;
				} else if (col<w+4) {		//zone2
					col_refl = col;
				} else { 					//zone3
					col_refl = 2*w + 7 - col;
				}
			} else if (fila<h+4) {
				fila_refl = fila;
				if (col<4) { 				//zone4
					col_refl = 7 - col;
				} else if (col<w+4) { 		//image
					col_refl = col;
				} else { 					//zone5
					col_refl =  2*w + 7 - col;
				}
			} else {
				fila_refl = 2*h + 7 - fila;
				if (col<4) { 				//zone6
					col_refl =	7 - col;
				} else if (col<w+4) {		//zone7
					col_refl = col;
				} else { 					//zone8
					col_refl =  2*w + 7 - col;
				}
			}
			aux[i] = im[(col_refl-4)+(w*(fila_refl-4))];
		}
	}
	
	/* Haralick algorithm 
	Haralick algorithm: coefficients k_1 to k_10 are computed in every pixel of the original image 
	(using the function get_neighbors_offset to get the index offsets of the neighbor pixels and 
	the function get_neighborhood to get the neighborhood of a pixel using those index offsets). 
	Once the coefficients are calculated, are computed
	C_2 = (k_2^2k_4 + k_2k_3k_5 + k_3^2k_6) / (k_2^2 + k_3^2)
	and
	C_3 = (k_2^3k_7 + k_2^2k_3k_8 + k_2k_3^2k_9 + k_3^3k_10) / (sqrt(k_2^2 + k_3^2)^3),
	and then the edge condition is evaluated in every pixel; $|C_2 / 2C_3| <= rhozero. */
	int i_zp, u, v, num_edges;
	num_edges = 0;
	double k[10];
	int *offsets = get_neighbors_offset(wx, 5);
	double acum;
	double C2, C3, denom, sintheta, costheta;
	for(fila=0;fila<h;fila++){
		for(col=0;col<w;col++){
			i = col + w*fila;				// original image & edges image index
			i_zp = (col+4) + wx*(fila+4);	// padded image index
			double *neighborhood = get_neighborhood(aux, i_zp, 5, offsets);
			// k1 to k10 (note: k1 (u=0) is not necessary)
			for(u=0;u<10;u++){
				acum = 0;
				for(v=0;v<25;v++){
					acum += neighborhood[v]*masks[u][v];
				}
				k[u] = acum;
			}
			// compute C2 and C3
			denom = sqrt( k[1]*k[1] + k[2]*k[2] );
			sintheta = - k[1] / denom;
			costheta = - k[2] / denom;
			C2 = k[3]*sintheta*sintheta + k[4]*sintheta*costheta + k[5]*costheta*costheta;
			C3 = k[6]*sintheta*sintheta*sintheta + k[7]*sintheta*sintheta*costheta +
				 k[8]*sintheta*costheta*costheta + k[9]*costheta*costheta*costheta;
			//if ((fabs(C2 / (3*C3))<=rhozero)&&(C3<=0)) {
			if ((fabs(C2 / (3*C3))<=rhozero)) {
				edges[i] = 255;
				num_edges += 1;
			}
			// free neighborhood
			free_neighborhood(neighborhood);
		}
	}
	
	/* Free memory */
	free(aux);
	free(edges);
		
	return edges;
}

// Marr-Hildreth edge detector with Gaussian and Laplacian kernels
float *edges_mh(double *im, int w, int h, float sigma, int n, float tzc, int padding_method){

	/* Generate Gaussian kernel */
	double *kernel = gaussian_kernel(n,sigma);
	
	/* Smooth input image with the Gaussian kernel */
	double *im_smoothed = conv2d(im, w, h, kernel, n, padding_method);
	
	/* Computation of the Laplacian of the smoothed image
	A 3x3 approximation of the laplacian operator is used */
	double operator[9] = {1, 1, 1, 1, -8, 1, 1, 1, 1};
	double *laplacian = conv2d(im_smoothed, w+n-1, h+n-1, operator, 3, padding_method);
	
	/* Calculate max absolute value of laplacian (required 
	for thresholding in zero-crossing (next) */
	double max_l = 0;
	int p;
	int pmax = (w+n+1)*(h+n+1);
	for (p=0;p<pmax;p++){
		if (abs(laplacian[p])>max_l){
			max_l = abs(laplacian[p]);
		}
	}
	
	/* Zero-crossing */
	float *zero_cross = calloc(w*h,sizeof(float));		
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
	
	/* Free memory */
	free(im_smoothed);
	free(laplacian);
	free_gaussian_kernel(kernel);
	
	/* Output image */
	return zero_cross;
}

// Marr-Hildreth edge detector with Laplacian of Gaussian kernel
float *edges_mh_log(double *im, int w, int h, float sigma, int n, float tzc, int padding_method){

	/* Generate LoG (Laplacian of Gaussian) kernel */
	double *kernel = LoG_kernel(n,sigma);
	
	/* Smooth input image with LoG kernel */
	double *im_smoothed = conv2d(im, w, h, kernel, n, padding_method);
	
	/* Calculate max absolute value of smoothed image 
	(required for thresholding in zero-crossing (next) */
		double max_l = 0;
	int p;
	int pmax = (w+n-1)*(h+n-1);
	for ( p=0; p<pmax; p++ ) {
		if ( abs(im_smoothed[p]) > max_l ){
			max_l = abs(im_smoothed[p]);
		}
	}
	
	/* Zero-crossing */
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
	
	/* Free memory */
	free(im_smoothed);
	free_gaussian_kernel(kernel);

	/* Output image */
	return zero_cross;
}
