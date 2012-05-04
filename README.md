REVIEW OF EDGE DETECTORS - IPOL 2012
====================================

About
-----

These are the source files implementing some classical algorithms of edge detection. 
	
The implemented algorithms are:
		
* First Derivative:
 * Roberts
 * Prewitt
 * Sobel
* Second Derivative:	
 * Marr-Hildreth
 * Haralick

The algorithms are implemented in C language and are compiled using cmake. 

2 - Filelist

	.:
	bin  			-	Binary files folder (empty until compilation).
	CMakeLists.txt  -	CMake configuration file.
	README  		- 	This read me file.
	src  			-	Source files folder.
	test			-	Test images folder.
	
	./bin: (empty until compilation)

	./example:
	lena_harlick.png	-	Output of Haralick algorithm using lena.png as input image.
	lena_mh.png			-	Output of Marr-Hildreth algorithm (Gaussian kernel) using lena.png as input image.
	lena_mhlog.png		-	Output of Marr-Hildreth algorithm (LoG kernel) using lena.png as input image.
	lena_prewitt.png	-	Output of Prewitt algorithm using lena.png as input image.
	lena_roberts.png	-	Output of Roberts algorithm using lena.png as input image.
	lena_sobel.png		-	Output of Sobel algorithm using lena.png as input image.
	
	./src:
	2dconvolution.c  	-	Implementation of the 2-D convolution.
	2dconvolution.h		-	(header file).
	CMakeLists.txt		-	CMake configuration file.
	gaussian_kernel.c	-	Implementation of the Gaussian and LoG kernel generation.
	gaussian_kernel.h	-	(header file).
	iio.c				-	IIO library C-file.
	iio.h     			-	IIO library (header file).
	test_fded.c      	-	First derivative edge detectors, main C file.
	test_haralick.c		-	Haralick algorithm, main C file.
	test_mh.c   		-	Marr-Hildreth algorithm (Gaussian kernel), main C-file.
	test_mh_log.c 		-	Marr-Hildreth algorithm (LoG kernel), main C-file.
	
	./test:
	lena.png  		-	Lena image (512x512px).
	molino.png  	-	Windmill image (1000x563px).
	oranges.png		-	Oranges image (536x480px).

3 - Compilation

	In the root directory:
		> cmake .
		> make

	Four executable files are created in the bin directory:

	bin/:
	test_fded  		-	First derivative edge detectors (Roberts, Prewitt & Sobel).
	test_haralick  	-	Haralick algorithm.
	test_mh  		-	Marr-Hildreth algorithm, using Gaussian kernel.
	test_mh_log		-	Marr-Hildreth algorithm, using LoG kernel.


4 - Execution

	First derivative edge detectors:
		Usage: test_fded input_image threshold
		Example:
					> bin/test_fded test/lena.png 0.1
					Input image loaded:	 512x512 image with 3 channel(s).
					images converted to grayscale
					Roberts image saved to roberts.png.
					Prewitt image saved to prewitt.png.
					Sobel image saved to sobel.png.
					Edge detection algorithms based on first derivative computation done.
					execution time: 0.280 s.

	Haralick algorithm:
		Usage: test_haralick input_image rhozero output
		Example:
					> bin/test_haralick test/lena.png 0.5 lena_harlick.png
					Input image loaded:	 512x512 image with 3 channel(s).
					images converted to grayscale
					60122 edge points found...
					Output Image saved in lena_harlick.png:	 512x512 image with 3 channel(s).
					haralick's edge detector computation done.
					execution time: 0.440 s.

	Marr-Hildreth algorithm (Gaussian Kernel):
		Usage: test_mh input_image sigma n tzc output_image
		Example:
					> bin/test_mh test/lena.png 2 13 0.1 lena_mh.png
					Input image loaded:	 512x512 image with 3 channel(s).
					images converted to grayscale
					Output Image saved in lena_mh.png:	 512x512 image with 3 channel(s).
					marr-hildreth edge detector computation done.
					execution time: 0.420 s.

	Marr-Hildreth algorithm (LoG kernel):
		Usage: test_mh_log input_image_1 sigma n tzc output
		Example:
					> bin/test_mh_log test/lena.png 2 17 0.1 lena_mhlog.png
					Input image loaded:	 512x512 image with 3 channel(s).
					images converted to grayscale
					Output Image saved in lena_mhlog.png:	 512x512 image with 3 channel(s).
					marr-hildreth edge detector computation done.
					execution time: 0.670 s.

5 - Test images

	Some test images are provided in the «test» folder.
	The output images using «lena.png» as input image, are also provided in the «example» directory.

6 - Copyrights

	Copyright (C) 2011-2012, Haldo Spontón <haldos@fing.edu.uy>
	Copyright (C) 2011-2012, Juan Cardelino <juanc@fing.edu.uy>

7 - Licence

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

