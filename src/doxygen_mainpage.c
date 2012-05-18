/*! \mainpage Review of edge detection algorithms
 *
 * \section intro_sec About

These are the source files implementing some classical algorithms of edge detection.

The implemented algorithms are:

- First Derivative:
 -  Roberts
 -  Prewitt
 -  Sobel
- Second Derivative:
 -  Marr-Hildreth
 -  Haralick

The algorithms are implemented in C language and are compiled using cmake. 

 * \section install_sec Compilation

In the root directory:

> cmake .

> make

 * \section git_sec Git Repository

https://github.com/haldos/edges

 */

/*!
 \defgroup fded First Derivative Edge Detection Algorithms
 \defgroup haralick Haralick Algorithm
 \defgroup mhg Marr-Hildreth (Gaussian) Algorithm
 \defgroup mhl Marr-Hildreth (LoG) Algorithm
*/
