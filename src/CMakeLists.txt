# Copyright (C) 2011-2012, Haldo Spontón <haldos@fing.edu.uy>
# Copyright (C) 2011-2012, Juan Cardelino <juanc@fing.edu.uy>

MESSAGE("libs inside: ${LIBS}")
SET(E_LIBS iio gaussian_kernel 2dconvolution edge_detectors)
ADD_LIBRARY(iio iio.c iio.h)
ADD_LIBRARY(gaussian_kernel SHARED gaussian_kernel.c gaussian_kernel.h)
ADD_LIBRARY(2dconv SHARED 2dconvolution.c 2dconvolution.h)
ADD_LIBRARY(edge_detectors SHARED edge_detectors.c edge_detectors.h)
