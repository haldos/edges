# Copyright (C) 2011-2013, Haldo Spontón <haldos@fing.edu.uy>
# Copyright (C) 2011-2013, Juan Cardelino <juanc@fing.edu.uy>

PROJECT(edge_detectors)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

SET(CMAKE_C_FLAGS "-std=c99")
SET(PROJECT_SOURCE_DIR src)
SET(PROJECT_BINARY_DIR bin)
SET( LIBS m)
 
FIND_PACKAGE(TIFF)
IF(TIFF_FOUND)
   include_directories(${TIFF_INCLUDE_DIR})
   SET( LIBS ${LIBS} ${TIFF_LIBRARIES})
ENDIF()

FIND_PACKAGE(JPEG)
if(JPEG_FOUND)
   include_directories(${JPEG_INCLUDE_DIR})
   SET( LIBS ${LIBS} ${JPEG_LIBRARIES})
ENDIF()

FIND_PACKAGE(PNG)
IF(PNG_FOUND)
   include_directories(${PNG_INCLUDE_DIR})
   SET( LIBS ${LIBS} ${PNG_LIBRARIES})
ENDIF()

MESSAGE("libraries: ${LIBS}")

ADD_SUBDIRECTORY(src)
MESSAGE("standard libraries: ${LIBS}")
MESSAGE("custom libraries: ${E_LIBS}")

ADD_EXECUTABLE(${PROJECT_BINARY_DIR}/edges_ipol ${PROJECT_SOURCE_DIR}/edges_ipol.c)
TARGET_LINK_LIBRARIES(${PROJECT_BINARY_DIR}/edges_ipol ${E_LIBS} ${LIBS})

INSTALL(TARGETS ${PROJECT_BINARY_DIR}/edges_ipol RUNTIME DESTINATION bin)
