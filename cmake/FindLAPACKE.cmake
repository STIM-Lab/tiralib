# - Try to find LAPACKE
#
# Once done this will define
#  LAPACKE_FOUND - System has LAPACKE
#  LAPACKE_INCLUDE_DIRS - The LAPACKE include directories
#  LAPACKE_LIBRARIES - The libraries needed to use LAPACKE
#  LAPACKE_DEFINITIONS - Compiler switches required for using LAPACKE
#
# Usually, LAPACKE requires LAPACK and the BLAS.  This module does
# not enforce anything about that.

find_path(LAPACKE_INCLUDE_DIR
          NAMES lapacke.h
          PATHS $ENV{LAPACK_ROOT} $ENV{LAPACK_PATH} ${INCLUDE_INSTALL_DIR}
          PATHS ENV INCLUDE)

find_library(LAPACKE_LIBRARY liblapacke lapacke
             PATHS $ENV{LAPACK_ROOT} $ENV{LAPACK_PATH} ${LIB_INSTALL_DIR}
             PATHS ENV LIBRARY_PATH
             PATHS ENV LD_LIBRARY_PATH)

if(MSVC)
	find_library(LAPACK_LIBRARY liblapack lapack
             PATHS $ENV{LAPACK_ROOT} $ENV{LAPACK_PATH} ${LIB_INSTALL_DIR}
             PATHS ENV LIBRARY_PATH
             PATHS ENV LD_LIBRARY_PATH)

	find_library(BLAS_LIBRARY libblas blas
             PATHS $ENV{LAPACK_ROOT} $ENV{LAPACK_PATH} ${LIB_INSTALL_DIR}
             PATHS ENV LIBRARY_PATH
             PATHS ENV LD_LIBRARY_PATH)
	
else()
	find_library(LAPACK REQUIRED)
	find_library(BLAS REQUIRED)
endif()
set(LAPACKE_LIBRARIES ${LAPACKE_LIBRARY} ${LAPACK_LIBRARY} ${BLAS_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LAPACKE DEFAULT_MSG
                                  LAPACKE_INCLUDE_DIR 
                                  LAPACKE_LIBRARIES)
mark_as_advanced(LAPACKE_INCLUDE_DIR LAPACKE_LIBRARIES)
