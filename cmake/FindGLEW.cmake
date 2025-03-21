# Copyright (c) 2012-2016 DreamWorks Animation LLC
#
# All rights reserved. This software is distributed under the
# Mozilla Public License 2.0 ( http://www.mozilla.org/MPL/2.0/ )
#
# Redistributions of source code must retain the above copyright
# and license notice and the following restrictions and disclaimer.
#
# *     Neither the name of DreamWorks Animation nor the names of
# its contributors may be used to endorse or promote products derived
# from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# IN NO EVENT SHALL THE COPYRIGHT HOLDERS' AND CONTRIBUTORS' AGGREGATE
# LIABILITY FOR ALL CLAIMS REGARDLESS OF THEIR BASIS EXCEED US$250.00.
#

#-*-cmake-*-
# - Find GLEW
#
# Author : Nicholas Yue yue.nicholas@gmail.com
#
# This auxiliary CMake file helps in find the GLEW headers and libraries
#
# GLEW_FOUND            set if Glew is found.
# GLEW_INCLUDE_DIR      GLEW's include directory
# GLEW_glew_LIBRARY        GLEW libraries
# GLEW_glewmx_LIBRARY      GLEWmx libraries (Mulitple Rendering Context)

FIND_PACKAGE ( PackageHandleStandardArgs )

FIND_PATH( GLEW_LOCATION include/GL/glew.h
  "$ENV{GLEW_ROOT}"
  NO_DEFAULT_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  )

FIND_PACKAGE_HANDLE_STANDARD_ARGS ( GLEW
  REQUIRED_VARS GLEW_LOCATION
  )

IF ( GLEW_LOCATION )

  SET( GLEW_INCLUDE_DIR "${GLEW_LOCATION}/include" CACHE STRING "GLEW include path")

  SET ( ORIGINAL_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})
  IF (GLEW_USE_STATIC_LIBS)
	IF (APPLE)
      SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
      FIND_LIBRARY ( GLEW_LIBRARY_PATH GLEW PATHS ${GLEW_LOCATION}/lib
		NO_DEFAULT_PATH
		NO_SYSTEM_ENVIRONMENT_PATH
		)
      FIND_LIBRARY ( GLEWmx_LIBRARY_PATH GLEWmx PATHS ${GLEW_LOCATION}/lib
		NO_DEFAULT_PATH
		NO_SYSTEM_ENVIRONMENT_PATH
		)
	  # MESSAGE ( "APPLE STATIC" )
	  # MESSAGE ( "GLEW_LIBRARY_PATH = " ${GLEW_LIBRARY_PATH} )
	ELSEIF (WIN32)
      # Link library
      SET(CMAKE_FIND_LIBRARY_SUFFIXES ".lib")
      FIND_LIBRARY ( GLEW_LIBRARY_PATH GLEW32S PATHS ${GLEW_LOCATION}/lib )
      FIND_LIBRARY ( GLEWmx_LIBRARY_PATH GLEW32MXS PATHS ${GLEW_LOCATION}/lib )
	ELSE (APPLE)
      SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
      FIND_LIBRARY ( GLEW_LIBRARY_PATH GLEW PATHS ${GLEW_LOCATION}/lib
		NO_DEFAULT_PATH
		NO_SYSTEM_ENVIRONMENT_PATH
		)
      FIND_LIBRARY ( GLEWmx_LIBRARY_PATH GLEWmx PATHS ${GLEW_LOCATION}/lib
		NO_DEFAULT_PATH
		NO_SYSTEM_ENVIRONMENT_PATH
		)
	  # MESSAGE ( "LINUX STATIC" )
	  # MESSAGE ( "GLEW_LIBRARY_PATH = " ${GLEW_LIBRARY_PATH} )
	ENDIF (APPLE)
  ELSE ()
	IF (APPLE)
      SET(CMAKE_FIND_LIBRARY_SUFFIXES ".dylib")
      FIND_LIBRARY ( GLEW_LIBRARY_PATH GLEW PATHS ${GLEW_LOCATION}/lib )
      FIND_LIBRARY ( GLEWmx_LIBRARY_PATH GLEWmx PATHS ${GLEW_LOCATION}/lib )
	ELSEIF (WIN32)
      # Link library
      SET(CMAKE_FIND_LIBRARY_SUFFIXES ".lib")
      FIND_LIBRARY ( GLEW_LIBRARY_PATH GLEW32 PATHS ${GLEW_LOCATION}/lib )
      FIND_LIBRARY ( GLEWmx_LIBRARY_PATH GLEW32mx PATHS ${GLEW_LOCATION}/lib )
      # Load library
      SET(CMAKE_FIND_LIBRARY_SUFFIXES ".dll")
      FIND_LIBRARY ( GLEW_DLL_PATH GLEW32 PATHS ${GLEW_LOCATION}/bin
		NO_DEFAULT_PATH
		NO_SYSTEM_ENVIRONMENT_PATH
		)
      FIND_LIBRARY ( GLEWmx_DLL_PATH GLEW32mx PATHS ${GLEW_LOCATION}/bin
		NO_DEFAULT_PATH
		NO_SYSTEM_ENVIRONMENT_PATH
		)
	ELSE (APPLE)
	  # Unices
      FIND_LIBRARY ( GLEW_LIBRARY_PATH GLEW PATHS ${GLEW_LOCATION}/lib
		NO_DEFAULT_PATH
		NO_SYSTEM_ENVIRONMENT_PATH
		)
      FIND_LIBRARY ( GLEWmx_LIBRARY_PATH GLEWmx PATHS ${GLEW_LOCATION}/lib
		NO_DEFAULT_PATH
		NO_SYSTEM_ENVIRONMENT_PATH
		)
	ENDIF (APPLE)
  ENDIF ()
  # MUST reset
  SET(CMAKE_FIND_LIBRARY_SUFFIXES ${ORIGINAL_CMAKE_FIND_LIBRARY_SUFFIXES})

  SET( GLEW_GLEW_LIBRARY ${GLEW_LIBRARY_PATH} CACHE STRING "GLEW library")
  SET( GLEW_GLEWmx_LIBRARY ${GLEWmx_LIBRARY_PATH} CACHE STRING "GLEWmx library")

ENDIF ()
