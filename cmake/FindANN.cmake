# - Try to find ANN
# Once done this will define
#  
#  ANN_FOUND        - system has ANN
#  ANN_INCLUDE_DIR  - the ANN include directory
#  ANN_LIBRARY      - Link these to use ANN
#   

IF (ANN_INCLUDE_DIRS)
  # Already in cache, be silent
  SET(ANN_FIND_QUIETLY TRUE)
ENDIF (ANN_INCLUDE_DIRS)

FIND_PATH( ANN_INCLUDE_DIR ANN/ANN.h
           PATHS "/usr/include" "C:/libs/ANN/include")

if( WIN32 )

	set(ANN_LIBRARY $ENV{ANN_PATH}\\ANN.lib)
	set(ANN_INCLUDE_DIR $ENV{ANN_PATH})

               
	 # Store the library dir. May be used for linking to dll!
	 # GET_FILENAME_COMPONENT( ANN_LIBRARY_DIR ${ANN_LIBRARY} PATH )

	 find_package_handle_standard_args(ANN DEFAULT_MSG ANN_INCLUDE_DIR)
 
else (WIN32)

 FIND_LIBRARY( ANN_LIBRARY
               NAMES ann ANN
               PATHS /lib /usr/lib /usr/lib64 /usr/local/lib )

endif( WIN32)


IF (ANN_INCLUDE_DIR AND ANN_LIBRARY)
  SET(ANN_FOUND TRUE)
ELSE (ANN_INCLUDE_DIR AND ANN_LIBRARY)
  SET( ANN_FOUND FALSE )
ENDIF (ANN_INCLUDE_DIR AND ANN_LIBRARY)

