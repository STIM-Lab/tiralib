# finds the STIM library (downloads it if it isn't present)
# set STIMLIB_PATH to the directory containing the stim subdirectory (the stim repository)

include(FindPackageHandleStandardArgs)



IF(NOT UNIX)
    set(CNPY_ROOT $ENV{PROGRAMFILES}/CNPY)
    find_path(CNPY_INCLUDE_DIR PATHS ${CNPY_ROOT} NAMES cnpy.h)
    find_library(CNPY_LIBRARY NAMES cnpy PATHS ${CNPY_ROOT})
ENDIF(NOT UNIX)


find_package_handle_standard_args(CNPY DEFAULT_MSG CNPY_LIBRARY CNPY_INCLUDE_DIR)

set(CNPY_LIBRARIES ${CNPY_LIBRARY})
set(CNPY_INCLUDE_DIRS ${CNPY_INCLUDE_DIR})