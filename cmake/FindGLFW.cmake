# Find GLFW 3
#
# GLFW_LIBRARIES
# GLFW_INCLUDE_DIRS.
# GLFW_FOUND
set(GLFW_ROOT $ENV{GLFW_ROOT})
IF(NOT UNIX)
    IF(NOT GLFW_ROOT)
        MESSAGE("ERROR: GLFW_ROOT must be set!")
    ENDIF(NOT GLFW_ROOT)

    FIND_PATH(GLFW_INCLUDE_DIRS DOC "Path to GLFW include directory."
            NAMES GLFW/glfw3.h
            PATHS ${GLFW_ROOT}/include)
    IF(MSVC15)
        FIND_LIBRARY(GLFW_LIBRARIES DOC "Absolute path to GLFW library."
                NAMES glfw3.lib
                PATHS ${GLFW_ROOT}/lib-vc2015)
    ELSEIF(MSVC13)
        FIND_LIBRARY(GLFW_LIBRARIES DOC "Absolute path to GLFW library."
                NAMES glfw3.lib
                PATHS ${GLFW_ROOT}/lib-vc2013)
    ELSEIF(MSVC12)
        FIND_LIBRARY(GLFW_LIBRARIES DOC "Absolute path to GLFW library."
                NAMES glfw3.lib
                PATHS ${GLFW_ROOT}/lib-vc2012)
    ELSEIF(MSVC10)
        FIND_LIBRARY(GLFW_LIBRARIES DOC "Absolute path to GLFW library."
                NAMES glfw3.lib
                PATHS ${GLFW_ROOT}/lib-vc2010)
    ELSEIF(MINGW)
        IF(CMAKE_CL_64)
            FIND_LIBRARY(GLFW_LIBRARIES DOC "Absolute path to GLFW library."
                    NAMES glfw3.dll
                    PATHS ${GLFW_ROOT}/lib-mingw-w64)
        ELSE(CMAKE_CL_64)
            FIND_LIBRARY(GLFW_LIBRARIES DOC "Absolute path to GLFW library."
                    NAMES glfw3.dll
                    PATHS ${GLFW_ROOT}/lib-mingw)
        ENDIF(CMAKE_CL_64)
    ELSE(MINGW)
        # Default to latest version of VC libs
        FIND_LIBRARY(GLFW_LIBRARIES DOC "Absolute path to GLFW library."
                NAMES glfw3.lib
                PATHS ${GLFW_ROOT}/lib-vc2015)
    ENDIF(MSVC15)
ELSE(NOT UNIX)
    FIND_PATH(GLFW_INCLUDE_DIRS DOC "Path to GLFW include directory."
            NAMES GLFW/glfw3.h
            PATHS
            /usr/include
            /usr/local/include
            /usr/target/include
            /sw/include
            /opt/local/include)

    FIND_LIBRARY(GLFW_LIBRARIES DOC "Absolute path to GLFW library."
            NAMES glfw3.dll
            PATHS
            /usr/local/lib
            /usr/lib
            /lib)
ENDIF(NOT UNIX)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GLFW DEFAULT_MSG GLFW_LIBRARIES GLFW_INCLUDE_DIRS)

mark_as_advanced(GLFW_INCLUDE_DIRS GLFW_LIBRARIES)
