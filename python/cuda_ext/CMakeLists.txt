cmake_minimum_required(VERSION 3.18)
project(tensorvote_cuda_ext LANGUAGES CXX CUDA)

# -----------------------------------------------------------------------
# Locate dependencies
# -----------------------------------------------------------------------

# Python interpreter and development headers (for pybind11)
find_package(Python3 REQUIRED COMPONENTS Interpreter Development)

# pybind11 — needs to be installed in the active venv:
#   pip install pybind11
# Then point CMAKE_PREFIX_PATH at the venv, or set pybind11_DIR explicitly.
find_package(pybind11 REQUIRED CONFIG)

# Boost — required because tensorvote.h unconditionally includes the Boost
# hypergeometric/beta headers (used by the 3D analytic plate path).
# Even though we only call the 2D TK/ATV paths at runtime, the headers are
# compiled and require Boost at build time.
# Typical locations: /usr, /usr/local, or set BOOST_ROOT env var.
find_package(Boost REQUIRED COMPONENTS headers)

# GLM — tensorvote.h uses glm::vec3, glm::mat3 in the 3D ATV path.
find_package(glm CONFIG REQUIRED)

# -----------------------------------------------------------------------
# Include paths
# -----------------------------------------------------------------------

# Root of tiralib (so that #include <tira/functions/tensorvote.h> resolves)
set(TIRALIB_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../")

# -----------------------------------------------------------------------
# CUDA settings
# -----------------------------------------------------------------------

# Allow __host__ __device__ functions in constexpr contexts (required by tira)
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-relaxed-constexpr")

# Generate code for the GPU present on this machine.
# Override with -DCMAKE_CUDA_ARCHITECTURES=<sm> if needed (e.g. 86 for Ampere).
if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES native)
endif()

# -----------------------------------------------------------------------
# Helper macro to build one pybind11 CUDA extension
# -----------------------------------------------------------------------
macro(add_tv_cuda_module target_name source_file)
    # pybind11_add_module does not understand CUDA sources natively, so we
    # build as a normal shared library and apply the same flags pybind11 would.
    add_library(${target_name} SHARED ${source_file})

    set_target_properties(${target_name} PROPERTIES
        # pybind11 modules are named without the "lib" prefix and with the
        # Python-specific suffix (.cpython-3xx-x86_64-linux-gnu.so).
        PREFIX ""
        SUFFIX ".${Python3_SOABI}${CMAKE_SHARED_MODULE_SUFFIX}"
        CUDA_SEPARABLE_COMPILATION ON
        POSITION_INDEPENDENT_CODE ON
    )

    target_compile_options(${target_name} PRIVATE
        $<$<COMPILE_LANGUAGE:CUDA>:
            --expt-relaxed-constexpr
            -Xcompiler -fPIC
        >
    )

    target_include_directories(${target_name} PRIVATE
        ${TIRALIB_INCLUDE_DIR}
        ${pybind11_INCLUDE_DIRS}
        ${Python3_INCLUDE_DIRS}
        ${Boost_INCLUDE_DIRS}
    )

    target_link_libraries(${target_name} PRIVATE
        pybind11::module
        glm::glm
    )

    # Install the .so next to the source (or set CMAKE_INSTALL_PREFIX)
    install(TARGETS ${target_name}
            LIBRARY DESTINATION "${CMAKE_CURRENT_SOURCE_DIR}")
endmacro()

# -----------------------------------------------------------------------
# Build both extensions
# -----------------------------------------------------------------------

add_tv_cuda_module(tensorvote_tk  tensorvote_tk.cu)
add_tv_cuda_module(tensorvote_atv tensorvote_atv.cu)
