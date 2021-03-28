cmake_minimum_required(VERSION 3.15)
project(GOOFIT CXX)

### Require out-of-source builds (needed to keep old build system) ###
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
message(STATUS ${LOC_PATH})
if(EXISTS ${LOC_PATH})
    message(FATAL_ERROR "You cannot build in the source directory (or any directory with a CMakeLists.txt file). Please make a build subdirectory.")
endif()

### Version number ###
set(GOOFIT_VERSION_MAJOR 2)
set(GOOFIT_VERSION_MINOR 0)
set(GOOFIT_VERSION_PATCH 0)
set(GOOFIT_VERSION "${GOOFIT_VERSION_MAJOR}.${GOOFIT_VERSION_MINOR}.${GOOFIT_VERSION_PATCH}")

configure_file (
    "${PROJECT_SOURCE_DIR}/include/goofit/Version.h.in"
    "${PROJECT_BINARY_DIR}/include/goofit/Version.h"
    )

add_subdirectory(extern/NVIDIA/thrust)

include(extern/GooStats/cmake-utilities/INITIALIZE_CUDA.cmake)
INITIALIZE_CUDA()

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING
        "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
        FORCE)
endif()

find_package(ROOT CONFIG REQUIRED COMPONENTS Minuit TreePlayer)

add_subdirectory(src)

add_library(goofit_lib INTERFACE)
target_link_libraries(goofit_lib INTERFACE goofit_base PDFs ROOT::Core ROOT::MathCore)

add_library(goofit_lib_NLL INTERFACE)
target_link_libraries(goofit_lib_NLL INTERFACE goofit_base PDFs_NLL ROOT::Core ROOT::MathCore)

# test
add_subdirectory(extern/google/googletest/googletest EXCLUDE_FROM_ALL)
enable_testing()
add_subdirectory(tests)
