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

include(FetchContent)

### C++ settings ###
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING
        "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
        FORCE)
endif()


### host & device ###
set(DEVICE_LISTING CUDA OMP CPP TBB Auto)
set(HOST_LISTING OMP CPP TBB Auto)
mark_as_advanced(DEVICE_LISTING HOST_LISTING)

set(GOOFIT_DEVICE Auto CACHE STRING "The compute device, options are ${DEVICE_LISTING}")
set(GOOFIT_HOST Auto CACHE STRING "The compute device, options are ${HOST_LISTING}")

if(NOT ${GOOFIT_DEVICE} IN_LIST DEVICE_LISTING)
    message(FATAL_ERROR "You must select a device from ${DEVICE_LISTING}, not ${GOOFIT_DEVICE}")
endif()
if(NOT ${GOOFIT_HOST} IN_LIST HOST_LISTING)
    message(FATAL_ERROR "You must select a host from ${HOST_LISTING}, not ${HOST_DEVICE}")
endif()

if(GOOFIT_DEVICE STREQUAL Auto)
    find_package(CUDA 10)
    if(CUDA_FOUND)
        set(GOOFIT_DEVICE CUDA CACHE STRING "The compute device, options are ${DEVICE_LISTING}" FORCE)
    else()
        find_package(OpenMP)
        if(OpenMP_CXX_FOUND OR OpenMP_FOUND)
            set(GOOFIT_DEVICE OMP CACHE STRING "The compute device, options are ${DEVICE_LISTING}" FORCE)
        else()
            set(GOOFIT_DEVICE CPP CACHE STRING "The compute device, options are ${DEVICE_LISTING}" FORCE)
        endif()
    endif()
    message(STATUS "Auto device selection: ${GOOFIT_DEVICE}")
endif()

if(GOOFIT_HOST STREQUAL Auto)
    if(GOOFIT_DEVICE STREQUAL OMP)
        set(GOOFIT_HOST OMP CACHE STRING "The host device, options are ${HOST_LISTING}" FORCE)
    elseif(GOOFIT_DEVICE STREQUAL TBB)
        set(GOOFIT_HOST TBB CACHE STRING "The host device, options are ${HOST_LISTING}" FORCE)
    else()
        set(GOOFIT_HOST CPP CACHE STRING "The host device, options are ${HOST_LISTING}" FORCE)
    endif()
    message(STATUS "Auto host selection: ${GOOFIT_HOST}")
endif()

set_property(CACHE GOOFIT_DEVICE PROPERTY STRINGS ${DEVICE_LISTING})
set_property(CACHE GOOFIT_HOST PROPERTY STRINGS ${HOST_LISTING})

# Thrust
set(THRUST_REQUIRED_SYSTEMS ${GOOFIT_HOST} ${GOOFIT_DEVICE})
FetchContent_Declare(
  thrust
  GIT_REPOSITORY https://github.com/NVIDIA/thrust
  GIT_TAG        1.12.0
  )
FetchContent_MakeAvailable(thrust)
thrust_create_target(MyThrust HOST ${GOOFIT_HOST} DEVICE ${GOOFIT_DEVICE})
target_compile_features(MyThrust INTERFACE cxx_std_14)

set(GOOFIT_ARCH Auto CACHE STRING "The GPU Archetecture, can be Auto, All, Common, a number, or a name")

if(GOOFIT_DEVICE STREQUAL CUDA)
    enable_language(CUDA)
    find_package(CUDA 10)

    set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)
    add_definitions("-DTARGET_SM35")
    cuda_select_nvcc_arch_flags(ARCH_FLAGS ${GOOFIT_ARCH})
    message(STATUS "Compiling for GPU arch: ${ARCH_FLAGS}")
endif()

find_package(ROOT CONFIG REQUIRED COMPONENTS Minuit TreePlayer)

function(goofit_add_library GNAME)
    add_library(${GNAME} STATIC ${ARGN})
    if(GOOFIT_DEVICE STREQUAL CUDA)
	target_compile_options(${GNAME} PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
		${ARCH_FLAGS}
		>)
    else()
        SET_SOURCE_FILES_PROPERTIES(${ARGN} PROPERTIES LANGUAGE CXX)
        target_compile_options(${GNAME} PUBLIC -x c++)
    endif()
endfunction()

add_subdirectory(src)

add_library(goofit_lib INTERFACE)
target_link_libraries(goofit_lib INTERFACE goofit_base PDFs ROOT::Core ROOT::MathCore)

# test
include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        release-1.8.0
)
set(googletest_BUILD_TESTS OFF)
FetchContent_MakeAvailable(googletest)
enable_testing()
add_subdirectory(tests)
