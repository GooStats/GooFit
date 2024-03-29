cmake_minimum_required(VERSION 3.15)
project(GOOFIT LANGUAGES CXX)
list(APPEND CMAKE_MODULE_PATH
        ${PROJECT_SOURCE_DIR}/extern/arsenm/sanitizers-cmake/cmake)
list(APPEND CMAKE_CXX_SOURCE_FILE_EXTENSIONS cu)
find_package(Sanitizers)

# Add clang-tidy if available
option(CLANG_TIDY "Add clang-tidy")
option(CLANG_TIDY_FIX "Perform fixes for Clang-Tidy - changes source inplace")
if (CLANG_TIDY OR CLANG_TIDY_FIX)
  find_program(
          CLANG_TIDY_EXE
          NAMES "clang-tidy"
          DOC "Path to clang-tidy executable")

  if (CLANG_TIDY_EXE)
    if (CLANG_TIDY_FIX)
      set(CMAKE_CXX_CLANG_TIDY "${CLANG_TIDY_EXE}" "-fix")
    else ()
      set(CMAKE_CXX_CLANG_TIDY "${CLANG_TIDY_EXE}")
    endif ()
  endif ()
  find_program(
          CLANG_FORMAT_EXE
          NAMES "clang-format"
          DOC "Path to clang-format executable")
  set(CMAKE_CXX_CLANG_FORMAT "${CLANG_FORMAT_EXE}" "-fix")
endif ()

# Add CCache if available and if supported
if (NOT CLANG_TIDY_FIX)
  find_program(CCACHE_PROGRAM ccache)
  if (CCACHE_PROGRAM)
    set(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
    set(CMAKE_CUDA_COMPILER_LAUNCHER "${CCACHE_PROGRAM}") # CMake 3.9+
  endif ()
endif ()

# Require out-of-source builds (needed to keep old build system) ###
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
message(STATUS ${LOC_PATH})
if (EXISTS ${LOC_PATH})
  message(
          FATAL_ERROR
          "You cannot build in the source directory "
          "(or any directory with a CMakeLists.txt file)."
          "Please make a build subdirectory.")
endif ()

# Version number ###
set(GOOFIT_VERSION_MAJOR 2)
set(GOOFIT_VERSION_MINOR 0)
set(GOOFIT_VERSION_PATCH 0)
set(GOOFIT_VERSION
        "${GOOFIT_VERSION_MAJOR}.${GOOFIT_VERSION_MINOR}.${GOOFIT_VERSION_PATCH}")

configure_file("${PROJECT_SOURCE_DIR}/include/goofit/Version.h.in"
        "${PROJECT_BINARY_DIR}/include/goofit/Version.h")

add_subdirectory(extern/NVIDIA/thrust)
include(extern/GooStats/cmake-utilities/INITIALIZE_CUDA.cmake)
initialize_cuda()

set(default_build_type "RelWithDebInfo")
if (NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
          STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
          "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif ()

find_package(ROOT CONFIG REQUIRED COMPONENTS Minuit TreePlayer)

add_subdirectory(src)

add_library(goofit_lib INTERFACE)
target_link_libraries(goofit_lib INTERFACE goofit_base PDFs ROOT::Core
        ROOT::MathCore)

add_library(goofit_lib_NLL INTERFACE)
target_link_libraries(goofit_lib_NLL INTERFACE goofit_base PDFs_NLL ROOT::Core
        ROOT::MathCore)

# test
enable_testing()
add_subdirectory(extern/google/googletest/googletest EXCLUDE_FROM_ALL)
add_subdirectory(tests)
