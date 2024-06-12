cmake_minimum_required(VERSION 3.0)

project(qci_interface VERSION 0.0)

if (CMAKE_CURRENT_SOURCE_DIR MATCHES "^${CMAKE_BINARY_DIRECTORY}$")
  message( FATAL_ERROR "CMake should not be run in the source directory" )
endif (CMAKE_CURRENT_SOURCE_DIR MATCHES "^${CMAKE_BINARY_DIRECTORY}$")

# git root directory should be one directory above the source directory
get_filename_component(GIT_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/.. ABSOLUTE)

# Cmake modules/macros are in a subdirectory to keep this file cleaner
# CMake modules implement some additional functionality to check
# for external dependencies
set(CMAKE_MODULE_PATH ${GIT_ROOT}/CMakeModules)

# if not specified by user, the standard build type is release
if(NOT CMAKE_BUILD_TYPE)
  # Release comes with -O3 by default
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

enable_language(CXX)

# the FORTRANCOMPILER check is in svn/CMakeModules
# it adds a compiler switch and sets flags
find_package(FORTRANCOMPILER)

# Add all the source files
file(GLOB QCI_SOURCES ${GIT_ROOT}/source/*.F90 ${GIT_ROOT}/source/*.f90 ${GIT_ROOT}/source/utils/*.F90 ${GIT_ROOT}/source/utils/*.f90)
# Remove interface
file(GLOB QCIINTERFACE ${GIT_ROOT}/source/QCI_interface.F90)

list(REMOVE_ITEM QCI_SOURCES ${QCIINTERFACE})

# Build a base library
add_library(qci_base ${QCI_SOURCES})
set_module_dir(qci_base)
    
# Build the interface library
add_library(qcilib ${QCIINTERFACE})
set_module_dir(qcilib)
set_module_depends(qcilib qci_base)
set_target_properties(qcilib PROPERTIES LINKER_LANGUAGE "Fortran")
set_target_properties(qcilib PROPERTIES COMPILE_DEFINITIONS "${COMPILE_DEFINITIONS};__QCI") 
target_link_libraries(qcilib qci_base)
