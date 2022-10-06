# - Try to find the Mindboggle cpp library
# Once done this will define
#
#  MINDBOGGLE_FOUND - system has MINDBOGGLE
#  MINDBOGGLE_INCLUDE_DIR - the MINDBOGGLE include directory
#  MINDBOGGLE_LIBRARY - the MINDBOGGLE library

if(MINDBOGGLE_FOUND)
    return()
endif()

find_path(MINDBOGGLE_INCLUDE_DIR MeshAnalyser.h
    PATHS
        $ENV{HOME}/Git/mindboggle
	    /media/DATOS/afernandez/libraries/mindboggle
    PATH_SUFFIXES vtk_cpp_tools
)

find_library(MINDBOGGLE_LIBRARY NAMES MeshAnalyser
    PATHS
        $ENV{HOME}/Git/mindboggle
        /media/DATOS/afernandez/libraries/mindboggle
    PATH_SUFFIXES vtk_cpp_tools/build
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MINDBOGGLE
    "\nMindboggle not found --- You can download it using:\n\tgit clone --recursive https://github.com/nipy/mindboggle.git ${CMAKE_SOURCE_DIR}/../mindboggle"
    MINDBOGGLE_INCLUDE_DIR)

find_package_handle_standard_args(libMeshAnalyser  DEFAULT_MSG
    MINDBOGGLE_LIBRARY MINDBOGGLE_INCLUDE_DIR)

mark_as_advanced(MINDBOGGLE_INCLUDE_DIR)

set(MINDBOGGLE_LIBRARY ${MINDBOGGLE_LIBRARY})
set(MINDBOGGLE_INCLUDE_DIR ${MINDBOGGLE_INCLUDE_DIR})
