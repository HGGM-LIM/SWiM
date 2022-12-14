# based in cmake from gptoolbox
cmake_minimum_required(VERSION 3.1)
set(CMAKE_BUILD_TYPE_INIT Release)
project(sulcalwidth_mex)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
# templated lambdas
set( CMAKE_CXX_STANDARD 14 )
set( CMAKE_CXX_STANDARD_REQUIRED ON )

# Move output to this dir (sulcalwidth/mex)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})

# https://stackoverflow.com/q/46724267/148668
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
if(COMPILER_SUPPORTS_MARCH_NATIVE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()

# From libigl/cmake/libigl.cmake libigl/cmake/LibiglWindows.cmake
if(MSVC)
  if("${MSVC_RUNTIME}" STREQUAL "")
    set(MSVC_RUNTIME "static")
  endif()
  if(${MSVC_RUNTIME} STREQUAL "static")
    message(STATUS "MSVC -> forcing use of statically-linked runtime.")
    foreach(config ${CMAKE_CONFIGURATION_TYPES})
      string(TOUPPER ${config} config)
      string(REPLACE /MD /MT CMAKE_C_FLAGS_${config} "${CMAKE_C_FLAGS_${config}}")
      string(REPLACE /MD /MT CMAKE_CXX_FLAGS_${config} "${CMAKE_CXX_FLAGS_${config}}")
    endforeach()
    string(REPLACE "/MDd" "/MTd" CMAKE_CXX_FLAGS_DEBUG ${CMAKE_CXX_FLAGS_DEBUG})
  else()
    message(STATUS "MSVC -> forcing use of dynamically-linked runtime.")
    foreach(config ${CMAKE_CONFIGURATION_TYPES})
      string(TOUPPER ${config} config)
      string(REPLACE /MT /MD CMAKE_C_FLAGS_${config} "${CMAKE_C_FLAGS_${config}}")
      string(REPLACE /MT /MD CMAKE_CXX_FLAGS_${config} "${CMAKE_CXX_FLAGS_${config}}")
    endforeach()
    string(REPLACE "/MTd" "/MDd" CMAKE_CXX_FLAGS_DEBUG ${CMAKE_CXX_FLAGS_DEBUG})
  endif()
endif()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
set(SULCALWIDTH_MEX_ROOT "${PROJECT_SOURCE_DIR}")
set(SULCALWIDTH_MEX_EXTERNAL "${SULCALWIDTH_MEX_ROOT}/external")

# Find matlab
if(MATLAB_PROXY)
  set(Matlab_ROOT_DIR "${SULCALWIDTH_MEX_ROOT}/external/matlab")
endif()

find_package(Matlab REQUIRED COMPONENTS MEX_COMPILER MX_LIBRARY ENG_LIBRARY)

# Find Mindboggle
find_package(Mindboggle REQUIRED)
include_directories(${MINDBOGGLE_INCLUDE_DIR})

# Find VTK
find_package(VTK REQUIRED NO_MODULE)
include(${VTK_USE_FILE})


set( CMAKE_SHARED_LIBRARY_SUFFIX ".${Matlab_MEX_EXTENSION}" )  # set suffix to .mexa64
set( CMAKE_SHARED_LIBRARY_PREFIX )      # remove the "lib" prefix
set( CMAKE_SHARED_LIBRARY_CXX_FLAGS )       # remove the -fPIC option. mex does not accept the "-fPIC" option
set( CMAKE_POSITION_INDEPENDENT_CODE ON)
set( CMAKE_INSTALL_RPATH "\$ORIGIN" )

function(compile_each I_MEX_SOURCES I_LIBS )
  foreach(MEX_SOURCE ${I_MEX_SOURCES})
    message(STATUS "MEX_SOURCE: ${MEX_SOURCE}")
    get_filename_component(MEX_OUTPUT ${MEX_SOURCE} NAME_WE)
    message(STATUS "MEX_OUTPUT: ${MEX_OUTPUT}")
    add_library(${MEX_OUTPUT} SHARED ${MEX_SOURCE})
    # Not sure if this should be (WIN32) or (MSVC OR MSYS OR MINGW)
    # https://stackoverflow.com/a/40152725/148668
    if(WIN32)
      # https://stackoverflow.com/a/11095194/148668
      set_target_properties(${MEX_OUTPUT} PROPERTIES LINK_FLAGS "/export:mexFunction")
    endif()
    target_link_libraries(${MEX_OUTPUT} ${I_LIBS} ${MINDBOGGLE_LIBRARY} ${VTK_LIBRARIES})
    target_include_directories(${MEX_OUTPUT} PUBLIC ${Matlab_INCLUDE_DIRS})
    list(APPEND MEX_OUTPUTS ${MEX_OUTPUT})
  endforeach()
endfunction()

#file( GLOB MEX_SOURCES *.cpp )
list(APPEND CORE_LIBS ${Matlab_LIBRARIES})

#list(APPEND MEX_SOURCES "impaste.cpp")

# LIBIGL
add_definitions("-DMEX")
compile_each("\
travelDepth.cpp;\
"
  "${CORE_LIBS}")

if(APPLE)
  find_library(FOUNDATION_LIBRARY Foundation)
  find_library(APPKIT_LIBRARY AppKit)
  mark_as_advanced (FOUNDATION_LIBRARY APPKIT_LIBRARY)
  message(STATUS "MEX_SOURCE: impaste.cpp")
  set(MEX_OUTPUT "impaste")
  add_library(${MEX_OUTPUT} SHARED "impaste.cpp" "paste.mm")
  target_link_libraries(${MEX_OUTPUT} 
    ${CORE_LIBS} ${FOUNDATION_LIBRARY} ${APPKIT_LIBRARY})
  target_include_directories(${MEX_OUTPUT} PUBLIC ${Matlab_INCLUDE_DIRS})
endif(APPLE)