cmake_minimum_required(VERSION 2.6)

project(dynamics CXX C)

# Set RPATH policy to suppress warning on Mac OS
# http://www.kitware.com/blog/home/post/510
set(CMAKE_MACOSX_RPATH 1)

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${CMAKE_CURRENT_SOURCE_DIR}/../gtsam/cmake")
include(GtsamMakeConfigFile)

# Load build type flags and default to Debug mode
include(GtsamBuildTypes)

# Use macros for creating tests/timing scripts
include(GtsamTesting)
include(GtsamPrinting)

#option(FEATURES_BUILD_SHARED_LIBRARY        "Enable/Disable building of a shared version of features" ON)
#option(FEATURES_BUILD_STATIC_LIBRARY        "Enable/Disable building of a static version of features" OFF)

# Sanity check building of libraries
#if (NOT FEATURES_BUILD_SHARED_LIBRARY AND NOT FEATURES_BUILD_STATIC_LIBRARY)
#    message(FATAL_ERROR "Both shared and static version of features library disabled - need to choose at least one!")
#endif()

# Find gtsam and boost
find_package(GTSAM REQUIRED NO_MODULE)

if(NOT Boost_USE_STATIC_LIBS)
	add_definitions(-DBOOST_ALL_DYN_LINK)
endif()
find_package(Boost 1.43 COMPONENTS serialization thread date_time REQUIRED)

# Only works with OpenCV 2.4.x
#find_package(OpenCV REQUIRED)
#message(STATUS "use opencv" ${OpenCV_VERSION})
#include_directories(${OpenCV_INCLUDE_DIRS})
#link_libraries(${OpenCV_LIBS_DEBUG})

#if ("${OpenCV_VERSION}" VERSION_EQUAL "3.0.0" OR "${OpenCV_VERSION}" VERSION_GREATER "3.0.0")
#  add_definitions(-DOPENCV_3_0_0)
#elseif ("${OpenCV_VERSION}" VERSION_GREATER "2.4.4")
#  add_definitions(-DOPENCV_2_4_5)
#elseif ("${OpenCV_VERSION}" VERSION_LESS "2.4.0")
#  message(FATAL_ERROR "cannnot work with opencv version lower than 2.4")
#endif()

# GCC 4.8+ complains about local typedefs which we use for shared_ptr etc.
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  if (NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.8)
    add_definitions(-Wno-unused-local-typedefs)
  endif()
endif()

# As of XCode 7, clang also complains about this
if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  if (NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 7.0)
    add_definitions(-Wno-unused-local-typedefs)
  endif()
endif()

#add_executable(dynamics cpp/dynamics.cpp cpp/dynamics.h)

add_subdirectory("cpp")


