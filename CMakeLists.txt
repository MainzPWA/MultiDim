cmake_minimum_required (VERSION 2.8)
project(MultiDim)
set(PROJ "MultiDim")

# check for in-source build, forbid it!
if(CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR AND NOT MSVC_IDE)
  message(FATAL_ERROR "\nIn-source build attempt detected!\n"
    "Please create a new directory (e.g. build) and run `cmake ..`.\n"
    "Also don't forget to delete the created CMakeCache.txt and CMakeFiles dir"
    )
endif()

message(STATUS "*** Build Type:      " ${CMAKE_BUILD_TYPE})
message(STATUS "*** Compiler Flags:  " ${DEFAULT_COMPILE_FLAGS})
message(STATUS "*** Install libs to: " ${LIBRARY_OUTPUT_PATH})
message(STATUS "*** Install bin to:  " ${EXECUTABLE_OUTPUT_PATH})

# use some sane default directories
if(NOT DEFINED EXECUTABLE_OUTPUT_PATH)
       set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")
endif()

# we check for empty string here, since the variable is indeed defined to an empty string
if(NOT CMAKE_BUILD_TYPE)
  # this also reflects this default in the GUI
        SET(CMAKE_BUILD_TYPE Debug CACHE STRING
    "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
    FORCE)
endif()

# really no optimization in debug mode
if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_DEBUG}   -O0")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native -s")
endif()

# require a fairly recent ROOT version
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
find_package(ROOT "5.28" REQUIRED)
link_directories(${ROOT_LIBRARY_DIR})

# set up various ROOT infrastructure things (dictionary, .rootmap)
set(DEPLIBS "")
ROOT_PREPARE_LIB("src" "src" "src" "src/${PROJ}Dict.cc" SRCS EXESRCS)

# construct the library from source files
include_directories("src" ${ROOT_INCLUDE_DIR})
add_library("lib${PROJ}" SHARED ${SRCS})
set_target_properties("lib${PROJ}" PROPERTIES LIBRARY_OUTPUT_NAME ${PROJ})

# construct the (stub) executable from source files and link with library
#add_executable(${PROJ} "src/${PROJ}.cc")
#target_link_libraries(${PROJ} "lib${PROJ}" ${ROOT_LIBRARIES}")
