# project
cmake_minimum_required (VERSION 2.6)
project (libviso2)

# gprof
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")
if(CMAKE_COMPILER_IS_GNUCXX)
  add_definitions(${CMAKE_CXX_FLAGS} "-Ofast")
  add_definitions(${CMAKE_CXX_FLAGS} "-fno-reorder-functions")
  add_definitions(${CMAKE_CXX_FLAGS} "-fno-inline")
endif()

# directories
set (LIBVISO2_SRC_DIR src)

# include directory
include_directories("${LIBVISO2_SRC_DIR}")

# use sse3 instruction set
SET(CMAKE_CXX_FLAGS "-msse3")

# sources
FILE(GLOB LIBVISO2_SRC_FILES "src/*.cpp")

# make release version
set(CMAKE_BUILD_TYPE Release)

# demo program
add_executable(viso2 ${LIBVISO2_SRC_FILES})
target_link_libraries (viso2 png)

