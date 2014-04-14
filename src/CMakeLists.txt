cmake_minimum_required(VERSION 2.8)

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
  ${CMAKE_SOURCE_DIR}/cmake
)

include(FindOpenCL)
include(FindDependencies)

option(DEBUG_MINING_AMD_OPENCL "Print performance counters")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")

if (WIN32)
  include_directories(${CMAKE_SOURCE_DIR}/msvc/include)
  set(OS_WINDOWS 1)
elseif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
  set(OS_LINUX 1)
  set(OS_COMMONUNIX 1)
elseif (APPLE)
  set(OS_DARWIN 1)
  set(OS_COMMONUNIX 1)
elseif(CMAKE_SYSTEM_NAME STREQUAL "FreeBSD")
  set(OS_FREEBSD 1) 
  set(OS_COMMONUNIX 1)
endif()

if(MSVC)
  add_definitions(
    -D_CRT_SECURE_NO_DEPRECATE
    -D_CRT_SECURE_NO_WARNINGS
    -D_CRT_NONSTDC_NO_DEPRECATE
  )
endif(MSVC)

include(TestBigEndian)
TEST_BIG_ENDIAN(IS_BIGENDIAN)
configure_file(
  ${CMAKE_SOURCE_DIR}/include/config.h.in
  ${CMAKE_BINARY_DIR}/include/config.h
)

include_directories(
  ${CMAKE_SOURCE_DIR}/include
  ${CMAKE_BINARY_DIR}/include
  ${CMAKE_SOURCE_DIR}/blkmaker
  ${OPENCL_INCLUDE_DIRS}
)

add_subdirectory(blkmaker)
add_subdirectory(common)
add_subdirectory(CPU)
add_subdirectory(OpenCL)
