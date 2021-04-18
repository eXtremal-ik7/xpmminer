cmake_minimum_required(VERSION 3.0)
project(xpmclient)

set (CMAKE_CXX_STANDARD 11)
option(STATIC_BUILD "Build with static libraries on Linux")
option(SANITIZER_ENABLED "Build with address sanitizer" OFF)
option(GPROF_ENABLED "Build with GNU profiler (use gprof ./exename -p > out.txt)" OFF)
option(EQUIHASH_ENABLED "Build equihash OpenCL GPU client" OFF)

if (SANITIZER_ENABLED)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=address")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address")
endif()

if (GPROF_ENABLED)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pg")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
endif()

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake
)
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
  ${CMAKE_SOURCE_DIR}/cmake
)

include(FindDependencies)
find_package(Protobuf REQUIRED)
find_package(OpenCL)

if (NOT CMAKE_CROSSCOMPILING)
  find_package(CUDA)
else()
  find_path(CUDA_TOOLKIT_INCLUDE
    device_functions.h
    PATHS ${CUDA_TOOLKIT_ROOT_DIR}
    PATH_SUFFIXES include
    NO_DEFAULT_PATH
  )
  if (NOT CUDA_TOOLKIT_INCLUDE STREQUAL "CUDA_TOOLKIT_INCLUDE-NOTFOUND")
    set(CUDA_FOUND 1)
  else()
    message("CUDA package not found")
  endif()
endif()

include_directories(
  ${GMP_INCLUDE_DIRECTORY}
  ${PROTOBUF_INCLUDE_DIRECTORY}
  ${OpenCL_INCLUDE_DIRS}
  ${OPENSSL_INCLUDE_DIR}
  ${CMAKE_SOURCE_DIR}/config4cpp/include
  ${CMAKE_BINARY_DIR}
)

if (CMAKE_COMPILER_IS_GNUCC)
  if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 6.0)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated -Wno-ignored-attributes")
  endif()
endif()

if (APPLE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DCL_SILENCE_DEPRECATION")
endif()

if (WIN32)
  message("Building for Win32")
  add_definitions(-D_WIN32 -D__WINDOWS__ -D__USE_MINGW_ANSI_STDIO=0)
else()
  message("Building for Linux")
  add_definitions(-DLINUX)
  if (STATIC_BUILD)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libgcc -static-libstdc++")
  endif()
endif()


set(LIBRARIES
  ${SODIUM_LIBRARY}
  ${GMP_LIBRARY}
  ${PROTOBUF_LIBRARY}
)

if (WIN32)
  set(LIBRARIES ${LIBRARIES} ws2_32 advapi32)
else()
  set(LIBRARIES ${LIBRARIES} dl)
endif()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
endif()

# XPM CUDA GPU Client

if (CUDA_FOUND)
  find_library(CUDA_driver_LIBRARY cuda
    ${CUDA_TOOLKIT_ROOT_DIR}/lib64
    ${CUDA_TOOLKIT_ROOT_DIR}/lib/x64
  )
  find_library(CUDA_nvrtc_LIBRARY nvrtc
    ${CUDA_TOOLKIT_ROOT_DIR}/lib64
    ${CUDA_TOOLKIT_ROOT_DIR}/lib/x64
  )

  add_executable(xpmclientnv
    xpmclient.cpp
    benchmarks.cpp
    cudautil.cpp
    sha256.cpp
    prime.cpp
    loguru.cpp
  )

  target_include_directories(xpmclientnv PUBLIC
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/xpm/cuda
    ${CUDA_INCLUDE_DIRS}
  )

  target_link_libraries(xpmclientnv ${LIBRARIES} ${CUDA_driver_LIBRARY} ${CUDA_nvrtc_LIBRARY})
  
  if (STATIC_BUILD)
    set_target_properties(xpmclientnv PROPERTIES LINK_SEARCH_END_STATIC 1)
  endif()
endif()
