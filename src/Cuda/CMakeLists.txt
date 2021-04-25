set(LIBRARIES
  ${GMP_LIBRARY}
)
if (WIN32)
  set(LIBRARIES ${LIBRARIES} ws2_32 advapi32)
else()
  set(LIBRARIES ${LIBRARIES} dl)
endif()

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

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
endif()

if (CUDA_FOUND)
  find_library(CUDA_driver_LIBRARY cuda
    ${CUDA_TOOLKIT_ROOT_DIR}/lib64
    ${CUDA_TOOLKIT_ROOT_DIR}/lib/x64
  )
  find_library(CUDA_nvrtc_LIBRARY nvrtc
    ${CUDA_TOOLKIT_ROOT_DIR}/lib64
    ${CUDA_TOOLKIT_ROOT_DIR}/lib/x64
  )

include_directories(
  ${GMP_INCLUDE_DIRECTORY}
  ${OPENSSL_INCLUDE_DIR}
  ${CMAKE_BINARY_DIR}
)

add_executable(xpmcuda
    xpmclient.cpp
    benchmarks.cpp
    cudautil.cpp
    sha256.cpp
    prime.cpp
    loguru.cpp
  )

target_include_directories(xpmcuda PUBLIC
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/xpm/cuda
    ${CUDA_INCLUDE_DIRS}
  )

target_link_libraries(xpmcuda
  common
  blkmaker
  ${LIBRARIES}
  ${OPENSSL_CRYPTO_LIBRARY}
  ${OPENSSL_SSL_LIBRARY}
  ${CUDA_driver_LIBRARY}
  ${CUDA_nvrtc_LIBRARY}
  )

endif()