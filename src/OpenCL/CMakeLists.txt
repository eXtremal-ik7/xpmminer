if (DEBUG_MINING_AMD_OPENCL)
  find_path(GPUPERFAPI_INCLUDE_DIRECTORY GPUPerfAPI.h)
  find_library(GPUPERFAPI_OPENCL_LIBRARY GPUPerfAPICL)
  find_library(GPUPERFAPI_OPENGL_LIBRARY GPUPerfAPIGL)

  if (OS_LINUX)
    add_definitions(-D_LINUX)
  endif()

  include_directories(${GPUPERFAPI_INCLUDE_DIRECTORY})
endif()


add_executable(xpmclminer
  xpmclminer.cpp
  benchmarks.cpp
  utils.cpp
)

target_link_libraries(xpmclminer
  common
  blkmaker
  ${GMP_LIBRARY}
  ${GMPXX_LIBRARY}
  ${SSL_LIBRARY}
  ${CRYPTO_LIBRARY}
  ${CURL_LIBRARY}
  ${JANSSON_LIBRARY}
  ${OPENCL_LIBRARIES}
  ${NCURSES_LIBRARY}  
)

if (OS_LINUX)
  target_link_libraries(xpmclminer pthread)
endif()

if (DEBUG_MINING_AMD_OPENCL)
  target_link_libraries(xpmclminer
    ${GPUPERFAPI_OPENCL_LIBRARY}
    ${GPUPERFAPI_OPENGL_LIBRARY}
  )
endif()

install(TARGETS xpmclminer DESTINATION bin)
install(FILES 
  ${CMAKE_SOURCE_DIR}/OpenCL/kernel.cl
  ${CMAKE_SOURCE_DIR}/OpenCL/sha256.h
  ${CMAKE_SOURCE_DIR}/OpenCL/fmt.h  
  DESTINATION share/xpmminer)
