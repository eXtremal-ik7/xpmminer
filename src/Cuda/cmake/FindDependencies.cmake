# GNU MP
find_path(GMP_INCLUDE_DIRECTORY gmp.h gmpxx.h)
find_library(GMP_LIBRARY gmp)
find_library(GMPXX_LIBRARY gmpxx)

# ZeroMQ
find_path(ZMQ_INCLUDE_DIRECTORY zmq.h)
if (WIN32)
  find_library(ZMQ_LIBRARY zmq.dll)
  find_library(SODIUM_LIBRARY sodium)
else()
  find_library(ZMQ_LIBRARY zmq)
  find_library(SODIUM_LIBRARY sodium)
endif()

# CLRadeonExtender
find_path(CLRX_INCLUDE_DIRECTORY CLRX/amdasm/Assembler.h)
find_path(CLRX_CONFIG_INCLUDE_DIRECTORY CLRX/Config.h)
if (WIN32)
  find_library(CLRX_AMDASM_LIBRARY CLRXAmdAsmStatic)
  find_library(CLRX_AMDBIN_LIBRARY CLRXAmdBinStatic)
  find_library(CLRX_UTILS_LIBRARY CLRXUtilsStatic)
else()
  find_library(CLRX_AMDASM_LIBRARY CLRXAmdAsm PATHS ${CMAKE_INSTALL_PREFIX}/lib64 STATIC)
  find_library(CLRX_AMDBIN_LIBRARY CLRXAmdBin ${CMAKE_INSTALL_PREFIX}/lib64 STATIC)
  find_library(CLRX_UTILS_LIBRARY CLRXUtils ${CMAKE_INSTALL_PREFIX}/lib64 STATIC)
endif()


set(CLRX_INCLUDE_DIRECTORIES ${CLRX_INCLUDE_DIRECTORY} ${CLRX_CONFIG_INCLUDE_DIRECTORY})
set(CLRX_LIBRARIES ${CLRX_AMDASM_LIBRARY} ${CLRX_AMDBIN_LIBRARY} ${CLRX_UTILS_LIBRARY})
