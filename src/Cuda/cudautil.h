#ifndef __CUDALIB_H_
#define __CUDALIB_H_

#include <cuda.h>
#include <nvrtc.h>
#include "loguru.hpp"
#include <string>
#include <vector>

#define NVRTC_SAFE_CALL(x) \
do { \
  nvrtcResult result = x; \
  if (result != NVRTC_SUCCESS) { \
    LOG_F(ERROR, "\nerror: %i\nfailed with error %s at %s:%d", static_cast<int>(result), nvrtcGetErrorString(result), __FILE__, __LINE__); \
    exit(1); \
  } \
} while(0)

#define CUDA_SAFE_CALL(x) \
do { \
  CUresult result = x; \
  if (result != CUDA_SUCCESS) { \
    const char *msg; \
    cuGetErrorName(result, &msg); \
    LOG_F(ERROR, "\nerror: %i\nfailed with error %s at %s:%d\n", static_cast<int>(result), msg, __FILE__, __LINE__); \
    exit(1); \
  } \
} while(0)


template<typename T>
class cudaBuffer {
public:
  size_t _size;
  T *_hostData;
  CUdeviceptr _deviceData;
  
public:
  cudaBuffer() : _size(0), _hostData(0), _deviceData(0) {}
  ~cudaBuffer() {
    delete[] _hostData;
    if (_deviceData) {
      CUresult result = cuMemFree(_deviceData);
      if (result != CUDA_SUCCESS) {
        const char *msg;
        cuGetErrorName(result, &msg);
        LOG_F(ERROR, "CUDA memory free failed with error %s code: %i\n", msg, static_cast<int>(result));
      }
    }
  }
  
  CUresult init(size_t size, bool hostNoAccess) {
    _size = size;
    if (!hostNoAccess)
      _hostData = new T[size];
    return cuMemAlloc(&_deviceData, sizeof(T)*size);
  }
  
  CUresult copyToDevice() {
    return cuMemcpyHtoD(_deviceData, _hostData, sizeof(T)*_size);
  }

  CUresult copyToDevice(CUstream stream) {
    return cuMemcpyHtoDAsync(_deviceData, _hostData, sizeof(T)*_size, stream);
  }  
  
  CUresult copyToDevice(T *hostData) {
    return cuMemcpyHtoD(_deviceData, hostData, sizeof(T)*_size);
  }
  
  CUresult copyToDevice(T *hostData, CUstream stream) {
    return cuMemcpyHtoDAsync(_deviceData, hostData, sizeof(T)*_size, stream);
  }  
  
  CUresult copyToHost() {
    return cuMemcpyDtoH(_hostData, _deviceData, sizeof(T)*_size);
  }
  
  CUresult copyToHost(CUstream stream) {
    return cuMemcpyDtoHAsync(_hostData, _deviceData, sizeof(T)*_size, stream);
  }  
  
  T& get(int index) {
    return _hostData[index];
  }
  
  T& operator[](int index) {
    return _hostData[index];
  }  
};


bool cudaCompileKernel(const char *kernelName,
                       const std::vector<const char*> &sources,
                       const char **arguments,
                       int argumentsNum,
                       CUmodule *module,
                       int majorComputeCapability,
                       int minorComputeCapability,
                       bool needRebuild);

#endif //__CUDALIB_H_
