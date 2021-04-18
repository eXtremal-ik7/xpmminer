#include "cudautil.h"
#include <string.h>
#include <fstream>
#include <iostream>
#include <memory>

bool cudaCompileKernel(const char *kernelName,
                       const std::vector<const char*> &sources,
                       const char **arguments,
                       int argumentsNum,
                       CUmodule *module,
                       int majorComputeCapability,
                       int,
                       bool needRebuild) 
{
  std::ifstream testfile(kernelName);
  if(needRebuild || !testfile) {
    LOG_F(INFO, "compiling ...");
    
    std::string sourceFile;
    for (auto &i: sources) {
      std::ifstream stream(i);
      std::string str((std::istreambuf_iterator<char>(stream)), std::istreambuf_iterator<char>());
      sourceFile.append(str);
    }
    
    LOG_F(INFO, "source: %u bytes", (unsigned)sourceFile.size());
    if(sourceFile.size() < 1){
      LOG_F(ERROR, "source files not found or empty");
      return false;
    }
    
    nvrtcProgram prog;
    NVRTC_SAFE_CALL(
      nvrtcCreateProgram(&prog,
                         sourceFile.c_str(),
                         "xpm.cu",
                         0,
                         NULL,
                         NULL));

    nvrtcResult compileResult = nvrtcCompileProgram(prog, argumentsNum, arguments);
    
    // Obtain compilation log from the program.
    size_t logSize;
    NVRTC_SAFE_CALL(nvrtcGetProgramLogSize(prog, &logSize));
    char *log = new char[logSize];
    NVRTC_SAFE_CALL(nvrtcGetProgramLog(prog, log));
    delete[] log;
    if (compileResult != NVRTC_SUCCESS) {
      LOG_F(ERROR, "nvrtcCompileProgram error: %s", nvrtcGetErrorString(compileResult));
      return false;
    }
    
    // Obtain PTX from the program.
    size_t ptxSize;
    NVRTC_SAFE_CALL(nvrtcGetPTXSize(prog, &ptxSize));
    char *ptx = new char[ptxSize];
    NVRTC_SAFE_CALL(nvrtcGetPTX(prog, ptx));
    
    // Destroy the program.
    NVRTC_SAFE_CALL(nvrtcDestroyProgram(&prog));
    
    {
      std::ofstream bin(kernelName, std::ofstream::binary | std::ofstream::trunc);
      bin.write(ptx, ptxSize);
      bin.close();      
    }
    
    delete[] ptx;
  }
  
  std::ifstream bfile(kernelName, std::ifstream::binary);
  if(!bfile) {
    return false;
  }  
  
  bfile.seekg(0, bfile.end);
  size_t binsize = bfile.tellg();
  bfile.seekg(0, bfile.beg);
  if(!binsize){
    LOG_F(ERROR, "%s empty", kernelName);
    return false;
  }
  
  std::unique_ptr<char[]> ptx(new char[binsize+1]);
  bfile.read(ptx.get(), binsize);
  bfile.close();
  
  CUresult result = cuModuleLoadDataEx(module, ptx.get(), 0, 0, 0);
  if (result != CUDA_SUCCESS) {
    if (result == CUDA_ERROR_INVALID_PTX || result == CUDA_ERROR_UNSUPPORTED_PTX_VERSION) {
      LOG_F(WARNING, "GPU Driver version too old, update recommended");
      LOG_F(WARNING, "Workaround: downgrade version in PTX to 6.0 ...");
      char *pv = strstr(ptx.get(), ".version ");
      if (pv) {
        pv[9] = '6';
        pv[11] = '0';
      }

      CUDA_SAFE_CALL(cuModuleLoadDataEx(module, ptx.get(), 0, 0, 0));
    } else {
      const char *msg;
      cuGetErrorName(result, &msg);
      LOG_F(ERROR, "Loading CUDA module failed with error %s", msg);
      return false;
    }
  }

  return true;
}
