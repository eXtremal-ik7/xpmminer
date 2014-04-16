#include "primecoin.h"
#include "system.h"
#include "utils.h"
#include <string.h>


static const char *gOpenCLKernelNames[] = {
  "searchNonce",
  "sieve",
  "FermatTestEnqueue",
  "FermatTestEnqueueBt",
  "sieveBenchmark",
  "multiplyBenchmark128",
  "multiplyBenchmark256",
  "multiplyBenchmark384",
  "multiplyBenchmark448",  
  "multiplyBenchmark512",  
  "modulo384to256test",
  "modulo512to384test",
  "modulo640to512test",
  "fermatTestBenchMark256",
  "fermatTestBenchMark384",
  "fermatTestBenchMark448"
};

uint32_t rand32()
{
  uint32_t result = rand();
  result = (result << 16) | rand();
  return result;
}

uint64_t rand64()
{
  uint64_t result = rand();
  result = (result << 16) | rand();
  result = (result << 16) | rand();
  result = (result << 16);
  return result;
}

static char *readFile(const char *name)
{
  std::unique_ptr<FILE, std::function<void(FILE*)>> sourceFile(
    fopen(name, "r"), [](FILE *file) { fclose(file); });
    
  if (!sourceFile)
    return logError<char*>(0, stderr, "Error: can't open file %s\n", name);
  
  fseek(sourceFile.get(), 0, SEEK_END);
  long sourceFileSize = ftell(sourceFile.get());
  if (sourceFileSize == 0)
    return logError<char*>(0, stderr, "Error: file %s is empty\n", name);
  
  std::unique_ptr<char[]> kernelSource(new char[sourceFileSize + 1]);
  fseek(sourceFile.get(), 0, SEEK_SET);
  if (fread(kernelSource.get(), 1, sourceFileSize, sourceFile.get()) != sourceFileSize)
    return logError<char*>(0, stderr, "Error: can't read from %s\n", name);
  
  kernelSource[sourceFileSize] = 0;
  return kernelSource.release();
}

int OpenCLInit(OpenCLPlatrormContext &ctx,
               const char *targetPlatformName,
               std::vector<unsigned> &deviceFilter,
               unsigned primorialIdx,
               unsigned sieveSize,
               unsigned sieveWindowSize,
               unsigned weaveDepth,
               unsigned extensionsNum,
               unsigned chainLength,
               bool useCPU,
               bool disableOpt)
{
  // Build compiler command line
  struct cmdLineArg {
    const char *name;
    unsigned value;
  };
  
  cmdLineArg arguments[] = {
    {"FixedPrimorial", primorialIdx},
    {"L1CacheSize", sieveWindowSize},
    {"WeaveDepth", weaveDepth},
    {"GroupSize", 256},
    {"ExtensionsNum", extensionsNum},
    {"ChainLength", chainLength},
    {"FixedRoundsNum", sieveSize / sieveWindowSize}
  };
  
  std::string cmdLine;
  
  cmdLine.append("-I");
  cmdLine.append(installPrefix());
  cmdLine.append("/share/xpmminer");
  cmdLine.push_back(' ');
  
  for (size_t i = 0; i < sizeof(arguments) / sizeof(cmdLineArg); i++) {
    char buffer[16];
    snprintf(buffer, sizeof(buffer), "%u", arguments[i].value);
    cmdLine.append(" -D");
    cmdLine.append(arguments[i].name);
    cmdLine.push_back('=');
    cmdLine.append(buffer);
  }
  
  cmdLine.push_back(' ');
  if (disableOpt)
    cmdLine.append("-cl-opt-disable -g ");
  
  // for isa dissasembly accesss
  if (strcmp(targetPlatformName, "Advanced Micro Devices, Inc.") == 0)
    cmdLine.append("-save-temps");
      
  const char *kernelFile = buildPath(PtData, "kernel.cl");
  const char *constantsFile = buildPath(PtData, "constants.h");
  makeConstantsHeader(constantsFile, weaveDepth+256);        
        
  // OpenCL preparing
  cl_int clResult;  
  cl_uint numPlatforms;
  if (clGetPlatformIDs(0, 0, &numPlatforms) != CL_SUCCESS ||
    !(numPlatforms > 0))
    return logError(1, stderr, "Error: no any OpenCL platforms found\n");
  
  // Find target OpenCL platform
  cl_platform_id targetPlatform = 0;
  {
    std::unique_ptr<cl_platform_id[]> platforms(new cl_platform_id[numPlatforms]);
    if (clGetPlatformIDs(numPlatforms, platforms.get(), 0) != CL_SUCCESS)
      return logError(1, stderr, "Error: no any OpenCL platforms found\n");
    
    for (cl_uint i = 0; i < numPlatforms; i++) {
      char platformName[128];
      if (clGetPlatformInfo(platforms[i], CL_PLATFORM_VENDOR,
                            sizeof(platformName), platformName, 0) != CL_SUCCESS)
        return logError(1, stderr, "Error: can't enumerate OpenCL platforms\n");

      printf("OpenCL platform available: %s\n", platformName);
      if (strcmp(platformName, targetPlatformName) == 0) {
        targetPlatform = platforms[i];
        break;
      }
    }
  }
  
  if (!targetPlatform)
    return logError(1, stderr, "Error: can't find plaform: %s\n", targetPlatformName);
  
  // Read kernel file
  const char *kernelSource = readFile(kernelFile);
  if (!kernelSource)
    return 1;
  
  // Retrieving OpenCL devices list
  ctx.devicesNum = 0;
  clGetDeviceIDs(targetPlatform,
                 useCPU ? CL_DEVICE_TYPE_CPU : CL_DEVICE_TYPE_GPU, 0, 0,
                 &ctx.devicesNum);
  if (!ctx.devicesNum)
    return logError(1, stderr, "Error: can't find '%s' GPUs\n", targetPlatformName);
  
  std::unique_ptr<cl_device_id[]> devices(new cl_device_id[ctx.devicesNum]);
  if (clGetDeviceIDs(targetPlatform,
                     useCPU ? CL_DEVICE_TYPE_CPU : CL_DEVICE_TYPE_GPU,
                     ctx.devicesNum,
                     devices.get(), 0) != CL_SUCCESS)
    return logError(1, stderr, "Error: can't query '%s' devices\n", targetPlatformName);
  
  ctx.devices.reset(new OpenCLDeviceContext[ctx.devicesNum]);
  for (cl_uint i = 0; i < ctx.devicesNum; i++) {
    OpenCLDeviceContext &deviceCtx = ctx.devices[i];
    deviceCtx.device = devices[i];
    
    // detect device name and number of compute units
    char deviceName[128] = {0};
    char boardName[128] = {0};
    cl_uint computeUnits;
    clGetDeviceInfo(deviceCtx.device, CL_DEVICE_NAME, sizeof(deviceName), deviceName, 0);
    clGetDeviceInfo(deviceCtx.device, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(computeUnits), &computeUnits, 0);
    fprintf(stderr, "[%u] %s; %u compute units\n", (unsigned)i, deviceName, computeUnits);
    
    deviceCtx.groupSize = 256;
    deviceCtx.groupsNum = computeUnits*4;
    
    // context creation
    cl_context_properties contextProperties[] = {
      CL_CONTEXT_PLATFORM, (cl_context_properties)targetPlatform, 0
    };      
    
    deviceCtx.context =
    clCreateContext(contextProperties, 1, &deviceCtx.device, 0, 0, &clResult);
    if (clResult != CL_SUCCESS || !deviceCtx.context)
      return logError(1, stderr, " * Error: can't create OpenCL context for this device");
    
    deviceCtx.queue =
      clCreateCommandQueue(deviceCtx.context, deviceCtx.device, 0, &clResult);
    if (clResult != CL_SUCCESS || !deviceCtx.queue)
      return logError(1, stderr, " * Error: can't create command queue for this device\n");

#ifdef DEBUG_MINING_AMD_OPENCL    
    if (GPA_OpenContext(deviceCtx.queue) != GPA_STATUS_OK ||
        GPA_EnableAllCounters() != GPA_STATUS_OK)
      return logError(1, stderr, "Error: GPA_OpenContext/GPA_EnableAllCounters failed\n");    
#endif
    
    deviceCtx.program = clCreateProgramWithSource(deviceCtx.context, 1, &kernelSource, 0, &clResult);
    if (clResult != CL_SUCCESS || deviceCtx.program == 0)
      return logError(1, stderr, "Error: can't compile OpenCL kernel %s\n", kernelFile);

    if (clBuildProgram(deviceCtx.program, 1, &deviceCtx.device, cmdLine.c_str(), 0, 0) != CL_SUCCESS) {    
      size_t logSize;
      clGetProgramBuildInfo(deviceCtx.program, deviceCtx.device, CL_PROGRAM_BUILD_LOG, 0, 0, &logSize);
      
      std::unique_ptr<char[]> log(new char[logSize]);
      clGetProgramBuildInfo(deviceCtx.program, deviceCtx.device, CL_PROGRAM_BUILD_LOG, logSize, log.get(), 0);
      fprintf(stderr, "%s\n", log.get());
      
      return logError(1, stderr, "Error: can't compile OpenCL kernel %s", kernelFile);
    }
    
    // retrieving handles of all used funtions
    deviceCtx.kernels.reset(new cl_kernel[CLKernelsNum]);
    for (size_t i = 0; i < CLKernelsNum; i++) {
      deviceCtx.kernels[i] = clCreateKernel(deviceCtx.program, gOpenCLKernelNames[i], &clResult);
      if (clResult != CL_SUCCESS)
        return logError(1, stderr, "Error: can't find function %s in kernel %s\n",
                        gOpenCLKernelNames[i],
                        kernelFile);
    }
  }

  return 0;
}

int OpenCLKernelsPrepare(OpenCLPlatrormContext &ctx,
                         mpz_class &primorial,
                         unsigned maxSieveSize,
                         unsigned maxWeaveDepth,
                         unsigned maxChainLength,
                         unsigned extensionsNum)
{
  uint32_t primorialBuffer[8];  
  {
    size_t limbsNumber;    
    memset(primorialBuffer, 0, sizeof(primorialBuffer));
    mpz_export(primorialBuffer, &limbsNumber, -1, 4, 0, 0, primorial.get_mpz_t());    
  }
  
  for (cl_uint i = 0; i < ctx.devicesNum; i++) {
    OpenCLDeviceContext &device = ctx.devices[i];
    size_t sieveBufferSize = maxSieveSize * (extensionsNum+1) / 8;
    size_t multipliersBufferSize = maxWeaveDepth*(maxChainLength+extensionsNum);
    if (clAlloc<uint8_t>(device, device.groupsNum*sieveBufferSize, &device.cunningham1, CL_MEM_READ_WRITE) != 0 ||
        clAlloc<uint8_t>(device, device.groupsNum*sieveBufferSize, &device.cunningham2, CL_MEM_READ_WRITE) != 0 ||
        clAlloc<uint8_t>(device, device.groupsNum*sieveBufferSize, &device.bitwin, CL_MEM_READ_WRITE) != 0 ||
        clAlloc<PrimecoinBlockHeader>(device, 1, &device.blockHeaderDevPtr, CL_MEM_READ_ONLY) != 0 ||
        clUpload<uint32_t>(device, primorialBuffer, 8, &device.primorialDevPtr, CL_MEM_READ_ONLY) ||
        clAlloc<GPUNonceAndHash>(device, device.groupsNum, &device.nonceAndHashDevPtr, CL_MEM_READ_WRITE) != 0 ||
        clAlloc<FermatQueue>(device, device.groupsNum, &device.c1QueueDevPtr, CL_MEM_READ_WRITE) != 0 ||
        clAlloc<FermatQueue>(device, device.groupsNum, &device.c2QueueDevPtr, CL_MEM_READ_WRITE) != 0 ||
        clAlloc<FermatQueue>(device, device.groupsNum, &device.btQueueDevPtr, CL_MEM_READ_WRITE) != 0 ||
        clAlloc<FermatTestResults>(device, device.groupsNum, &device.resultsDevPtr,CL_MEM_READ_WRITE) != 0) 
      return 1;
    
    clSetKernelArg(device.kernels[CLKernelSearchNonce], 0, sizeof(device.blockHeaderDevPtr), &device.blockHeaderDevPtr);
    clSetKernelArg(device.kernels[CLKernelSearchNonce], 1, sizeof(device.nonceAndHashDevPtr), &device.nonceAndHashDevPtr);
    
    clSetKernelArg(device.kernels[CLKernelSieve], 0, sizeof(device.cunningham1), &device.cunningham1);
    clSetKernelArg(device.kernels[CLKernelSieve], 1, sizeof(device.cunningham2), &device.cunningham2);
    clSetKernelArg(device.kernels[CLKernelSieve], 2, sizeof(device.bitwin), &device.bitwin);
    clSetKernelArg(device.kernels[CLKernelSieve], 3, sizeof(device.primorialDevPtr), &device.primorialDevPtr);
    clSetKernelArg(device.kernels[CLKernelSieve], 4, sizeof(device.nonceAndHashDevPtr), &device.nonceAndHashDevPtr);
    
    clSetKernelArg(device.kernels[CLFermatTestEnqueue], 0, sizeof(device.nonceAndHashDevPtr), &device.nonceAndHashDevPtr);
    clSetKernelArg(device.kernels[CLFermatTestEnqueue], 1, sizeof(device.primorialDevPtr), &device.primorialDevPtr);
    clSetKernelArg(device.kernels[CLFermatTestEnqueue], 2, sizeof(device.cunningham1), &device.cunningham1);
    clSetKernelArg(device.kernels[CLFermatTestEnqueue], 3, sizeof(device.cunningham2), &device.cunningham2);
    clSetKernelArg(device.kernels[CLFermatTestEnqueue], 4, sizeof(device.bitwin), &device.bitwin);
    clSetKernelArg(device.kernels[CLFermatTestEnqueue], 5, sizeof(device.c1QueueDevPtr), &device.c1QueueDevPtr);
    clSetKernelArg(device.kernels[CLFermatTestEnqueue], 6, sizeof(device.c2QueueDevPtr), &device.c2QueueDevPtr);
    clSetKernelArg(device.kernels[CLFermatTestEnqueue], 7, sizeof(device.btQueueDevPtr), &device.btQueueDevPtr);    
    clSetKernelArg(device.kernels[CLFermatTestEnqueue], 8, sizeof(device.resultsDevPtr), &device.resultsDevPtr);
    
    clSetKernelArg(device.kernels[CLFermatTestEnqueueBt], 0, sizeof(device.nonceAndHashDevPtr), &device.nonceAndHashDevPtr);
    clSetKernelArg(device.kernels[CLFermatTestEnqueueBt], 1, sizeof(device.primorialDevPtr), &device.primorialDevPtr);
    clSetKernelArg(device.kernels[CLFermatTestEnqueueBt], 2, sizeof(device.cunningham1), &device.cunningham1);
    clSetKernelArg(device.kernels[CLFermatTestEnqueueBt], 3, sizeof(device.cunningham2), &device.cunningham2);
    clSetKernelArg(device.kernels[CLFermatTestEnqueueBt], 4, sizeof(device.bitwin), &device.bitwin);
    clSetKernelArg(device.kernels[CLFermatTestEnqueueBt], 5, sizeof(device.c1QueueDevPtr), &device.c1QueueDevPtr);
    clSetKernelArg(device.kernels[CLFermatTestEnqueueBt], 6, sizeof(device.c2QueueDevPtr), &device.c2QueueDevPtr);
    clSetKernelArg(device.kernels[CLFermatTestEnqueueBt], 7, sizeof(device.btQueueDevPtr), &device.btQueueDevPtr);    
    clSetKernelArg(device.kernels[CLFermatTestEnqueueBt], 8, sizeof(device.resultsDevPtr), &device.resultsDevPtr);    
  }
  
  return 0;
}


bool OpenCLNewBlockPrepare(OpenCLDeviceContext &device,
                           unsigned groupsNum,
                           PrimecoinBlockHeader &header,
                           GPUNonceAndHash *nonceZeroBuffer,
                           FermatQueue *queueZeroBuffer)
                           
{
  cl_int results[5];
  cl_event events[5];

  results[0] = clEnqueueWriteBuffer(device.queue, device.blockHeaderDevPtr, CL_TRUE, 0,
                                    sizeof(header), &header,
                                   0, NULL, &events[0]);
  results[1] = clEnqueueWriteBuffer(device.queue, device.nonceAndHashDevPtr, CL_TRUE, 0,
                                   sizeof(GPUNonceAndHash)*groupsNum, nonceZeroBuffer,
                                   0, NULL, &events[1]);
  results[2] = clEnqueueWriteBuffer(device.queue, device.c1QueueDevPtr, CL_TRUE, 0,
                                    sizeof(FermatQueue)*groupsNum, queueZeroBuffer,
                                    0, NULL, &events[2]);
  results[3] = clEnqueueWriteBuffer(device.queue, device.c2QueueDevPtr, CL_TRUE, 0,
                                    sizeof(FermatQueue)*groupsNum, queueZeroBuffer,
                                    0, NULL, &events[3]);
  results[4] = clEnqueueWriteBuffer(device.queue, device.btQueueDevPtr, CL_TRUE, 0,
                                    sizeof(FermatQueue)*groupsNum, queueZeroBuffer,
                                    0, NULL, &events[4]);  
  
  if (results[0] != CL_SUCCESS ||
      results[1] != CL_SUCCESS ||
      results[2] != CL_SUCCESS ||
      results[3] != CL_SUCCESS ||
      results[4] != CL_SUCCESS) {
    fprintf(stderr, " * Error: cannot enqueue write buffer operation\n");
    return false;
  }
  
  if (clWaitForEvents(sizeof(events)/sizeof(cl_event), events) != CL_SUCCESS) {
    fprintf(stderr, " * Error: cannot write data to OpenCL device\n");
    return false;
  }
  
  return true;
}


bool OpenCLMiningRound(OpenCLDeviceContext &device,
                       unsigned groupsNum,
                       FermatTestResults *results)
{
  size_t globalThreads[1] = { groupsNum*device.groupSize};
  size_t localThreads[1] = { device.groupSize };
  cl_event events[5];
  cl_int result;
  if ((result = clEnqueueNDRangeKernel(device.queue, device.kernels[CLKernelSearchNonce],
                                       1, 0, globalThreads, localThreads, 0, 0, &events[0])) != CL_SUCCESS) {
     fprintf(stderr, "CLKernelSearchNonce launch error!\n");
     return false;
  }
  
  if ((result = clEnqueueNDRangeKernel(device.queue, device.kernels[CLKernelSieve],
                                       1, 0, globalThreads, localThreads, 1, &events[0], &events[1])) != CL_SUCCESS) {
     fprintf(stderr, "CLKernelSieve launch error!\n");
     return false;
  }
  
  if ((result = clEnqueueNDRangeKernel(device.queue, device.kernels[CLFermatTestEnqueue],
                                       1, 0, globalThreads, localThreads, 1, &events[1], &events[2])) != CL_SUCCESS) {
     fprintf(stderr, "CLFermatTestEnqueue launch error!\n");
     return false;
  }  
  
  if ((result = clEnqueueNDRangeKernel(device.queue, device.kernels[CLFermatTestEnqueueBt],
                                       1, 0, globalThreads, localThreads, 1, &events[1], &events[3])) != CL_SUCCESS) {
     fprintf(stderr, "CLFermatTestEnqueueBt launch error!\n");
     return false;
  }    
  
  if ((result = clEnqueueReadBuffer(device.queue, device.resultsDevPtr, CL_TRUE,
                                    0, sizeof(FermatTestResults)*groupsNum,
                                    results, 2, &events[2], &events[4])) != CL_SUCCESS) {
    fprintf(stderr, "nonceAndHash query error!\n");
    return false;
  }

  if (clWaitForEvents(1, &events[4]) != CL_SUCCESS) {
    fprintf(stderr, "clWaitForEvents error!\n");
    return false;
  }
      
  clReleaseEvent(events[0]);
  clReleaseEvent(events[1]);
  clReleaseEvent(events[2]);
  clReleaseEvent(events[3]);
  clReleaseEvent(events[4]);  
  return true;
}



void dumpArrayAsC(FILE *hFile, const char *constantName, uint32_t *array, unsigned size)
{
  fprintf(hFile, "__constant uint %s[%u] = {\n  ", constantName, size);
  for (unsigned i = 0; i < size; i++) {
    if (i) {
      fprintf(hFile, ", ");
      if ((i % 8) == 0)
        fprintf(hFile, "\n  ");
    }
    
    fprintf(hFile, "0x%08X", array[i]);
  }
  fprintf(hFile, "\n};\n\n");  
}

void dumpArrayAsC(FILE *hFile, const char *constantName, uint64_t *array, unsigned size)
{
  fprintf(hFile, "__constant ulong %s[%u] = {\n  ", constantName, size);
  for (unsigned i = 0; i < size; i++) {
    if (i) {
      fprintf(hFile, ", ");
      if ((i % 8) == 0)
        fprintf(hFile, "\n  ");
    }
    
    fprintf(hFile, "0x%016llX", (unsigned long long)array[i]);
  }
  fprintf(hFile, "\n};\n\n");  
}

int makeConstantsHeader(const char *fileName, unsigned weaveDepth)
{
  PrimeSource primeSource(1000000, weaveDepth);
  std::unique_ptr<FILE, std::function<void(FILE*)>> headerFile(
    fopen(fileName, "w+"), [](FILE *file) { fclose(file); });
  
  if (!headerFile)
    return logError<int>(0, stderr, "Error: can't open file %s\n", fileName);
  
  dumpArrayAsC(headerFile.get(), "primes", primeSource.primesPtr(), weaveDepth);
  dumpArrayAsC(headerFile.get(), "multipliers64", primeSource.multiplier64Ptr(), weaveDepth);
  dumpArrayAsC(headerFile.get(), "offsets64", primeSource.offsets64Ptr(), weaveDepth);
}

#ifdef DEBUG_MINING_AMD_OPENCL
void printAMDCounters(gpa_uint32 sessionId, gpa_uint32 sampleId)
{
  gpa_uint32 count;
  GPA_GetNumCounters(&count);
  for (gpa_uint32 index = 0; index < count; index++) {
    const char *name;
    const char *descrption;
    GPA_GetCounterName(index, &name);
    GPA_GetCounterDescription(index, &descrption);

    GPA_Type type;
    char counterText[256];
    
    union {
      gpa_uint32 ResultUInt32;
      gpa_uint64 ResultUInt64;
      gpa_float32 ResultFloat32;
      gpa_float64 ResultFloat64;      
    };
    
    if (GPA_GetCounterDataType(index, &type) != GPA_STATUS_OK)
      continue;
    
    strcpy(counterText, "can't get counter");
    switch (type) {
      case GPA_TYPE_INT32 :
      case GPA_TYPE_UINT32 :
        if (GPA_GetSampleUInt32(sessionId, sampleId, index, &ResultUInt32) == GPA_STATUS_OK)
          snprintf(counterText, sizeof(counterText), "%u", ResultUInt32);
        break;
      case GPA_TYPE_INT64 :
      case GPA_TYPE_UINT64 :
        if (GPA_GetSampleUInt64(sessionId, sampleId, index, &ResultUInt64) == GPA_STATUS_OK)
        snprintf(counterText, sizeof(counterText), "%llu", ResultUInt64);
        break;
      case GPA_TYPE_FLOAT32 :
        if (GPA_GetSampleFloat32(sessionId, sampleId, index, &ResultFloat32) == GPA_STATUS_OK)
          snprintf(counterText, sizeof(counterText), "%.3f", ResultFloat32);
        break;
      case GPA_TYPE_FLOAT64 :
        if (GPA_GetSampleFloat64(sessionId, sampleId, index, &ResultFloat64) == GPA_STATUS_OK)
          snprintf(counterText, sizeof(counterText), "%.3lf", ResultFloat64);        
        break; 
    }
    
    fprintf(stderr, "%s = %s\n", name, counterText); 
  }
}
#endif
