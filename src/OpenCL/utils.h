#include "Debug.h"
#include <CL/cl.h>
#include <memory>
#include <vector>

enum OpenCLKernels {
  CLKernelSearchNonce = 0,
  CLKernelSieve,
  CLFermatTestEnqueue,
  CLFermatTestEnqueueBt,
  CLKernelSieveBenchmark,
  CLKernelMultiplyBenchmark128,
  CLKernelMultiplyBenchmark256,
  CLKernelMultiplyBenchmark384,
  CLKernelMultiplyBenchmark448,
  CLKernelMultiplyBenchmark512,
  CLKernelModulo384to256test,
  CLKernelModulo512to384test,
  CLKernelModulo640to512test,
  CLKernelFermatTestBenchmark256,
  CLKernelFermatTestBenchmark384,
  CLKernelFermatTestBenchmark448,
  CLKernelEmpty,
  CLKernelsNum
}; 

const unsigned GPUSieveWindowSize = 16384/2;
const unsigned GPUSieveMaxRoundsNum = 128*2;
const unsigned GPUMultiprecisionLimbs = 12;
const unsigned GPUMaxSieveSize = GPUSieveWindowSize*GPUSieveMaxRoundsNum;

#define FermatQueueChunks 16
#define FermatQueueBufferSize ((FermatQueueChunks*256) + (FermatQueueChunks*256)/8)

#pragma pack(push, 1)
struct GPUNonceAndHash {
  uint32_t hash[8*256];
  uint32_t nonce[256];
  uint32_t currentNonce;
  uint32_t totalNonces;  
  uint32_t align[2];
};

struct GPUResult {
  union {
    uint32_t multiplier;
    uint32_t multipliersNum;
  };
  uint32_t chainLength;
};

struct FermatQueue {
  uint32_t size;
  uint32_t _align1[3];
  
  uint32_t chainElements[3*4*FermatQueueBufferSize];
  uint32_t chainMultipliers[FermatQueueBufferSize];
  uint32_t chainLengths[FermatQueueBufferSize];
  uint32_t nonces[FermatQueueBufferSize];
};

struct FermatTestResults {
  uint32_t size;
  uint32_t _align1[3];
  uint32_t resultTypes[256*FermatQueueChunks];  
  uint32_t resultMultipliers[256*FermatQueueChunks];
  uint32_t resultChainLength[256*FermatQueueChunks];
  uint32_t resultNonces[256*FermatQueueChunks];
};

#pragma pack(pop)

struct OpenCLDeviceContext {
  cl_context context;
  cl_device_id device;
  cl_command_queue queue;
  std::unique_ptr<cl_kernel[]> kernels;
  
  unsigned groupSize;
  unsigned groupsNum;
  
  cl_mem primesDevPtr;
  cl_mem multipliers64DevPtr;
  cl_mem offsets64DevPtr;
  
  // bitfields, results of weaving
  cl_mem cunningham1;
  cl_mem cunningham2;
  cl_mem bitwin;
  cl_mem blockHeaderDevPtr;
  cl_mem primorialDevPtr;
  cl_mem nonceAndHashDevPtr;
  cl_mem c1QueueDevPtr;
  cl_mem c2QueueDevPtr;
  cl_mem btQueueDevPtr;
  cl_mem resultsDevPtr;
};

struct OpenCLPlatrormContext {
  cl_uint devicesNum;
  cl_program program;
  std::unique_ptr<OpenCLDeviceContext[]> devices;
};

uint32_t rand32();
uint64_t rand64();

template<typename T>
static T logError(T error, FILE *stream, const char *fmt, ...)
{
  va_list arguments;
  va_start(arguments, fmt);
  vfprintf(stream, fmt, arguments);
  return error;
}

template<typename T>
int clUpload(OpenCLDeviceContext &device,
             T *data,
             size_t N,
             cl_mem *buffer,
             cl_mem_flags flags)
{
  cl_int result;
  cl_event event;
  *buffer = clCreateBuffer(device.context, flags, sizeof(T)*N, 0, &result);
  if (!*buffer || result != CL_SUCCESS)
    return logError(1, stderr,
                    " * Error(clAlloc): cannot allocate %u*%u elements at device ?",
                    (unsigned)sizeof(T), (unsigned)N);
  
  if (clEnqueueWriteBuffer(device.queue, *buffer, CL_TRUE, 0,
                           sizeof(T)*N, data, 0, NULL, &event) != CL_SUCCESS)
    return logError(1, stderr,
                    " * Error(clAlloc): cannot allocate %u*%u elements at device ?",
                    (unsigned)sizeof(T), (unsigned)N);
    
  if (clWaitForEvents(1, &event) != CL_SUCCESS)
    logError(1, stderr, " * Error(clAlloc): clWaitForEvents error!\n");
  
  clReleaseEvent(event);
  dbgPrint(" * info: %u*%u bytes allocated\n", (unsigned)sizeof(T), (unsigned)N);  
  return 0;
}

template<typename T>
int clAlloc(OpenCLDeviceContext &device,
            size_t N,
            cl_mem *buffer,
            cl_mem_flags flags)
{
  cl_int result;
  cl_event event;
  *buffer = clCreateBuffer(device.context, flags, sizeof(T)*N, 0, &result);
  if (!*buffer || result != CL_SUCCESS)
    return logError(1, stderr,
                    " * Error(clAlloc): cannot allocate %u*%u elements at device ?",
                    (unsigned)sizeof(T), (unsigned)N);
  
  dbgPrint(" * info: %u*%u bytes allocated\n", (unsigned)sizeof(T), (unsigned)N);  
  return 0;
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
               bool disableOpt);

int OpenCLKernelsPrepare(OpenCLPlatrormContext &platform,
                         OpenCLDeviceContext &device,
                         PrimeSource &primeSource,
                         mpz_class &primorial,
                         unsigned maxSieveSize,
                         unsigned weaveDepth,
                         unsigned maxChainLength,
                         unsigned extensionsNum);

bool OpenCLNewBlockPrepare(OpenCLDeviceContext &device,
                           unsigned groupsNum,
                           PrimecoinBlockHeader &header,
                           GPUNonceAndHash *nonceZeroBuffer,
                           FermatQueue *queueZeroBuffer);

bool OpenCLMiningRound(OpenCLDeviceContext &device,
                       unsigned groupsNum,
                       FermatTestResults *results);
