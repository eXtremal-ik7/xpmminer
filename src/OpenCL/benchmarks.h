#include <stdint.h>
#include <gmpxx.h>

class OpenCLDeviceContext;
class PrimeSource;

void sieveBenchmark(PrimeSource &primeSource,
                    mpz_class &primorial,
                    OpenCLDeviceContext &device,
                    double difficulty,
                    bool checkResults); 

void gpuMinerBenchmark(OpenCLDeviceContext &device,
                       double difficulty,
                       unsigned roundsNum,
                       bool checkResults);

void multiplyBenchmark(OpenCLDeviceContext &device,
                       unsigned mulOperandSize,
                       uint32_t elementsNum);

void moduloBenchmark(OpenCLDeviceContext &device,
                     unsigned dividendOperandSize,
                     unsigned divisorOperandSize,
                     unsigned elementsNum);

void fermatTestBenchmark(OpenCLDeviceContext &device,
                         unsigned operandSize,
                         unsigned elementsNum);