#include "CSieveOfEratosthenesL1Ext.h"
#include "benchmarks.h"
#include "primecoin.h"
#include "system.h"
#include "utils.h"

extern unsigned gMultipliersSize;
extern unsigned gWeaveDepth;
extern unsigned gMaxSieveSize;
extern unsigned gExtensionsNum;
extern unsigned gL1CacheSize;
extern unsigned gPrimorial;
extern unsigned gResultsBufferSize;
extern unsigned gSieveSize;


bool sieveResultsTest(const PrimeSource &primeSource,
                      mpz_class &fixedMultiplier,
                      const uint8_t *bytefield,
                      unsigned sieveSize,
                      unsigned chainLength,
                      unsigned depth,
                      unsigned extensionsNum,
                      unsigned type,
                      unsigned *candidatesCount)
{
  *candidatesCount = 0;
  unsigned sieveWords = sieveSize/32;
  unsigned i = 0;
  for (unsigned extNum = 0; extNum <= extensionsNum; extNum++) {
    const uint32_t *ptr = ((const uint32_t*)bytefield) + extNum*sieveWords;
    while (i < sieveWords) {
      uint32_t sieveWord = ptr[i];
      for (unsigned j = 0; j < 32; j++, sieveWord >>= 1) {    
        if (!(sieveWord & 0x1))
          continue;        

        unsigned M = i*32 + j/8 + 4*(j&0x7);
        if (M == 0)
          continue;
        
        ++*candidatesCount;
        unsigned candidateMultiplier = M * (1 << extNum);
        
        mpz_class chainOrigin = fixedMultiplier*candidateMultiplier;
        if (type == PRIME_CHAIN_CUNNINGHAM1) {
          if (!trialDivisionChainTest(primeSource, chainOrigin, true, chainLength, depth)) {
            fprintf(stderr, "multiplier %u can't pass trial division test\n", candidateMultiplier);
            return false;
          }
        } else if (type == PRIME_CHAIN_CUNNINGHAM2) {        
          if (!trialDivisionChainTest(primeSource, chainOrigin, false, chainLength, depth)) {
            fprintf(stderr, "multiplier %u can't pass trial division test\n", candidateMultiplier);
            return false;
          }
        } else {
          mpz_class chainOriginExtra = chainOrigin;                
          if (!trialDivisionChainTest(primeSource, chainOrigin, true, (chainLength+1)/2, depth) ||
              !trialDivisionChainTest(primeSource, chainOriginExtra, false, chainLength/2, depth)) {
            fprintf(stderr, "multiplier %u can't pass trial division test\n", candidateMultiplier);            
            return false;
          }
        }
      }
      
      i++;
      
    }
    
    i = sieveWords / 2;
  }
  
  return true;
}


// Function measure OpenCL kernel working time
// Source data building on CPU
void sieveBenchmark(PrimeSource &primeSource,
                    mpz_class &primorial,
                    OpenCLDeviceContext &device,
                    double difficulty,
                    bool checkResults)
{ 
  unsigned groupsNum = checkResults ? 2 : device.groupsNum;
  
  CPrimalityTestParams testParams(bitsFromDifficulty(difficulty));
  PrimecoinBlockHeader header;
  generateRandomHeader(&header, difficulty);
  
  typedef CSieveOfEratosthenesL1Ext* CSieveOfEratosthenesL1ExtPtr;
  std::unique_ptr<CSieveOfEratosthenesL1ExtPtr[]> sieves(new CSieveOfEratosthenesL1Ext*[groupsNum]);
  for (unsigned i = 0; i < groupsNum; i++)
    sieves[i] = new CSieveOfEratosthenesL1Ext(primeSource);

  mpz_class fixedMultipliers[groupsNum];
  uint32_t fixedMultipliersBufferGPU[GPUMultiprecisionLimbs*groupsNum];
  memset(fixedMultipliersBufferGPU, 0, sizeof(fixedMultipliersBufferGPU));
  for (unsigned i = 0; i < groupsNum; i++) {
    size_t limbsNumber;
    mpz_class blockHeaderHash;    
    if (!updateBlock(&header, blockHeaderHash, primeSource, testParams)) {
      fprintf(stderr, "sieveBenchmark error: cannot find block with correct SHA256 hash\n");
      exit(1);
    }

    fixedMultipliers[i] = primorial * blockHeaderHash;
    mpz_export(&fixedMultipliersBufferGPU[GPUMultiprecisionLimbs*i], &limbsNumber,
               -1, 4, 0, 0, fixedMultipliers[i].get_mpz_t());
  }
  
  cl_mem fixedMultipliersDevPtr;
  if (clUpload<uint32_t>(device, fixedMultipliersBufferGPU,
                         GPUMultiprecisionLimbs*groupsNum,
                         &fixedMultipliersDevPtr, CL_MEM_READ_WRITE)) {
    fprintf(stderr, "sieveBenchmark error: cannot allocate memory for sieve bitfields\n");
    exit(1);
  }
    
  cl_kernel kernel = device.kernels[CLKernelSieveBenchmark];
  clSetKernelArg(kernel, 0, sizeof(fixedMultipliersDevPtr), &fixedMultipliersDevPtr);
  clSetKernelArg(kernel, 1, sizeof(device.cunningham1), &device.cunningham1);
  clSetKernelArg(kernel, 2, sizeof(device.cunningham2), &device.cunningham2);
  clSetKernelArg(kernel, 3, sizeof(device.bitwin), &device.bitwin);
  clSetKernelArg(kernel, 4, sizeof(device.primesDevPtr), &device.primesDevPtr);
  clSetKernelArg(kernel, 5, sizeof(device.multipliers64DevPtr), &device.multipliers64DevPtr);
  clSetKernelArg(kernel, 6, sizeof(device.offsets64DevPtr), &device.offsets64DevPtr);    
  
  size_t sieveBufferSize = GPUMaxSieveSize * (gExtensionsNum+1) / 8;
  std::unique_ptr<uint8_t[]> cunningham1Bytefield(new uint8_t[sieveBufferSize]);
  std::unique_ptr<uint8_t[]> cunningham2Bytefield(new uint8_t[sieveBufferSize]);
  std::unique_ptr<uint8_t[]> bitwinBytefield(new uint8_t[sieveBufferSize]);
  
  for (unsigned i = 2; i <= 16; i += 2) {
    unsigned sieveSize = CSieveOfEratosthenesL1Ext::L1CacheSize*i;
    unsigned realSieveSize =
      sieveSize + sieveSize/2*CSieveOfEratosthenesL1Ext::ExtensionsNum;
    
    timeMark beginPoint = getTimeMark();    

#ifdef DEBUG_MINING_AMD_OPENCL    
    gpa_uint32 passCount = 1;       
    gpa_uint32 sessionId;
#else
    unsigned passCount = 1;
#endif

#ifdef DEBUG_MINING_AMD_OPENCL      
      if (GPA_BeginSession(&sessionId) != GPA_STATUS_OK) {
        fprintf(stderr, "GPA_BeginSession error!\n");
        return;
      }

      if (GPA_GetPassCount(&passCount) != GPA_STATUS_OK) {
        fprintf(stderr, "GPA_GetPassCount error!\n");
        return;
      }
      fprintf(stderr, " * Running %u passes\n", passCount);
#endif      
      
      for (unsigned pass = 0; pass < passCount; pass++) {
#ifdef DEBUG_MINING_AMD_OPENCL
        GPA_BeginPass();
        GPA_BeginSample(1);
#endif
        cl_uint roundsNum = sieveSize / GPUSieveWindowSize;
        size_t globalThreads[1] = { groupsNum*device.groupSize};
        size_t localThreads[1] = { device.groupSize };
        cl_event event;
        cl_int result;
        
        cl_kernel kernel = device.kernels[CLKernelSieveBenchmark];      
        clSetKernelArg(kernel, 7, sizeof(roundsNum), &roundsNum);
        if ((result = clEnqueueNDRangeKernel(device.queue, kernel, 1, 0,
                                             globalThreads, localThreads, 0, 0, &event)) != CL_SUCCESS) {
          fprintf(stderr, "clEnqueueNDRangeKernel error!\n");
          return;
        }
          
        if (clWaitForEvents(1, &event) != CL_SUCCESS) {
          fprintf(stderr, "clWaitForEvents error!\n");
          return;
        }
          
          clReleaseEvent(event);      
#ifdef DEBUG_MINING_AMD_OPENCL          
          GPA_EndSample();
          GPA_EndPass();      
#endif
      }
      
#ifdef DEBUG_MINING_AMD_OPENCL
      GPA_EndSession();      
      bool isReady = false;
      while (!isReady)
        GPA_IsSessionReady(&isReady, sessionId);
      
      printAMDCounters(sessionId, 1);
#endif      

    timeMark middlePoint = getTimeMark();
    
    // Calculate CPU sieve time and compare it to GPU time
    unsigned candidatesCount[groupsNum];

    { 
      for (unsigned groupIdx = 0; groupIdx < groupsNum; groupIdx++) {
        sieves[groupIdx]->reset(sieveSize,
                                (unsigned)difficulty,
                                gWeaveDepth + gPrimorial,
                                fixedMultipliers[groupIdx]);
        sieves[groupIdx]->Weave();
        candidatesCount[groupIdx] = 0;
      }
    }
    
    timeMark endPoint = getTimeMark();
    for (unsigned groupIdx = 0; groupIdx < groupsNum; groupIdx++) {
      sieves[groupIdx]->resetCandidateIterator();
      
      unsigned multiplier, type;
      while (sieves[groupIdx]->GetNextCandidateMultiplier(multiplier, type))
        candidatesCount[groupIdx]++;
    }
    
    // GPU results checking
    for (unsigned groupIdx = 0; groupIdx < groupsNum; groupIdx++) {
      cl_event events[3];
      cl_int results[3];
        
      results[0] = clEnqueueReadBuffer(device.queue, device.cunningham1, CL_TRUE,
                                       groupIdx*sieveBufferSize, sieveBufferSize,
                                       cunningham1Bytefield.get(), 0, NULL, &events[0]);

      results[1] = clEnqueueReadBuffer(device.queue, device.cunningham2, CL_TRUE,
                                       groupIdx*sieveBufferSize, sieveBufferSize,
                                       cunningham2Bytefield.get(), 0, NULL, &events[1]);        
        
      results[2] = clEnqueueReadBuffer(device.queue, device.bitwin, CL_TRUE,
                                       groupIdx*sieveBufferSize, sieveBufferSize,
                                       bitwinBytefield.get(), 0, NULL, &events[2]);
      if (results[0] != CL_SUCCESS ||
          results[1] != CL_SUCCESS ||
          results[2] != CL_SUCCESS ||
          clWaitForEvents(3, events) != CL_SUCCESS) {
        fprintf(stderr, " * Error: cannot copy memory from device to host\n");
        exit(1);
      }
        
      clReleaseEvent(events[0]);
      clReleaseEvent(events[1]);
      clReleaseEvent(events[2]);        
        
      if (!checkResults)
        continue;
        
      unsigned allCandidatesCount = 0, count;
        
      if (!sieveResultsTest(primeSource,
                            fixedMultipliers[groupIdx],
                            cunningham1Bytefield.get(),
                            sieveSize,
                            (unsigned)difficulty,
                            gWeaveDepth+gPrimorial,
                            gExtensionsNum,
                            PRIME_CHAIN_CUNNINGHAM1,
                            &count)) {
        fprintf(stderr, " * Error: OpenCL device produces wrong results! "
                        "(cunningham 1 chain, groupNum=%u)\n", groupIdx);
        exit(1);
      }
        
      allCandidatesCount += count;
        
      if (!sieveResultsTest(primeSource,
                            fixedMultipliers[groupIdx],
                            cunningham2Bytefield.get(),
                            sieveSize,
                            (unsigned)difficulty,
                            gWeaveDepth+gPrimorial,
                            gExtensionsNum,
                            PRIME_CHAIN_CUNNINGHAM2,
                            &count)) {
        fprintf(stderr, " * Error: OpenCL device produces wrong results! "
                        "(cunningham 2 chain, groupNum=%u)\n", groupIdx);
        exit(1);
      }
        
      allCandidatesCount += count;        
       
      if (!sieveResultsTest(primeSource,
                            fixedMultipliers[groupIdx],
                            bitwinBytefield.get(),
                            sieveSize,
                            (unsigned)difficulty,
                            gWeaveDepth+gPrimorial,
                            gExtensionsNum,
                            PRIME_CHAIN_BI_TWIN,
                            &count)) {
        fprintf(stderr, " * Error: OpenCL device produces wrong results! "
                        "(bitwin chain, groupNum=%u)\n", groupIdx);
        exit(1);
      }
        
      allCandidatesCount += count;
      
      if (allCandidatesCount != candidatesCount[groupIdx]) {
        fprintf(stderr,
                " * Error: different multipliers number detected, CPU: %u, GPU: %u",
                candidatesCount[groupIdx],
                allCandidatesCount);
        exit(1);
      }
    }
      
    
    double tm = usDiff(beginPoint, middlePoint) / 1000.0 / passCount;
    double hpTm = usDiff(middlePoint, endPoint) / 1000.0;
    printf("(%u, real %u) "
           "GPU(time %.3lfms speed %.3lf) "
           "L1Ext/Fastest CPU(time %.3lfms speed %.3lf) %.3lf times faster\n",
           sieveSize, realSieveSize,
           tm, realSieveSize/tm,
           hpTm, realSieveSize/hpTm,
           hpTm/tm);
    fflush(stdout);
  }
}


void gpuMinerBenchmark(OpenCLDeviceContext &device, double difficulty, unsigned roundsNum, bool compareWithCPU)
{
  mpz_class primorial;
  PrimeSource primeSource(1000000, gWeaveDepth+256);
  CPrimalityTestParams testParams(bitsFromDifficulty(difficulty));
  PrimorialFast(gPrimorial, primorial, primeSource);

  unsigned groupsNum = compareWithCPU ? 1 : device.groupsNum;
  PrimecoinBlockHeader header;
  mpz_class blockHeaderHash;

  std::unique_ptr<GPUNonceAndHash[]> nonceAndHash(new GPUNonceAndHash[groupsNum]); 
  std::unique_ptr<FermatQueue[]> queue(new FermatQueue[groupsNum]);
  std::unique_ptr<FermatTestResults[]> results_(new FermatTestResults[groupsNum]);
  memset(nonceAndHash.get(), 0, sizeof(GPUNonceAndHash)*groupsNum);
  memset(queue.get(), 0, sizeof(FermatQueue)*groupsNum);
  
  generateRandomHeader(&header, difficulty); 
  
  OpenCLNewBlockPrepare(device, groupsNum, header,
                        nonceAndHash.get(), queue.get());
  
  for (unsigned round = 0; round < roundsNum; round++) {
    timeMark gpuBegin = getTimeMark();
    OpenCLMiningRound(device, groupsNum, results_.get());
    timeMark gpuEnd = getTimeMark();
    
    if (compareWithCPU) {
      for (unsigned groupIdx = 0; groupIdx < groupsNum; groupIdx++) {
        FermatTestResults &results = results_[groupIdx];
        printf("[group=%u] %u results\n", groupIdx, results.size);
        for (unsigned mIdx = 0; mIdx < results.size; mIdx++) {
          if (results.resultMultipliers[mIdx] == 0)
            continue;
          
          unsigned nonce = results.resultNonces[mIdx];
          unsigned triedMultiplier = results.resultMultipliers[mIdx];
          unsigned type = results.resultTypes[mIdx];          
          
          uint8_t hash1[32];
          uint8_t hashData[32];    
          header.nonce = nonce;        
          sha256(hash1, &header, 80);
          sha256(hashData, hash1, 32);
                    
          mpz_class blockHeaderHash;
          mpz_import(blockHeaderHash.get_mpz_t(),
                     32 / sizeof(unsigned long),
                     -1,
                     sizeof(unsigned long),
                     -1,
                     0,
                     hashData);
          if (!ProbablePrimalityTestWithTrialDivisionFast(blockHeaderHash, 1000, primeSource, testParams)) {
            fprintf(stderr, " * Error: invalid nonce found on OpenCL device\n");
            
            for (unsigned i = 0; i < 1000; i++) {
              if (blockHeaderHash % primeSource.prime(i) == 0) {
                fprintf(stderr, "   * block header hash not prime (can be divided by %u)\n", primeSource.prime(i));
                return;
              }
            }      
            
            unsigned nLength;
            if (!FermatProbablePrimalityTestFast(blockHeaderHash, nLength, testParams, true))
              fprintf(stderr, "   * Fermat test failed\n");
            return;
          }

          mpz_class fixedMultiplier = blockHeaderHash*primorial;     
          mpz_class chainOrigin = fixedMultiplier * triedMultiplier;          

          {
            mpz_class chainOriginExtra = chainOrigin;
            unsigned chainLength = (unsigned)difficulty;
            if (type == PRIME_CHAIN_CUNNINGHAM1) {
              if (!trialDivisionChainTest(primeSource, chainOriginExtra, true, chainLength, gWeaveDepth+gPrimorial)) {
                fprintf(stderr, " * Error: C1-multiplier %u can't pass trial division test\n", triedMultiplier);
                continue;
              }
            } else if (type == PRIME_CHAIN_CUNNINGHAM2) {
              if (!trialDivisionChainTest(primeSource, chainOriginExtra, false, chainLength, gWeaveDepth+gPrimorial)) {
                fprintf(stderr, " * Error: C2-multiplier %u can't pass trial division test\n", triedMultiplier);
                continue;
              }
            } else {
              mpz_class chainOriginExtra2 = chainOrigin;
              if (!trialDivisionChainTest(primeSource, chainOriginExtra, true, (chainLength+1)/2, gWeaveDepth+gPrimorial) ||
                  !trialDivisionChainTest(primeSource, chainOriginExtra2, false, chainLength/2, gWeaveDepth+gPrimorial)) {
                fprintf(stderr, " * Error: BT-multiplier %u can't pass trial division test\n", triedMultiplier);
                continue;
              }
            }
          }
          
          testParams.candidateType = type;
          ProbablePrimeChainTestFast(chainOrigin, testParams);          
          if (results.resultChainLength[mIdx] != chainLengthFromBits(testParams.chainLength)) {
            fprintf(stderr, "multiplier [%u]%u difference, CPU: %u, GPU: %u\n",
                    nonce, triedMultiplier,
                    (unsigned)chainLengthFromBits(testParams.chainLength),
                    (unsigned)results.resultChainLength[mIdx]);
          }
        }
      }
    }
  
  double tm = usDiff(gpuBegin, gpuEnd);
  unsigned realSieveSize = gSieveSize + gSieveSize/2*gExtensionsNum;  
  printf("  * [%u] round time: %.3lf milliseconds, speed: %.3lf M\n", round,
         tm / 1000.0,
         realSieveSize*groupsNum / tm);
  fflush(stdout);
  }
}


void multiplyBenchmark(OpenCLDeviceContext &device,
                       unsigned mulOperandSize,
                       uint32_t elementsNum)
{
  unsigned limbsNum = elementsNum*mulOperandSize;
  std::unique_ptr<uint32_t[]> m1(new uint32_t[limbsNum]);
  std::unique_ptr<uint32_t[]> m2(new uint32_t[limbsNum]);
  std::unique_ptr<uint32_t[]> mR(new uint32_t[limbsNum*2]);
  std::unique_ptr<uint32_t[]> cpuR(new uint32_t[limbsNum*2]);

  for (unsigned i = 0; i < elementsNum; i++) {
    for (unsigned j = 0; j < mulOperandSize; j++) {
      m1[i*mulOperandSize + j] = rand32();
      m2[i*mulOperandSize + j] = rand32();
    }
  }
  
  cl_mem bufferM1;
  cl_mem bufferM2;
  cl_mem bufferResult;
  clUpload(device, m1.get(), limbsNum, &bufferM1, CL_MEM_READ_WRITE);
  clUpload(device, m2.get(), limbsNum, &bufferM2, CL_MEM_READ_WRITE);
  clAlloc<uint32_t>(device, limbsNum*2, &bufferResult, CL_MEM_READ_WRITE);

  cl_kernel kernel;
  if (mulOperandSize == 128/32) {
    kernel = device.kernels[CLKernelMultiplyBenchmark128];
  } else if (mulOperandSize == 256/32) {
    kernel = device.kernels[CLKernelMultiplyBenchmark256];
  } else if (mulOperandSize == 384/32) {
    kernel = device.kernels[CLKernelMultiplyBenchmark384];
  } else if (mulOperandSize == 448/32) {
    kernel = device.kernels[CLKernelMultiplyBenchmark448];
  }  else if (mulOperandSize == 512/32) {
    kernel = device.kernels[CLKernelMultiplyBenchmark512];  
  } else {
    fprintf(stderr, "Can't multiply %u-size operands on OpenCL device\n", mulOperandSize*32);
    return;
  }
  
  clSetKernelArg(kernel, 0, sizeof(bufferM1), &bufferM1);
  clSetKernelArg(kernel, 1, sizeof(bufferM2), &bufferM2);
  clSetKernelArg(kernel, 2, sizeof(bufferResult), &bufferResult);
  clSetKernelArg(kernel, 3, sizeof(elementsNum), &elementsNum);
  
  std::unique_ptr<mpz_class[]> cpuM1(new mpz_class[elementsNum]);
  std::unique_ptr<mpz_class[]> cpuM2(new mpz_class[elementsNum]);
  std::unique_ptr<mpz_class[]> cpuResult(new mpz_class[elementsNum]);
  
  for (unsigned i = 0; i < elementsNum; i++) {
    mpz_import(cpuM1[i].get_mpz_t(), mulOperandSize, -1, 4, 0, 0, &m1[i*mulOperandSize]);
    mpz_import(cpuM2[i].get_mpz_t(), mulOperandSize, -1, 4, 0, 0, &m2[i*mulOperandSize]);
    mpz_import(cpuResult[i].get_mpz_t(), mulOperandSize*2, -1, 4, 0, 0, &mR[i*mulOperandSize*2]);
  }
  
  timeMark gpuBegin = getTimeMark();
  
  {
    size_t globalThreads[1] = { device.groupsNum*device.groupSize };
    size_t localThreads[1] = { device.groupSize };
    cl_event event;
    cl_int result;
    if ((result = clEnqueueNDRangeKernel(device.queue,
                                         kernel,
                                         1,
                                         0,
                                         globalThreads,
                                         localThreads,
                                         0, 0, &event)) != CL_SUCCESS) {
      fprintf(stderr, "clEnqueueNDRangeKernel error!\n");
      return;
    }
    
    if (clWaitForEvents(1, &event) != CL_SUCCESS) {
      fprintf(stderr, "clWaitForEvents error!\n");
      return;
    }
    
    clReleaseEvent(event);
  }
  
  timeMark gpuEnd = getTimeMark();
  
  // Запуск умножения на CPU
  for (unsigned i = 0; i < elementsNum; i++) {
    unsigned gmpLimbsNum = cpuM1[i].get_mpz_t()->_mp_size;
    mp_limb_t *Operand1 = cpuM1[i].get_mpz_t()->_mp_d;
    mp_limb_t *Operand2 = cpuM2[i].get_mpz_t()->_mp_d;
    mp_limb_t *target = (mp_limb_t*)&cpuR[i*mulOperandSize*2];
    for (unsigned j = 0; j < 512; j++) {
      mpn_mul_n(target, Operand1, Operand2, gmpLimbsNum);
      memcpy(Operand1, target, gmpLimbsNum*sizeof(mp_limb_t));
    }
  }

  timeMark cpuEnd = getTimeMark();

  {
    cl_event event;
    cl_int result;
  
    result = clEnqueueReadBuffer(device.queue, bufferResult, CL_TRUE,
                                 0, limbsNum*2*4,
                                 mR.get(), 0, NULL, &event);
    if (clWaitForEvents(1, &event) != CL_SUCCESS) {
      fprintf(stderr, "clWaitForEvents error!\n");
      return;
    }
    
    clReleaseEvent(event);
  }

  for (unsigned i = 0; i < elementsNum; i++) {
    if (memcmp(&mR[i*mulOperandSize*2], &cpuR[i*mulOperandSize*2], 4*mulOperandSize*2) != 0) {
      fprintf(stderr, "element index: %u\n", i);
      fprintf(stderr, "gmp: ");
      for (unsigned j = 0; j < mulOperandSize*2; j++)
        fprintf(stderr, "%08X ", cpuR[i*mulOperandSize*2 + j]);
      fprintf(stderr, "\ngpu: ");
      for (unsigned j = 0; j < mulOperandSize*2; j++)
        fprintf(stderr, "%08X ", mR[i*mulOperandSize*2 + j]);
      fprintf(stderr, "\n");
      fprintf(stderr, "results differ!\n");
      break;
    }
  }
  
  double cpuTime = usDiff(gpuEnd, cpuEnd) / 1000000.0;
  double gpuTime = usDiff(gpuBegin, gpuEnd) / 1000000.0;
  printf("Multiply %u bits CPU time: %.3lf, GPU time: %.3lf (%.3lf times faster)\n",
         mulOperandSize*32, cpuTime, gpuTime, cpuTime/gpuTime);
          
}


void moduloBenchmark(OpenCLDeviceContext &device,
                     unsigned dividendOperandSize,
                     unsigned divisorOperandSize,
                     unsigned elementsNum)
  
{
  unsigned dividendLimbsNum = elementsNum*dividendOperandSize;
  unsigned divisorLimbsNum = elementsNum*divisorOperandSize;
  std::unique_ptr<uint32_t[]> dividends(new uint32_t[dividendLimbsNum]);
  std::unique_ptr<uint32_t[]> divisors(new uint32_t[divisorLimbsNum]);
  std::unique_ptr<uint32_t[]> gpuModulos(new uint32_t[divisorLimbsNum]);
  std::unique_ptr<uint32_t[]> cpuModulos(new uint32_t[divisorLimbsNum]);
  
  for (unsigned i = 0; i < elementsNum; i++) {
    for (unsigned j = 0; j < dividendOperandSize; j++)
      dividends[i*dividendOperandSize + j] = (j == (dividendOperandSize-4)) ? 1 + rand() % 2 : 0;
    switch (i%4) {
      case 0 :
        for (unsigned j = 0; j < divisorOperandSize; j++)
          divisors[i*divisorOperandSize + j] = rand32() << (rand() % 4);
        break;
      case 1 :
        for (unsigned j = 0; j < divisorOperandSize; j++)
          divisors[i*divisorOperandSize + j] = (j >= (divisorOperandSize-1)) ? 0 : rand32() << (rand() % 4);
        break;
      case 2 :
        for (unsigned j = 0; j < divisorOperandSize; j++)
          divisors[i*divisorOperandSize + j] = (j >= (divisorOperandSize-2)) ? 0 : rand32() << (rand() % 4);
        break;
      case 3 :
        for (unsigned j = 0; j < divisorOperandSize; j++)
          divisors[i*divisorOperandSize + j] = (j >= (divisorOperandSize-3)) ? 0 : rand32() << (rand() % 4);
        break;        
    }
  }
  
  cl_mem gpuDividendsBuffer;
  cl_mem gpuDivisorsBuffer;
  cl_mem gpuModulosBuffer;
  cl_kernel kernel;
  clUpload(device, dividends.get(), dividendLimbsNum, &gpuDividendsBuffer, CL_MEM_READ_ONLY);
  clUpload(device, divisors.get(), divisorLimbsNum, &gpuDivisorsBuffer, CL_MEM_READ_ONLY);
  clAlloc<uint32_t>(device, divisorLimbsNum, &gpuModulosBuffer, CL_MEM_READ_WRITE);
  if (dividendOperandSize == 384/32 && divisorOperandSize == 256/32) {
    kernel = device.kernels[CLKernelModulo384to256test];
  } else if (dividendOperandSize == 512/32 && divisorOperandSize == 384/32) {
    kernel = device.kernels[CLKernelModulo512to384test];
  } else if (dividendOperandSize == 640/32 && divisorOperandSize == 512/32) {
    kernel = device.kernels[CLKernelModulo640to512test];
  } else {
    fprintf(stderr, "Can't divide %ubits to %ubits in OpenCL (no such implementation)\n",
            dividendOperandSize*32,
            divisorOperandSize*32);
    return;
  }
  clSetKernelArg(kernel, 0, sizeof(gpuDividendsBuffer), &gpuDividendsBuffer);
  clSetKernelArg(kernel, 1, sizeof(gpuDivisorsBuffer), &gpuDivisorsBuffer);
  clSetKernelArg(kernel, 2, sizeof(gpuModulosBuffer), &gpuModulosBuffer);
  clSetKernelArg(kernel, 3, sizeof(elementsNum), &elementsNum);

  std::unique_ptr<mpz_class[]> cpuDividendsBuffer(new mpz_class[elementsNum]);
  std::unique_ptr<mpz_class[]> cpuDivisorsBuffer(new mpz_class[elementsNum]);
  std::unique_ptr<mpz_class[]> cpuModulosBuffer(new mpz_class[elementsNum]);
  for (unsigned i = 0; i < elementsNum; i++) {
    mpz_import(cpuDividendsBuffer[i].get_mpz_t(), dividendOperandSize, -1, 4, 0, 0, &dividends[i*dividendOperandSize]);
    mpz_import(cpuDivisorsBuffer[i].get_mpz_t(), divisorOperandSize, -1, 4, 0, 0, &divisors[i*divisorOperandSize]);
    mpz_import(cpuModulosBuffer[i].get_mpz_t(), divisorOperandSize, -1, 4, 0, 0, &cpuModulos[i*divisorOperandSize]);
  }

  timeMark gpuBegin = getTimeMark();
  
  {
    size_t globalThreads[1] = { device.groupsNum*device.groupSize};
    size_t localThreads[1] = { device.groupSize };
    cl_event event;
    cl_int result;
    if ((result = clEnqueueNDRangeKernel(device.queue,
                                         kernel,
                                         1,
                                         0,
                                         globalThreads,
                                         localThreads,
                                         0, 0, &event)) != CL_SUCCESS) {
      fprintf(stderr, "clEnqueueNDRangeKernel error!\n");
      return;
    }
      
    if (clWaitForEvents(1, &event) != CL_SUCCESS) {
      fprintf(stderr, "clWaitForEvents error!\n");
      return;
    }
      
    clReleaseEvent(event);
  }
  
  timeMark gpuEnd = getTimeMark();

  for (unsigned i = 0; i < elementsNum; i++)
    cpuModulosBuffer[i] = cpuDividendsBuffer[i] % cpuDivisorsBuffer[i];
  
  timeMark cpuEnd = getTimeMark();

  {
    cl_event event;
    cl_int result;
    
    result = clEnqueueReadBuffer(device.queue, gpuModulosBuffer, CL_TRUE,
                                 0, divisorLimbsNum*4,
                                 gpuModulos.get(), 0, NULL, &event);
    if (clWaitForEvents(1, &event) != CL_SUCCESS) {
      fprintf(stderr, "clWaitForEvents error!\n");
      return;
    }
    
    clReleaseEvent(event);
  }

  memset(cpuModulos.get(), 0, 4*divisorOperandSize*elementsNum);
  for (unsigned i = 0; i < elementsNum; i++) {
    size_t exportedLimbs;
    mpz_export(&cpuModulos[i*divisorOperandSize], &exportedLimbs, -1, 4, 0, 0, cpuModulosBuffer[i].get_mpz_t());
    if (memcmp(&gpuModulos[i*divisorOperandSize], &cpuModulos[i*divisorOperandSize], 4*divisorOperandSize) != 0) {
      fprintf(stderr, "element index: %u\n", i);
      fprintf(stderr, "gmp: ");
      for (unsigned j = 0; j < divisorOperandSize; j++)
        fprintf(stderr, "%08X ", cpuModulos[i*dividendOperandSize + j]);
      fprintf(stderr, "\ngpu: ");
      for (unsigned j = 0; j < divisorOperandSize; j++)
        fprintf(stderr, "%08X ", gpuModulos[i*dividendOperandSize + j]);
      fprintf(stderr, "\n");
      fprintf(stderr, "results differ!\n");
      break;
    }
  }
  
  double cpuTime = usDiff(gpuEnd, cpuEnd) / 1000000.0;
  double gpuTime = usDiff(gpuBegin, gpuEnd) / 1000000.0 / 32;
  printf("Modulo (%u %% %u) CPU time: %.3lf, GPU time: %.3lf (%.3lf times faster)\n",
         dividendOperandSize*32, divisorOperandSize*32,
         cpuTime, gpuTime, cpuTime/gpuTime);
}

void fermatTestBenchmark(OpenCLDeviceContext &device,
                         unsigned operandSize,
                         unsigned elementsNum)
{ 
  unsigned numberLimbsNum = elementsNum*operandSize;
  std::unique_ptr<uint32_t[]> numbers(new uint32_t[numberLimbsNum]);
  std::unique_ptr<uint32_t[]> gpuResults(new uint32_t[numberLimbsNum]);
  std::unique_ptr<uint32_t[]> cpuResults(new uint32_t[numberLimbsNum]);
  
  
  for (unsigned i = 0; i < elementsNum; i++) {
    for (unsigned j = 0; j < operandSize; j++)
      numbers[i*operandSize + j] = (j == operandSize-1) ? (1 << (i % 32)) : rand32();
    numbers[i*operandSize] |= 0x1; 
  }

  cl_mem gpuNumbersBuffer;
  cl_mem gpuResultsBuffer;
  cl_kernel kernel;
  clUpload(device, numbers.get(), numberLimbsNum, &gpuNumbersBuffer, CL_MEM_READ_ONLY);
  clAlloc<uint32_t>(device, numberLimbsNum, &gpuResultsBuffer, CL_MEM_READ_WRITE);
  if (operandSize == 256/32) {
    kernel = device.kernels[CLKernelFermatTestBenchmark256];
  } else if (operandSize == 384/32) {
    kernel = device.kernels[CLKernelFermatTestBenchmark384];
  } else if (operandSize == 448/32) {
    kernel = device.kernels[CLKernelFermatTestBenchmark448];
  } else {
    fprintf(stderr, "Can't do Fermat test on %ubit operand\n", operandSize*32);
    return;
  }
  
  clSetKernelArg(kernel, 0, sizeof(gpuNumbersBuffer), &gpuNumbersBuffer);
  clSetKernelArg(kernel, 1, sizeof(gpuResultsBuffer), &gpuResultsBuffer);
  clSetKernelArg(kernel, 2, sizeof(elementsNum), &elementsNum);

  std::unique_ptr<mpz_t[]> cpuNumbersBuffer(new mpz_t[elementsNum]);
  std::unique_ptr<mpz_t[]> cpuResultsBuffer(new mpz_t[elementsNum]);
  mpz_class mpzTwo = 2;
  mpz_class mpzE;
  mpz_import(mpzE.get_mpz_t(), operandSize, -1, 4, 0, 0, numbers.get());
  for (unsigned i = 0; i < elementsNum; i++) {
    mpz_init(cpuNumbersBuffer[i]);
    mpz_init(cpuResultsBuffer[i]);
    mpz_import(cpuNumbersBuffer[i], operandSize, -1, 4, 0, 0, &numbers[i*operandSize]);
    mpz_import(cpuResultsBuffer[i], operandSize, -1, 4, 0, 0, &cpuResults[i*operandSize]);
  }

  timeMark gpuBegin = getTimeMark();
  
#ifdef DEBUG_MINING_AMD_OPENCL    
  gpa_uint32 passCount = 1;       
  gpa_uint32 sessionId;
#else
  unsigned passCount = 1;
#endif  
  
#ifdef DEBUG_MINING_AMD_OPENCL      
  if (GPA_BeginSession(&sessionId) != GPA_STATUS_OK) {
    fprintf(stderr, "GPA_BeginSession error!\n");
    return;
  }
  
  if (GPA_GetPassCount(&passCount) != GPA_STATUS_OK) {
    fprintf(stderr, "GPA_GetPassCount error!\n");
    return;
  }
  fprintf(stderr, " * Running %u passes\n", passCount);
#endif   
  
  for (unsigned i = 0; i < passCount; i++) {
#ifdef DEBUG_MINING_AMD_OPENCL
    GPA_BeginPass();
    GPA_BeginSample(1);
#endif    
    size_t globalThreads[1] = { device.groupsNum*device.groupSize };
    size_t localThreads[1] = { device.groupSize };
    cl_event event;
    cl_int result;
    if ((result = clEnqueueNDRangeKernel(device.queue,
                                         kernel,
                                         1,
                                         0,
                                         globalThreads,
                                         localThreads,
                                         0, 0, &event)) != CL_SUCCESS) {
      fprintf(stderr, "clEnqueueNDRangeKernel error!\n");
      return;
    }
      
    cl_int error;
    if ((error = clWaitForEvents(1, &event)) != CL_SUCCESS) {
      fprintf(stderr, "clWaitForEvents error %i!\n", error);
      return;
    }
      
    clReleaseEvent(event);
#ifdef DEBUG_MINING_AMD_OPENCL          
    GPA_EndSample();
    GPA_EndPass();      
#endif    
  }
  timeMark gpuEnd = getTimeMark();
  
#ifdef DEBUG_MINING_AMD_OPENCL
  GPA_EndSession();      
  bool isReady = false;
  while (!isReady)
    GPA_IsSessionReady(&isReady, sessionId);
  
  printAMDCounters(sessionId, 1);
#endif        

  for (unsigned i = 0; i < elementsNum; i++) {
    mpz_sub_ui(mpzE.get_mpz_t(), cpuNumbersBuffer[i], 1);
    mpz_powm(cpuResultsBuffer[i], mpzTwo.get_mpz_t(), mpzE.get_mpz_t(), cpuNumbersBuffer[i]);
  }
  timeMark cpuEnd = getTimeMark();    

  {
    cl_event event;
    cl_int result;
    
    result = clEnqueueReadBuffer(device.queue, gpuResultsBuffer, CL_TRUE,
                                 0, numberLimbsNum*4,
                                 gpuResults.get(), 0, NULL, &event);
    if (clWaitForEvents(1, &event) != CL_SUCCESS) {
      fprintf(stderr, "clWaitForEvents error!\n");
      return;
    }
    
    clReleaseEvent(event);
  }
  
  memset(cpuResults.get(), 0, 4*operandSize*elementsNum);
  for (unsigned i = 0; i < elementsNum; i++) {
    size_t exportedLimbs;
    mpz_export(&cpuResults[i*operandSize], &exportedLimbs, -1, 4, 0, 0, cpuResultsBuffer[i]);
    if (memcmp(&gpuResults[i*operandSize], &cpuResults[i*operandSize], 4*operandSize) != 0) {
      fprintf(stderr, "element index: %u\n", i);
      fprintf(stderr, "gmp: ");
      for (unsigned j = 0; j < operandSize; j++)
        fprintf(stderr, "%08X ", cpuResults[i*operandSize + j]);
      fprintf(stderr, "\ngpu: ");
      for (unsigned j = 0; j < operandSize; j++)
        fprintf(stderr, "%08X ", gpuResults[i*operandSize + j]);
      fprintf(stderr, "\n");
      fprintf(stderr, "results differ!\n");
      break;
    }
  }
  
  double cpuTime = usDiff(gpuEnd, cpuEnd) / 1000000.0;
  double gpuTime = usDiff(gpuBegin, gpuEnd) / 1000000.0 ;
  printf("Fermat test %u bits CPU time: %.3lf, GPU time: %.3lf (%.3lf times faster)\n",
         operandSize*32, cpuTime, gpuTime, cpuTime/gpuTime);  
}
