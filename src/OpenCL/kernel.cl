typedef unsigned char uint8_t;
typedef unsigned int uint32_t;
typedef unsigned long uint64_t;

#include "sha256.h"
#include "fmt.h"
#pragma OPENCL EXTENSION cl_amd_printf : enable

// You must define this parameters through compiler command line
// #define FixedPrimorial 19
// #define L1CacheSize (16384)
// #define WeaveDepth (6144+256)
// #define GroupSize (256)
// #define ExtensionsNum (9)
// #define ChainLength (10)
// #define FixedRoundsNum (80)

#define GPUMultiprecisionLimbs (12)
#define MaxRoundsNum (256)
#define MaxChainLength (20)
//#define MaxWeaveDepth (6144+256)
#define MaxSieveSize (L1CacheSize*MaxRoundsNum)
#define MaxSieveBufferSize (MaxSieveSize*(ExtensionsNum+1))
//#define MaxMultipliersBufferSize ((MaxChainLength+ExtensionsNum)*MaxWeaveDepth)
#define FixedSieveSize (FixedRoundsNum*L1CacheSize)

#define FirstLargePrimeIndex (2048)

__constant uint32_t binvert_limb_table[128] = {
  0x01, 0xAB, 0xCD, 0xB7, 0x39, 0xA3, 0xC5, 0xEF,
  0xF1, 0x1B, 0x3D, 0xA7, 0x29, 0x13, 0x35, 0xDF,
  0xE1, 0x8B, 0xAD, 0x97, 0x19, 0x83, 0xA5, 0xCF,
  0xD1, 0xFB, 0x1D, 0x87, 0x09, 0xF3, 0x15, 0xBF,
  0xC1, 0x6B, 0x8D, 0x77, 0xF9, 0x63, 0x85, 0xAF,
  0xB1, 0xDB, 0xFD, 0x67, 0xE9, 0xD3, 0xF5, 0x9F,
  0xA1, 0x4B, 0x6D, 0x57, 0xD9, 0x43, 0x65, 0x8F,
  0x91, 0xBB, 0xDD, 0x47, 0xC9, 0xB3, 0xD5, 0x7F,
  0x81, 0x2B, 0x4D, 0x37, 0xB9, 0x23, 0x45, 0x6F,
  0x71, 0x9B, 0xBD, 0x27, 0xA9, 0x93, 0xB5, 0x5F,
  0x61, 0x0B, 0x2D, 0x17, 0x99, 0x03, 0x25, 0x4F,
  0x51, 0x7B, 0x9D, 0x07, 0x89, 0x73, 0x95, 0x3F,
  0x41, 0xEB, 0x0D, 0xF7, 0x79, 0xE3, 0x05, 0x2F,
  0x31, 0x5B, 0x7D, 0xE7, 0x69, 0x53, 0x75, 0x1F,
  0x21, 0xCB, 0xED, 0xD7, 0x59, 0xC3, 0xE5, 0x0F,
  0x11, 0x3B, 0x5D, 0xC7, 0x49, 0x33, 0x55, 0xFF
};


void dbgPrint_v4(uint4 Reg)
{
  printf("%08X %08X %08X %08X ", Reg.x, Reg.y, Reg.z, Reg.w);
}

// Taking a modulo of long integer to 32-bit integer

#define longModuloByMulRound(mod, number, inversedMultiplier, shift) {\
  uint64_t dividend = (mod << 32) + number;\
  uint64_t quotient = mul_hi(dividend, inversedMultiplier) >> shift;\
  mod = dividend - quotient*divisor;\
}

uint32_t longModuloByMul256(uint4 nl0, uint4 nl1,
                            uint32_t divisor,
                            uint64_t inversedMultiplier,
                            uint32_t shift)
{
  uint64_t mod = 0;
  longModuloByMulRound(mod, nl1.w, inversedMultiplier, shift);
  longModuloByMulRound(mod, nl1.z, inversedMultiplier, shift);
  longModuloByMulRound(mod, nl1.y, inversedMultiplier, shift);
  longModuloByMulRound(mod, nl1.x, inversedMultiplier, shift);
  longModuloByMulRound(mod, nl0.w, inversedMultiplier, shift);
  longModuloByMulRound(mod, nl0.z, inversedMultiplier, shift);
  longModuloByMulRound(mod, nl0.y, inversedMultiplier, shift);
  longModuloByMulRound(mod, nl0.x, inversedMultiplier, shift);
  return mod;
}

uint32_t longModuloByMul384(uint4 nl0, uint4 nl1, uint4 nl2,
                            uint32_t divisor,
                            uint64_t inversedMultiplier,
                            uint32_t shift)
{
  uint64_t mod = 0;
  longModuloByMulRound(mod, nl2.w, inversedMultiplier, shift);
  longModuloByMulRound(mod, nl2.z, inversedMultiplier, shift);
  longModuloByMulRound(mod, nl2.y, inversedMultiplier, shift);
  longModuloByMulRound(mod, nl2.x, inversedMultiplier, shift);  
  longModuloByMulRound(mod, nl1.w, inversedMultiplier, shift);
  longModuloByMulRound(mod, nl1.z, inversedMultiplier, shift);
  longModuloByMulRound(mod, nl1.y, inversedMultiplier, shift);
  longModuloByMulRound(mod, nl1.x, inversedMultiplier, shift);  
  longModuloByMulRound(mod, nl0.w, inversedMultiplier, shift);
  longModuloByMulRound(mod, nl0.z, inversedMultiplier, shift);
  longModuloByMulRound(mod, nl0.y, inversedMultiplier, shift);
  longModuloByMulRound(mod, nl0.x, inversedMultiplier, shift);    
  return mod;
}

// Extended Euclidian algorithm
uint32_t intInvert(uint32_t a, uint32_t mod)
{
  uint32_t rem0 = mod, rem1 = a % mod, rem2;
  uint32_t aux0 = 0, aux1 = 1, aux2;
  uint32_t quotient, inverse;
  
  while (1) {
    if (rem1 <= 1) {
      inverse = aux1;
      break;
    }
    
    rem2 = rem0 % rem1;
    quotient = rem0 / rem1;
    aux2 = -quotient * aux1 + aux0;
    
    if (rem2 <= 1) {
      inverse = aux2;
      break;
    }
    
    rem0 = rem1 % rem2;
    quotient = rem1 / rem2;
    aux0 = -quotient * aux2 + aux1;
    
    if (rem0 <= 1) {
      inverse = aux0;
      break;
    }
    
    rem1 = rem2 % rem0;
    quotient = rem2 / rem0;
    aux1 = -quotient * aux0 + aux2;
  }
  
  return (inverse + mod) % mod;
}

unsigned int calculateOffset(unsigned int currentPrime,
                             unsigned int currentPrimeMod,
                             unsigned int offset)
{
  return offset >= currentPrimeMod ?
    offset - currentPrimeMod : currentPrime + offset - currentPrimeMod;
}


void clFillMemoryByGroup(__global void *buffer, unsigned size)
{
  __global uint32_t *ptr32 = (__global uint32_t*)buffer;
  for (unsigned i = get_local_id(0); i < size/4; i += GroupSize)
    ptr32[i] = 0xFFFFFFFF;
}

void flushLayer(unsigned layer,
                __global uint8 *cunningham1,
                __global uint8 *cunningham2,
                __global uint8 *bitwin,
                uint8 CR1,
                uint8 CR2,
                unsigned threadId)
{
  if (layer < ChainLength/2) {
    cunningham1[threadId] &= CR1;
    cunningham2[threadId] &= CR2;
    
    CR1 &= CR2;    
    bitwin[threadId] &= CR1;
  } else if (layer < (ChainLength + 1)/2) {
    cunningham1[threadId] &= CR1;
    cunningham2[threadId] &= CR2;
    bitwin[threadId] &= CR1;
  } else if (layer < ChainLength) {
    cunningham1[threadId] &= CR1;
    cunningham2[threadId] &= CR2;
  }
}

// GPU weave procedure. For understanding, look first to CSieveOfEratosthenesL1Ext
void weave(uint4 M0, uint4 M1, uint4 M2,
           __global uint32_t *cunningham1Bitfield,
           __global uint32_t *cunningham2Bitfield,
           __global uint32_t *bitwinBitfield,
           __local uint8_t *localCunningham1,
           __local uint8_t *localCunningham2,
           __constant uint32_t *primes,
           __global uint64_t *multipliers64,
           __global uint32_t *offsets64,
           unsigned roundsNum)
{
  const unsigned primesPerThread = WeaveDepth / GroupSize;  
  const unsigned layersNum = ChainLength + ExtensionsNum;
  const unsigned threadId = get_local_id(0);
  unsigned sieveBytes = L1CacheSize * roundsNum / 8;
  unsigned sieveWords = sieveBytes / 4;

  uint32_t inverseModulos[256]; // <--- large buffer for global memory placing
  uint32_t inverseModulosCurrent[WeaveDepth / GroupSize];          
  
  // Set all bits of output buffers to 1
  clFillMemoryByGroup(cunningham1Bitfield, sieveBytes);
  clFillMemoryByGroup(cunningham2Bitfield, sieveBytes);
  clFillMemoryByGroup(bitwinBitfield, sieveBytes);
  for (unsigned extNum = 1; extNum <= ExtensionsNum; extNum++) {
    clFillMemoryByGroup(cunningham1Bitfield + extNum*sieveWords + sieveWords/2, sieveBytes/2);
    clFillMemoryByGroup(cunningham2Bitfield + extNum*sieveWords + sieveWords/2, sieveBytes/2);
    clFillMemoryByGroup(bitwinBitfield + extNum*sieveWords + sieveWords/2, sieveBytes/2);
  }

  unsigned primeIdx = FixedPrimorial+threadId;
  for (unsigned j = 0; j < WeaveDepth/get_local_size(0); j++, primeIdx += GroupSize) {
    unsigned currentPrime = primes[primeIdx];
    unsigned mod = longModuloByMul384(M0, M1, M2,
                                  currentPrime,
                                  multipliers64[primeIdx],
                                  offsets64[primeIdx]);          

    inverseModulos[j] = intInvert(mod, currentPrime);
  }

  const unsigned L1CacheWords = L1CacheSize/4;
  __local uint32_t *localCunningham1_32 = (__local uint32_t*)localCunningham1;
  __local uint32_t *localCunningham2_32 = (__local uint32_t*)localCunningham2;

  for (unsigned round = 0; round < roundsNum; round += 8) {
    unsigned lowIdx = L1CacheSize * round;
    for (unsigned i = 0; i < primesPerThread; i++)
      inverseModulosCurrent[i] = inverseModulos[i];
    
    for (unsigned layer = 0; layer < layersNum; layer++) {
      if (layer >= ChainLength && round < roundsNum/2)
        break;

      barrier(CLK_LOCAL_MEM_FENCE);      
      for (unsigned i = 0, index = get_local_id(0); i < L1CacheWords/GroupSize; i++, index += GroupSize) {
        uint32_t X = 0xFFFFFFFF;
        localCunningham1_32[index] = X;
        localCunningham2_32[index] = X;        
      }
     
      unsigned primeIdx;
      for (unsigned j = 0, primeIdx = FixedPrimorial+threadId; j < FirstLargePrimeIndex/get_local_size(0); j++, primeIdx += GroupSize) {
        unsigned offset;
        unsigned offset2;
        const uint32_t currentPrime = primes[primeIdx];
        const uint32_t inverseModulo = inverseModulosCurrent[j];
        const uint32_t currentPrimeMod = lowIdx % currentPrime;

        offset = calculateOffset(currentPrime, currentPrimeMod, inverseModulo);
        offset2 = calculateOffset(currentPrime, currentPrimeMod, currentPrime - inverseModulo);
        
        for (unsigned iter = 0; iter < 8; iter++) {        
          barrier(CLK_LOCAL_MEM_FENCE);            
          uint8_t fill = ~(1 << iter);

          unsigned maxOffset = max(offset, offset2);
          while (maxOffset < L1CacheSize) {
            localCunningham1[offset] &= fill;
            localCunningham2[offset2] &= fill;
            offset += currentPrime;
            offset2 += currentPrime;
            maxOffset += currentPrime;
          }
          
          if (offset < L1CacheSize) {
            localCunningham1[offset] &= fill;
            offset += currentPrime;
          }
          
          if (offset2 < L1CacheSize) {
            localCunningham2[offset2] &= fill;
            offset2 += currentPrime;
          }          

          offset -= L1CacheSize;
          offset2 -= L1CacheSize;        
        }
        
        inverseModulosCurrent[j] = (inverseModulo & 0x1) ?
          (inverseModulo + currentPrime) / 2 : inverseModulo / 2;
      }


      for (unsigned j = FirstLargePrimeIndex/GroupSize,
             primeIdx = FirstLargePrimeIndex+FixedPrimorial+threadId;
           j < primesPerThread;
           j++, primeIdx += GroupSize) {
        const uint32_t currentPrime = primes[primeIdx];
        const uint32_t currentPrimeMod = lowIdx % currentPrime;        
        const uint32_t inverseModulo = inverseModulosCurrent[j];
       
        unsigned offset = calculateOffset(currentPrime, currentPrimeMod, inverseModulo);
        unsigned offset2 = calculateOffset(currentPrime, currentPrimeMod, currentPrime - inverseModulo);
      
        for (unsigned iter = 0; iter < 8; iter++) {
          uint8_t fill = ~(1 << iter);          
          if (offset < L1CacheSize) {
            localCunningham1[offset] &= fill;
            offset += currentPrime;
          }
          if (offset2 < L1CacheSize) {
            localCunningham2[offset2] &= fill;
            offset2 += currentPrime;
          }
          
          offset -= L1CacheSize;
          offset2 -= L1CacheSize;
        }
        
        inverseModulosCurrent[j] = (inverseModulo & 0x1) ?
          (inverseModulo + currentPrime) / 2 : inverseModulo / 2;
      } 

      barrier(CLK_LOCAL_MEM_FENCE);

      // Flush window to "main" range
      __global uint32_t *cunningham1 = cunningham1Bitfield + lowIdx/32;
      __global uint32_t *cunningham2 = cunningham2Bitfield + lowIdx/32;
      __global uint32_t *bitwin = bitwinBitfield + lowIdx/32;

      uint8 CR1 = vload8(threadId, localCunningham1_32);
      uint8 CR2 = vload8(threadId, localCunningham2_32);
      
      flushLayer(layer, cunningham1, cunningham2, bitwin, CR1, CR2, threadId);
      
      // Flush window to extensions
      if (round < roundsNum/2)
        continue;
      
      unsigned extNumMin = layer < ChainLength ? 1 : layer - ChainLength + 1;
      unsigned extNumMax = layer < ExtensionsNum ? layer : ExtensionsNum;
      for (unsigned extNum = extNumMin; extNum <= extNumMax; extNum++) {
        __global uint32_t *extCunningham1 = cunningham1Bitfield + extNum*sieveWords + lowIdx/32;
        __global uint32_t *extCunningham2 = cunningham2Bitfield + extNum*sieveWords + lowIdx/32;
        __global uint32_t *extBitwin = bitwinBitfield + extNum*sieveWords + lowIdx/32;
        flushLayer(layer - extNum, extCunningham1, extCunningham2, extBitwin, CR1, CR2, threadId);
      }
    }
  }
}           

                     
__kernel void sieveBenchmark(__global uint32_t *fixedMultipliers,
                             __global uint32_t *cunningham1Bitfield,
                             __global uint32_t *cunningham2Bitfield,
                             __global uint32_t *bitwinBitfield,
                             __constant uint32_t *primes,
                             __global uint64_t *multipliers64,
                             __global uint32_t *offsets64,                             
                             unsigned roundsNum)
{
  __local uint8_t localCunningham1[L1CacheSize];
  __local uint8_t localCunningham2[L1CacheSize];

  __global uint4 *M = (__global uint4*)(fixedMultipliers + get_group_id(0) * GPUMultiprecisionLimbs);
  weave(M[0], M[1], M[2],
        cunningham1Bitfield + get_group_id(0) * MaxSieveBufferSize/32,
        cunningham2Bitfield + get_group_id(0) * MaxSieveBufferSize/32,
        bitwinBitfield + get_group_id(0) * MaxSieveBufferSize/32,
        localCunningham1,
        localCunningham2,
        primes,
        multipliers64,
        offsets64,
        roundsNum);
        
  return;
}

uint32_t add128(uint4 *A, uint4 B)
{
  *A += B; 
  uint4 carry = -convert_uint4((*A) < B);
  
  (*A).y += carry.x; carry.y += ((*A).y < carry.x);
  (*A).z += carry.y; carry.z += ((*A).z < carry.y);
  (*A).w += carry.z;
  return carry.w + ((*A).w < carry.z); 
}

uint32_t add128Carry(uint4 *A, uint4 B, uint32_t externalCarry)
{
  *A += B;
  uint4 carry = -convert_uint4((*A) < B);
 
  (*A).x += externalCarry; carry.x += ((*A).x < externalCarry);
  (*A).y += carry.x; carry.y += ((*A).y < carry.x);
  (*A).z += carry.y; carry.z += ((*A).z < carry.y);
  (*A).w += carry.z;
  return carry.w + ((*A).w < carry.z); 
}

uint32_t add256(uint4 *a0, uint4 *a1, uint4 b0, uint4 b1)
{
  return add128Carry(a1, b1, add128(a0, b0));
}

uint32_t add384(uint4 *a0, uint4 *a1, uint4 *a2, uint4 b0, uint4 b1, uint4 b2)
{
  return add128Carry(a2, b2, add128Carry(a1, b1, add128(a0, b0)));
}

uint32_t add512(uint4 *a0, uint4 *a1, uint4 *a2, uint4 *a3, uint4 b0, uint4 b1, uint4 b2, uint4 b3)
{
  return add128Carry(a3, b3, add128Carry(a2, b2, add128Carry(a1, b1, add128(a0, b0))));
}

uint32_t sub64Borrow(uint2 *A, uint2 B, uint32_t externalBorrow)
{
  uint2 borrow = -convert_uint2((*A) < B);
  *A -= B;
  
  borrow.x += (*A).x < externalBorrow; (*A).x -= externalBorrow;
  borrow.y += (*A).y < borrow.x; (*A).y -= borrow.x;
  return borrow.y;
}

uint32_t sub128(uint4 *A, uint4 B)
{
  uint4 borrow = -convert_uint4((*A) < B);
  *A -= B;
 
  borrow.y += (*A).y < borrow.x; (*A).y -= borrow.x;
  borrow.z += (*A).z < borrow.y; (*A).z -= borrow.y;
  borrow.w += (*A).w < borrow.z; (*A).w -= borrow.z;
  return borrow.w;
}

uint32_t sub128Borrow(uint4 *A, uint4 B, uint32_t externalBorrow)
{
  uint4 borrow = -convert_uint4((*A) < B);
  *A -= B;
  
  borrow.x += (*A).x < externalBorrow; (*A).x -= externalBorrow;
  borrow.y += (*A).y < borrow.x; (*A).y -= borrow.x;
  borrow.z += (*A).z < borrow.y; (*A).z -= borrow.y;
  borrow.w += (*A).w < borrow.z; (*A).w -= borrow.z;
  return borrow.w;
}

uint32_t sub256(uint4 *a0, uint4 *a1, uint4 b0, uint4 b1)
{
  return sub128Borrow(a1, b1, sub128(a0, b0));
}

uint32_t sub384(uint4 *a0, uint4 *a1, uint4 *a2, uint4 b0, uint4 b1, uint4 b2)
{
  return sub128Borrow(a2, b2, sub128Borrow(a1, b1, sub128(a0, b0)));
}

uint32_t sub448(uint4 *a0, uint4 *a1, uint4 *a2, uint2 *a3, uint4 b0, uint4 b1, uint4 b2, uint2 b3)
{
  return sub64Borrow(a3, b3, sub128Borrow(a2, b2, sub128Borrow(a1, b1, sub128(a0, b0))));
}

void mul128round(uint4 op1l0, uint32_t m1, uint32_t m2,
                 uint64_t *R0, uint64_t *R1, uint64_t *R2, uint64_t *R3)
{
  uint4 l0 = mul_hi(op1l0, m1);
  uint4 l1 = op1l0 * m2;
  
  *R0 += l0.x; *R0 += l1.x;
  *R1 += l0.y; *R1 += l1.y; *R1 += (*R0 >> 32);
  *R2 += l0.z; *R2 += l1.z;
  *R3 += l0.w; *R3 += l1.w;
}

void mul128schoolBook_v3(uint4 op1l0, uint4 op2l0, uint4 *rl0, uint4 *rl1)
{
#define b1 op2l0.w
#define b2 op2l0.z
#define b3 op2l0.y
#define b4 op2l0.x

  ulong R0x = op1l0.x * b4;
  ulong R0y = op1l0.y * b4;
  ulong R0z = op1l0.z * b4;
  ulong R0w = op1l0.w * b4;
  ulong R1x = mul_hi(op1l0.x, b1);
  ulong R1y = mul_hi(op1l0.y, b1);
  ulong R1z = mul_hi(op1l0.z, b1);
  ulong R1w = mul_hi(op1l0.w, b1);
  
  mul128round(op1l0, b4, b3, &R0y, &R0z, &R0w, &R1x);
  mul128round(op1l0, b3, b2, &R0z, &R0w, &R1x, &R1y);
  mul128round(op1l0, b2, b1, &R0w, &R1x, &R1y, &R1z);
  R1y += (R1x >> 32);
  R1z += (R1y >> 32);
  R1w += (R1z >> 32);
  *rl0 = (uint4){R0x, R0y, R0z, R0w};
  *rl1 = (uint4){R1x, R1y, R1z, R1w};
#undef b1
#undef b2
#undef b3
#undef b4
}


void mul256round_v3(uint4 op1l0, uint4 op1l1, uint32_t m1, uint32_t m2,
                    uint64_t *R0, uint64_t *R1, uint64_t *R2, uint64_t *R3,
                    uint64_t *R4, uint64_t *R5, uint64_t *R6, uint64_t *R7)
{
  uint4 m1l0 = mul_hi(op1l0, m1);
  uint4 m1l1 = mul_hi(op1l1, m1);
  uint4 m2l0 = op1l0 * m2;
  uint4 m2l1 = op1l1 * m2;
  
  union {
    uint2 v32;
    ulong v64;
  } Int;  
  *R0 += m1l0.x; *R0 += m2l0.x;
  *R1 += m1l0.y; *R1 += m2l0.y; Int.v64 = *R0; *R1 += Int.v32.y;
  *R2 += m1l0.z; *R2 += m2l0.z;
  *R3 += m1l0.w; *R3 += m2l0.w;
  *R4 += m1l1.x; *R4 += m2l1.x;
  *R5 += m1l1.y; *R5 += m2l1.y;
  *R6 += m1l1.z; *R6 += m2l1.z;
  *R7 += m1l1.w; *R7 += m2l1.w;  
}


void mul384round_v3(uint4 op1l0, uint4 op1l1, uint4 op1l2, uint32_t m1, uint32_t m2,
                    uint64_t *R0, uint64_t *R1, uint64_t *R2, uint64_t *R3,
                    uint64_t *R4, uint64_t *R5, uint64_t *R6, uint64_t *R7,
                    uint64_t *R8, uint64_t *R9, uint64_t *R10, uint64_t *R11)
{
  uint4 m1l0 = mul_hi(op1l0, m1);
  uint4 m1l1 = mul_hi(op1l1, m1);
  uint4 m1l2 = mul_hi(op1l2, m1);
  uint4 m2l0 = op1l0 * m2;
  uint4 m2l1 = op1l1 * m2;
  uint4 m2l2 = op1l2 * m2;
  
  union {
    uint2 v32;
    ulong v64;
  } Int;
  *R0 += m1l0.x; *R0 += m2l0.x;
  *R1 += m1l0.y; *R1 += m2l0.y; Int.v64 = *R0; *R1 += Int.v32.y;
  *R2 += m1l0.z; *R2 += m2l0.z;
  *R3 += m1l0.w; *R3 += m2l0.w;
  *R4 += m1l1.x; *R4 += m2l1.x;
  *R5 += m1l1.y; *R5 += m2l1.y;
  *R6 += m1l1.z; *R6 += m2l1.z;
  *R7 += m1l1.w; *R7 += m2l1.w;  
  *R8 += m1l2.x; *R8 += m2l2.x;
  *R9 += m1l2.y; *R9 += m2l2.y;
  *R10 += m1l2.z; *R10 += m2l2.z;
  *R11 += m1l2.w; *R11 += m2l2.w;
}


void mul448round_v3(uint4 op1l0, uint4 op1l1, uint4 op1l2, uint2 op1l3, uint32_t m1, uint32_t m2,
                    uint64_t *R0, uint64_t *R1, uint64_t *R2, uint64_t *R3,
                    uint64_t *R4, uint64_t *R5, uint64_t *R6, uint64_t *R7,
                    uint64_t *R8, uint64_t *R9, uint64_t *R10, uint64_t *R11,
                    uint64_t *R12, uint64_t *R13)
{
  uint4 m1l0 = mul_hi(op1l0, m1);
  uint4 m1l1 = mul_hi(op1l1, m1);
  uint4 m1l2 = mul_hi(op1l2, m1);
  uint2 m1l3 = mul_hi(op1l3, m1);
  uint4 m2l0 = op1l0 * m2;
  uint4 m2l1 = op1l1 * m2;
  uint4 m2l2 = op1l2 * m2;
  uint2 m2l3 = op1l3 * m2;
  
  union {
    uint2 v32;
    ulong v64;
  } Int;
  *R0 += m1l0.x; *R0 += m2l0.x;
  *R1 += m1l0.y; *R1 += m2l0.y; Int.v64 = *R0; *R1 += Int.v32.y;
  *R2 += m1l0.z; *R2 += m2l0.z;
  *R3 += m1l0.w; *R3 += m2l0.w;
  *R4 += m1l1.x; *R4 += m2l1.x;
  *R5 += m1l1.y; *R5 += m2l1.y;
  *R6 += m1l1.z; *R6 += m2l1.z;
  *R7 += m1l1.w; *R7 += m2l1.w;  
  *R8 += m1l2.x; *R8 += m2l2.x;
  *R9 += m1l2.y; *R9 += m2l2.y;
  *R10 += m1l2.z; *R10 += m2l2.z;
  *R11 += m1l2.w; *R11 += m2l2.w;
  *R12 += m1l3.x; *R12 += m2l3.x;
  *R13 += m1l3.y; *R13 += m2l3.y;
}


void mul256schoolBook_v3(uint4 op1l0, uint4 op1l1,
                         uint4 op2l0, uint4 op2l1,
                         uint4 *rl0, uint4 *rl1, uint4 *rl2, uint4 *rl3)
{
#define b1 op2l1.w
#define b2 op2l1.z
#define b3 op2l1.y
#define b4 op2l1.x
#define b5 op2l0.w
#define b6 op2l0.z
#define b7 op2l0.y
#define b8 op2l0.x
  ulong R0x = op1l0.x * b8;
  ulong R0y = op1l0.y * b8;
  ulong R0z = op1l0.z * b8;
  ulong R0w = op1l0.w * b8;
  ulong R1x = op1l1.x * b8;
  ulong R1y = op1l1.y * b8;
  ulong R1z = op1l1.z * b8;
  ulong R1w = op1l1.w * b8;  
  ulong R2x = mul_hi(op1l0.x, b1);
  ulong R2y = mul_hi(op1l0.y, b1);
  ulong R2z = mul_hi(op1l0.z, b1);
  ulong R2w = mul_hi(op1l0.w, b1);
  ulong R3x = mul_hi(op1l1.x, b1);
  ulong R3y = mul_hi(op1l1.y, b1);
  ulong R3z = mul_hi(op1l1.z, b1);
  ulong R3w = mul_hi(op1l1.w, b1);  
  
  mul256round_v3(op1l0, op1l1, b8, b7, &R0y, &R0z, &R0w, &R1x, &R1y, &R1z, &R1w, &R2x); 
  mul256round_v3(op1l0, op1l1, b7, b6, &R0z, &R0w, &R1x, &R1y, &R1z, &R1w, &R2x, &R2y); 
  mul256round_v3(op1l0, op1l1, b6, b5, &R0w, &R1x, &R1y, &R1z, &R1w, &R2x, &R2y, &R2z); 
  mul256round_v3(op1l0, op1l1, b5, b4, &R1x, &R1y, &R1z, &R1w, &R2x, &R2y, &R2z, &R2w); 
  mul256round_v3(op1l0, op1l1, b4, b3, &R1y, &R1z, &R1w, &R2x, &R2y, &R2z, &R2w, &R3x); 
  mul256round_v3(op1l0, op1l1, b3, b2, &R1z, &R1w, &R2x, &R2y, &R2z, &R2w, &R3x, &R3y); 
  mul256round_v3(op1l0, op1l1, b2, b1, &R1w, &R2x, &R2y, &R2z, &R2w, &R3x, &R3y, &R3z);   
  
  union {
    uint2 v32;
    ulong v64;
  } Int;

  Int.v64 = R2x; R2y += Int.v32.y;
  Int.v64 = R2y; R2z += Int.v32.y;
  Int.v64 = R2z; R2w += Int.v32.y;
  Int.v64 = R2w; R3x += Int.v32.y;
  Int.v64 = R3x; R3y += Int.v32.y;
  Int.v64 = R3y; R3z += Int.v32.y;
  Int.v64 = R3z; R3w += Int.v32.y;
  
  *rl0 = (uint4){R0x, R0y, R0z, R0w};
  *rl1 = (uint4){R1x, R1y, R1z, R1w};
  *rl2 = (uint4){R2x, R2y, R2z, R2w};
  *rl3 = (uint4){R3x, R3y, R3z, R3w};
#undef b1
#undef b2
#undef b3
#undef b4
#undef b5
#undef b6
#undef b7
#undef b8  
}


void mul384schoolBook_v3(uint4 op1l0, uint4 op1l1, uint4 op1l2, // low --> hi
                         uint4 op2l0, uint4 op2l1, uint4 op2l2, // low --> hi
                         uint4 *rl0, uint4 *rl1, uint4 *rl2, uint4 *rl3, uint4 *rl4, uint4 *rl5)
{
#define b1 op2l2.w
#define b2 op2l2.z
#define b3 op2l2.y
#define b4 op2l2.x
#define b5 op2l1.w
#define b6 op2l1.z
#define b7 op2l1.y
#define b8 op2l1.x
#define b9 op2l0.w
#define b10 op2l0.z
#define b11 op2l0.y
#define b12 op2l0.x  
  ulong R0x = op1l0.x * b12;
  ulong R0y = op1l0.y * b12;
  ulong R0z = op1l0.z * b12;
  ulong R0w = op1l0.w * b12;
  ulong R1x = op1l1.x * b12;
  ulong R1y = op1l1.y * b12;
  ulong R1z = op1l1.z * b12;
  ulong R1w = op1l1.w * b12;  
  ulong R2x = op1l2.x * b12;
  ulong R2y = op1l2.y * b12;
  ulong R2z = op1l2.z * b12;
  ulong R2w = op1l2.w * b12;
  ulong R3x = mul_hi(op1l0.x, b1);
  ulong R3y = mul_hi(op1l0.y, b1);
  ulong R3z = mul_hi(op1l0.z, b1);
  ulong R3w = mul_hi(op1l0.w, b1);
  ulong R4x = mul_hi(op1l1.x, b1);
  ulong R4y = mul_hi(op1l1.y, b1);
  ulong R4z = mul_hi(op1l1.z, b1);
  ulong R4w = mul_hi(op1l1.w, b1);
  ulong R5x = mul_hi(op1l2.x, b1);
  ulong R5y = mul_hi(op1l2.y, b1);
  ulong R5z = mul_hi(op1l2.z, b1);
  ulong R5w = mul_hi(op1l2.w, b1);
  
  mul384round_v3(op1l0, op1l1, op1l2, b12, b11, &R0y, &R0z, &R0w, &R1x, &R1y, &R1z, &R1w, &R2x, &R2y, &R2z, &R2w, &R3x);
  mul384round_v3(op1l0, op1l1, op1l2, b11, b10, &R0z, &R0w, &R1x, &R1y, &R1z, &R1w, &R2x, &R2y, &R2z, &R2w, &R3x, &R3y);
  mul384round_v3(op1l0, op1l1, op1l2, b10,  b9, &R0w, &R1x, &R1y, &R1z, &R1w, &R2x, &R2y, &R2z, &R2w, &R3x, &R3y, &R3z);
  mul384round_v3(op1l0, op1l1, op1l2,  b9,  b8, &R1x, &R1y, &R1z, &R1w, &R2x, &R2y, &R2z, &R2w, &R3x, &R3y, &R3z, &R3w);
  mul384round_v3(op1l0, op1l1, op1l2,  b8,  b7, &R1y, &R1z, &R1w, &R2x, &R2y, &R2z, &R2w, &R3x, &R3y, &R3z, &R3w, &R4x);
  mul384round_v3(op1l0, op1l1, op1l2,  b7,  b6, &R1z, &R1w, &R2x, &R2y, &R2z, &R2w, &R3x, &R3y, &R3z, &R3w, &R4x, &R4y);
  mul384round_v3(op1l0, op1l1, op1l2,  b6,  b5, &R1w, &R2x, &R2y, &R2z, &R2w, &R3x, &R3y, &R3z, &R3w, &R4x, &R4y, &R4z);
  mul384round_v3(op1l0, op1l1, op1l2,  b5,  b4, &R2x, &R2y, &R2z, &R2w, &R3x, &R3y, &R3z, &R3w, &R4x, &R4y, &R4z, &R4w);
  mul384round_v3(op1l0, op1l1, op1l2,  b4,  b3, &R2y, &R2z, &R2w, &R3x, &R3y, &R3z, &R3w, &R4x, &R4y, &R4z, &R4w, &R5x);
  mul384round_v3(op1l0, op1l1, op1l2,  b3,  b2, &R2z, &R2w, &R3x, &R3y, &R3z, &R3w, &R4x, &R4y, &R4z, &R4w, &R5x, &R5y);
  mul384round_v3(op1l0, op1l1, op1l2,  b2,  b1, &R2w, &R3x, &R3y, &R3z, &R3w, &R4x, &R4y, &R4z, &R4w, &R5x, &R5y, &R5z);

  union {
    uint2 v32;
    ulong v64;
  } Int;

  Int.v64 = R3x; R3y += Int.v32.y;
  Int.v64 = R3y; R3z += Int.v32.y;
  Int.v64 = R3z; R3w += Int.v32.y;
  Int.v64 = R3w; R4x += Int.v32.y;
  Int.v64 = R4x; R4y += Int.v32.y;
  Int.v64 = R4y; R4z += Int.v32.y;
  Int.v64 = R4z; R4w += Int.v32.y;
  Int.v64 = R4w; R5x += Int.v32.y;
  Int.v64 = R5x; R5y += Int.v32.y;
  Int.v64 = R5y; R5z += Int.v32.y;
  Int.v64 = R5z; R5w += Int.v32.y;
  
  *rl0 = (uint4){R0x, R0y, R0z, R0w};
  *rl1 = (uint4){R1x, R1y, R1z, R1w};
  *rl2 = (uint4){R2x, R2y, R2z, R2w};
  *rl3 = (uint4){R3x, R3y, R3z, R3w};
  *rl4 = (uint4){R4x, R4y, R4z, R4w};
  *rl5 = (uint4){R5x, R5y, R5z, R5w};
  
#undef b1
#undef b2
#undef b3
#undef b4
#undef b5
#undef b6
#undef b7
#undef b8
#undef b9
#undef b10
#undef b11
#undef b12
}

void mul448schoolBook_v3(uint4 op1l0, uint4 op1l1, uint4 op1l2, uint2 op1l3, // low --> hi
                         uint4 op2l0, uint4 op2l1, uint4 op2l2, uint2 op2l3, // low --> hi
                         uint4 *rl0, uint4 *rl1, uint4 *rl2, uint4 *rl3, uint4 *rl4, uint4 *rl5, uint4 *rl6)
{
#define b1 op2l3.y
#define b2 op2l3.x
#define b3 op2l2.w
#define b4 op2l2.z
#define b5 op2l2.y
#define b6 op2l2.x
#define b7 op2l1.w
#define b8 op2l1.z
#define b9 op2l1.y
#define b10 op2l1.x
#define b11 op2l0.w
#define b12 op2l0.z
#define b13 op2l0.y
#define b14 op2l0.x  
  
  ulong R0x = op1l0.x * b14;
  ulong R0y = op1l0.y * b14;
  ulong R0z = op1l0.z * b14;
  ulong R0w = op1l0.w * b14;
  ulong R1x = op1l1.x * b14;
  ulong R1y = op1l1.y * b14;
  ulong R1z = op1l1.z * b14;
  ulong R1w = op1l1.w * b14;  
  ulong R2x = op1l2.x * b14;
  ulong R2y = op1l2.y * b14;
  ulong R2z = op1l2.z * b14;
  ulong R2w = op1l2.w * b14;
  ulong R3x = op1l3.x * b14;
  ulong R3y = op1l3.y * b14;
  ulong R3z = mul_hi(op1l0.x, b1);
  ulong R3w = mul_hi(op1l0.y, b1);
  ulong R4x = mul_hi(op1l0.z, b1);
  ulong R4y = mul_hi(op1l0.w, b1);
  ulong R4z = mul_hi(op1l1.x, b1);
  ulong R4w = mul_hi(op1l1.y, b1);
  ulong R5x = mul_hi(op1l1.z, b1);
  ulong R5y = mul_hi(op1l1.w, b1);
  ulong R5z = mul_hi(op1l2.x, b1);
  ulong R5w = mul_hi(op1l2.y, b1);
  ulong R6x = mul_hi(op1l2.z, b1);
  ulong R6y = mul_hi(op1l2.w, b1);
  ulong R6z = mul_hi(op1l3.x, b1);
  ulong R6w = mul_hi(op1l3.y, b1);
  
  mul448round_v3(op1l0, op1l1, op1l2, op1l3, b14, b13, &R0y, &R0z, &R0w, &R1x, &R1y, &R1z, &R1w, &R2x, &R2y, &R2z, &R2w, &R3x, &R3y, &R3z);
  mul448round_v3(op1l0, op1l1, op1l2, op1l3, b13, b12, &R0z, &R0w, &R1x, &R1y, &R1z, &R1w, &R2x, &R2y, &R2z, &R2w, &R3x, &R3y, &R3z, &R3w);
  mul448round_v3(op1l0, op1l1, op1l2, op1l3, b12, b11, &R0w, &R1x, &R1y, &R1z, &R1w, &R2x, &R2y, &R2z, &R2w, &R3x, &R3y, &R3z, &R3w, &R4x);
  mul448round_v3(op1l0, op1l1, op1l2, op1l3, b11, b10, &R1x, &R1y, &R1z, &R1w, &R2x, &R2y, &R2z, &R2w, &R3x, &R3y, &R3z, &R3w, &R4x, &R4y);
  mul448round_v3(op1l0, op1l1, op1l2, op1l3, b10,  b9, &R1y, &R1z, &R1w, &R2x, &R2y, &R2z, &R2w, &R3x, &R3y, &R3z, &R3w, &R4x, &R4y, &R4z);
  mul448round_v3(op1l0, op1l1, op1l2, op1l3,  b9,  b8, &R1z, &R1w, &R2x, &R2y, &R2z, &R2w, &R3x, &R3y, &R3z, &R3w, &R4x, &R4y, &R4z, &R4w);
  mul448round_v3(op1l0, op1l1, op1l2, op1l3,  b8,  b7, &R1w, &R2x, &R2y, &R2z, &R2w, &R3x, &R3y, &R3z, &R3w, &R4x, &R4y, &R4z, &R4w, &R5x);
  mul448round_v3(op1l0, op1l1, op1l2, op1l3,  b7,  b6, &R2x, &R2y, &R2z, &R2w, &R3x, &R3y, &R3z, &R3w, &R4x, &R4y, &R4z, &R4w, &R5x, &R5y);
  mul448round_v3(op1l0, op1l1, op1l2, op1l3,  b6,  b5, &R2y, &R2z, &R2w, &R3x, &R3y, &R3z, &R3w, &R4x, &R4y, &R4z, &R4w, &R5x, &R5y, &R5z);
  mul448round_v3(op1l0, op1l1, op1l2, op1l3,  b5,  b4, &R2z, &R2w, &R3x, &R3y, &R3z, &R3w, &R4x, &R4y, &R4z, &R4w, &R5x, &R5y, &R5z, &R5w);
  mul448round_v3(op1l0, op1l1, op1l2, op1l3,  b4,  b3, &R2w, &R3x, &R3y, &R3z, &R3w, &R4x, &R4y, &R4z, &R4w, &R5x, &R5y, &R5z, &R5w, &R6x);
  mul448round_v3(op1l0, op1l1, op1l2, op1l3,  b3,  b2, &R3x, &R3y, &R3z, &R3w, &R4x, &R4y, &R4z, &R4w, &R5x, &R5y, &R5z, &R5w, &R6x, &R6y);  
  mul448round_v3(op1l0, op1l1, op1l2, op1l3,  b2,  b1, &R3y, &R3z, &R3w, &R4x, &R4y, &R4z, &R4w, &R5x, &R5y, &R5z, &R5w, &R6x, &R6y, &R6z);
  
  union {
    uint2 v32;
    ulong v64;
  } Int;
  
  Int.v64 = R3z; R3w += Int.v32.y;
  Int.v64 = R3w; R4x += Int.v32.y;
  Int.v64 = R4x; R4y += Int.v32.y;
  Int.v64 = R4y; R4z += Int.v32.y;
  Int.v64 = R4z; R4w += Int.v32.y;
  Int.v64 = R4w; R5x += Int.v32.y;
  Int.v64 = R5x; R5y += Int.v32.y;
  Int.v64 = R5y; R5z += Int.v32.y;
  Int.v64 = R5z; R5w += Int.v32.y;
  Int.v64 = R5w; R6x += Int.v32.y;
  Int.v64 = R6x; R6y += Int.v32.y;
  Int.v64 = R6y; R6z += Int.v32.y;
  Int.v64 = R6z; R6w += Int.v32.y;  
  
  *rl0 = (uint4){R0x, R0y, R0z, R0w};
  *rl1 = (uint4){R1x, R1y, R1z, R1w};
  *rl2 = (uint4){R2x, R2y, R2z, R2w};
  *rl3 = (uint4){R3x, R3y, R3z, R3w};
  *rl4 = (uint4){R4x, R4y, R4z, R4w};
  *rl5 = (uint4){R5x, R5y, R5z, R5w};
  *rl6 = (uint4){R6x, R6y, R6z, R6w};
  
#undef b1
#undef b2
#undef b3
#undef b4
#undef b5
#undef b6
#undef b7
#undef b8
#undef b9
#undef b10
#undef b11
#undef b12
#undef b13
#undef b14
}


void mul512schoolBook_v3(uint4 op1l0, uint4 op1l1, uint4 op1l2, uint4 op1l3, // low --> hi
                         uint4 op2l0, uint4 op2l1, uint4 op2l2, uint4 op2l3, // low --> hi
                         uint4 *rl0, uint4 *rl1, uint4 *rl2, uint4 *rl3, uint4 *rl4, uint4 *rl5, uint4 *rl6, uint4 *rl7)
{
}

__kernel void multiplyBenchmark128(__global uint32_t *m1,
                                   __global uint32_t *m2,
                                   __global uint32_t *out,
                                   unsigned elementsNum)
{
  __global uint4 *M1 = (__global uint4*)m1;
  __global uint4 *M2 = (__global uint4*)m2;
  __global uint4 *OUT = (__global uint4*)out;
  
  unsigned globalSize = get_global_size(0);
  for (unsigned i = get_global_id(0); i < elementsNum; i += globalSize) {
    uint4 op0l0 = M1[i];
    uint4 op1l0 = M2[i];
    uint4 ResultLimb1;
    uint4 ResultLimb2;
    
    for (unsigned repeatNum = 0; repeatNum < 512; repeatNum++) {
      mul128schoolBook_v3(op0l0, op1l0, &ResultLimb1, &ResultLimb2);
      op0l0 = ResultLimb1;
    }
    
    OUT[i*2] = ResultLimb1;
    OUT[i*2+1] = ResultLimb2;
  }
}

__kernel void multiplyBenchmark256(__global uint32_t *m1,
                                   __global uint32_t *m2,
                                   __global uint32_t *out,
                                   unsigned elementsNum)
{
  __global uint4 *M1 = (__global uint4*)m1;
  __global uint4 *M2 = (__global uint4*)m2;
  __global uint4 *OUT = (__global uint4*)out;
  
  unsigned globalSize = get_global_size(0);
  for (unsigned i = get_global_id(0); i < elementsNum; i += globalSize) {
    uint4 op1l0 = M1[i*2];
    uint4 op1l1 = M1[i*2+1];
    uint4 op2l0 = M2[i*2];
    uint4 op2l1 = M2[i*2+1];
    uint4 ResultLimb1;
    uint4 ResultLimb2;
    uint4 ResultLimb3;
    uint4 ResultLimb4;    
    
    for (unsigned repeatNum = 0; repeatNum < 512; repeatNum++) {
      mul256schoolBook_v3(op1l0, op1l1, op2l0, op2l1,
                          &ResultLimb1, &ResultLimb2, &ResultLimb3, &ResultLimb4);
      op1l0 = ResultLimb1;
      op1l1 = ResultLimb2;
    }
    
    OUT[i*4] = ResultLimb1;
    OUT[i*4+1] = ResultLimb2;
    OUT[i*4+2] = ResultLimb3;
    OUT[i*4+3] = ResultLimb4;    
  }
}

__kernel void multiplyBenchmark384(__global uint32_t *m1,
                                   __global uint32_t *m2,
                                   __global uint32_t *out,
                                   unsigned elementsNum)
{
  __global uint4 *M1 = (__global uint4*)m1;
  __global uint4 *M2 = (__global uint4*)m2;
  __global uint4 *OUT = (__global uint4*)out;
  
  unsigned globalSize = get_global_size(0);
  for (unsigned i = get_global_id(0); i < elementsNum; i += globalSize) {
    uint4 op1l0 = M1[i*3];
    uint4 op1l1 = M1[i*3+1];
    uint4 op1l2 = M1[i*3+2];
    uint4 op2l0 = M2[i*3];
    uint4 op2l1 = M2[i*3+1];
    uint4 op2l2 = M2[i*3+2];
    uint4 ResultLimb1;
    uint4 ResultLimb2;
    uint4 ResultLimb3;
    uint4 ResultLimb4;    
    uint4 ResultLimb5;  
    uint4 ResultLimb6;      
    
    for (unsigned repeatNum = 0; repeatNum < 512; repeatNum++) {
      mul384schoolBook_v3(op1l0, op1l1, op1l2, op2l0, op2l1, op2l2,
                          &ResultLimb1, &ResultLimb2, &ResultLimb3, &ResultLimb4, &ResultLimb5, &ResultLimb6);
      op1l0 = ResultLimb1;
      op1l1 = ResultLimb2;
      op1l2 = ResultLimb3;
    }
    
    OUT[i*6] = ResultLimb1;
    OUT[i*6+1] = ResultLimb2;
    OUT[i*6+2] = ResultLimb3;
    OUT[i*6+3] = ResultLimb4;    
    OUT[i*6+4] = ResultLimb5;
    OUT[i*6+5] = ResultLimb6;        
  }
}


__kernel void multiplyBenchmark448(__global uint32_t *m1,
                                   __global uint32_t *m2,
                                   __global uint32_t *out,
                                   unsigned elementsNum)
{
  __global uint4 *M1 = (__global uint4*)m1;
  __global uint4 *M2 = (__global uint4*)m2;
  __global uint4 *OUT = (__global uint4*)out;
  
  unsigned globalSize = get_global_size(0);
  for (unsigned i = get_global_id(0); i < elementsNum; i += globalSize) {
    uint4 op1l0 = *(__global uint4*)(m1 + i*14); //M1[i*3];
    uint4 op1l1 = *(__global uint4*)(m1 + i*14 + 4);
    uint4 op1l2 = *(__global uint4*)(m1 + i*14 + 8);
    uint2 op1l3 = *(__global uint2*)(m1 + i*14 + 12);
    uint4 op2l0 = *(__global uint4*)(m2 + i*14);
    uint4 op2l1 = *(__global uint4*)(m2 + i*14 + 4);
    uint4 op2l2 = *(__global uint4*)(m2 + i*14 + 8);
    uint2 op2l3 = *(__global uint2*)(m2 + i*14 + 12);
    uint4 ResultLimb1;
    uint4 ResultLimb2;
    uint4 ResultLimb3;
    uint4 ResultLimb4;    
    uint4 ResultLimb5;  
    uint4 ResultLimb6;
    uint4 ResultLimb7;
    
    for (unsigned repeatNum = 0; repeatNum < 512; repeatNum++) {
      mul448schoolBook_v3(op1l0, op1l1, op1l2, op1l3, op2l0, op2l1, op2l2, op2l3,
                          &ResultLimb1, &ResultLimb2, &ResultLimb3, &ResultLimb4, &ResultLimb5, &ResultLimb6, &ResultLimb7);
      op1l0 = ResultLimb1;
      op1l1 = ResultLimb2;
      op1l2 = ResultLimb3;
      op1l3 = ResultLimb4.xy;
    }
    
    OUT[i*7] = ResultLimb1;
    OUT[i*7+1] = ResultLimb2;
    OUT[i*7+2] = ResultLimb3;
    OUT[i*7+3] = ResultLimb4;
    OUT[i*7+4] = ResultLimb5;
    OUT[i*7+5] = ResultLimb6;
    OUT[i*7+6] = ResultLimb7;
  }
}


__kernel void multiplyBenchmark512(__global uint32_t *m1,
                                   __global uint32_t *m2,
                                   __global uint32_t *out,
                                   unsigned elementsNum)
{
  __global uint4 *M1 = (__global uint4*)m1;
  __global uint4 *M2 = (__global uint4*)m2;
  __global uint4 *OUT = (__global uint4*)out;
  
  unsigned globalSize = get_global_size(0);
  for (unsigned i = get_global_id(0); i < elementsNum; i += globalSize) {
    uint4 op1l0 = M1[i*4];
    uint4 op1l1 = M1[i*4+1];
    uint4 op1l2 = M1[i*4+2];
    uint4 op1l3 = M1[i*4+3];
    uint4 op2l0 = M2[i*4];
    uint4 op2l1 = M2[i*4+1];
    uint4 op2l2 = M2[i*4+2];
    uint4 op2l3 = M2[i*4+3];
    uint4 ResultLimb1;
    uint4 ResultLimb2;
    uint4 ResultLimb3;
    uint4 ResultLimb4;    
    uint4 ResultLimb5;  
    uint4 ResultLimb6;      
    uint4 ResultLimb7;  
    uint4 ResultLimb8;     
    
    for (unsigned repeatNum = 0; repeatNum < 512; repeatNum++) {
      mul512schoolBook_v3(op1l0, op1l1, op1l2, op1l3, op2l0, op2l1, op2l2, op2l3,
                          &ResultLimb1, &ResultLimb2, &ResultLimb3, &ResultLimb4, &ResultLimb5, &ResultLimb6, &ResultLimb7, &ResultLimb8);
      op1l0 = ResultLimb1;
      op1l1 = ResultLimb2;
      op1l2 = ResultLimb3;
      op1l3 = ResultLimb4;
    }
    
    OUT[i*8] = ResultLimb1;
    OUT[i*8+1] = ResultLimb2;
    OUT[i*8+2] = ResultLimb3;
    OUT[i*8+3] = ResultLimb4;    
    OUT[i*8+4] = ResultLimb5;
    OUT[i*8+5] = ResultLimb6;        
    OUT[i*8+6] = ResultLimb7;
    OUT[i*8+7] = ResultLimb8;       
  }
}

void lshiftByLimb2(uint4 *limbs1,
                   uint4 *limbs2)
{
  (*limbs2).yzw = (*limbs2).xyz; (*limbs2).x = (*limbs1).w;
  (*limbs1).yzw = (*limbs1).xyz; (*limbs1).x = 0;
}

void rshiftByLimb2(uint4 *limbs1,
                   uint4 *limbs2)
{
  (*limbs1).xyz = (*limbs1).yzw; (*limbs1).w = (*limbs2).x;
  (*limbs2).xyz = (*limbs2).yzw; (*limbs2).w = 0;
}

void lshiftByLimb3(uint4 *limbs1,
                   uint4 *limbs2,
                   uint4 *limbs3)
{
  (*limbs3).yzw = (*limbs3).xyz; (*limbs3).x = (*limbs2).w;
  (*limbs2).yzw = (*limbs2).xyz; (*limbs2).x = (*limbs1).w;
  (*limbs1).yzw = (*limbs1).xyz; (*limbs1).x = 0;
}

void rshiftByLimb3(uint4 *limbs1,
                   uint4 *limbs2,
                   uint4 *limbs3)
{
  (*limbs1).xyz = (*limbs1).yzw; (*limbs1).w = (*limbs2).x;
  (*limbs2).xyz = (*limbs2).yzw; (*limbs2).w = (*limbs3).x;
  (*limbs3).xyz = (*limbs3).yzw; (*limbs3).w = 0;
}

void lshiftByLimb4(uint4 *limbs1,
                   uint4 *limbs2,
                   uint4 *limbs3,
                   uint4 *limbs4)
{
  (*limbs4).yzw = (*limbs4).xyz; (*limbs4).x = (*limbs3).w;  
  (*limbs3).yzw = (*limbs3).xyz; (*limbs3).x = (*limbs2).w;
  (*limbs2).yzw = (*limbs2).xyz; (*limbs2).x = (*limbs1).w;
  (*limbs1).yzw = (*limbs1).xyz; (*limbs1).x = 0;
}

void rshiftByLimb4(uint4 *limbs1,
                   uint4 *limbs2,
                   uint4 *limbs3,
                   uint4 *limbs4)
{
  (*limbs1).xyz = (*limbs1).yzw; (*limbs1).w = (*limbs2).x;
  (*limbs2).xyz = (*limbs2).yzw; (*limbs2).w = (*limbs3).x;
  (*limbs3).xyz = (*limbs3).yzw; (*limbs3).w = (*limbs4).x;
  (*limbs4).xyz = (*limbs4).yzw; (*limbs4).w = 0;
}

void lshiftByLimb5(uint4 *limbs1,
                   uint4 *limbs2,
                   uint4 *limbs3,
                   uint4 *limbs4,
                   uint4 *limbs5)
{
  (*limbs5).yzw = (*limbs5).xyz; (*limbs5).x = (*limbs4).w;
  (*limbs4).yzw = (*limbs4).xyz; (*limbs4).x = (*limbs3).w;  
  (*limbs3).yzw = (*limbs3).xyz; (*limbs3).x = (*limbs2).w;
  (*limbs2).yzw = (*limbs2).xyz; (*limbs2).x = (*limbs1).w;
  (*limbs1).yzw = (*limbs1).xyz; (*limbs1).x = 0;
}

void rshiftByLimb5(uint4 *limbs1,
                   uint4 *limbs2,
                   uint4 *limbs3,
                   uint4 *limbs4,
                   uint4 *limbs5)
{
  (*limbs1).xyz = (*limbs1).yzw; (*limbs1).w = (*limbs2).x;
  (*limbs2).xyz = (*limbs2).yzw; (*limbs2).w = (*limbs3).x;
  (*limbs3).xyz = (*limbs3).yzw; (*limbs3).w = (*limbs4).x;
  (*limbs4).xyz = (*limbs4).yzw; (*limbs4).w = (*limbs5).x;
  (*limbs5).xyz = (*limbs5).yzw; (*limbs5).w = 0;
}

void rshiftByLimb6(uint4 *limbs1,
                   uint4 *limbs2,
                   uint4 *limbs3,
                   uint4 *limbs4,
                   uint4 *limbs5,
                   uint4 *limbs6)
{
  (*limbs1).xyz = (*limbs1).yzw; (*limbs1).w = (*limbs2).x;
  (*limbs2).xyz = (*limbs2).yzw; (*limbs2).w = (*limbs3).x;
  (*limbs3).xyz = (*limbs3).yzw; (*limbs3).w = (*limbs4).x;
  (*limbs4).xyz = (*limbs4).yzw; (*limbs4).w = (*limbs5).x;
  (*limbs5).xyz = (*limbs5).yzw; (*limbs5).w = (*limbs6).x;
  (*limbs6).xyz = (*limbs6).yzw; (*limbs6).w = 0;
}

void lshift2(uint4 *limbs1, uint4 *limbs2, unsigned count)
{
  if (!count)
    return;  
  unsigned lowBitsCount = 32 - count;  
  
  {
    uint4 lowBits = {
      (*limbs1).w >> lowBitsCount,
      (*limbs2).x >> lowBitsCount,
      (*limbs2).y >> lowBitsCount,
      (*limbs2).z >> lowBitsCount
    };
    (*limbs2) = ((*limbs2) << count) | lowBits;
  }  

  {
    uint4 lowBits = {
      0,
      (*limbs1).x >> lowBitsCount,
      (*limbs1).y >> lowBitsCount,
      (*limbs1).z >> lowBitsCount
    };
    (*limbs1) = ((*limbs1) << count) | lowBits;
  }
}

void rshift2(uint4 *limbs1, uint4 *limbs2, unsigned count)
{
  if (!count)
    return;  
  unsigned lowBitsCount = 32 - count;  
  
  {
    uint4 lowBits = {
      (*limbs1).y << lowBitsCount,
      (*limbs1).z << lowBitsCount,
      (*limbs1).w << lowBitsCount,
      (*limbs2).x << lowBitsCount
    };
    (*limbs1) = ((*limbs1) >> count) | lowBits;
  }  
  
  {
    uint4 lowBits = {
      (*limbs2).y << lowBitsCount,
      (*limbs2).z << lowBitsCount,
      (*limbs2).w << lowBitsCount,
      0
    };
    (*limbs2) = ((*limbs2) >> count) | lowBits;
  }  
}

void lshift3(uint4 *limbs1, uint4 *limbs2, uint4 *limbs3, unsigned count)
{
  if (!count)
    return;  
  unsigned lowBitsCount = 32 - count;
  
  {
    uint4 lowBits = {
      (*limbs2).w >> lowBitsCount,
      (*limbs3).x >> lowBitsCount,
      (*limbs3).y >> lowBitsCount,
      (*limbs3).z >> lowBitsCount
    };
    (*limbs3) = ((*limbs3) << count) | lowBits;
  }    
  
  {
    uint4 lowBits = {
      (*limbs1).w >> lowBitsCount,
      (*limbs2).x >> lowBitsCount,
      (*limbs2).y >> lowBitsCount,
      (*limbs2).z >> lowBitsCount
    };
    (*limbs2) = ((*limbs2) << count) | lowBits;
  }  
  
  {
    uint4 lowBits = {
      0,
      (*limbs1).x >> lowBitsCount,
      (*limbs1).y >> lowBitsCount,
      (*limbs1).z >> lowBitsCount
    };
    (*limbs1) = ((*limbs1) << count) | lowBits;
  }
}

void rshift3(uint4 *limbs1, uint4 *limbs2, uint4 *limbs3, unsigned count)
{
  if (!count)
    return;
  unsigned lowBitsCount = 32 - count;  
  
  {
    uint4 lowBits = {
      (*limbs1).y << lowBitsCount,
      (*limbs1).z << lowBitsCount,
      (*limbs1).w << lowBitsCount,
      (*limbs2).x << lowBitsCount
    };
    (*limbs1) = ((*limbs1) >> count) | lowBits;
  }    
  
  {
    uint4 lowBits = {
      (*limbs2).y << lowBitsCount,
      (*limbs2).z << lowBitsCount,
      (*limbs2).w << lowBitsCount,
      (*limbs3).x << lowBitsCount
    };
    (*limbs2) = ((*limbs2) >> count) | lowBits;
  }  
  
  {
    uint4 lowBits = {
      (*limbs3).y << lowBitsCount,
      (*limbs3).z << lowBitsCount,
      (*limbs3).w << lowBitsCount,
      0
    };
    (*limbs3) = ((*limbs3) >> count) | lowBits;
  }  
}

void lshift4(uint4 *limbs1, uint4 *limbs2, uint4 *limbs3, uint4 *limbs4, unsigned count)
{
  if (!count)
    return;  
  unsigned lowBitsCount = 32 - count;
  
  {
    uint4 lowBits = {
      (*limbs3).w >> lowBitsCount,
      (*limbs4).x >> lowBitsCount,
      (*limbs4).y >> lowBitsCount,
      (*limbs4).z >> lowBitsCount
    };
    (*limbs4) = ((*limbs4) << count) | lowBits;
  }      
  
  {
    uint4 lowBits = {
      (*limbs2).w >> lowBitsCount,
      (*limbs3).x >> lowBitsCount,
      (*limbs3).y >> lowBitsCount,
      (*limbs3).z >> lowBitsCount
    };
    (*limbs3) = ((*limbs3) << count) | lowBits;
  }    
  
  {
    uint4 lowBits = {
      (*limbs1).w >> lowBitsCount,
      (*limbs2).x >> lowBitsCount,
      (*limbs2).y >> lowBitsCount,
      (*limbs2).z >> lowBitsCount
    };
    (*limbs2) = ((*limbs2) << count) | lowBits;
  }  
  
  {
    uint4 lowBits = {
      0,
      (*limbs1).x >> lowBitsCount,
      (*limbs1).y >> lowBitsCount,
      (*limbs1).z >> lowBitsCount
    };
    (*limbs1) = ((*limbs1) << count) | lowBits;
  }
}

void lshift5(uint4 *limbs1, uint4 *limbs2, uint4 *limbs3, uint4 *limbs4, uint4 *limbs5, unsigned count)
{
  if (!count)
    return;  
  unsigned lowBitsCount = 32 - count;

  {
    uint4 lowBits = {
      (*limbs4).w >> lowBitsCount,
      (*limbs5).x >> lowBitsCount,
      (*limbs5).y >> lowBitsCount,
      (*limbs5).z >> lowBitsCount
    };
    (*limbs5) = ((*limbs5) << count) | lowBits;
  }      
  
  {
    uint4 lowBits = {
      (*limbs3).w >> lowBitsCount,
      (*limbs4).x >> lowBitsCount,
      (*limbs4).y >> lowBitsCount,
      (*limbs4).z >> lowBitsCount
    };
    (*limbs4) = ((*limbs4) << count) | lowBits;
  }      
  
  {
    uint4 lowBits = {
      (*limbs2).w >> lowBitsCount,
      (*limbs3).x >> lowBitsCount,
      (*limbs3).y >> lowBitsCount,
      (*limbs3).z >> lowBitsCount
    };
    (*limbs3) = ((*limbs3) << count) | lowBits;
  }    
  
  {
    uint4 lowBits = {
      (*limbs1).w >> lowBitsCount,
      (*limbs2).x >> lowBitsCount,
      (*limbs2).y >> lowBitsCount,
      (*limbs2).z >> lowBitsCount
    };
    (*limbs2) = ((*limbs2) << count) | lowBits;
  }  
  
  {
    uint4 lowBits = {
      0,
      (*limbs1).x >> lowBitsCount,
      (*limbs1).y >> lowBitsCount,
      (*limbs1).z >> lowBitsCount
    };
    (*limbs1) = ((*limbs1) << count) | lowBits;
  }
}


void rshift4(uint4 *limbs1, uint4 *limbs2, uint4 *limbs3, uint4 *limbs4, unsigned count)
{
  if (!count)
    return;
  unsigned lowBitsCount = 32 - count;  
  
  {
    uint4 lowBits = {
      (*limbs1).y << lowBitsCount,
      (*limbs1).z << lowBitsCount,
      (*limbs1).w << lowBitsCount,
      (*limbs2).x << lowBitsCount
    };
    (*limbs1) = ((*limbs1) >> count) | lowBits;
  }    
  
  {
    uint4 lowBits = {
      (*limbs2).y << lowBitsCount,
      (*limbs2).z << lowBitsCount,
      (*limbs2).w << lowBitsCount,
      (*limbs3).x << lowBitsCount
    };
    (*limbs2) = ((*limbs2) >> count) | lowBits;
  }  
  
  {
    uint4 lowBits = {
      (*limbs3).y << lowBitsCount,
      (*limbs3).z << lowBitsCount,
      (*limbs3).w << lowBitsCount,
      (*limbs4).x << lowBitsCount,
    };
    (*limbs3) = ((*limbs3) >> count) | lowBits;
  }
  
  {
    uint4 lowBits = {
      (*limbs4).y << lowBitsCount,
      (*limbs4).z << lowBitsCount,
      (*limbs4).w << lowBitsCount,
      0
    };
    (*limbs4) = ((*limbs4) >> count) | lowBits;
  }    
}


void rshift5(uint4 *limbs1, uint4 *limbs2, uint4 *limbs3, uint4 *limbs4, uint4 *limbs5, unsigned count)
{
  if (!count)
    return;
  unsigned lowBitsCount = 32 - count;  
  
  {
    uint4 lowBits = {
      (*limbs1).y << lowBitsCount,
      (*limbs1).z << lowBitsCount,
      (*limbs1).w << lowBitsCount,
      (*limbs2).x << lowBitsCount
    };
    (*limbs1) = ((*limbs1) >> count) | lowBits;
  }    
  
  {
    uint4 lowBits = {
      (*limbs2).y << lowBitsCount,
      (*limbs2).z << lowBitsCount,
      (*limbs2).w << lowBitsCount,
      (*limbs3).x << lowBitsCount
    };
    (*limbs2) = ((*limbs2) >> count) | lowBits;
  }  
  
  {
    uint4 lowBits = {
      (*limbs3).y << lowBitsCount,
      (*limbs3).z << lowBitsCount,
      (*limbs3).w << lowBitsCount,
      (*limbs4).x << lowBitsCount,
    };
    (*limbs3) = ((*limbs3) >> count) | lowBits;
  }
  
  {
    uint4 lowBits = {
      (*limbs4).y << lowBitsCount,
      (*limbs4).z << lowBitsCount,
      (*limbs4).w << lowBitsCount,
      (*limbs5).x << lowBitsCount,
    };
    (*limbs4) = ((*limbs4) >> count) | lowBits;
  }
  
  {
    uint4 lowBits = {
      (*limbs5).y << lowBitsCount,
      (*limbs5).z << lowBitsCount,
      (*limbs5).w << lowBitsCount,
      0
    };
    (*limbs5) = ((*limbs5) >> count) | lowBits;
  }    
}


// Calculate (256bit){b1, b0} -= (256bit)({a1, a0} >> 32) * (32bit)M, returns highest product limb
void subMul1_v2(uint4 *b0, uint4 *b1, uint4 *b2,
                uint4 a0, uint4 a1,
                uint32_t M)
{
#define bringBorrow(Data, Borrow, NextBorrow) NextBorrow += (Data < Borrow); Data -= Borrow;
  
  uint4 Mv4 = {M, M, M, M};
  uint32_t clow;
  uint4 c1 = {0, 0, 0, 0};
  uint4 c2 = {0, 0, 0, 0};
  
  {
    uint4 a0M = a0*Mv4;
    uint4 a0Mhi = mul_hi(a0, Mv4);
  
    clow = (*b0).w < a0M.x;
    (*b0).w -= a0M.x;
  
    c1.xyz -= convert_uint3((*b1).xyz < a0M.yzw);
    (*b1).xyz -= a0M.yzw;

    c1 -= convert_uint4((*b1) < a0Mhi);
    (*b1) -= a0Mhi;
  }
  
  {
    uint4 a1M = a1*Mv4;
    uint4 a1Mhi = mul_hi(a1, Mv4);
    
    c1.w += ((*b1).w < a1M.x);
    (*b1).w -= a1M.x;

    c2.xyz -= convert_uint3((*b2).xyz < a1M.yzw);
    (*b2).xyz -= a1M.yzw;

    c2 -= convert_uint4((*b2) < a1Mhi);
    (*b2) -= a1Mhi;
  }
  
  bringBorrow((*b1).x, clow, c1.x);
  bringBorrow((*b1).y, c1.x, c1.y);
  bringBorrow((*b1).z, c1.y, c1.z);
  bringBorrow((*b1).w, c1.z, c1.w);
  bringBorrow((*b2).x, c1.w, c2.x);
  bringBorrow((*b2).y, c2.x, c2.y);
  bringBorrow((*b2).z, c2.y, c2.z);
  bringBorrow((*b2).w, c2.z, c2.w);
#undef bringBorrow
}

// // Calculate (384bit){b2, b1, b0} -= (384bit)({a2, a1, a0} >> 32) * (32bit)M, returns highest product limb
void subMul1_v3(uint4 *b0, uint4 *b1, uint4 *b2, uint4 *b3,
                uint4 a0, uint4 a1, uint4 a2,
                uint32_t M)
{
#define bringBorrow(Data, Borrow, NextBorrow) NextBorrow += (Data < Borrow); Data -= Borrow;
  
  uint4 Mv4 = {M, M, M, M};
  uint32_t clow;
  uint4 c1 = {0, 0, 0, 0};
  uint4 c2 = {0, 0, 0, 0};
  uint4 c3 = {0, 0, 0, 0};
  
  {
    uint4 a0M = a0*Mv4;
    uint4 a0Mhi = mul_hi(a0, Mv4);
    
    clow = (*b0).w < a0M.x;
    (*b0).w -= a0M.x;
    
    c1.xyz -= convert_uint3((*b1).xyz < a0M.yzw);
    (*b1).xyz -= a0M.yzw;
    
    c1 -= convert_uint4((*b1) < a0Mhi);
    (*b1) -= a0Mhi;
  }
  
  {
    uint4 a1M = a1*Mv4;
    uint4 a1Mhi = mul_hi(a1, Mv4);
    
    c1.w += ((*b1).w < a1M.x);
    (*b1).w -= a1M.x;
    
    c2.xyz -= convert_uint3((*b2).xyz < a1M.yzw);
    (*b2).xyz -= a1M.yzw;
    
    c2 -= convert_uint4((*b2) < a1Mhi);
    (*b2) -= a1Mhi;
  }
  
  {
    uint4 a2M = a2*Mv4;
    uint4 a2Mhi = mul_hi(a2, Mv4);
    
    c2.w += ((*b2).w < a2M.x);
    (*b2).w -= a2M.x;
    
    c3.xyz -= convert_uint3((*b3).xyz < a2M.yzw);
    (*b3).xyz -= a2M.yzw;
    c3 -= convert_uint4((*b3) < a2Mhi);
    (*b3) -= a2Mhi;
  }
  
  bringBorrow((*b1).x, clow, c1.x);
  bringBorrow((*b1).y, c1.x, c1.y);
  bringBorrow((*b1).z, c1.y, c1.z);
  bringBorrow((*b1).w, c1.z, c1.w);
  bringBorrow((*b2).x, c1.w, c2.x);
  bringBorrow((*b2).y, c2.x, c2.y);
  bringBorrow((*b2).z, c2.y, c2.z);
  bringBorrow((*b2).w, c2.z, c2.w);
  bringBorrow((*b3).x, c2.w, c3.x);
  bringBorrow((*b3).y, c3.x, c3.y);
  bringBorrow((*b3).z, c3.y, c3.z);
  bringBorrow((*b3).w, c3.z, c3.w);
#undef bringBorrow
}

// // Calculate (512bit){b3, b2, b1, b0} -= (512bit)({a3, a2, a1, a0} >> 32) * (32bit)M, returns highest product limb
void subMul1_512(uint4 *b0, uint4 *b1, uint4 *b2, uint4 *b3, uint4 *b4,
                uint4 a0, uint4 a1, uint4 a2, uint4 a3,
                uint32_t M)
{
#define bringBorrow(Data, Borrow, NextBorrow) NextBorrow += (Data < Borrow); Data -= Borrow;
  
  uint4 Mv4 = {M, M, M, M};
  uint32_t clow;
  uint4 c1 = {0, 0, 0, 0};
  uint4 c2 = {0, 0, 0, 0};
  uint4 c3 = {0, 0, 0, 0};
  uint4 c4 = {0, 0, 0, 0};
  
  {
    uint4 a0M = a0*Mv4;
    uint4 a0Mhi = mul_hi(a0, Mv4);
    
    clow = (*b0).w < a0M.x;
    (*b0).w -= a0M.x;
    
    c1.xyz -= convert_uint3((*b1).xyz < a0M.yzw);
    (*b1).xyz -= a0M.yzw;
    
    c1 -= convert_uint4((*b1) < a0Mhi);
    (*b1) -= a0Mhi;
  }
  
  {
    uint4 a1M = a1*Mv4;
    uint4 a1Mhi = mul_hi(a1, Mv4);
    
    c1.w += ((*b1).w < a1M.x);
    (*b1).w -= a1M.x;
    
    c2.xyz -= convert_uint3((*b2).xyz < a1M.yzw);
    (*b2).xyz -= a1M.yzw;
    
    c2 -= convert_uint4((*b2) < a1Mhi);
    (*b2) -= a1Mhi;
  }
  
  {
    uint4 a2M = a2*Mv4;
    uint4 a2Mhi = mul_hi(a2, Mv4);
    
    c2.w += ((*b2).w < a2M.x);
    (*b2).w -= a2M.x;
    
    c3.xyz -= convert_uint3((*b3).xyz < a2M.yzw);
    (*b3).xyz -= a2M.yzw;
    c3 -= convert_uint4((*b3) < a2Mhi);
    (*b3) -= a2Mhi;
  }
  
  {
    uint4 a3M = a3*Mv4;
    uint4 a3Mhi = mul_hi(a3, Mv4);
    
    c3.w += ((*b3).w < a3M.x);
    (*b3).w -= a3M.x;
    
    c4.xyz -= convert_uint3((*b4).xyz < a3M.yzw);
    (*b4).xyz -= a3M.yzw;
    c4 -= convert_uint4((*b4) < a3Mhi);
    (*b4) -= a3Mhi;
  }  
  
  bringBorrow((*b1).x, clow, c1.x);
  bringBorrow((*b1).y, c1.x, c1.y);
  bringBorrow((*b1).z, c1.y, c1.z);
  bringBorrow((*b1).w, c1.z, c1.w);
  bringBorrow((*b2).x, c1.w, c2.x);
  bringBorrow((*b2).y, c2.x, c2.y);
  bringBorrow((*b2).z, c2.y, c2.z);
  bringBorrow((*b2).w, c2.z, c2.w);
  bringBorrow((*b3).x, c2.w, c3.x);
  bringBorrow((*b3).y, c3.x, c3.y);
  bringBorrow((*b3).z, c3.y, c3.z);
  bringBorrow((*b3).w, c3.z, c3.w);
  bringBorrow((*b4).x, c3.w, c4.x);
  bringBorrow((*b4).y, c4.x, c4.y);
  bringBorrow((*b4).z, c4.y, c4.z);
  bringBorrow((*b4).w, c4.z, c4.w);  
#undef bringBorrow
}


uint2 modulo384to256(uint4 dividendLimbs0,
                     uint4 dividendLimbs1,
                     uint4 dividendLimbs2,
                     uint4 divisorLimbs0,
                     uint4 divisorLimbs1,
                     uint4 *moduloLimbs0,
                     uint4 *moduloLimbs1)
{
  // Detect dividend and divisor limbs count (remove trailing zero limbs)
  unsigned dividendLimbs = 12;
  unsigned divisorLimbs = 8;

  while (divisorLimbs && !divisorLimbs1.w) {
    lshiftByLimb2(&divisorLimbs0, &divisorLimbs1);
    divisorLimbs--;
  }
  
  // Normalize dividend and divisor (high bit of divisor must be set to 1)
  unsigned normalizeShiftCount = 0;  
  uint32_t bit = 0x80000000;
  while (!(divisorLimbs1.w & bit)) {
    normalizeShiftCount++;
    bit >>= 1;
  }
  
  lshift3(&dividendLimbs0, &dividendLimbs1, &dividendLimbs2, normalizeShiftCount);
  lshift2(&divisorLimbs0, &divisorLimbs1, normalizeShiftCount);
  
  
  while (dividendLimbs && !dividendLimbs2.w) {
    lshiftByLimb3(&dividendLimbs0, &dividendLimbs1, &dividendLimbs2);
    dividendLimbs--;
  }    

  for (unsigned i = 0; i < (dividendLimbs - divisorLimbs); i++) {

    uint32_t i32quotient;
    if (dividendLimbs2.w == divisorLimbs1.w) {
      i32quotient = 0xFFFFFFFF;
    } else {
      uint64_t i64dividend = (((uint64_t)dividendLimbs2.w) << 32) | dividendLimbs2.z;
      i32quotient = i64dividend / divisorLimbs1.w;
    }

    subMul1_v2(&dividendLimbs0, &dividendLimbs1, &dividendLimbs2,
               divisorLimbs0, divisorLimbs1,
               i32quotient);     
    uint32_t borrow = dividendLimbs2.w;    
    lshiftByLimb3(&dividendLimbs0, &dividendLimbs1, &dividendLimbs2);
    if (borrow) {
      add256(&dividendLimbs1, &dividendLimbs2, divisorLimbs0, divisorLimbs1);
      if (dividendLimbs2.w > divisorLimbs1.w)
        add256(&dividendLimbs1, &dividendLimbs2, divisorLimbs0, divisorLimbs1);
    }     
  }
 
  rshift3(&dividendLimbs0, &dividendLimbs1, &dividendLimbs2, normalizeShiftCount);
  for (unsigned i = 0; i < (8-divisorLimbs); i++)
    rshiftByLimb3(&dividendLimbs0, &dividendLimbs1, &dividendLimbs2);
  
  *moduloLimbs0 = dividendLimbs1;
  *moduloLimbs1 = dividendLimbs2;
  return (uint2){divisorLimbs, 32-normalizeShiftCount};
}


uint2 modulo512to384(uint4 dividendLimbs0,
                     uint4 dividendLimbs1,
                     uint4 dividendLimbs2,
                     uint4 dividendLimbs3,
                     uint4 divisorLimbs0,
                     uint4 divisorLimbs1,
                     uint4 divisorLimbs2,
                     uint4 *moduloLimbs0,
                     uint4 *moduloLimbs1,
                     uint4 *moduloLimbs2)
{
  // Detect dividend and divisor limbs count (remove trailing zero limbs)
  unsigned dividendLimbs = 16;
  unsigned divisorLimbs = 12;
  
  while (divisorLimbs && !divisorLimbs2.w) {
    lshiftByLimb3(&divisorLimbs0, &divisorLimbs1, &divisorLimbs2);
    divisorLimbs--;
  }  
  
  // Normalize dividend and divisor (high bit of divisor must be set to 1)
  unsigned normalizeShiftCount = 0;  
  uint32_t bit = 0x80000000;
  while (!(divisorLimbs2.w & bit)) {
    normalizeShiftCount++;
    bit >>= 1;  
  }    
  
  lshift4(&dividendLimbs0, &dividendLimbs1, &dividendLimbs2, &dividendLimbs3, normalizeShiftCount);
  lshift3(&divisorLimbs0, &divisorLimbs1, &divisorLimbs2, normalizeShiftCount);    

    
  while (dividendLimbs && !dividendLimbs3.w) {
    lshiftByLimb4(&dividendLimbs0, &dividendLimbs1, &dividendLimbs2, &dividendLimbs3);
    dividendLimbs--;
  }  

  for (unsigned i = 0; i < (dividendLimbs - divisorLimbs); i++) {
    uint32_t i32quotient;
    if (dividendLimbs3.w == divisorLimbs2.w) {
      i32quotient = 0xFFFFFFFF;
    } else {
      uint64_t i64dividend = (((uint64_t)dividendLimbs3.w) << 32) | dividendLimbs3.z;
      i32quotient = i64dividend / divisorLimbs2.w;
    }

    subMul1_v3(&dividendLimbs0, &dividendLimbs1, &dividendLimbs2, &dividendLimbs3,
               divisorLimbs0, divisorLimbs1, divisorLimbs2,
               i32quotient);    
    uint32_t borrow = dividendLimbs3.w;    
    lshiftByLimb4(&dividendLimbs0, &dividendLimbs1, &dividendLimbs2, &dividendLimbs3);
    if (borrow) {
      add384(&dividendLimbs1, &dividendLimbs2, &dividendLimbs3, divisorLimbs0, divisorLimbs1, divisorLimbs2);
      if (dividendLimbs3.w > divisorLimbs2.w)
        add384(&dividendLimbs1, &dividendLimbs2, &dividendLimbs3, divisorLimbs0, divisorLimbs1, divisorLimbs2);
    }
  }
 
  rshift4(&dividendLimbs0, &dividendLimbs1, &dividendLimbs2, &dividendLimbs3, normalizeShiftCount);
  for (unsigned i = 0; i < (12-divisorLimbs); i++)
    rshiftByLimb4(&dividendLimbs0, &dividendLimbs1, &dividendLimbs2, &dividendLimbs3);
    
  *moduloLimbs0 = dividendLimbs1;
  *moduloLimbs1 = dividendLimbs2;
  *moduloLimbs2 = dividendLimbs3;
  return (uint2){divisorLimbs, 32-normalizeShiftCount};
}


uint2 modulo640to512(uint4 dividendLimbs0,
                     uint4 dividendLimbs1,
                     uint4 dividendLimbs2,
                     uint4 dividendLimbs3,
                     uint4 dividendLimbs4,
                     uint4 divisorLimbs0,
                     uint4 divisorLimbs1,
                     uint4 divisorLimbs2,
                     uint4 divisorLimbs3,
                     uint4 *moduloLimbs0,
                     uint4 *moduloLimbs1,
                     uint4 *moduloLimbs2,
                     uint4 *moduloLimbs3)
{
  // Detect dividend and divisor limbs count (remove trailing zero limbs)
  unsigned dividendLimbs = 20;
  unsigned divisorLimbs = 16;
  
  while (divisorLimbs && !divisorLimbs3.w) {
    lshiftByLimb4(&divisorLimbs0, &divisorLimbs1, &divisorLimbs2, &divisorLimbs3);
    divisorLimbs--;
  }  
  
  // Normalize dividend and divisor (high bit of divisor must be set to 1)
  unsigned normalizeShiftCount = 0;  
  uint32_t bit = 0x80000000;
  while (!(divisorLimbs3.w & bit)) {
    normalizeShiftCount++;
    bit >>= 1;  
  }    
  
  lshift5(&dividendLimbs0, &dividendLimbs1, &dividendLimbs2, &dividendLimbs3, &dividendLimbs4, normalizeShiftCount);
  lshift4(&divisorLimbs0, &divisorLimbs1, &divisorLimbs2, &divisorLimbs3, normalizeShiftCount);    
  
  
  while (dividendLimbs && !dividendLimbs4.w) {
    lshiftByLimb5(&dividendLimbs0, &dividendLimbs1, &dividendLimbs2, &dividendLimbs3, &dividendLimbs4);
    dividendLimbs--;
  }  
  
  for (unsigned i = 0; i < (dividendLimbs - divisorLimbs); i++) {
    uint32_t i32quotient;
    if (dividendLimbs4.w == divisorLimbs3.w) {
      i32quotient = 0xFFFFFFFF;
    } else {
      uint64_t i64dividend = (((uint64_t)dividendLimbs4.w) << 32) | dividendLimbs4.z;
      i32quotient = i64dividend / divisorLimbs3.w;
    }
    
    subMul1_512(&dividendLimbs0, &dividendLimbs1, &dividendLimbs2, &dividendLimbs3, &dividendLimbs4,
                 divisorLimbs0, divisorLimbs1, divisorLimbs2, divisorLimbs3,
                 i32quotient);    
    uint32_t borrow = dividendLimbs4.w;
    lshiftByLimb5(&dividendLimbs0, &dividendLimbs1, &dividendLimbs2, &dividendLimbs3, &dividendLimbs4);
    if (borrow) {
      add512(&dividendLimbs1, &dividendLimbs2, &dividendLimbs3, &dividendLimbs4, divisorLimbs0, divisorLimbs1, divisorLimbs2, divisorLimbs3);
      if (dividendLimbs4.w > divisorLimbs3.w)
        add512(&dividendLimbs1, &dividendLimbs2, &dividendLimbs3, &dividendLimbs4, divisorLimbs0, divisorLimbs1, divisorLimbs2, divisorLimbs3);
    }
  }
  
  rshift5(&dividendLimbs0, &dividendLimbs1, &dividendLimbs2, &dividendLimbs3, &dividendLimbs4, normalizeShiftCount);
  for (unsigned i = 0; i < (16-divisorLimbs); i++)
    rshiftByLimb5(&dividendLimbs0, &dividendLimbs1, &dividendLimbs2, &dividendLimbs3, &dividendLimbs4);
  
  *moduloLimbs0 = dividendLimbs1;
  *moduloLimbs1 = dividendLimbs2;
  *moduloLimbs2 = dividendLimbs3;
  *moduloLimbs3 = dividendLimbs4;
  return (uint2){divisorLimbs, 32-normalizeShiftCount};
}


uint32_t invert_limb(uint32_t limb)
{
  uint32_t inv = binvert_limb_table[(limb/2) & 0x7F];
  inv = 2*inv - inv*inv*limb;
  inv = 2*inv - inv*inv*limb;
  return -inv;
}

void redc256_round_v3(uint4 op1l0, uint4 op1l1, uint32_t m1, uint32_t *m2, uint32_t invm,
                      uint64_t *R0, uint64_t *R1, uint64_t *R2, uint64_t *R3,
                      uint64_t *R4, uint64_t *R5, uint64_t *R6, uint64_t *R7)
{
  uint4 m1l0 = mul_hi(op1l0, m1);
  uint4 m1l1 = mul_hi(op1l1, m1);
  
  *m2 = invm * ((uint32_t)*R0 + m1l0.x);
  uint4 m2l0 = op1l0 * (*m2);
  uint4 m2l1 = op1l1 * (*m2);
  
  union {
    uint2 v32;
    ulong v64;
  } Int;  
  
  *R0 += m1l0.x; *R0 += m2l0.x;
  *R1 += m1l0.y; *R1 += m2l0.y; Int.v64 = *R0; *R1 += Int.v32.y;
  *R2 += m1l0.z; *R2 += m2l0.z;
  *R3 += m1l0.w; *R3 += m2l0.w;
  *R4 += m1l1.x; *R4 += m2l1.x;
  *R5 += m1l1.y; *R5 += m2l1.y;
  *R6 += m1l1.z; *R6 += m2l1.z;
  *R7 += m1l1.w; *R7 += m2l1.w;  
}


void redc384_round_v3(uint4 op1l0, uint4 op1l1, uint4 op1l2, uint32_t m1, uint32_t *m2, uint32_t invm,
                      uint64_t *R0, uint64_t *R1, uint64_t *R2, uint64_t *R3,
                      uint64_t *R4, uint64_t *R5, uint64_t *R6, uint64_t *R7,
                      uint64_t *R8, uint64_t *R9, uint64_t *R10, uint64_t *R11)
{
  uint4 m1l0 = mul_hi(op1l0, m1);
  uint4 m1l1 = mul_hi(op1l1, m1);
  uint4 m1l2 = mul_hi(op1l2, m1);
  
  *m2 = invm * ((uint32_t)*R0 + m1l0.x);
  uint4 m2l0 = op1l0 * (*m2);
  uint4 m2l1 = op1l1 * (*m2);
  uint4 m2l2 = op1l2 * (*m2);
  
  union {
    uint2 v32;
    ulong v64;
  } Int; 
  
  *R0 += m1l0.x; *R0 += m2l0.x;
  *R1 += m1l0.y; *R1 += m2l0.y; Int.v64 = *R0; *R1 += Int.v32.y;
  *R2 += m1l0.z; *R2 += m2l0.z;
  *R3 += m1l0.w; *R3 += m2l0.w;
  *R4 += m1l1.x; *R4 += m2l1.x;
  *R5 += m1l1.y; *R5 += m2l1.y;
  *R6 += m1l1.z; *R6 += m2l1.z;
  *R7 += m1l1.w; *R7 += m2l1.w;  
  *R8 += m1l2.x; *R8 += m2l2.x;
  *R9 += m1l2.y; *R9 += m2l2.y;
  *R10 += m1l2.z; *R10 += m2l2.z;
  *R11 += m1l2.w; *R11 += m2l2.w;
}


void redc448_round_v3(uint4 op1l0, uint4 op1l1, uint4 op1l2, uint2 op1l3, uint32_t m1, uint32_t *m2, uint32_t invm,
                      uint64_t *R0, uint64_t *R1, uint64_t *R2, uint64_t *R3,
                      uint64_t *R4, uint64_t *R5, uint64_t *R6, uint64_t *R7,
                      uint64_t *R8, uint64_t *R9, uint64_t *R10, uint64_t *R11,
                      uint64_t *R12, uint64_t *R13)
{
  uint4 m1l0 = mul_hi(op1l0, m1);
  uint4 m1l1 = mul_hi(op1l1, m1);
  uint4 m1l2 = mul_hi(op1l2, m1);
  uint2 m1l3 = mul_hi(op1l3, m1);
  
  *m2 = invm * ((uint32_t)*R0 + m1l0.x);
  uint4 m2l0 = op1l0 * (*m2);
  uint4 m2l1 = op1l1 * (*m2);
  uint4 m2l2 = op1l2 * (*m2);
  uint2 m2l3 = op1l3 * (*m2);
  
  union {
    uint2 v32;
    ulong v64;
  } Int; 
  
  *R0 += m1l0.x; *R0 += m2l0.x;
  *R1 += m1l0.y; *R1 += m2l0.y; Int.v64 = *R0; *R1 += Int.v32.y;
  *R2 += m1l0.z; *R2 += m2l0.z;
  *R3 += m1l0.w; *R3 += m2l0.w;
  *R4 += m1l1.x; *R4 += m2l1.x;
  *R5 += m1l1.y; *R5 += m2l1.y;
  *R6 += m1l1.z; *R6 += m2l1.z;
  *R7 += m1l1.w; *R7 += m2l1.w;  
  *R8 += m1l2.x; *R8 += m2l2.x;
  *R9 += m1l2.y; *R9 += m2l2.y;
  *R10 += m1l2.z; *R10 += m2l2.z;
  *R11 += m1l2.w; *R11 += m2l2.w;
  *R12 += m1l3.x; *R12 += m2l3.x;
  *R13 += m1l3.y; *R13 += m2l3.y;
}


void redc1_256_v3(uint4 limbs0, uint4 limbs1, uint4 limbs2, uint4 limbs3,
                  uint4 moduloLimb0, uint4 moduloLimb1,
                  uint32_t invm,
                  uint4 *ResultLimb0, uint4 *ResultLimb1)
{
  ulong R0x = limbs0.x;
  ulong R0y = limbs0.y;
  ulong R0z = limbs0.z;
  ulong R0w = limbs0.w;
  ulong R1x = limbs1.x;
  ulong R1y = limbs1.y;
  ulong R1z = limbs1.z;
  ulong R1w = limbs1.w;
  ulong R2x = limbs2.x;
  ulong R2y = limbs2.y;
  ulong R2z = limbs2.z;
  ulong R2w = limbs2.w;
  ulong R3x = limbs3.x;
  ulong R3y = limbs3.y;
  ulong R3z = limbs3.z;
  ulong R3w = limbs3.w;
  
  uint32_t i0, i1, i2, i3, i4, i5, i6, i7;
  union {
    uint2 v32;
    ulong v64;
  } Int;     
  
  i0 = limbs0.x * invm;
  {
    uint4 M1l0 = moduloLimb0 * i0;
    uint4 M1l1 = moduloLimb1 * i0;
    R0x += M1l0.x;
    R0y += M1l0.y; Int.v64 = R0x; R0y += Int.v32.y;
    R0z += M1l0.z;
    R0w += M1l0.w;
    R1x += M1l1.x;
    R1y += M1l1.y;
    R1z += M1l1.z;
    R1w += M1l1.w;
  }
  
  redc256_round_v3(moduloLimb0, moduloLimb1, i0, &i1, invm, &R0y, &R0z, &R0w, &R1x, &R1y, &R1z, &R1w, &R2x);
  redc256_round_v3(moduloLimb0, moduloLimb1, i1, &i2, invm, &R0z, &R0w, &R1x, &R1y, &R1z, &R1w, &R2x, &R2y);
  redc256_round_v3(moduloLimb0, moduloLimb1, i2, &i3, invm, &R0w, &R1x, &R1y, &R1z, &R1w, &R2x, &R2y, &R2z);
  redc256_round_v3(moduloLimb0, moduloLimb1, i3, &i4, invm, &R1x, &R1y, &R1z, &R1w, &R2x, &R2y, &R2z, &R2w);
  redc256_round_v3(moduloLimb0, moduloLimb1, i4, &i5, invm, &R1y, &R1z, &R1w, &R2x, &R2y, &R2z, &R2w, &R3x);
  redc256_round_v3(moduloLimb0, moduloLimb1, i5, &i6, invm, &R1z, &R1w, &R2x, &R2y, &R2z, &R2w, &R3x, &R3y);
  redc256_round_v3(moduloLimb0, moduloLimb1, i6, &i7, invm, &R1w, &R2x, &R2y, &R2z, &R2w, &R3x, &R3y, &R3z);
  
  {
    uint4 M1l0 = mul_hi(moduloLimb0, i7);
    uint4 M1l1 = mul_hi(moduloLimb1, i7);
    R2x += M1l0.x;
    R2y += M1l0.y; Int.v64 = R2x; R2y += Int.v32.y;
    R2z += M1l0.z; Int.v64 = R2y; R2z += Int.v32.y;
    R2w += M1l0.w; Int.v64 = R2z; R2w += Int.v32.y;
    R3x += M1l1.x; Int.v64 = R2w; R3x += Int.v32.y;
    R3y += M1l1.y; Int.v64 = R3x; R3y += Int.v32.y;
    R3z += M1l1.z; Int.v64 = R3y; R3z += Int.v32.y;
    R3w += M1l1.w; Int.v64 = R3z; R3w += Int.v32.y;
  }

  *ResultLimb0 = (uint4){R2x, R2y, R2z, R2w};
  *ResultLimb1 = (uint4){R3x, R3y, R3z, R3w};  
  
  {
    Int.v64 = R3w;
    uint4 l0 = Int.v32.y ? moduloLimb0 : (uint4){0, 0, 0, 0};
    uint4 l1 = Int.v32.y ? moduloLimb1 : (uint4){0, 0, 0, 0};   
    sub256(ResultLimb0, ResultLimb1, l0, l1);
  }  
}

void redc1_384_v3(uint4 limbs0, uint4 limbs1, uint4 limbs2, uint4 limbs3, uint4 limbs4, uint4 limbs5,
                  uint4 moduloLimb0, uint4 moduloLimb1, uint4 moduloLimb2,
                  uint32_t invm,
                  uint4 *ResultLimb0, uint4 *ResultLimb1, uint4 *ResultLimb2)
{
  ulong R0x = limbs0.x;
  ulong R0y = limbs0.y;
  ulong R0z = limbs0.z;
  ulong R0w = limbs0.w;
  ulong R1x = limbs1.x;
  ulong R1y = limbs1.y;
  ulong R1z = limbs1.z;
  ulong R1w = limbs1.w;
  ulong R2x = limbs2.x;
  ulong R2y = limbs2.y;
  ulong R2z = limbs2.z;
  ulong R2w = limbs2.w;
  ulong R3x = limbs3.x;
  ulong R3y = limbs3.y;
  ulong R3z = limbs3.z;
  ulong R3w = limbs3.w;
  ulong R4x = limbs4.x;
  ulong R4y = limbs4.y;
  ulong R4z = limbs4.z;
  ulong R4w = limbs4.w;
  ulong R5x = limbs5.x;
  ulong R5y = limbs5.y;
  ulong R5z = limbs5.z;
  ulong R5w = limbs5.w;  
  
  uint32_t i0, i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11;
  
  union {
    uint2 v32;
    ulong v64;
  } Int;   
  
  i0 = limbs0.x * invm;
  {
    uint4 M1l0 = moduloLimb0 * i0;
    uint4 M1l1 = moduloLimb1 * i0;
    uint4 M1l2 = moduloLimb2 * i0;
    R0x += M1l0.x;
    R0y += M1l0.y; Int.v64 = R0x; R0y += Int.v32.y;
    R0z += M1l0.z;
    R0w += M1l0.w;
    R1x += M1l1.x;
    R1y += M1l1.y;
    R1z += M1l1.z;
    R1w += M1l1.w;
    R2x += M1l2.x;
    R2y += M1l2.y;
    R2z += M1l2.z;
    R2w += M1l2.w;
  }
  
  redc384_round_v3(moduloLimb0, moduloLimb1, moduloLimb2, i0, &i1, invm, &R0y, &R0z, &R0w, &R1x, &R1y, &R1z, &R1w, &R2x, &R2y, &R2z, &R2w, &R3x);
  redc384_round_v3(moduloLimb0, moduloLimb1, moduloLimb2, i1, &i2, invm, &R0z, &R0w, &R1x, &R1y, &R1z, &R1w, &R2x, &R2y, &R2z, &R2w, &R3x, &R3y);
  redc384_round_v3(moduloLimb0, moduloLimb1, moduloLimb2, i2, &i3, invm, &R0w, &R1x, &R1y, &R1z, &R1w, &R2x, &R2y, &R2z, &R2w, &R3x, &R3y, &R3z);
  redc384_round_v3(moduloLimb0, moduloLimb1, moduloLimb2, i3, &i4, invm, &R1x, &R1y, &R1z, &R1w, &R2x, &R2y, &R2z, &R2w, &R3x, &R3y, &R3z, &R3w);
  redc384_round_v3(moduloLimb0, moduloLimb1, moduloLimb2, i4, &i5, invm, &R1y, &R1z, &R1w, &R2x, &R2y, &R2z, &R2w, &R3x, &R3y, &R3z, &R3w, &R4x);
  redc384_round_v3(moduloLimb0, moduloLimb1, moduloLimb2, i5, &i6, invm, &R1z, &R1w, &R2x, &R2y, &R2z, &R2w, &R3x, &R3y, &R3z, &R3w, &R4x, &R4y);
  redc384_round_v3(moduloLimb0, moduloLimb1, moduloLimb2, i6, &i7, invm, &R1w, &R2x, &R2y, &R2z, &R2w, &R3x, &R3y, &R3z, &R3w, &R4x, &R4y, &R4z);
  redc384_round_v3(moduloLimb0, moduloLimb1, moduloLimb2, i7, &i8, invm, &R2x, &R2y, &R2z, &R2w, &R3x, &R3y, &R3z, &R3w, &R4x, &R4y, &R4z, &R4w);
  redc384_round_v3(moduloLimb0, moduloLimb1, moduloLimb2, i8, &i9, invm, &R2y, &R2z, &R2w, &R3x, &R3y, &R3z, &R3w, &R4x, &R4y, &R4z, &R4w, &R5x);
  redc384_round_v3(moduloLimb0, moduloLimb1, moduloLimb2, i9, &i10, invm, &R2z, &R2w, &R3x, &R3y, &R3z, &R3w, &R4x, &R4y, &R4z, &R4w, &R5x, &R5y);
  redc384_round_v3(moduloLimb0, moduloLimb1, moduloLimb2, i10, &i11, invm, &R2w, &R3x, &R3y, &R3z, &R3w, &R4x, &R4y, &R4z, &R4w, &R5x, &R5y, &R5z);
    
  {
    uint4 M1l0 = mul_hi(moduloLimb0, i11);
    uint4 M1l1 = mul_hi(moduloLimb1, i11);    
    uint4 M1l2 = mul_hi(moduloLimb2, i11);        
    R3x += M1l0.x;
    R3y += M1l0.y; Int.v64 = R3x; R3y += Int.v32.y;
    R3z += M1l0.z; Int.v64 = R3y; R3z += Int.v32.y;
    R3w += M1l0.w; Int.v64 = R3z; R3w += Int.v32.y;
    R4x += M1l1.x; Int.v64 = R3w; R4x += Int.v32.y;
    R4y += M1l1.y; Int.v64 = R4x; R4y += Int.v32.y;
    R4z += M1l1.z; Int.v64 = R4y; R4z += Int.v32.y;
    R4w += M1l1.w; Int.v64 = R4z; R4w += Int.v32.y;
    R5x += M1l2.x; Int.v64 = R4w; R5x += Int.v32.y;
    R5y += M1l2.y; Int.v64 = R5x; R5y += Int.v32.y;
    R5z += M1l2.z; Int.v64 = R5y; R5z += Int.v32.y;
    R5w += M1l2.w; Int.v64 = R5z; R5w += Int.v32.y;
  }
  
  *ResultLimb0 = (uint4){R3x, R3y, R3z, R3w};
  *ResultLimb1 = (uint4){R4x, R4y, R4z, R4w};  
  *ResultLimb2 = (uint4){R5x, R5y, R5z, R5w};    

 
  {
    Int.v64 = R5w;
    uint4 l0 = Int.v32.y ? moduloLimb0 : (uint4){0, 0, 0, 0};
    uint4 l1 = Int.v32.y ? moduloLimb1 : (uint4){0, 0, 0, 0};
    uint4 l2 = Int.v32.y ? moduloLimb2 : (uint4){0, 0, 0, 0};    
    sub384(ResultLimb0, ResultLimb1, ResultLimb2, l0, l1, l2);
  }
}




void redc1_448_v3(uint4 limbs0, uint4 limbs1, uint4 limbs2, uint4 limbs3, uint4 limbs4, uint4 limbs5, uint4 limbs6,
                  uint4 moduloLimb0, uint4 moduloLimb1, uint4 moduloLimb2, uint2 moduloLimb3,
                  uint32_t invm,
                  uint4 *ResultLimb0, uint4 *ResultLimb1, uint4 *ResultLimb2, uint2 *ResultLimb3)
{
  ulong R0x = limbs0.x;
  ulong R0y = limbs0.y;
  ulong R0z = limbs0.z;
  ulong R0w = limbs0.w;
  ulong R1x = limbs1.x;
  ulong R1y = limbs1.y;
  ulong R1z = limbs1.z;
  ulong R1w = limbs1.w;
  ulong R2x = limbs2.x;
  ulong R2y = limbs2.y;
  ulong R2z = limbs2.z;
  ulong R2w = limbs2.w;
  ulong R3x = limbs3.x;
  ulong R3y = limbs3.y;
  ulong R3z = limbs3.z;
  ulong R3w = limbs3.w;
  ulong R4x = limbs4.x;
  ulong R4y = limbs4.y;
  ulong R4z = limbs4.z;
  ulong R4w = limbs4.w;
  ulong R5x = limbs5.x;
  ulong R5y = limbs5.y;
  ulong R5z = limbs5.z;
  ulong R5w = limbs5.w;
  ulong R6x = limbs6.x;
  ulong R6y = limbs6.y;
  ulong R6z = limbs6.z;
  ulong R6w = limbs6.w;
  
  uint32_t i0, i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13;
  
  union {
    uint2 v32;
    ulong v64;
  } Int;   
  
  i0 = limbs0.x * invm;
  {
    uint4 M1l0 = moduloLimb0 * i0;
    uint4 M1l1 = moduloLimb1 * i0;
    uint4 M1l2 = moduloLimb2 * i0;
    uint2 M1l3 = moduloLimb3 * i0;
    R0x += M1l0.x;
    R0y += M1l0.y; Int.v64 = R0x; R0y += Int.v32.y;
    R0z += M1l0.z;
    R0w += M1l0.w;
    R1x += M1l1.x;
    R1y += M1l1.y;
    R1z += M1l1.z;
    R1w += M1l1.w;
    R2x += M1l2.x;
    R2y += M1l2.y;
    R2z += M1l2.z;
    R2w += M1l2.w;
    R3x += M1l3.x;
    R3y += M1l3.y;
  }
  
  redc448_round_v3(moduloLimb0, moduloLimb1, moduloLimb2, moduloLimb3, i0,  &i1,  invm, &R0y, &R0z, &R0w, &R1x, &R1y, &R1z, &R1w, &R2x, &R2y, &R2z, &R2w, &R3x, &R3y, &R3z);
  redc448_round_v3(moduloLimb0, moduloLimb1, moduloLimb2, moduloLimb3, i1,  &i2,  invm, &R0z, &R0w, &R1x, &R1y, &R1z, &R1w, &R2x, &R2y, &R2z, &R2w, &R3x, &R3y, &R3z, &R3w);  
  redc448_round_v3(moduloLimb0, moduloLimb1, moduloLimb2, moduloLimb3, i2,  &i3,  invm, &R0w, &R1x, &R1y, &R1z, &R1w, &R2x, &R2y, &R2z, &R2w, &R3x, &R3y, &R3z, &R3w, &R4x);    
  redc448_round_v3(moduloLimb0, moduloLimb1, moduloLimb2, moduloLimb3, i3,  &i4,  invm, &R1x, &R1y, &R1z, &R1w, &R2x, &R2y, &R2z, &R2w, &R3x, &R3y, &R3z, &R3w, &R4x, &R4y);
  redc448_round_v3(moduloLimb0, moduloLimb1, moduloLimb2, moduloLimb3, i4,  &i5,  invm, &R1y, &R1z, &R1w, &R2x, &R2y, &R2z, &R2w, &R3x, &R3y, &R3z, &R3w, &R4x, &R4y, &R4z);
  redc448_round_v3(moduloLimb0, moduloLimb1, moduloLimb2, moduloLimb3, i5,  &i6,  invm, &R1z, &R1w, &R2x, &R2y, &R2z, &R2w, &R3x, &R3y, &R3z, &R3w, &R4x, &R4y, &R4z, &R4w);  
  redc448_round_v3(moduloLimb0, moduloLimb1, moduloLimb2, moduloLimb3, i6,  &i7,  invm, &R1w, &R2x, &R2y, &R2z, &R2w, &R3x, &R3y, &R3z, &R3w, &R4x, &R4y, &R4z, &R4w, &R5x);
  redc448_round_v3(moduloLimb0, moduloLimb1, moduloLimb2, moduloLimb3, i7,  &i8,  invm, &R2x, &R2y, &R2z, &R2w, &R3x, &R3y, &R3z, &R3w, &R4x, &R4y, &R4z, &R4w, &R5x, &R5y);  
  redc448_round_v3(moduloLimb0, moduloLimb1, moduloLimb2, moduloLimb3, i8,  &i9,  invm, &R2y, &R2z, &R2w, &R3x, &R3y, &R3z, &R3w, &R4x, &R4y, &R4z, &R4w, &R5x, &R5y, &R5z);
  redc448_round_v3(moduloLimb0, moduloLimb1, moduloLimb2, moduloLimb3, i9 , &i10, invm, &R2z, &R2w, &R3x, &R3y, &R3z, &R3w, &R4x, &R4y, &R4z, &R4w, &R5x, &R5y, &R5z, &R5w);  
  redc448_round_v3(moduloLimb0, moduloLimb1, moduloLimb2, moduloLimb3, i10, &i11, invm, &R2w, &R3x, &R3y, &R3z, &R3w, &R4x, &R4y, &R4z, &R4w, &R5x, &R5y, &R5z, &R5w, &R6x);    
  redc448_round_v3(moduloLimb0, moduloLimb1, moduloLimb2, moduloLimb3, i11, &i12, invm, &R3x, &R3y, &R3z, &R3w, &R4x, &R4y, &R4z, &R4w, &R5x, &R5y, &R5z, &R5w, &R6x, &R6y);      
  redc448_round_v3(moduloLimb0, moduloLimb1, moduloLimb2, moduloLimb3, i12, &i13, invm, &R3y, &R3z, &R3w, &R4x, &R4y, &R4z, &R4w, &R5x, &R5y, &R5z, &R5w, &R6x, &R6y, &R6z);        
  
  {
    uint4 M1l0 = mul_hi(moduloLimb0, i13);
    uint4 M1l1 = mul_hi(moduloLimb1, i13);
    uint4 M1l2 = mul_hi(moduloLimb2, i13);
    uint2 M1l3 = mul_hi(moduloLimb3, i13);
    
    R3z += M1l0.x;
    R3w += M1l0.y; Int.v64 = R3z; R3w += Int.v32.y;
    R4x += M1l0.z; Int.v64 = R3w; R4x += Int.v32.y;
    R4y += M1l0.w; Int.v64 = R4x; R4y += Int.v32.y;
    R4z += M1l1.x; Int.v64 = R4y; R4z += Int.v32.y;
    R4w += M1l1.y; Int.v64 = R4z; R4w += Int.v32.y;
    R5x += M1l1.z; Int.v64 = R4w; R5x += Int.v32.y;
    R5y += M1l1.w; Int.v64 = R5x; R5y += Int.v32.y;
    R5z += M1l2.x; Int.v64 = R5y; R5z += Int.v32.y;
    R5w += M1l2.y; Int.v64 = R5z; R5w += Int.v32.y;
    R6x += M1l2.z; Int.v64 = R5w; R6x += Int.v32.y;
    R6y += M1l2.w; Int.v64 = R6x; R6y += Int.v32.y;
    R6z += M1l3.x; Int.v64 = R6y; R6z += Int.v32.y;
    R6w += M1l3.y; Int.v64 = R6z; R6w += Int.v32.y;
  }
  
  *ResultLimb0 = (uint4){R3z, R3w, R4x, R4y};
  *ResultLimb1 = (uint4){R4z, R4w, R5x, R5y};  
  *ResultLimb2 = (uint4){R5z, R5w, R6x, R6y};    
  *ResultLimb3 = (uint2){R6z, R6w};
  
  
  {
    Int.v64 = R6w;
    uint4 l0 = Int.v32.y ? moduloLimb0 : (uint4){0, 0, 0, 0};
    uint4 l1 = Int.v32.y ? moduloLimb1 : (uint4){0, 0, 0, 0};
    uint4 l2 = Int.v32.y ? moduloLimb2 : (uint4){0, 0, 0, 0};    
    uint2 l3 = Int.v32.y ? moduloLimb3 : (uint2){0, 0};
    sub448(ResultLimb0, ResultLimb1, ResultLimb2, ResultLimb3, l0, l1, l2, l3);
  }
}




uint32_t dec128(uint4 *a0)
{
  --(*a0).x;
  (*a0).y -= ((*a0).x == 0xFFFFFFFF);
  (*a0).z -= ((*a0).y == 0xFFFFFFFF);
  (*a0).w -= ((*a0).z == 0xFFFFFFFF);
  return (*a0).w == 0xFFFFFFFF;
}

uint32_t dec256(uint4 *a0, uint4 *a1)
{
  --(*a0).x;
  (*a0).y -= ((*a0).x == 0xFFFFFFFF);
  (*a0).z -= ((*a0).y == 0xFFFFFFFF);
  (*a0).w -= ((*a0).z == 0xFFFFFFFF);
  (*a1).x -= ((*a0).w == 0xFFFFFFFF);
  (*a1).y -= ((*a1).x == 0xFFFFFFFF);
  (*a1).z -= ((*a1).y == 0xFFFFFFFF);
  (*a1).w -= ((*a1).w == 0xFFFFFFFF);
  return (*a1).w == 0xFFFFFFFF;
}

void inc384(uint4 *a0, uint4 *a1, uint4 *a2)
{
  ++(*a0).x;
  (*a0).y += ((*a0).x == 0);
  (*a0).z += ((*a0).y == 0);
  (*a0).w += ((*a0).z == 0);
  
  (*a1).x += ((*a0).w == 0);
  (*a1).y += ((*a1).x == 0);
  (*a1).z += ((*a1).y == 0);
  (*a1).w += ((*a1).z == 0);
  
  (*a2).x += ((*a1).w == 0);
  (*a2).y += ((*a2).x == 0);
  (*a2).z += ((*a2).y == 0);
  (*a2).w += ((*a2).z == 0);
}

void dec384(uint4 *a0, uint4 *a1, uint4 *a2)
{
  --(*a0).x;
  (*a0).y -= ((*a0).x == 0xFFFFFFFF);
  (*a0).z -= ((*a0).y == 0xFFFFFFFF);
  (*a0).w -= ((*a0).z == 0xFFFFFFFF);
  (*a1).x -= ((*a0).w == 0xFFFFFFFF);
  (*a1).y -= ((*a1).x == 0xFFFFFFFF);
  (*a1).z -= ((*a1).y == 0xFFFFFFFF);
  (*a1).w -= ((*a1).z == 0xFFFFFFFF);
  (*a2).x -= ((*a1).w == 0xFFFFFFFF);
  (*a2).y -= ((*a2).x == 0xFFFFFFFF);
  (*a2).z -= ((*a2).y == 0xFFFFFFFF);
  (*a2).w -= ((*a2).z == 0xFFFFFFFF);  
}

void montgomeryMul256(uint4 *rl0, uint4 *rl1,
                      uint4 ml0, uint4 ml1, 
                      uint4 modl0, uint4 modl1,
                      uint32_t inverted)
{
  uint4 m0, m1, m2, m3;
  mul256schoolBook_v3(*rl0, *rl1, ml0, ml1, &m0, &m1, &m2, &m3);
  redc1_256_v3(m0, m1, m2, m3, modl0, modl1, inverted, rl0, rl1);
}

void FermatTest256(uint4 limbs0, uint4 limbs1,
                   uint4 *resultLimbs0, uint4 *resultLimbs1)
{
  // Invert lowest modulo limb
  __private uint4 precomputed[2*16*2];
  uint32_t inverted = invert_limb(limbs0.x);
  
  uint2 bitSize;
  uint4 BmLimbs0, BmLimbs1;
  {
    uint4 dl2 = {1, 0, 0, 0};
    uint4 dl1 = {0, 0, 0, 0};
    uint4 dl0 = {0, 0, 0, 0};
  
    // Retrieve of "1" in Montgomery representation    
    modulo384to256(dl0, dl1, dl2, limbs0, limbs1, &BmLimbs0, &BmLimbs1);
    precomputed[2*0 + 0] = BmLimbs0;
    precomputed[2*0 + 1] = BmLimbs1;
    
    // Retrieve of "2" in Montgomery representation    
    dl2.x = 2;
    bitSize = modulo384to256(dl0, dl1, dl2, limbs0, limbs1, &BmLimbs0, &BmLimbs1);
    precomputed[2*1 + 0] = BmLimbs0;
    precomputed[2*1 + 1] = BmLimbs1;
    --bitSize.y;
    if (bitSize.y == 0) {
      --bitSize.x;
      bitSize.y = 32;
    }
  }
  
  // Calcutate powers range 2-16 and cache in global memory
  uint4 m0, m1, m2, m3;  
  uint4 redcl0 = BmLimbs0,
        redcl1 = BmLimbs1;
  for (unsigned i = 2; i < 16; i++) {
    montgomeryMul256(&redcl0, &redcl1, BmLimbs0, BmLimbs1, limbs0, limbs1, inverted);
    precomputed[2*i + 0] = redcl0;
    precomputed[2*i + 1] = redcl1;
  }
  
  uint4 exp0 = limbs0;
  uint4 exp1 = limbs1;
  --exp0.x;
  
  for (unsigned i = 0; i < 8-bitSize.x; i++)
    lshiftByLimb2(&exp0, &exp1);
  
  barrier(CLK_LOCAL_MEM_FENCE);  
  exp1.w <<= (bitSize.y ? 32-bitSize.y : 0);  
  redcl0 = BmLimbs0;
  redcl1 = BmLimbs1;
  
  unsigned shiftCount = 0;      
  unsigned square = 1;
  unsigned groupSize = (bitSize.y % 4) ? (bitSize.y % 4) : 4;
  unsigned bitcount = bitSize.y;
  while (bitSize.x) {
    while (bitcount) {
      uint4 mult0, mult1;
      if (!square) {
        unsigned index = 2 * (exp1.w >> (32 - shiftCount));
        mult0 = precomputed[index];
        mult1 = precomputed[index+1];
      }

      groupSize = square ? groupSize : 1;
      exp1.w <<= shiftCount;
      bitcount -= shiftCount;
      shiftCount = square ? groupSize : 0;            
      while (groupSize) {
        mult0 = square ? redcl0 : mult0;
        mult1 = square ? redcl1 : mult1;
        montgomeryMul256(&redcl0, &redcl1, mult0, mult1, limbs0, limbs1, inverted);
        groupSize--;
      }
      
      square ^= 1;
      groupSize = 4;
    }
    
    lshiftByLimb2(&exp0, &exp1);
    --bitSize.x;
    bitcount = 32;
  }
  
  barrier(CLK_LOCAL_MEM_FENCE);
  redc1_256_v3(redcl0, redcl1, 0, 0, limbs0, limbs1, inverted, resultLimbs0, resultLimbs1); 
  return;
}



__kernel void fermatTestBenchMark256(__global uint32_t *numbers,
                                     __global uint32_t *result,
                                     unsigned elementsNum)
{
  __global uint4 *numbersPtr = (__global uint4*)numbers;
  __global uint4 *resultPtr = (__global uint4*)result;
  
  unsigned globalSize = get_global_size(0);
  for (unsigned repeatNum = 0; repeatNum < 1; repeatNum++) {
    for (unsigned i = get_global_id(0); i < elementsNum; i += globalSize) {
      uint4 numbersLimbs0 = numbersPtr[i*2];
      uint4 numbersLimbs1 = numbersPtr[i*2+1];
      
      uint4 resultLimbs0;
      uint4 resultLimbs1;
      FermatTest256(numbersLimbs0, numbersLimbs1, &resultLimbs0, &resultLimbs1);
      
      resultPtr[i*2] = resultLimbs0;
      resultPtr[i*2+1] = resultLimbs1;
    }
  }
}

void montgomeryMul384(uint4 *rl0, uint4 *rl1, uint4 *rl2,
                      uint4 ml0, uint4 ml1, uint4 ml2, 
                      uint4 modl0, uint4 modl1, uint4 modl2,
                      uint32_t inverted)
{
  uint4 m0, m1, m2, m3, m4, m5;
  mul384schoolBook_v3(*rl0, *rl1, *rl2, ml0, ml1, ml2, &m0, &m1, &m2, &m3, &m4, &m5);
  redc1_384_v3(m0, m1, m2, m3, m4, m5, modl0, modl1, modl2, inverted, rl0, rl1, rl2);
}

void FermatTest384(uint4 limbs0, uint4 limbs1, uint4 limbs2,
                   uint4 *resultLimbs0, uint4 *resultLimbs1, uint4 *resultLimbs2)
{
  // Invert lowest modulo limb
  __private uint4 precomputed[3*16];
  uint32_t inverted = invert_limb(limbs0.x);

  uint2 bitSize;
  uint4 BmLimbs0, BmLimbs1, BmLimbs2;
  {
    uint4 dl3 = {1, 0, 0, 0};
    uint4 dl2 = {0, 0, 0, 0};
    uint4 dl1 = {0, 0, 0, 0};
    uint4 dl0 = {0, 0, 0, 0};
   
    // Retrieve of "1" in Montgomery representation        
    modulo512to384(dl0, dl1, dl2, dl3, limbs0, limbs1, limbs2, &BmLimbs0, &BmLimbs1, &BmLimbs2);
    precomputed[3*0 + 0] = BmLimbs0;
    precomputed[3*0 + 1] = BmLimbs1;
    precomputed[3*0 + 2] = BmLimbs2;    
    
    // Retrieve of "2" in Montgomery representation        
    dl3.x = 2;
    bitSize = modulo512to384(dl0, dl1, dl2, dl3, limbs0, limbs1, limbs2, &BmLimbs0, &BmLimbs1, &BmLimbs2);
    precomputed[3*1 + 0] = BmLimbs0;
    precomputed[3*1 + 1] = BmLimbs1;
    precomputed[3*1 + 2] = BmLimbs2;        
    --bitSize.y;
    if (bitSize.y == 0) {
      --bitSize.x;
      bitSize.y = 32;
    }
  }

  // Calcutate powers range 2-16 and cache in global memory  
  uint4 m0, m1, m2, m3, m4, m5;  
  uint4 redcl0 = BmLimbs0,
        redcl1 = BmLimbs1,
        redcl2 = BmLimbs2;
  for (unsigned i = 2; i < 16; i++) {
    montgomeryMul384(&redcl0, &redcl1, &redcl2, BmLimbs0, BmLimbs1, BmLimbs2, limbs0, limbs1, limbs2, inverted);
    precomputed[3*i + 0] = redcl0;
    precomputed[3*i + 1] = redcl1;
    precomputed[3*i + 2] = redcl2;
  }

  uint4 exp0 = limbs0;
  uint4 exp1 = limbs1;
  uint4 exp2 = limbs2;
  --exp0.x;
 
  for (unsigned i = 0; i < 12-bitSize.x; i++)
    lshiftByLimb3(&exp0, &exp1, &exp2);
  
  barrier(CLK_LOCAL_MEM_FENCE);  
  exp2.w <<= (bitSize.y ? 32-bitSize.y : 0);  
  redcl0 = BmLimbs0;
  redcl1 = BmLimbs1;
  redcl2 = BmLimbs2;
  
  unsigned shiftCount = 0;
  unsigned square = 1;
  unsigned groupSize = (bitSize.y % 4) ? (bitSize.y % 4) : 4;
  unsigned bitcount = bitSize.y;
  while (bitSize.x) {
    while (bitcount) {
      uint4 mult0, mult1, mult2;
      if (!square) {
        unsigned index = 3 * (exp2.w >> (32 - shiftCount));
        mult0 = precomputed[index];
        mult1 = precomputed[index+1];
        mult2 = precomputed[index+2];
      }

      groupSize = square ? groupSize : 1;
      exp2.w <<= shiftCount;
      bitcount -= shiftCount;
      shiftCount = square ? groupSize : 0;    
      while (groupSize) {
        mult0 = square ? redcl0 : mult0;
        mult1 = square ? redcl1 : mult1;
        mult2 = square ? redcl2 : mult2;
        montgomeryMul384(&redcl0, &redcl1, &redcl2, mult0, mult1, mult2, limbs0, limbs1, limbs2, inverted);
        groupSize--;
      }

      square ^= 1;
      groupSize = 4;
    }

    lshiftByLimb3(&exp0, &exp1, &exp2);
    --bitSize.x;
    bitcount = 32;
  }

  barrier(CLK_LOCAL_MEM_FENCE);
  redc1_384_v3(redcl0, redcl1, redcl2, 0, 0, 0, limbs0, limbs1, limbs2, inverted, resultLimbs0, resultLimbs1, resultLimbs2); 
  return;
}


__kernel void fermatTestBenchMark384(__global uint32_t *numbers,
                                     __global uint32_t *result,
                                     unsigned elementsNum)
{
  __global uint4 *numbersPtr = (__global uint4*)numbers;
  __global uint4 *resultPtr = (__global uint4*)result;
  
  unsigned globalSize = get_global_size(0);
  for (unsigned repeatNum = 0; repeatNum < 1; repeatNum++) {
    for (unsigned i = get_global_id(0); i < elementsNum; i += globalSize) {
      uint4 numbersLimbs0 = numbersPtr[i*3];
      uint4 numbersLimbs1 = numbersPtr[i*3+1];
      uint4 numbersLimbs2 = numbersPtr[i*3+2];
      
      uint4 resultLimbs0;
      uint4 resultLimbs1;
      uint4 resultLimbs2;
      FermatTest384(numbersLimbs0, numbersLimbs1, numbersLimbs2, &resultLimbs0, &resultLimbs1, &resultLimbs2);
      
      resultPtr[i*3] = resultLimbs0;
      resultPtr[i*3+1] = resultLimbs1;
      resultPtr[i*3+2] = resultLimbs2;
    }
  }
}

void montgomeryMul448(uint4 *rl0, uint4 *rl1, uint4 *rl2, uint2 *rl3,
                      uint4 ml0, uint4 ml1, uint4 ml2, uint2 ml3,
                      uint4 modl0, uint4 modl1, uint4 modl2, uint2 modl3,
                      uint32_t inverted)
{
  uint4 m0, m1, m2, m3, m4, m5, m6;
  mul448schoolBook_v3(*rl0, *rl1, *rl2, *rl3, ml0, ml1, ml2, ml3, &m0, &m1, &m2, &m3, &m4, &m5, &m6);
  redc1_448_v3(m0, m1, m2, m3, m4, m5, m6, modl0, modl1, modl2, modl3, inverted, rl0, rl1, rl2, rl3);
}

void FermatTest448(uint4 limbs0, uint4 limbs1, uint4 limbs2, uint2 limbs3,
                   uint4 *resultLimbs0, uint4 *resultLimbs1, uint4 *resultLimbs2, uint2 *resultLimbs3)
{
  // Invert lowest modulo limb
  __private uint4 precomputed[4*16];
  uint32_t inverted = invert_limb(limbs0.x);
  
  uint2 bitSize;
  uint4 BmLimbs0, BmLimbs1, BmLimbs2, BmLimbs3;
  {
    uint4 dl4 = {0, 0, 0, 0};
    uint4 dl3 = {0, 0, 1, 0};
    uint4 dl2 = {0, 0, 0, 0};
    uint4 dl1 = {0, 0, 0, 0};
    uint4 dl0 = {0, 0, 0, 0};    
        
    // Retrieve of "1" in Montgomery representation        
    modulo640to512(dl0, dl1, dl2, dl3, dl4,
                   limbs0, limbs1, limbs2, (uint4){limbs3.x, limbs3.y, 0, 0},
                   &BmLimbs0, &BmLimbs1, &BmLimbs2, &BmLimbs3);
    precomputed[4*0 + 0] = BmLimbs0;
    precomputed[4*0 + 1] = BmLimbs1;
    precomputed[4*0 + 2] = BmLimbs2;    
    precomputed[4*0 + 3] = BmLimbs3;
    
    // Retrieve of "2" in Montgomery representation        
    dl3.z = 2;
    bitSize = modulo640to512(dl0, dl1, dl2, dl3, dl4,
                             limbs0, limbs1, limbs2, (uint4){limbs3.x, limbs3.y, 0, 0},
                             &BmLimbs0, &BmLimbs1, &BmLimbs2, &BmLimbs3);
    precomputed[4*1 + 0] = BmLimbs0;
    precomputed[4*1 + 1] = BmLimbs1;
    precomputed[4*1 + 2] = BmLimbs2;
    precomputed[4*1 + 3] = BmLimbs3;
    --bitSize.y;
    if (bitSize.y == 0) {
      --bitSize.x;
      bitSize.y = 32;
    }
  }

  // Calcutate powers range 2-16 and cache in global memory  
  uint4 redcl0 = BmLimbs0,
        redcl1 = BmLimbs1,
        redcl2 = BmLimbs2;
  uint2 redcl3 = BmLimbs3.xy;
  
  for (unsigned i = 2; i < 16; i++) {
    montgomeryMul448(&redcl0, &redcl1, &redcl2, &redcl3,
                     BmLimbs0, BmLimbs1, BmLimbs2, BmLimbs3.xy,
                     limbs0, limbs1, limbs2, limbs3,
                     inverted);
    precomputed[4*i + 0] = redcl0;
    precomputed[4*i + 1] = redcl1;
    precomputed[4*i + 2] = redcl2;
    precomputed[4*i + 3] = (uint4){redcl3.x, redcl3.y, 0, 0};
  }
  
  uint4 exp0 = limbs0;
  uint4 exp1 = limbs1;
  uint4 exp2 = limbs2;
  uint4 exp3 = {limbs3.x, limbs3.y, 0, 0};
  --exp0.x;
  
  for (unsigned i = 0; i < 16-bitSize.x; i++)
    lshiftByLimb4(&exp0, &exp1, &exp2, &exp3);
  
  barrier(CLK_LOCAL_MEM_FENCE);  
  exp3.w <<= (bitSize.y ? 32-bitSize.y : 0);  
  redcl0 = BmLimbs0;
  redcl1 = BmLimbs1;
  redcl2 = BmLimbs2;
  redcl3 = BmLimbs3.xy;
  
  unsigned shiftCount = 0;
  unsigned square = 1;
  unsigned groupSize = (bitSize.y % 4) ? (bitSize.y % 4) : 4;
  unsigned bitcount = bitSize.y;
  while (bitSize.x) {
    while (bitcount) {
      uint4 mult0, mult1, mult2;
      uint2 mult3;
      if (!square) {
        unsigned index = 4 * (exp3.w >> (32 - shiftCount));
        mult0 = precomputed[index];
        mult1 = precomputed[index+1];
        mult2 = precomputed[index+2];
        mult3 = precomputed[index+3].xy;
      }
      
      groupSize = square ? groupSize : 1;
      exp3.w <<= shiftCount;
      bitcount -= shiftCount;
      shiftCount = square ? groupSize : 0;    
      while (groupSize) {
        mult0 = square ? redcl0 : mult0;
        mult1 = square ? redcl1 : mult1;
        mult2 = square ? redcl2 : mult2;
        mult3 = square ? redcl3 : mult3;
        montgomeryMul448(&redcl0, &redcl1, &redcl2, &redcl3,
                         mult0, mult1, mult2, mult3,
                         limbs0, limbs1, limbs2, limbs3,
                         inverted);
        groupSize--;
      }
      
      square ^= 1;
      groupSize = 4;
    }
    
    lshiftByLimb4(&exp0, &exp1, &exp2, &exp3);
    --bitSize.x;
    bitcount = 32;
  }
  
  barrier(CLK_LOCAL_MEM_FENCE);
  redc1_448_v3(redcl0, redcl1, redcl2, (uint4){redcl3.x, redcl3.y, 0, 0}, 0, 0, 0,
               limbs0, limbs1, limbs2, limbs3,
               inverted,
               resultLimbs0, resultLimbs1, resultLimbs2, resultLimbs3); 
}

__kernel void fermatTestBenchMark448(__global uint32_t *numbers,
                                     __global uint32_t *result,
                                     unsigned elementsNum)
{
  unsigned globalSize = get_global_size(0);
  for (unsigned repeatNum = 0; repeatNum < 1; repeatNum++) {
    for (unsigned i = get_global_id(0); i < elementsNum; i += globalSize) {
      uint4 numbersLimbs0 = {numbers[i*14], numbers[i*14+1], numbers[i*14+2], numbers[i*14+3]};
      uint4 numbersLimbs1 = {numbers[i*14+4], numbers[i*14+5], numbers[i*14+6], numbers[i*14+7]};
      uint4 numbersLimbs2 = {numbers[i*14+8], numbers[i*14+9], numbers[i*14+10], numbers[i*14+11]};
      uint2 numbersLimbs3 = {numbers[i*14+12], numbers[i*14+13]};      
      
      uint4 resultLimbs0;
      uint4 resultLimbs1;
      uint4 resultLimbs2;
      uint2 resultLimbs3;
      FermatTest448(numbersLimbs0, numbersLimbs1, numbersLimbs2, numbersLimbs3,
                    &resultLimbs0, &resultLimbs1, &resultLimbs2, &resultLimbs3);
      
      result[i*14] = resultLimbs0.x;
      result[i*14+1] = resultLimbs0.y;
      result[i*14+2] = resultLimbs0.z;
      result[i*14+3] = resultLimbs0.w;
      result[i*14+4] = resultLimbs1.x;
      result[i*14+5] = resultLimbs1.y;
      result[i*14+6] = resultLimbs1.z;
      result[i*14+7] = resultLimbs1.w;      
      result[i*14+8] = resultLimbs2.x;
      result[i*14+9] = resultLimbs2.y;
      result[i*14+10] = resultLimbs2.z;
      result[i*14+11] = resultLimbs2.w;      
      result[i*14+12] = resultLimbs3.x;
      result[i*14+13] = resultLimbs3.y;
    }
  }
}


__kernel void modulo384to256test(__global uint32_t *dividends,
                                 __global uint32_t *divisors,
                                 __global uint32_t *modulos,
                                 unsigned elementsNum)
{
  __global uint4 *dividendPtr = (__global uint4*)dividends;
  __global uint4 *divisorPtr = (__global uint4*)divisors;
  __global uint4 *modulosPtr = (__global uint4*)modulos;
  
  unsigned globalSize = get_global_size(0);
  
  for (unsigned repeatNum = 0; repeatNum < 32; repeatNum++) {
    for (unsigned i = get_global_id(0); i < elementsNum; i += globalSize) {
      uint4 dividendLimbs0 = dividendPtr[i*3];
      uint4 dividendLimbs1 = dividendPtr[i*3+1];
      uint4 dividendLimbs2 = dividendPtr[i*3+2];
      
      uint4 divisorLimbs0 = divisorPtr[i*2];
      uint4 divisorLimbs1 = divisorPtr[i*2+1];
      
      uint4 moduloLimb0;
      uint4 moduloLimb1;
      
      modulo384to256(dividendLimbs0, dividendLimbs1, dividendLimbs2,
                     divisorLimbs0, divisorLimbs1,
                     &moduloLimb0, &moduloLimb1);
      
      modulosPtr[i*2] = moduloLimb0;
      modulosPtr[i*2+1] = moduloLimb1;
    }
  }
}

__kernel void modulo512to384test(__global uint32_t *dividends,
                                 __global uint32_t *divisors,
                                 __global uint32_t *modulos,
                                 unsigned elementsNum)
{
  __global uint4 *dividendPtr = (__global uint4*)dividends;
  __global uint4 *divisorPtr = (__global uint4*)divisors;
  __global uint4 *modulosPtr = (__global uint4*)modulos;
  
  unsigned globalSize = get_global_size(0);
  
  for (unsigned repeatNum = 0; repeatNum < 32; repeatNum++) {
    for (unsigned i = get_global_id(0); i < elementsNum; i += globalSize) {
      uint4 dividendLimbs0 = dividendPtr[i*4];
      uint4 dividendLimbs1 = dividendPtr[i*4+1];
      uint4 dividendLimbs2 = dividendPtr[i*4+2];
      uint4 dividendLimbs3 = dividendPtr[i*4+3];
      
      uint4 divisorLimbs0 = divisorPtr[i*3];
      uint4 divisorLimbs1 = divisorPtr[i*3+1];
      uint4 divisorLimbs2 = divisorPtr[i*3+2];
      
      uint4 moduloLimb0;
      uint4 moduloLimb1;
      uint4 moduloLimb2;
      
      modulo512to384(dividendLimbs0, dividendLimbs1, dividendLimbs2, dividendLimbs3,
                     divisorLimbs0, divisorLimbs1, divisorLimbs2,
                     &moduloLimb0, &moduloLimb1, &moduloLimb2);
      
      modulosPtr[i*3] = moduloLimb0;
      modulosPtr[i*3+1] = moduloLimb1;
      modulosPtr[i*3+2] = moduloLimb2;
    }
  }
}

__kernel void modulo640to512test(__global uint32_t *dividends,
                                 __global uint32_t *divisors,
                                 __global uint32_t *modulos,
                                 unsigned elementsNum)
{
  __global uint4 *dividendPtr = (__global uint4*)dividends;
  __global uint4 *divisorPtr = (__global uint4*)divisors;
  __global uint4 *modulosPtr = (__global uint4*)modulos;
  
  unsigned globalSize = get_global_size(0);
  
  for (unsigned repeatNum = 0; repeatNum < 32; repeatNum++) {
    for (unsigned i = get_global_id(0); i < elementsNum; i += globalSize) {
      uint4 dividendLimbs0 = dividendPtr[i*5];
      uint4 dividendLimbs1 = dividendPtr[i*5+1];
      uint4 dividendLimbs2 = dividendPtr[i*5+2];
      uint4 dividendLimbs3 = dividendPtr[i*5+3];
      uint4 dividendLimbs4 = dividendPtr[i*5+4];
      
      uint4 divisorLimbs0 = divisorPtr[i*4];
      uint4 divisorLimbs1 = divisorPtr[i*4+1];
      uint4 divisorLimbs2 = divisorPtr[i*4+2];
      uint4 divisorLimbs3 = divisorPtr[i*4+3];
      
      uint4 moduloLimb0;
      uint4 moduloLimb1;
      uint4 moduloLimb2;
      uint4 moduloLimb3;
      
      modulo640to512(dividendLimbs0, dividendLimbs1, dividendLimbs2, dividendLimbs3, dividendLimbs4,
                     divisorLimbs0, divisorLimbs1, divisorLimbs2, divisorLimbs3,
                     &moduloLimb0, &moduloLimb1, &moduloLimb2, &moduloLimb3);
      
      modulosPtr[i*4] = moduloLimb0;
      modulosPtr[i*4+1] = moduloLimb1;
      modulosPtr[i*4+2] = moduloLimb2;
      modulosPtr[i*4+3] = moduloLimb3;
    }
  }
}

__kernel void searchNonce(__constant uint4 *block,
                          __global struct GPUNonceAndHash *nonceAndHash,
                          __constant uint32_t *primes,
                          __global uint64_t *multipliers64,
                          __global uint32_t *offsets64)
{
  __local uint4 hash[2*GroupSize];
  __local uint8_t hasModulo[GroupSize];
  __global struct GPUNonceAndHash *context = &nonceAndHash[get_group_id(0)];
  uint32_t nonceIdx = context->currentNonce + 1;

  barrier(CLK_GLOBAL_MEM_FENCE);
  if (nonceIdx >= context->totalNonces) {
    uint4 m0 = block[0];
    uint4 m1 = block[1];
    uint4 m2 = block[2];
    uint4 m3 = block[3];
    uint4 m4 = block[4];
    
    sha256SwapByteOrder(&m0);
    sha256SwapByteOrder(&m1);
    sha256SwapByteOrder(&m2);
    sha256SwapByteOrder(&m3);
    sha256SwapByteOrder(&m4);
    
    uint4 hash1l0, hash1l1;
    uint4 hash2l0, hash2l1;
    uint4 targetHashl0, targetHashl1;
    uint32_t targetNonce;

    uint32_t currentNonce = context->nonce[nonceIdx-1];
    currentNonce = (currentNonce == 0) ?
      get_global_id(0) : currentNonce - (currentNonce % GroupSize) + get_global_size(0) + get_local_id(0);
    unsigned trialDivisionPassedCounter = 0;
    
    while (trialDivisionPassedCounter < GroupSize) {
      m4.w = EndianSwap(currentNonce);    
      SHA256_fresh(&hash1l0, &hash1l1, m0, m1, m2, m3);
      sha256(&hash1l0, &hash1l1,
             m4, (uint4){0x80000000, 0, 0, 0}, (uint4){0, 0, 0, 0}, (uint4){0, 0, 0, 0x00000280});        
      SHA256_fresh(&hash2l0, &hash2l1,
                   hash1l0, hash1l1, (uint4){0x80000000, 0, 0, 0}, (uint4){0, 0, 0, 0x00000100});
      sha256SwapByteOrder(&hash2l0);
      sha256SwapByteOrder(&hash2l1);
      barrier(CLK_LOCAL_MEM_FENCE);      
      hash[get_local_id(0)*2 + 0] = hash2l0;
      hash[get_local_id(0)*2 + 1] = hash2l1;    
      barrier(CLK_LOCAL_MEM_FENCE);        
      
      for (unsigned i = 0; (i < get_local_size(0)) && (trialDivisionPassedCounter < GroupSize); i++) {
        hash1l0 = hash[i*2 + 0];
        hash1l1 = hash[i*2 + 1];
        if (!((hash1l1.w >> 31) & hash1l0.x))
          continue;
        
        unsigned passed = 1;    
        for (unsigned j = 0, divIdx = get_local_id(0)+1; divIdx < 1024+1; j++, divIdx += GroupSize) {
          barrier(CLK_LOCAL_MEM_FENCE);      
          hasModulo[get_local_id(0)] = (longModuloByMul256(hash1l0,
                                                           hash1l1,
                                                           primes[divIdx],
                                                           multipliers64[divIdx],
                                                           offsets64[divIdx]) == 0);
          barrier(CLK_LOCAL_MEM_FENCE);              
          uint4 M = {0, 0, 0, 0};
          for (unsigned k = 0; k < GroupSize/16; k++)
            M |= ((__local uint4*)hasModulo)[k];
          barrier(CLK_LOCAL_MEM_FENCE);               
          M.xy |= M.zw;
          M.x |= M.y;
          if (M.x) {
            passed = 0;
            break;
          }
        }
        
        if (passed) {
          if (get_local_id(0) == trialDivisionPassedCounter) {
            targetHashl0 = hash1l0;
            targetHashl1 = hash1l1;
            targetNonce = currentNonce - (currentNonce % GroupSize) + i;
          }
          trialDivisionPassedCounter++;
        }
      }
      
      currentNonce += get_global_size(0);
    }
    barrier(CLK_LOCAL_MEM_FENCE);               
    
    unsigned isPrime;
    {
      uint4 modPowL0, modPowL1;
      FermatTest256(targetHashl0, targetHashl1, &modPowL0, &modPowL1);    
      --modPowL0.x;
      modPowL0 |= modPowL1;
      modPowL0.xy |= modPowL1.zw;
      modPowL0.x |= modPowL0.y;    
      isPrime = modPowL0.x == 0;
    }
    
    hasModulo[get_local_id(0)] = isPrime;
    barrier(CLK_LOCAL_MEM_FENCE);

    uint32_t index, sum = 0;
    unsigned threadGrId = get_local_id(0) / 4;
    for (unsigned i = 0; i < GroupSize / 4; i++)
      sum += (i < threadGrId) ? ((__local uint32_t*)hasModulo)[i] : 0;
    
    index = (sum & 0xFF) + ((sum >> 8) & 0xFF) + ((sum >> 16) & 0xFF) + (sum >> 24);
    index += (get_local_id(0) % 4 >= 1) ? hasModulo[get_local_id(0) - 1] : 0;
    index += (get_local_id(0) % 4 >= 2) ? hasModulo[get_local_id(0) - 2] : 0;
    index += (get_local_id(0) % 4 >= 3) ? hasModulo[get_local_id(0) - 3] : 0;
    
    if (isPrime) {
      context->hash[index*2] = targetHashl0;
      context->hash[index*2 + 1] = targetHashl1;
      context->nonce[index] = targetNonce;
    }
    
    nonceIdx = 0;
    context->currentNonce = 0;
    if (get_local_id(0) == GroupSize - 1)
      context->totalNonces = index;
    barrier(CLK_GLOBAL_MEM_FENCE);           
  } else {
    context->currentNonce = nonceIdx;
  }
}


__kernel void sieve(__global uint32_t *cunningham1Bitfield,
                    __global uint32_t *cunningham2Bitfield,
                    __global uint32_t *bitwinBitfield,
                    __constant uint4 *primorial,
                    __global struct GPUNonceAndHash *nonceAndHash,
                    __constant uint32_t *primes,
                    __global uint64_t *multipliers64,
                    __global uint32_t *offsets64)
{
  __local uint8_t localCunningham1[L1CacheSize];
  __local uint8_t localCunningham2[L1CacheSize];  
  __global struct GPUNonceAndHash *context = &nonceAndHash[get_group_id(0)];
  uint32_t nonceIdx = context->currentNonce;
    
  uint4 M0, M1, M2, M3;
  mul256schoolBook_v3(context->hash[nonceIdx*2], context->hash[nonceIdx*2+1],
                      primorial[0], primorial[1],
                      &M0, &M1, &M2, &M3);
  
  weave(M0, M1, M2,
        cunningham1Bitfield + get_group_id(0) * MaxSieveBufferSize/32,
        cunningham2Bitfield + get_group_id(0) * MaxSieveBufferSize/32,
        bitwinBitfield + get_group_id(0) * MaxSieveBufferSize/32,
        localCunningham1,
        localCunningham2,
        primes,
        multipliers64,
        offsets64,
        get_local_size(0) - GroupSize + FixedRoundsNum);
}

void mul384_1(uint4 l0, uint4 l1, uint4 l2, uint32_t m,
              uint4 *r0, uint4 *r1, uint4 *r2)
{
  *r0 = l0 * m;
  *r1 = l1 * m;
  *r2 = l2 * m;
  
  uint4 h0 = mul_hi(l0, m);
  uint4 h1 = mul_hi(l1, m);
  uint4 h2 = mul_hi(l2, m);
  
  add384(r0, r1, r2,
         (uint4){0, h0.x, h0.y, h0.z},
         (uint4){h0.w, h1.x, h1.y, h1.z},
         (uint4){h1.w, h2.x, h2.y, h2.z});
}


unsigned extractMultipliers2(__global struct GPUNonceAndHash *sieve,
                             __constant uint4 *primorial,
                             __global uint32_t *ptr,
                             __local uint32_t *multipliersPerThread,
                             __local uint32_t *multipliersNum,
                             __global struct FermatQueue *queue,
                             unsigned modifier,
                             unsigned *newSize)
{
  uint32_t localMultipliers[32];
  uint32_t localIndex = 0;
  
  unsigned c32 = get_local_size(0)/8;
  unsigned cExt = c32 - (32 - ExtensionsNum);  
  
  unsigned i = get_local_id(0);
  unsigned sieveWords = FixedSieveSize/32;
  for (unsigned extNum = 0; extNum <= cExt; extNum++) {
    __global uint32_t *lPtr = ptr + extNum*sieveWords;
    while (i < sieveWords) {
      uint32_t word = lPtr[i];
      if (word == 0) {
        i += GroupSize;
        continue;
      }
      
      for (unsigned j = 0; j < c32; j++, word >>= 1) {
        if (word & 0x1) {
          unsigned M = 8*L1CacheSize*(i*4/L1CacheSize) + (i*4 + j/8) % L1CacheSize + (j&0x7)*L1CacheSize;
          localMultipliers[localIndex++] = M << extNum;
        }
      }
      
      i += GroupSize;
    }
    
    i = get_local_id(0) + sieveWords / 2;    
  }
  
  multipliersPerThread[get_local_id(0)] = localIndex;
  barrier(CLK_LOCAL_MEM_FENCE);
  
  unsigned globalIndex = 0;
  for (unsigned i = 0; i < GroupSize; i++)
    globalIndex += ((i < get_local_id(0)) ? multipliersPerThread[i] : 0);

  if (get_local_id(0) == GroupSize-1)
    *multipliersNum = globalIndex + multipliersPerThread[GroupSize-1];
  barrier(CLK_LOCAL_MEM_FENCE);
  
  unsigned nonceIdx = sieve->currentNonce;
  uint4 hashl0 = sieve->hash[2*nonceIdx];
  uint4 hashl1 = sieve->hash[2*nonceIdx+1];
  uint32_t nonce = sieve->nonce[nonceIdx];

  uint32_t position = queue->position;
  uint32_t size = queue->size;
  
  uint4 m0, m1, m2, m3;
  globalIndex += position + size;
  mul256schoolBook_v3(hashl0, hashl1, primorial[0], primorial[1], &m0, &m1, &m2, &m3);
  for (unsigned i = 0; i < localIndex; i++) {
    uint4 chOrl0, chOrl1, chOrl2;
    unsigned bufferIdx = (globalIndex+i) % FermatQueueBufferSize;

    mul384_1(m0, m1, m2, localMultipliers[i], &chOrl0, &chOrl1, &chOrl2);
    chOrl0.x += modifier;
    queue->chainOrigins[3*bufferIdx] = chOrl0;
    queue->chainOrigins[3*bufferIdx + 1] = chOrl1;
    queue->chainOrigins[3*bufferIdx + 2] = chOrl2;
    
    queue->multipliers[bufferIdx] = localMultipliers[i];
    queue->chainLengths[bufferIdx] = 0;
    queue->nonces[bufferIdx] = nonce;
  }

  *newSize = size + *multipliersNum;
  return position;
}





void doFermatTestC12(__global struct GPUNonceAndHash *context,
                     __global struct FermatQueue *groupQueue,
                     __global struct FermatTestResults *groupResults,
                     __global uint32_t *bitfield,
                     __local uint32_t *multipliersPerThread,
                     __local uint8_t *isPrime,   
                     __local uint8_t *isResult,
                     __local uint32_t *multipliersNum,
                     __local uint32_t *primesNum,
                     __local uint32_t *resultsNum,
                     __constant uint4 *primorial,
                     unsigned type,
                     unsigned *outputSize)
  
{
  uint32_t queueSize;
  uint32_t position = extractMultipliers2(context, primorial, bitfield,
                                          multipliersPerThread, multipliersNum,
                                          groupQueue, (type == 1 ? -1 : 1), &queueSize);
  barrier(CLK_GLOBAL_MEM_FENCE);
  while (queueSize >= GroupSize) {
    unsigned bufferIdx = (position + get_local_id(0)) % FermatQueueBufferSize ;
    uint4 chl0 = groupQueue->chainOrigins[3*bufferIdx];
    uint4 chl1 = groupQueue->chainOrigins[3*bufferIdx + 1];
    uint4 chl2 = groupQueue->chainOrigins[3*bufferIdx + 2];
    uint32_t chainLength = groupQueue->chainLengths[bufferIdx];
    uint32_t multiplier = groupQueue->multipliers[bufferIdx];    
    uint32_t nonce = groupQueue->nonces[bufferIdx];
    unsigned chainContinue;
    
    uint4 modpowl0, modpowl1, modpowl2;
    FermatTest384(chl0, chl1, chl2, &modpowl0, &modpowl1, &modpowl2);
    --modpowl0.x;
    modpowl0 |= modpowl1;
    modpowl0 |= modpowl2;
    modpowl0.xy |= modpowl0.zw;
    modpowl0.x |= modpowl0.y;
    chainContinue = !modpowl0.x;
    
    isPrime[get_local_id(0)] = chainContinue;
    isResult[get_local_id(0)] = !chainContinue && chainLength >= 1;
    barrier(CLK_LOCAL_MEM_FENCE);
    
    lshift3(&chl0, &chl1, &chl2, 1);
    chl0.x += (type == 1 ? 1 : -1);
    
    unsigned primeSum = 0, primeIndex;
    unsigned resultSum = 0, resultIndex;
    unsigned threadGrId = get_local_id(0) / 4;
    for (unsigned i = 0; i < GroupSize / 4; i++) {
      primeSum += ((i < threadGrId) ? ((__local uint32_t*)isPrime)[i] : 0);
      resultSum += ((i < threadGrId) ? ((__local uint32_t*)isResult)[i] : 0);
    }
    
    primeIndex = (primeSum & 0xFF) + ((primeSum >> 8) & 0xFF) + ((primeSum >> 16) & 0xFF) + (primeSum >> 24);
    primeIndex += (get_local_id(0) % 4 >= 1) ? isPrime[get_local_id(0) - 1] : 0;
    primeIndex += (get_local_id(0) % 4 >= 2) ? isPrime[get_local_id(0) - 2] : 0;
    primeIndex += (get_local_id(0) % 4 >= 3) ? isPrime[get_local_id(0) - 3] : 0;
    
    resultIndex = (resultSum & 0xFF) + ((resultSum >> 8) & 0xFF) + ((resultSum >> 16) & 0xFF) + (resultSum >> 24);
    resultIndex += (get_local_id(0) % 4 >= 1) ? isResult[get_local_id(0) - 1] : 0;
    resultIndex += (get_local_id(0) % 4 >= 2) ? isResult[get_local_id(0) - 2] : 0;
    resultIndex += (get_local_id(0) % 4 >= 3) ? isResult[get_local_id(0) - 3] : 0;
    
    if (get_local_id(0) == GroupSize-1) {
      *primesNum = primeIndex + isPrime[GroupSize-1];
      *resultsNum = resultIndex + isResult[GroupSize-1];
    }
    barrier(CLK_LOCAL_MEM_FENCE);
    
    if (chainContinue) {
      chainLength++;
      primeIndex = (primeIndex + position + queueSize) % FermatQueueBufferSize;
      groupQueue->chainOrigins[3*primeIndex] = chl0;
      groupQueue->chainOrigins[3*primeIndex+1] = chl1;
      groupQueue->chainOrigins[3*primeIndex+2] = chl2;
      groupQueue->chainLengths[primeIndex] = chainLength;
      groupQueue->nonces[primeIndex] = nonce;
      groupQueue->multipliers[primeIndex] = multiplier;
    } else if (chainLength >= 1) {
      resultIndex += *outputSize;
      groupResults->resultTypes[resultIndex] = type;
      groupResults->resultMultipliers[resultIndex] = multiplier;
      groupResults->resultChainLength[resultIndex] = chainLength;
      groupResults->resultNonces[resultIndex] = nonce;
    }    
    
    *outputSize += *resultsNum;
    position += GroupSize;
    queueSize -= (GroupSize - *primesNum);
  }

  groupQueue->position = position % FermatQueueBufferSize; 
  groupQueue->size = queueSize;
}

unsigned bitwinChainLength(uint32_t C1, uint32_t C2)
{
  return (C1 > C2) ? 2*C2 + 1 : 2*C1;
}

void doFermatTestBt(__global struct GPUNonceAndHash *context,
                    __global struct FermatQueue *groupQueue,
                    __global struct FermatTestResults *groupResults,
                    __global uint32_t *bitfield,
                    __local uint32_t *multipliersPerThread,
                    __local uint8_t *isPrime,   
                    __local uint8_t *isResult,
                    __local uint32_t *multipliersNum,
                    __local uint32_t *primesNum,
                    __local uint32_t *resultsNum,
                    __constant uint4 *primorial,
                    unsigned *outputSize)

{
  uint32_t queueSize;
  uint32_t position = extractMultipliers2(context, primorial, bitfield,
                                          multipliersPerThread, multipliersNum,
                                          groupQueue, -1, &queueSize);
  
  barrier(CLK_GLOBAL_MEM_FENCE);
  while (queueSize >= GroupSize) {
    unsigned bufferIdx = (position + get_local_id(0)) % FermatQueueBufferSize;
    uint4 chl0 = groupQueue->chainOrigins[3*bufferIdx];
    uint4 chl1 = groupQueue->chainOrigins[3*bufferIdx + 1];
    uint4 chl2 = groupQueue->chainOrigins[3*bufferIdx + 2];
    uint32_t chainLength = groupQueue->chainLengths[bufferIdx];
    uint32_t multiplier = groupQueue->multipliers[bufferIdx];    
    uint32_t nonce = groupQueue->nonces[bufferIdx];
    uint32_t c1ChainLength = chainLength >> 16;
    
    unsigned chainContinue;
    unsigned usedChainLength;
    unsigned switchedType;
    
    uint4 modpowl0, modpowl1, modpowl2;
    FermatTest384(chl0, chl1, chl2, &modpowl0, &modpowl1, &modpowl2);
    --modpowl0.x;
    modpowl0 |= modpowl1;
    modpowl0 |= modpowl2;
    modpowl0.xy |= modpowl0.zw;
    modpowl0.x |= modpowl0.y;
    
    chainContinue = !modpowl0.x ? 1 : (!c1ChainLength && chainLength >= 2);
    switchedType = modpowl0.x ? chainContinue : 0;
    usedChainLength = chainContinue ? 0 : bitwinChainLength(c1ChainLength, chainLength & 0xFFFF);

    isPrime[get_local_id(0)] = chainContinue;
    isResult[get_local_id(0)] = usedChainLength >= 1;
    barrier(CLK_LOCAL_MEM_FENCE);
    
    if (switchedType) {
      // 384-bit right shift more faster than fetch from global memory (?)
      rshift3(&chl0, &chl1, &chl2, chainLength & 0xFF);
      chl0.x += 2;
      chainLength <<= 16;
    } else {
      lshift3(&chl0, &chl1, &chl2, 1);
      chl0.x += !c1ChainLength ? 1 : -1;
      chainLength++;
    }
    
    unsigned primeSum = 0, primeIndex;
    unsigned resultSum = 0, resultIndex;
    unsigned threadGrId = get_local_id(0) / 4;
    for (unsigned i = 0; i < GroupSize / 4; i++) {
      primeSum += ((i < threadGrId) ? ((__local uint32_t*)isPrime)[i] : 0);
      resultSum += ((i < threadGrId) ? ((__local uint32_t*)isResult)[i] : 0);
    }
    
    primeIndex = (primeSum & 0xFF) + ((primeSum >> 8) & 0xFF) + ((primeSum >> 16) & 0xFF) + (primeSum >> 24);
    primeIndex += (get_local_id(0) % 4 >= 1) ? isPrime[get_local_id(0) - 1] : 0;
    primeIndex += (get_local_id(0) % 4 >= 2) ? isPrime[get_local_id(0) - 2] : 0;
    primeIndex += (get_local_id(0) % 4 >= 3) ? isPrime[get_local_id(0) - 3] : 0;
    
    resultIndex = (resultSum & 0xFF) + ((resultSum >> 8) & 0xFF) + ((resultSum >> 16) & 0xFF) + (resultSum >> 24);
    resultIndex += (get_local_id(0) % 4 >= 1) ? isResult[get_local_id(0) - 1] : 0;
    resultIndex += (get_local_id(0) % 4 >= 2) ? isResult[get_local_id(0) - 2] : 0;
    resultIndex += (get_local_id(0) % 4 >= 3) ? isResult[get_local_id(0) - 3] : 0;
    
    if (get_local_id(0) == GroupSize-1) {
      *primesNum = primeIndex + isPrime[GroupSize-1];
      *resultsNum = resultIndex + isResult[GroupSize-1];
    }
    barrier(CLK_LOCAL_MEM_FENCE);
    
    if (chainContinue) {
      primeIndex = (primeIndex + position + queueSize) % FermatQueueBufferSize;
      groupQueue->chainOrigins[3*primeIndex] = chl0;
      groupQueue->chainOrigins[3*primeIndex+1] = chl1;
      groupQueue->chainOrigins[3*primeIndex+2] = chl2;
      groupQueue->chainLengths[primeIndex] = chainLength;
      groupQueue->nonces[primeIndex] = nonce;
      groupQueue->multipliers[primeIndex] = multiplier;
    } else if (usedChainLength >= 1) {
      resultIndex += *outputSize;
      groupResults->resultTypes[resultIndex] = 3;
      groupResults->resultMultipliers[resultIndex] = multiplier;
      groupResults->resultChainLength[resultIndex] = usedChainLength;
      groupResults->resultNonces[resultIndex] = nonce;
    }    
    
    *outputSize += *resultsNum;
    position += GroupSize;
    queueSize -= (GroupSize - *primesNum);
  }
  
  groupQueue->position = position % FermatQueueBufferSize; 
  groupQueue->size = queueSize;
}

__kernel void FermatTestEnqueue(__global struct GPUNonceAndHash *nonceAndHash,
                                __constant uint4 *primorial,
                                __global uint32_t *cunningham1Bitfield,
                                __global uint32_t *cunningham2Bitfield,
                                __global uint32_t *bitwinBitfield,
                                __global struct FermatQueue *c1Queue,
                                __global struct FermatQueue *c2Queue,
                                __global struct FermatQueue *btQueue,
                                __global struct FermatTestResults *results)
{
  __global struct GPUNonceAndHash *context = &nonceAndHash[get_group_id(0)];  
  __global struct FermatQueue *c1GroupQueue = &c1Queue[get_group_id(0)];
  __global struct FermatQueue *c2GroupQueue = &c2Queue[get_group_id(0)];
  __global struct FermatTestResults *groupResults = &results[get_group_id(0)];
  __global uint32_t *cunningham1Ptr = cunningham1Bitfield + get_group_id(0) * MaxSieveBufferSize/32;
  __global uint32_t *cunningham2Ptr = cunningham2Bitfield + get_group_id(0) * MaxSieveBufferSize/32;

  __local uint32_t multipliersPerThread[GroupSize];
  __local uint8_t isPrime[GroupSize];
  __local uint8_t isResult[GroupSize];  
  __local uint32_t multipliersNum;
  __local uint32_t primesNum;
  __local uint32_t resultsNum;
  
  unsigned outputSize = 0;
  doFermatTestC12(context, c1GroupQueue, groupResults, cunningham1Ptr,
                  multipliersPerThread, isPrime, isResult,
                  &multipliersNum, &primesNum, &resultsNum,
                  primorial, 1, &outputSize);
  
  doFermatTestC12(context, c2GroupQueue, groupResults, cunningham2Ptr,
                  multipliersPerThread, isPrime, isResult,
                  &multipliersNum, &primesNum, &resultsNum,
                  primorial, 2, &outputSize);

  groupResults->size = outputSize;
}

__kernel void FermatTestEnqueueBt(__global struct GPUNonceAndHash *nonceAndHash,
                                  __constant uint4 *primorial,
                                  __global uint32_t *cunningham1Bitfield,
                                  __global uint32_t *cunningham2Bitfield,
                                  __global uint32_t *bitwinBitfield,
                                  __global struct FermatQueue *c1Queue,
                                  __global struct FermatQueue *c2Queue,
                                  __global struct FermatQueue *btQueue,
                                  __global struct FermatTestResults *results)
{
  __global struct GPUNonceAndHash *context = &nonceAndHash[get_group_id(0)];  
  __global struct FermatQueue *btGroupQueue = &btQueue[get_group_id(0)];
  __global struct FermatTestResults *groupResults = &results[get_group_id(0)];
  __global uint32_t *bitwinPtr = bitwinBitfield + get_group_id(0) * MaxSieveBufferSize/32;
  
  __local uint32_t multipliersPerThread[GroupSize];
  __local uint8_t isPrime[GroupSize];
  __local uint8_t isResult[GroupSize];  
  __local uint32_t multipliersNum;
  __local uint32_t primesNum;
  __local uint32_t resultsNum;
  
  unsigned outputSize = groupResults->size;
  
  doFermatTestBt(context, btGroupQueue, groupResults, bitwinPtr,
                 multipliersPerThread, isPrime, isResult,
                 &multipliersNum, &primesNum, &resultsNum,
                 primorial, &outputSize);
  
  groupResults->size = outputSize;
}

__kernel void empty()
{
}