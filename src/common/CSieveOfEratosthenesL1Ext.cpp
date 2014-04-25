#include "CSieveOfEratosthenesL1Ext.h"
#include "config.h"
#include "system.h"
#include <limits>
#include <vector>

static inline unsigned calculateOffset(unsigned lowIdx,
                                       unsigned currentPrime,
                                       unsigned currentPrimeMod,
                                       unsigned offset)
{
  return offset >= currentPrimeMod ?
    lowIdx - currentPrimeMod + offset :
    lowIdx + (currentPrime - currentPrimeMod) + offset;
}

// Extended Euclidian algorithm, get inverse modulo
static unsigned intInvert(unsigned a, unsigned mod)
{
  int rem0 = mod, rem1 = a % mod, rem2;
  int aux0 = 0, aux1 = 1, aux2;
  int quotient, inverse;
  
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

void CSieveOfEratosthenesL1Ext::reset(unsigned sieveSize,
                                      unsigned chainLength,
                                      unsigned depth,
                                      const uint32_t *fixedMultiplier)
{
  this->nSieveSize = sieveSize;
  this->chainLength = chainLength;
  this->depth = depth;
  memcpy(fixedFactor, fixedMultiplier, fixedMultiplier[0]*4 + 4);
  sieveWords = nSieveSize / 64;
  layersNum = chainLength + ExtensionsNum;
  
  memset(cunningham1Bitfield, 0, sieveWords*8);
  memset(cunningham2Bitfield, 0, sieveWords*8);
  memset(bitwinBitfield, 0, sieveWords*8);
  for (unsigned i = 1; i <= ExtensionsNum; i++) {
    memset(cunningham1Bitfield + i*sieveWords + sieveWords/2, 0, (sieveWords/2)*8);
    memset(cunningham2Bitfield + i*sieveWords + sieveWords/2, 0, (sieveWords/2)*8);
    memset(bitwinBitfield + i*sieveWords + sieveWords/2, 0, (sieveWords/2)*8);
  }
}

void CSieveOfEratosthenesL1Ext::reset(unsigned sieveSize,
                                      unsigned chainLength,
                                      unsigned depth,
                                      const mpz_class &fixedMultiplier)
{
  size_t limbsNumber;
  this->nSieveSize = sieveSize;
  this->chainLength = chainLength;
  this->depth = depth;
  mpzExport(fixedMultiplier, fixedFactor);
  sieveWords = nSieveSize / 64;
  layersNum = chainLength + ExtensionsNum;
  
  memset(cunningham1Bitfield, 0, sieveWords*8);
  memset(cunningham2Bitfield, 0, sieveWords*8);
  memset(bitwinBitfield, 0, sieveWords*8);
  for (unsigned i = 1; i <= ExtensionsNum; i++) {
    memset(cunningham1Bitfield + i*sieveWords + sieveWords/2, 0, (sieveWords/2)*8);
    memset(cunningham2Bitfield + i*sieveWords + sieveWords/2, 0, (sieveWords/2)*8);
    memset(bitwinBitfield + i*sieveWords + sieveWords/2, 0, (sieveWords/2)*8);
  }
}


unsigned int CSieveOfEratosthenesL1Ext::GetCandidateCount()
{
  unsigned candidateCount = 0;
  
  // Main range
  for (unsigned i = 0; i < sieveWords; i++) {
    uint64_t value =
      ~(cunningham1Bitfield[i] & cunningham2Bitfield[i] & bitwinBitfield[i]);
    if (!value)
      continue;
    for (unsigned j = 0; j < 64; j++, value >>= 1)
      candidateCount += (value & 0x1);
  }  
  
  // Extensions
  uint64_t *extCunningham1 = cunningham1Bitfield + sieveWords + sieveWords/2;
  uint64_t *extCunningham2 = cunningham2Bitfield + sieveWords + sieveWords/2;
  uint64_t *extBitwin = bitwinBitfield + sieveWords + sieveWords/2;    
  for (unsigned extNum = 0; extNum < ExtensionsNum; extNum++) {
    for (unsigned i = 0; i < sieveWords/2; i++) {
      uint64_t value =
        ~(extCunningham1[i]) | (~extCunningham2[i]) | (~extBitwin[i]);
      if (!value)
        continue;
      for (unsigned j = 0; j < 64; j++, value >>= 1)
        candidateCount += (value & 0x1);
    }
    
    extCunningham1 += sieveWords;
    extCunningham2 += sieveWords;
    extBitwin += sieveWords;
  }

  return candidateCount;
}

bool CSieveOfEratosthenesL1Ext::GetNextCandidateMultiplier(unsigned &nVariableMultiplier,
                                                          unsigned &nCandidateType)
{
  do {
    if (word) {
      for (; bitOffset < 64; bitOffset++, word >>= 1) {
        if (word & 0x1) {
          nCandidateType = candidateType;
          nVariableMultiplier = (wordIdx*64 + bitOffset) * (1 << extensionNum);
          bitOffset++;
          word >>= 1;
          return true;
        }
      }
    }
  } while (increment());
  
  return false;
}


static uint32_t longModuloByMul(uint32_t *number,
                                uint32_t divider,
                                uint64_t multiplier,
                                unsigned shift)
{
  uint64_t mod = 0;
  size_t limbsNum = number[0];
  for (uint32_t i = limbsNum; i > 0; i--) {
#if (INT128_SIZE == 16)
    __int128 dividend = (mod << 32) + (uint64_t)number[i];
    uint64_t quote = (dividend*multiplier) >> shift;    
#else
    uint64_t dividend = (mod << 32) + (uint64_t)number[i];
    uint64_t quote = dividend / (uint64_t)divider;    
#endif
    mod = (uint64_t)dividend - quote*divider;
  }
  
  return mod;  
}

void CSieveOfEratosthenesL1Ext::Weave()
{
#ifdef DEBUG_MINING  
  timeMark beginPoint = getTimeMark();
#endif
  
  unsigned roundsNum = nSieveSize / L1CacheSize;  
  uint64_t a_mod_xy;
  for (unsigned i = 0; i < depth; i++) {
    unsigned currentPrime = _primeSource[i];
    if (i == 0 || _primeSource.isNewMultiplier(i)) {
      // a_mod_xy = bnFixedFactor % __primeSource.primesCombined(i)
      a_mod_xy = longModuloByMul(fixedFactor,
                                 _primeSource.primesCombined(i),
                                 _primeSource.combinedMultiplier(i),
                                 _primeSource.combinedOffset(i));
    }
    
    // unsigned a_mod_xy_mod_x = a_mod_xy % currentPrime;
    unsigned qtmp = (a_mod_xy * _primeSource.multiplier(i)) >> _primeSource.offset(i);
    unsigned a_mod_xy_mod_x = a_mod_xy - (qtmp*currentPrime);
    unsigned inverseModulo = a_mod_xy_mod_x ? 
      intInvert(a_mod_xy_mod_x, currentPrime) : 0;
    
    if (inverseModulo) {
      for (unsigned layer = 0; layer < chainLength; layer++) {
        cunningham1Multipliers[layersNum * i + layer] = inverseModulo;
        cunningham2Multipliers[layersNum * i + layer] = currentPrime - inverseModulo;
        inverseModulo = (inverseModulo & 0x1) ?
          (inverseModulo + currentPrime) / 2 : inverseModulo / 2;
      }
      
      unsigned extensionLowIdx = L1CacheSize*(roundsNum/2);
      for (unsigned layer = chainLength; layer < layersNum; layer++) {
        // unsigned currentPrimeMod = extensionLowIdx % currentPrime;
        qtmp = ((uint64_t)extensionLowIdx * _primeSource.multiplier(i)) >> _primeSource.offset(i);  
        unsigned currentPrimeMod = extensionLowIdx - (qtmp*currentPrime);

        cunningham1Multipliers[layersNum * i + layer] =
          calculateOffset(extensionLowIdx, currentPrime, currentPrimeMod, inverseModulo);
        cunningham2Multipliers[layersNum * i + layer] =
          calculateOffset(extensionLowIdx, currentPrime, currentPrimeMod, currentPrime - inverseModulo);
        inverseModulo = (inverseModulo & 0x1) ?
          (inverseModulo + currentPrime) / 2 : inverseModulo / 2;        
      }
    } else {
      memset(&cunningham1Multipliers[layersNum * i], 0, 4*layersNum);
      memset(&cunningham2Multipliers[layersNum * i], 0, 4*layersNum); 
    }
  }
  
#ifdef DEBUG_MINING  
  timeMark prepareEndPoint = getTimeMark();
#endif

  for (unsigned round = 0; round < roundsNum; round++) {
    unsigned lowIdx = L1CacheSize * round;
    for (unsigned layer = 0; layer < layersNum; layer++) {
      if (layer >= chainLength && round < roundsNum/2)
        break;

      // Build layer for all primes from 2 to _primeSource[depth]
      memset(&local, 0, sizeof(local));
      for (unsigned primeIdx = 0; primeIdx < depth; primeIdx++) {
        unsigned currentPrime = _primeSource[primeIdx];
        unsigned offset;

        offset = cunningham1Multipliers[layersNum * primeIdx + layer];
        if (!offset)
          continue;

        offset -= lowIdx;
        while (offset < L1CacheSize) {
          setValue(local.cunningham1, offset);
          offset += currentPrime;          
        }
        cunningham1Multipliers[layersNum * primeIdx + layer] = offset + lowIdx;
        
        offset = cunningham2Multipliers[layersNum * primeIdx + layer] - lowIdx;
        while (offset < L1CacheSize) {
          setValue(local.cunningham2, offset);
          offset += currentPrime;          
        }
        cunningham2Multipliers[layersNum * primeIdx + layer] = offset + lowIdx;
      }

      // Map layer to main bitfield [N..xN]
      uint64_t *cunningham1 = cunningham1Bitfield + lowIdx/64;
      uint64_t *cunningham2 = cunningham2Bitfield + lowIdx/64;
      uint64_t *bitwin = bitwinBitfield + lowIdx/64;
      
      if (layer < chainLength/2) {
        for (unsigned i = 0; i < L1CacheSize/64; i++) {
          cunningham1[i] |= local.cunningham1[i];
          cunningham2[i] |= local.cunningham2[i];
          bitwin[i] |= local.cunningham1[i] | local.cunningham2[i];
        }
      } else if (layer < (chainLength + 1)/2) {
        for (unsigned i = 0; i < L1CacheSize/64; i++) {
          cunningham1[i] |= local.cunningham1[i];
          cunningham2[i] |= local.cunningham2[i];
          bitwin[i] |= local.cunningham1[i];
        }
      } else if (layer < chainLength) {
        for (unsigned i = 0; i < L1CacheSize/64; i++) {
          cunningham1[i] |= local.cunningham1[i];
          cunningham2[i] |= local.cunningham2[i];
        }
      }

      // Map layer to extensions [2N..2xN], [4N..4xN], ...
      // For extensions used range [N/2..N] instead of [1..N]
      if (round < roundsNum/2)
        continue;
      
      // ext[0] <- layers [0..сhainLength]
      // ext[1] <- layers [1..chainLength+1]
      // ...
      // ext[M] <- layers [M..chainLength+M]
      //
      // Detect extensions, affected by current layer
      for (unsigned extNum =
             layer < chainLength ? 1 : layer - chainLength + 1;
           extNum <= std::min(layer, (unsigned)ExtensionsNum);
           extNum++) {
        uint64_t *extCunningham1 = cunningham1Bitfield + extNum*sieveWords + lowIdx/64;
        uint64_t *extCunningham2 = cunningham2Bitfield + extNum*sieveWords + lowIdx/64;
        uint64_t *extBitwin = bitwinBitfield + extNum*sieveWords + lowIdx/64;
      
        if (layer - extNum < chainLength/2) {
          for (unsigned i = 0; i < L1CacheSize/64; i++) {
            extCunningham1[i] |= local.cunningham1[i];
            extCunningham2[i] |= local.cunningham2[i];
            extBitwin[i] |= local.cunningham1[i] | local.cunningham2[i];
          }
        } else if (layer - extNum < (chainLength+1)/2) {
          for (unsigned i = 0; i < L1CacheSize/64; i++) {
            extCunningham1[i] |= local.cunningham1[i];
            extCunningham2[i] |= local.cunningham2[i];
            extBitwin[i] |= local.cunningham1[i];
          }          
        } else {
          for (unsigned i = 0; i < L1CacheSize/64; i++) {
            extCunningham1[i] |= local.cunningham1[i];
            extCunningham2[i] |= local.cunningham2[i];
          }          
        }
      }
    }
  }

#ifdef DEBUG_MINING
  timeMark endPoint = getTimeMark();
  fprintf(stderr,
          "L1Ext: prepare %.3lfmsec, sieve %.3lfmsec\n",
          usDiff(beginPoint, prepareEndPoint) / 1000.0,
          usDiff(prepareEndPoint, endPoint) / 1000.0);
#endif  
}

bool CSieveOfEratosthenesL1Ext::fastSelfTest()
{
  mpz_class fixedMultiplier;
  mpzImport(fixedFactor, fixedMultiplier);

  resetCandidateIterator();
  unsigned multiplier;
  unsigned type;
  unsigned fakePrimes = 0;
  while (GetNextCandidateMultiplier(multiplier, type)) {
    mpz_class chainOrigin = fixedMultiplier*multiplier;
    if (type == PRIME_CHAIN_CUNNINGHAM1) {
      if (!trialDivisionChainTest(_primeSource, chainOrigin, true, chainLength, depth))
        return false;
    } else if (type == PRIME_CHAIN_CUNNINGHAM2) {
      if (!trialDivisionChainTest(_primeSource, chainOrigin, false, chainLength, depth))
        return false;
    } else {
      mpz_class chainOriginExtra = chainOrigin;
      if (!trialDivisionChainTest(_primeSource, chainOrigin, true, (chainLength+1)/2, depth) ||
          !trialDivisionChainTest(_primeSource, chainOriginExtra, false, chainLength/2, depth))
        return false;
    }
  }
  
  return true;
}
