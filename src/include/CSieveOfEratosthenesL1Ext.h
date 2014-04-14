#include "primecoin.h"
#include <gmpxx.h>

class CSieveOfEratosthenesL1Ext { 
public:
  enum {
    L1CacheSize = 131072,
    ExtensionsNum = 9,
    MaxSieveSize = L1CacheSize*32,
    MaxDepth = 32768,
    MaxChainLength = 20
  };
  
private:
  const PrimeSource &_primeSource;
  unsigned nSieveSize; // size of the sieve
  unsigned sieveWords;
  unsigned chainLength; // target of the prime chain to search for
  unsigned depth;
  unsigned layersNum;
  uint32_t fixedFactor[256];
  
  // bit maps of the sieve, index represents the variable part of multiplier
  uint64_t cunningham1Bitfield[(ExtensionsNum+1)*(MaxSieveSize/64)];
  uint64_t cunningham2Bitfield[(ExtensionsNum+1)*(MaxSieveSize/64)];
  uint64_t bitwinBitfield[(ExtensionsNum+1)*(MaxSieveSize/64)];
  uint32_t cunningham1Multipliers[MaxDepth*(MaxChainLength+ExtensionsNum)];
  uint32_t cunningham2Multipliers[MaxDepth*(MaxChainLength+ExtensionsNum)];
  
  struct {
    uint64_t cunningham1[L1CacheSize/64];
    uint64_t cunningham2[L1CacheSize/64];
  } local;
  
  // iteration context
  uint64_t word;
  uint64_t *ptr;
  uint64_t *end;
  size_t wordIdx;
  unsigned bitOffset;
  unsigned extensionNum;
  unsigned candidateType;
  
  inline void setValue(uint64_t *bitfield, unsigned offset) {
    *(bitfield + (offset >> 6)) |= (1ULL << (offset&0x3F));
  }
  
  void setContext(uint64_t *bitfield) {
    ptr = bitfield;
    end = ptr + sieveWords;
    word = (~*ptr) >> 1;
    bitOffset = 1;
    wordIdx = 0;
    extensionNum = 0;
  }
  
  bool increment() {
    if (++ptr == end) {
      if (++extensionNum > ExtensionsNum) {
        if (++candidateType > PRIME_CHAIN_BI_TWIN) {
          return false; 
        } else {
          switch (candidateType) {
            case PRIME_CHAIN_CUNNINGHAM2 :
              setContext(cunningham2Bitfield);
              break;
            case PRIME_CHAIN_BI_TWIN :
              setContext(bitwinBitfield);
              break;
          }
        }
      } else {
        ptr += sieveWords/2;
        end += sieveWords;
        word = ~(*ptr);
        wordIdx = sieveWords/2;
        bitOffset = 0;    
      }
    } else {
      word = ~(*ptr);
      wordIdx++;
      bitOffset = 0;      
    }
    
    return true;
  }
  
public:
  CSieveOfEratosthenesL1Ext(const PrimeSource &primeSource) :
    _primeSource(primeSource) {}
  
  void reset(unsigned sieveSize,
             unsigned chainLength,
             unsigned depth,
             const uint32_t *fixedMultiplier);
  
  void reset(unsigned sieveSize,
             unsigned chainLength,
             unsigned depth,
             const mpz_class &fixedMultiplier);
  
  unsigned int GetCandidateCount();
  bool GetNextCandidateMultiplier(unsigned &nVariableMultiplier, unsigned &nCandidateType);
  
  void Weave();
 
  void resetCandidateIterator() {
    candidateType = PRIME_CHAIN_CUNNINGHAM1;
    setContext(cunningham1Bitfield);
  }
  
  bool fastSelfTest();
};

