#ifndef __PRIMECOIN_H_
#define __PRIMECOIN_H_

#include <stddef.h>
#include <stdint.h>
#include <vector>

#include "gmpxx.h"

static const uint32_t DifficultyFractionalBits = 24;
static const uint32_t DifficultyFractionalMask = (1U << DifficultyFractionalBits) - 1;
static const uint32_t DifficultyChainLengthMask = ~DifficultyFractionalMask;

const unsigned MaxChainLength = 20;

enum {
  PRIME_CHAIN_CUNNINGHAM1 = 1u,
  PRIME_CHAIN_CUNNINGHAM2 = 2u,
  PRIME_CHAIN_BI_TWIN     = 3u,
};

class PrimeSource {
private:
  uint32_t _primesNum;
  uint32_t *_primes;
  
  uint32_t *_isNewMultiplier;
  uint32_t *_primesCombined;  
  uint32_t *_multipliers;
  uint64_t *_multipliers64;
  uint64_t *_combinedMultipliers;
  uint32_t *_offsets;
  uint32_t *_offsets64;
  uint32_t *_combinedOffsets;

public:
  PrimeSource(uint32_t primesNum, unsigned inversedMultipliersNum);
  ~PrimeSource() { delete[] _primes; }
 
  uint32_t primesNum() const { return _primesNum; }
  uint32_t prime(uint32_t N) const { return _primes[N]; }
  uint32_t operator[](uint32_t N) const { return _primes[N]; }
  
  uint32_t *primesPtr() { return _primes; }
  uint32_t *isNewMultipliersPtr() { return _isNewMultiplier; }
  uint32_t *primesCombinedPtr() { return _primesCombined; }
  uint32_t *multipliersPtr() { return _multipliers; }
  uint64_t *multiplier64Ptr() { return _multipliers64; }
  uint64_t *combinedMultipliersPtr() { return _combinedMultipliers; }
  uint32_t *offsetsPtr() { return _offsets; }
  uint32_t *offsets64Ptr() { return _offsets64; }
  uint32_t *combinedOffsetsPtr() { return _combinedOffsets; }
  
  bool isNewMultiplier(unsigned N) const { return _isNewMultiplier[N]; }
  unsigned primesCombined(unsigned N) const { return _primesCombined[N]; }
  uint32_t multiplier(unsigned N) const { return _multipliers[N]; }
  uint64_t multiplier64(unsigned N) const { return _multipliers64[N]; }
  unsigned offset(unsigned N) const { return _offsets[N]; }
  unsigned offset64(unsigned N) const { return _offsets64[N]; }
  uint64_t combinedMultiplier(unsigned N) const { return _combinedMultipliers[N]; }
  unsigned combinedOffset(unsigned N) const { return _combinedOffsets[N]; }
};

#pragma pack(push, 1)
struct PrimecoinBlockHeader {
  int32_t version;
  uint8_t hashPrevBlock[32];
  uint8_t hashMerkleRoot[32];
  uint32_t time;
  uint32_t bits;
  uint32_t nonce;

  uint8_t multiplier[256];
};
#pragma pack(pop)


uint32_t bitsFromDifficulty(double difficulty);
uint32_t bitsFromChainLengthWithoutFractional(uint32_t chainLength);
uint32_t chainLengthFromBits(uint32_t bits);
double difficultyFromBits(uint32_t bits);

void incrementChainLengthInBits(uint32_t *bits);

void generateRandomHeader(PrimecoinBlockHeader *header, double difficulty);

void mpzExport(const mpz_class &N, uint32_t *limbs);
void mpzImport(const uint32_t *limbs, mpz_class &N);

void generatePrimes(uint32_t *out, unsigned primesNum);

void generateMultipliers(uint32_t *primesCombined,
                         uint32_t *isNewMultiplier,
                         uint32_t *multipliers,
                         uint32_t *offsets,
                         uint64_t *multipliers64,
                         uint32_t *offsets64,
                         uint64_t *combinedMultipliers,
                         uint32_t *combinedOffsets,
                         uint32_t *primes,
                         unsigned primesNum,
                         unsigned multipliersNum);

bool trialDivisionChainTest(const PrimeSource &primeSource,
                            mpz_class &N,
                            bool fSophieGermain,
                            unsigned chainLength,
                            unsigned depth);

bool sha256(void *out, const void *data, size_t size);

class CSieveOfEratosthenesL1Ext;

struct CPrimalityTestParams {
  mpz_t mpzE;
  mpz_t mpzR;
  mpz_t mpzRplusOne;
  
  mpz_class mpzOriginMinusOne;
  mpz_class mpzOriginPlusOne;
  mpz_class N;
  
  unsigned int bits;
  unsigned int candidateType;
  unsigned int chainLength;
  
  CPrimalityTestParams(unsigned int difficultyBits) {
    bits = difficultyBits;
    chainLength = 0;
    mpz_init(mpzE);
    mpz_init(mpzR);
    mpz_init(mpzRplusOne);
  }
  
  ~CPrimalityTestParams() {
    mpz_clear(mpzE);
    mpz_clear(mpzR);
    mpz_clear(mpzRplusOne);
  }
};

bool ProbablePrimeChainTestFast(const mpz_class& mpzPrimeChainOrigin,
                                CPrimalityTestParams& testParams);

void PrimorialFast(unsigned int length,
                   mpz_class& bnPrimorial,
                   const PrimeSource &primeSource);

bool FermatProbablePrimalityTestFast(const mpz_class &n,
                                     unsigned int& nLength, 
                                     CPrimalityTestParams &testParams,
                                     bool fFastFail = false);

bool ProbablePrimalityTestWithTrialDivisionFast(const mpz_class &candidate,
                                                unsigned trialDivisionLimit,
                                                const PrimeSource &primeSource,
                                                CPrimalityTestParams &testParams);

bool MineProbablePrimeChainFast(PrimecoinBlockHeader &header,
                                CSieveOfEratosthenesL1Ext *sieve,                                
                                mpz_class &blockHeaderHash,
                                mpz_class &primorial,
                                unsigned int& nProbableChainLength,
                                unsigned int& nTests,
                                unsigned int& nPrimesHit,
                                CPrimalityTestParams &testParams,
                                const PrimeSource &primeSource,
                                uint64_t *foundChains);

bool updateBlock(PrimecoinBlockHeader *header,
                 mpz_class &blockHeaderHash,
                 const PrimeSource &primeSource,
                 CPrimalityTestParams &testParams,
                 unsigned nonceIncrement = 1);

#endif //__PRIMECOIN_H_
