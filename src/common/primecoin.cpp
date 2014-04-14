#include "CSieveOfEratosthenesL1Ext.h"
#include "primecoin.h"
#include "system.h"

#include <stdlib.h>
#include <time.h>
#include <limits>

#include <openssl/bn.h>
#include <openssl/sha.h>

#include <algorithm>

extern unsigned gDebug;
extern unsigned gSieveSize;
extern unsigned gWeaveDepth;

uint32_t bitsFromDifficulty(double difficulty)
{
  uint32_t chainLength = (uint32_t)difficulty;
  uint32_t fractionalPart = (uint32_t)
    ((difficulty - (double)chainLength) * (1U << DifficultyFractionalBits));
    
  return (chainLength << DifficultyFractionalBits) | fractionalPart;
}

uint32_t bitsFromChainLengthWithoutFractional(uint32_t chainLength)
{
  return chainLength << DifficultyFractionalBits;
}

uint32_t chainLengthFromBits(uint32_t bits)
{
  return bits >> DifficultyFractionalBits;
}

double difficultyFromBits(uint32_t bits)
{
  return (double)(bits >> DifficultyFractionalBits) +
         (bits & DifficultyFractionalMask) * (1.0 / 16777216);
}

void incrementChainLengthInBits(uint32_t *bits)
{
  *bits += (1 << DifficultyFractionalBits);
}

void generateRandomHeader(PrimecoinBlockHeader *header, double difficulty)
{
  header->version = 2;
  for (unsigned i = 0; i < sizeof(header->hashPrevBlock); i++) {
    header->hashPrevBlock[i] = rand() % 0xFF;
    header->hashMerkleRoot[i] = rand() % 0xFF;
  }
  
  header->time = time(0);
  header->bits = bitsFromDifficulty(difficulty);
  header->nonce = 0;
}

// "Magic number" search function (for replace division to multiplication & shift)
static bool findMultiplierForConstantDivision(unsigned maxDividendBits,
                                              const mpz_class &divisor,
                                              mpz_class *multiplier,
                                              unsigned *offset)
{
  mpz_class two = 2;
  
  mpz_class maxDividend = 1;
  maxDividend <<= maxDividendBits;
  maxDividend--;
  
  mpz_class nc = (maxDividend / divisor) * divisor - 1;
  for (unsigned i = 0; i < 2*maxDividendBits + 1; i++) {
    mpz_class powOf2 = 1;
    powOf2 <<= i;
    if (powOf2 > nc*(divisor - 1 - ((powOf2 - 1) % divisor))) {
      *multiplier = (powOf2 + divisor - 1 - ((powOf2 - 1) % divisor)) / divisor;
      *offset = i;
      return true;
    }
  }
  
  return false;
}

void generatePrimes(uint32_t *out, unsigned primesNum)
{
  size_t index = 0;
  
  std::vector<bool> sieve(primesNum, false);
  for (uint32_t F = 2; F*F < primesNum; F++) {
    if (sieve[F])
      continue;
    for (uint32_t C = F*F; C < primesNum; C += F)
      sieve[C] = true;
  }
  
  for (unsigned int n = 2; n < primesNum; n++) {
    if (!sieve[n])
      out[index++] = n;
  }
}

void mpzExport(const mpz_class &N, uint32_t *limbs)
{
  size_t limbsNumber;
  mpz_export(limbs+1, &limbsNumber, -1, 4, 0, 0, N.get_mpz_t());
  limbs[0] = limbsNumber;
}

void mpzImport(const uint32_t *limbs, mpz_class &N)
{
  mpz_import(N.get_mpz_t(), limbs[0], -1, 4, 0, 0, limbs+1);
}

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
                         unsigned multipliersNum)
{
  unsigned primeCombinedIdx = 0;
  mpz_class divisor;
  mpz_class inversedMultiplier;
  mpz_class inversedMultiplier64;
  mpz_class inversedCombinedMultiplier;
  unsigned offset;
  unsigned offset64;
  unsigned combinedOffset;
  unsigned primeCombined;
  for (unsigned i = 0; i < multipliersNum; i++) {
    findMultiplierForConstantDivision(31, primes[i], &inversedMultiplier, &offset);
    if (offset < 32)
      findMultiplierForConstantDivision(32, primes[i], &inversedMultiplier, &offset);
    
    findMultiplierForConstantDivision(63, primes[i], &inversedMultiplier64, &offset64);    
    if (offset64 < 64)
      findMultiplierForConstantDivision(64, primes[i], &inversedMultiplier64, &offset64);          
    if (primeCombinedIdx <= i) {
      primeCombined = 1;
      while (primeCombined < (1U << 31) / primes[primeCombinedIdx])
        primeCombined *= primes[primeCombinedIdx++];
      
      divisor = primeCombined;
      findMultiplierForConstantDivision(63, divisor, &inversedCombinedMultiplier, &combinedOffset);
      isNewMultiplier[i] = 1;
    } else {
      isNewMultiplier[i] = 0;
    }
    
    size_t size;
    primesCombined[i] = primeCombined;    
    mpz_export(&multipliers[i], &size, -1, 4, 0, 0, inversedMultiplier.get_mpz_t());
    mpz_export(&multipliers64[i], &size, -1, 4, 0, 0, inversedMultiplier64.get_mpz_t());    
    mpz_export(&combinedMultipliers[i], &size, -1, 4, 0, 0, inversedCombinedMultiplier.get_mpz_t());
    offsets[i] = offset;
    offsets64[i] = offset64;
    combinedOffsets[i] = combinedOffset;
  }
}

bool trialDivisionChainTest(const PrimeSource &primeSource,
                            mpz_class &N,
                            bool fSophieGermain,
                            unsigned chainLength,
                            unsigned depth)
{
  N += (fSophieGermain ? -1 : 1);
  for (unsigned i = 0; i < chainLength; i++) {
    for (unsigned divIdx = 0; divIdx < depth; divIdx++) { 
      if (mpz_tdiv_ui(N.get_mpz_t(), primeSource[divIdx]) == 0) {
        fprintf(stderr, " * divisor: [%u]%u\n", divIdx, primeSource[divIdx]);
        return false;
      }
    }
    
    N <<= 1;
    N += (fSophieGermain ? 1 : -1);
  }
  
  return true;
}


PrimeSource::PrimeSource(uint32_t primesNum, unsigned inversedMultipliersNum) :
  _primesNum(primesNum)
{
  // Use Eratosthenes sieve for search first N primes
  _primes = new uint32_t[primesNum];

  _primesCombined = new uint32_t[inversedMultipliersNum];
  _isNewMultiplier = new uint32_t[inversedMultipliersNum];
  _multipliers = new uint32_t[inversedMultipliersNum];
  _offsets = new uint32_t[inversedMultipliersNum];
  _multipliers64 = new uint64_t[inversedMultipliersNum];
  _offsets64 = new uint32_t[inversedMultipliersNum];
  _combinedMultipliers = new uint64_t[inversedMultipliersNum];
  _combinedOffsets = new uint32_t[inversedMultipliersNum];
   
  generatePrimes(_primes, primesNum);
  generateMultipliers(_primesCombined, _isNewMultiplier,
                      _multipliers, _offsets,
                      _multipliers64, _offsets64,
                      _combinedMultipliers, _combinedOffsets,
                      _primes, primesNum, inversedMultipliersNum);
}

bool sha256(void *out, const void *data, size_t size)
{
  SHA256_CTX ctx;
  SHA256_Init(&ctx);
  SHA256_Update(&ctx, data, size);
  SHA256_Final((unsigned char*)out, &ctx);
  return true;
}


bool updateBlock(PrimecoinBlockHeader *header,
                 mpz_class &blockHeaderHash,
                 const PrimeSource &primeSource,
                 CPrimalityTestParams &testParams,
                 unsigned nonceIncrement)
{
  uint8_t hash1[32];
  uint8_t hashData[32];
  while (header->nonce < 0xFFFF0000) {
    header->nonce += nonceIncrement;
    sha256(hash1, header, 80);
    sha256(hashData, hash1, 32);
    // sha256-hash must be greater than 2^255
    if (!(hashData[31] & 0x80))
      continue;    
    
    mpz_import(blockHeaderHash.get_mpz_t(),
               32 / sizeof(unsigned long),
               -1,
               sizeof(unsigned long),
               -1,
               0,
               hashData);
    
    if (ProbablePrimalityTestWithTrialDivisionFast(blockHeaderHash, 1000, primeSource, testParams))
      break;
  }
  
  return header->nonce < 0xFFFF0000;
}

void PrimorialFast(unsigned int length,
                   mpz_class& bnPrimorial,
                   const PrimeSource &primeSource)
{
  bnPrimorial = 1;
  for (size_t i = 0; i < length; i++)
    bnPrimorial *= primeSource.prime(i);
}

bool FermatProbablePrimalityTestFast(const mpz_class &n,
                                     unsigned int& nLength, 
                                     CPrimalityTestParams &testParams,
                                     bool fFastFail)
{
  static const mpz_class mpzTwo = 2;
  mpz_t& mpzE = testParams.mpzE;
  mpz_t& mpzR = testParams.mpzR;
  
  mpz_sub_ui(mpzE, n.get_mpz_t(), 1);
  mpz_powm(mpzR, mpzTwo.get_mpz_t(), mpzE, n.get_mpz_t());
  if (mpz_cmp_ui(mpzR, 1) == 0)
    return true;
  if (fFastFail)
    return false;
  
  // Fermat test failed, calculate chain length (integer and fractional part)
  mpz_sub(mpzE, n.get_mpz_t(), mpzR);
  mpz_mul_2exp(mpzR, mpzE, DifficultyFractionalBits);
  mpz_tdiv_q(mpzE, mpzR, n.get_mpz_t());
  unsigned int fractionalLength = mpz_get_ui(mpzE);
  
  nLength = (nLength & DifficultyChainLengthMask) | fractionalLength;
  return false;
}

static bool EulerLagrangeLifchitzPrimalityTestFast(const mpz_class& n,
                                                   bool fSophieGermain,
                                                   unsigned int& nLength,
                                                   CPrimalityTestParams& testParams)
{
  static const mpz_class mpzTwo = 2;  
  
  // Faster GMP version
  mpz_t& mpzE = testParams.mpzE;
  mpz_t& mpzR = testParams.mpzR;
  mpz_t& mpzRplusOne = testParams.mpzRplusOne;
  
  mpz_sub_ui(mpzE, n.get_mpz_t(), 1);
  mpz_tdiv_q_2exp(mpzE, mpzE, 1);
  mpz_powm(mpzR, mpzTwo.get_mpz_t(), mpzE, n.get_mpz_t());
  unsigned int nMod8 = mpz_get_ui(n.get_mpz_t()) % 8;
  bool fPassedTest = false;
  if (fSophieGermain && (nMod8 == 7)) // Euler & Lagrange
    fPassedTest = !mpz_cmp_ui(mpzR, 1);
  else if (fSophieGermain && (nMod8 == 3)) // Lifchitz
  {
    mpz_add_ui(mpzRplusOne, mpzR, 1);
    fPassedTest = !mpz_cmp(mpzRplusOne, n.get_mpz_t());
  }
  else if ((!fSophieGermain) && (nMod8 == 5)) // Lifchitz
  {
    mpz_add_ui(mpzRplusOne, mpzR, 1);
    fPassedTest = !mpz_cmp(mpzRplusOne, n.get_mpz_t());
  }
  else if ((!fSophieGermain) && (nMod8 == 1)) // LifChitz
    fPassedTest = !mpz_cmp_ui(mpzR, 1);
  
  if (fPassedTest)
  {
    return true;
  }
  
  // Failed test, calculate fractional length
  mpz_mul(mpzE, mpzR, mpzR);
  mpz_tdiv_r(mpzR, mpzE, n.get_mpz_t()); // derive Fermat test remainder
  
  mpz_sub(mpzE, n.get_mpz_t(), mpzR);
  mpz_mul_2exp(mpzR, mpzE, DifficultyFractionalBits);
  mpz_tdiv_q(mpzE, mpzR, n.get_mpz_t());
  unsigned int nFractionalLength = mpz_get_ui(mpzE);
  
  nLength = (nLength & DifficultyChainLengthMask) | nFractionalLength;
  return false;
}

bool ProbablePrimalityTestWithTrialDivisionFast(const mpz_class &candidate,
                                                unsigned trialDivisionLimit,
                                                const PrimeSource &primeSource,
                                                CPrimalityTestParams &testParams)
{
  for (unsigned i = 0; i < trialDivisionLimit; i++) {
    if (candidate % primeSource.prime(i) == 0)
      return false;
  }
  unsigned nLength = 0;
  return FermatProbablePrimalityTestFast(candidate, nLength, testParams, true);
}


// Test Probable Cunningham Chain for: n
// fSophieGermain:
//   true - Test for Cunningham Chain of first kind (n, 2n+1, 4n+3, ...)
//   false - Test for Cunningham Chain of second kind (n, 2n-1, 4n-3, ...)
// Return value:
//   true - Probable Cunningham Chain found (length at least 2)
//   false - Not Cunningham Chain
static bool ProbableCunninghamChainTestFast(const mpz_class& n,
                                            bool fSophieGermain,
                                            bool fFermatTest,
                                            unsigned int& nProbableChainLength,
                                            CPrimalityTestParams& testParams)
{
  nProbableChainLength = 0;
  
  // Fermat test for n first
  if (!FermatProbablePrimalityTestFast(n, nProbableChainLength, testParams, true))
    return false;
  
  // Euler-Lagrange-Lifchitz test for the following numbers in chain
  mpz_class &N = testParams.N;
  N = n;
  while (true) {
    incrementChainLengthInBits(&nProbableChainLength);
    N <<= 1;
    N += (fSophieGermain? 1 : (-1));
    if (fFermatTest) {
      if (!FermatProbablePrimalityTestFast(N, nProbableChainLength, testParams))
        break;
    } else {
      if (!EulerLagrangeLifchitzPrimalityTestFast(N, fSophieGermain, nProbableChainLength, testParams))
        break;
    }
  }
  
  return (chainLengthFromBits(nProbableChainLength) >= 2);
}

// Test probable prime chain for: nOrigin
// Return value:
//   true - Probable prime chain found (one of nChainLength meeting target)
//   false - prime chain too short (none of nChainLength meeting target)
bool ProbablePrimeChainTestFast(const mpz_class& mpzPrimeChainOrigin,
                                CPrimalityTestParams& testParams)
{
  const unsigned int nBits = testParams.bits;
  const unsigned int nCandidateType = testParams.candidateType;
  unsigned int& nChainLength = testParams.chainLength;
  mpz_class& mpzOriginMinusOne = testParams.mpzOriginMinusOne;
  mpz_class& mpzOriginPlusOne = testParams.mpzOriginPlusOne;
  nChainLength = 0;
  
  // Test for Cunningham Chain of first kind
  if (nCandidateType == PRIME_CHAIN_CUNNINGHAM1) {
    mpzOriginMinusOne = mpzPrimeChainOrigin - 1;
    ProbableCunninghamChainTestFast(mpzOriginMinusOne, true, false, nChainLength, testParams);
  } else if (nCandidateType == PRIME_CHAIN_CUNNINGHAM2) {
    // Test for Cunningham Chain of second kind
    mpzOriginPlusOne = mpzPrimeChainOrigin + 1;
    ProbableCunninghamChainTestFast(mpzOriginPlusOne, false, false, nChainLength, testParams);
  } else {
    unsigned int nChainLengthCunningham1 = 0;
    unsigned int nChainLengthCunningham2 = 0;
    mpzOriginMinusOne = mpzPrimeChainOrigin - 1;
    if (ProbableCunninghamChainTestFast(mpzOriginMinusOne, true, false, nChainLengthCunningham1, testParams)) {
      mpzOriginPlusOne = mpzPrimeChainOrigin + 1;
      ProbableCunninghamChainTestFast(mpzOriginPlusOne, false, false, nChainLengthCunningham2, testParams);
      // Figure out BiTwin Chain length
      // BiTwin Chain allows a single prime at the end for odd length chain
      nChainLength =
      (chainLengthFromBits(nChainLengthCunningham1) > chainLengthFromBits(nChainLengthCunningham2))?
      (nChainLengthCunningham2 + bitsFromChainLengthWithoutFractional(chainLengthFromBits(nChainLengthCunningham2)+1)) :
      (nChainLengthCunningham1 + bitsFromChainLengthWithoutFractional(chainLengthFromBits(nChainLengthCunningham1)));
    }
  }
  
  return (nChainLength >= nBits);
}


bool MineProbablePrimeChainFast(PrimecoinBlockHeader &header,
                                CSieveOfEratosthenesL1Ext *sieve,
                                mpz_class &blockHeaderHash,
                                mpz_class &primorial,
                                unsigned int& nProbableChainLength,
                                unsigned int& nTests,
                                unsigned int& nPrimesHit,
                                CPrimalityTestParams &testParams,
                                const PrimeSource &primeSource,
                                uint64_t *foundChains)
{
  timeMark sieveBegin = getTimeMark();
  mpz_class hashMultiplier = blockHeaderHash*primorial;
  sieve->reset(gSieveSize, chainLengthFromBits(header.bits), gWeaveDepth, hashMultiplier);
  sieve->Weave();
  timeMark sieveEnd = getTimeMark();  
  if (gDebug) {
    fprintf(stderr,
            " * sieve %.3lfmsec: %u@%u/%u ",
            usDiff(sieveBegin, sieveEnd) / 1000.0,
            sieve->GetCandidateCount(),
            gSieveSize,
            gWeaveDepth);
  }  
  
  nTests = 0;
  nPrimesHit = 0;
  unsigned nTriedMultiplier;
  mpz_class bnChainOrigin;
  
  unsigned int &nChainLength = testParams.chainLength;
  unsigned int &nCandidateType = testParams.candidateType;      
  sieve->resetCandidateIterator();
  while (true) {
    nTests++;
    if (!sieve->GetNextCandidateMultiplier(nTriedMultiplier, nCandidateType)) {
      timeMark primalityTestEnd = getTimeMark();
      if (gDebug) {
        fprintf(stderr,
                " primality Fermat test %.3lfmsec\n",
                usDiff(sieveEnd, primalityTestEnd) / 1000.0);
      }
      
      return false;
    }
    
    bnChainOrigin = hashMultiplier;
    bnChainOrigin *= nTriedMultiplier;
    nChainLength = 0;
    if (ProbablePrimeChainTestFast(bnChainOrigin, testParams)) {
      uint8_t buffer[256];
      BIGNUM *xxx = 0;
      mpz_class targetMultiplier = primorial*nTriedMultiplier;
      BN_dec2bn(&xxx, targetMultiplier.get_str().c_str());
      BN_bn2mpi(xxx, buffer);
      header.multiplier[0] = buffer[3];
      std::reverse_copy(buffer+4, buffer+4+buffer[3], header.multiplier+1);
      fprintf(stderr, "targetMultiplier=%s\n", targetMultiplier.get_str().c_str());
      return true;
    }
    
    nProbableChainLength = nChainLength;
    if (chainLengthFromBits(nProbableChainLength) >= 1) {
      foundChains[chainLengthFromBits(nProbableChainLength)]++;
      nPrimesHit++;
    }
  }
  
  return false;
}