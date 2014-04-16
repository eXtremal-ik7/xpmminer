#ifndef __RIPPED_H_
#define __RIPPED_H_

typedef unsigned long sieve_word_t;
typedef uint64_t uint64;
#define ULONG_MAX_VALUE std::numeric_limits<unsigned long>::max()
#define UINT_MAX_VALUE std::numeric_limits<unsigned int>::max()

class CBlockIndex;

const unsigned int nFractionalBits = 24;
const unsigned int TARGET_FRACTIONAL_MASK = (1u<<nFractionalBits) - 1;
const unsigned int TARGET_LENGTH_MASK = ~TARGET_FRACTIONAL_MASK;
const uint64_t nFractionalDifficultyMax = (1llu << (nFractionalBits + 32));
const uint64_t nFractionalDifficultyMin = (1llu << 32);
const uint64_t nFractionalDifficultyThreshold = (1llu << (8 + 32));

static int nSieveTargetLength = -1;


static unsigned int TargetGetLength(unsigned int nBits)
{
  return ((nBits & TARGET_LENGTH_MASK) >> nFractionalBits);
}

#endif //__RIPPED_H_
