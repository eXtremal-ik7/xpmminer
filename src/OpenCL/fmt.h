#ifdef __HOST__
#include <stdint.h>
struct uint4 {
  uint32_t x;
  uint32_t y;
  uint32_t z;
  uint32_t w;
};
#endif

#define FermatQueueChunks 16
#define FermatQueueBufferSize (FermatQueueChunks*256)

#pragma pack(push, 1)
struct GPUNonceAndHash {
  uint4 hash[2*256];
  uint32_t nonce[256];
  uint32_t currentNonce;
  uint32_t totalNonces;
  uint32_t align[2];
};

struct FermatQueue {
  uint32_t position;
  uint32_t size;
  uint32_t _align1[2];
  
  uint4 chainOrigins[3*FermatQueueBufferSize];
  uint32_t multipliers[FermatQueueBufferSize];
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
