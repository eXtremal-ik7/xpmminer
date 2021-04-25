/*
 * xpmclient.h
 *
 *  Created on: 01.05.2014
 *      Author: mad
 */

#ifndef XPMCLIENT_H_
#define XPMCLIENT_H_


#include <gmp.h>
#include <gmpxx.h>

#include "cudautil.h"
#include "uint256.h"
#include "sha256.h"

#define FERMAT_PIPELINES 2

#define PW 512        // Pipeline width (number of hashes to store)
#define SW 16         // maximum number of sieves in one iteration
#define MSO 128*1024  // max sieve output
#define MFS 2*SW*MSO  // max fermat size

const unsigned maxHashPrimorial = 16;

extern unsigned gPrimes[96*1024];
extern std::vector<unsigned> gPrimes2;

struct stats_t {
	
	unsigned id;
	unsigned errors;
	unsigned fps;
	double primeprob;
	double cpd;
	
	stats_t(){
		id = 0;
		errors = 0;
		fps = 0;
		primeprob = 0;
		cpd = 0;
	}
	
};


struct config_t {
	
	uint32_t N;
	uint32_t SIZE;
	uint32_t STRIPES;
	uint32_t WIDTH;
	uint32_t PCOUNT;
	uint32_t TARGET;
	uint32_t LIMIT13;
	uint32_t LIMIT14;
	uint32_t LIMIT15;
};


struct CUDADeviceInfo {
  int index;
  CUdevice device;
  CUcontext context;
  int majorComputeCapability;
  int minorComputeCapability;
};


template<typename T> class lifoBuffer {
private:
  T *_data;
  size_t _size;
  size_t _readPos;
  size_t _writePos;
  
  size_t nextPos(size_t pos) { return (pos+1) % _size; }
  
public:
  lifoBuffer(size_t size) : _size(size), _readPos(0), _writePos(0) {
    _data = new T[size];
  }
  
  size_t readPos() const { return _readPos; }
  size_t writePos() const { return _writePos; }
  T *data() const { return _data; }
  T& get(size_t index) const { return _data[index]; }
  
  bool empty() const {
    return _readPos == _writePos;
  }  
  
  size_t remaining() const {
    return _writePos >= _readPos ?
    _writePos - _readPos :
    _size - (_readPos - _writePos);
  }
  
  void clear() {
    _readPos = _writePos;
  }
  
  size_t push(const T& element) {
    size_t oldWritePos = _writePos;
    size_t nextWritePos = nextPos(_writePos);
    if (nextWritePos != _readPos) {
      _data[_writePos] = element;
      _writePos = nextWritePos;
    }
    
    return oldWritePos;
  }
  
  size_t pop() {
    size_t oldReadPos = _readPos;
    if (!empty())
      _readPos = nextPos(_readPos);
    return oldReadPos;
  }
};

class PrimeMiner {
public:
	
	struct block_t {
		
		static const int CURRENT_VERSION = 2;
		
		int version;
		uint256 hashPrevBlock;
		uint256 hashMerkleRoot;
		unsigned int time;
		unsigned int bits;
		unsigned int nonce;
		
	};
  
	struct search_t {
		
		cudaBuffer<uint32_t> midstate;
		cudaBuffer<uint32_t> found;
    cudaBuffer<uint32_t> primorialBitField;
		cudaBuffer<uint32_t> count;   
    
		
	};  
	
	struct hash_t {
		
		unsigned iter;
		unsigned nonce;
		unsigned time;
		uint256 hash;
		mpz_class shash;
    mpz_class primorial;
    unsigned primorialIdx;
	};
	
	struct fermat_t {
		uint32_t index;
    uint32_t hashid;
		uint8_t origin;
		uint8_t chainpos;
		uint8_t type;
		uint8_t reserved;
	};

  struct info_t {
    cudaBuffer<fermat_t> info;
    cudaBuffer<uint32_t> count;
  };
  
  struct pipeline_t {
    unsigned current;
    unsigned bsize;
    cudaBuffer<uint32_t> input;
    cudaBuffer<uint8_t> output;
    info_t buffer[2];
  };
	
  PrimeMiner(unsigned id, unsigned threads, unsigned sievePerRound, unsigned depth, unsigned LSize);
	~PrimeMiner();
	
  bool Initialize(CUcontext context, CUdevice device, CUmodule module);
	
  config_t getConfig() { return mConfig; }
	
	bool MakeExit;
	void Mining();
	
private:
  void FermatInit(pipeline_t &fermat, unsigned mfs);  
  
  void FermatDispatch(pipeline_t &fermat,
                      cudaBuffer<fermat_t>  sieveBuffers[SW][FERMAT_PIPELINES][2],
                      cudaBuffer<uint32_t> candidatesCountBuffers[SW][2],
                      unsigned pipelineIdx,
                      int ridx,
                      int widx,
                      uint64_t &testCount,
                      uint64_t &fermatCount,
                      CUfunction fermatKernel,
                      unsigned sievePerRound);  
  
	
	unsigned mID;
	unsigned mThreads;
	
	config_t mConfig;
  unsigned mSievePerRound;
	unsigned mBlockSize;
	uint32_t mDepth;
  unsigned mLSize;  

  CUcontext _context;
  CUstream mSieveStream;
	CUstream mHMFermatStream;

	CUfunction mHashMod;
	CUfunction mSieveSetup;
	CUfunction mSieve;
	CUfunction mSieveSearch;
	CUfunction mFermatSetup;
	CUfunction mFermatKernel352;
  CUfunction mFermatKernel320;  
	CUfunction mFermatCheck;
  info_t final;
  cudaBuffer<uint32_t> hashBuf;
};

#endif /* XPMCLIENT_H_ */
