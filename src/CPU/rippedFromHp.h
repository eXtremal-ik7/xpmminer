#include <stdlib.h>
#include <stdint.h>
#include <gmpxx.h>
#include <vector>
#include "primecoin.h"
#include "ripped.h"

class CSieveOfEratosthenesHp
{
    unsigned int nSieveSize; // size of the sieve
    unsigned int nSieveFilterPrimes; // filter a certain number of primes
    unsigned int nSieveExtensions; // extend the sieve a given number of times
    unsigned int nBits; // target of the prime chain to search for
    mpz_class mpzHash; // hash of the block header
    mpz_class mpzFixedMultiplier; // fixed round multiplier

    // final set of candidates for probable primality checking
    sieve_word_t *vfCandidates;
    sieve_word_t *vfCompositeBiTwin;
    sieve_word_t *vfCompositeCunningham1;
    sieve_word_t *vfCompositeCunningham2;
    sieve_word_t *vfCompositeLayerCC1;
    sieve_word_t *vfCompositeLayerCC2;

    // extended sets
    sieve_word_t *vfExtendedCandidates;
    sieve_word_t *vfExtendedCompositeBiTwin;
    sieve_word_t *vfExtendedCompositeCunningham1;
    sieve_word_t *vfExtendedCompositeCunningham2;

    // divisible multipliers
    unsigned int *vCunningham1Multipliers;
    unsigned int *vCunningham2Multipliers;

    static const unsigned int nWordBits = 8 * sizeof(sieve_word_t);
    unsigned int nCandidatesWords;
    unsigned int nCandidatesBytes;

    unsigned int nPrimeSeq; // prime sequence number currently being processed
    unsigned int nCandidateCount; // cached total count of candidates
    unsigned int nCandidateMultiplier; // current candidate for power test
    unsigned int nCandidateIndex; // internal candidate index
    bool fCandidateIsExtended; // is the current candidate in the extended part
    unsigned int nCandidateActiveExtension; // which extension is active

    unsigned int nChainLength; // target chain length
    unsigned int nSieveLayers; // sieve layers
    unsigned int nPrimes; // number of times to weave the sieve
    unsigned int nL1CacheElements; // number of bits that can be stored in L1 cache

    CBlockIndex* pindexPrev;

    // previous parameters
    unsigned int nCandidatesBytesPrev;
    unsigned int nSieveExtensionsPrev;
    unsigned int nMultiplierBytesPrev;

    bool fIsReady;
    bool fIsDepleted;
    
    const PrimeSource &vPrimes;

    unsigned int GetWordNum(unsigned int nBitNum) {
        return nBitNum / nWordBits;
    }

    sieve_word_t GetBitMask(unsigned int nBitNum) {
        return (sieve_word_t)1 << (nBitNum % nWordBits);
    }

    void ProcessMultiplier(sieve_word_t *vfComposites, const unsigned int nMinMultiplier, const unsigned int nMaxMultiplier, const PrimeSource& vPrimes, unsigned int *vMultipliers, unsigned int nLayerSeq);

    void freeArrays()
    {
        if (vfCandidates)
            free(vfCandidates);
        if (vfCompositeBiTwin)
            free(vfCompositeBiTwin);
        if (vfCompositeCunningham1)
            free(vfCompositeCunningham1);
        if (vfCompositeCunningham2)
            free(vfCompositeCunningham2);
        if (vfCompositeLayerCC1)
            free(vfCompositeLayerCC1);
        if (vfCompositeLayerCC2)
            free(vfCompositeLayerCC2);
        if (vfExtendedCandidates)
            free(vfExtendedCandidates);
        if (vfExtendedCompositeBiTwin)
            free(vfExtendedCompositeBiTwin);
        if (vfExtendedCompositeCunningham1)
            free(vfExtendedCompositeCunningham1);
        if (vfExtendedCompositeCunningham2)
            free(vfExtendedCompositeCunningham2);
        if (vCunningham1Multipliers)
            free(vCunningham1Multipliers);
        if (vCunningham2Multipliers)
            free(vCunningham2Multipliers);
        vfCandidates = NULL;
        vfCompositeBiTwin = NULL;
        vfCompositeCunningham1 = NULL;
        vfCompositeCunningham2 = NULL;
        vfCompositeLayerCC1 = NULL;
        vfCompositeLayerCC2 = NULL;
        vfExtendedCandidates = NULL;
        vfExtendedCompositeBiTwin = NULL;
        vfExtendedCompositeCunningham1 = NULL;
        vfExtendedCompositeCunningham2 = NULL;
    }

public:
    CSieveOfEratosthenesHp(const PrimeSource &primeSource) : vPrimes(primeSource)
    {
        nSieveSize = 0;
        nSieveFilterPrimes = 0;
        nSieveExtensions = 0;
        nBits = 0;
        mpzHash = 0;
        mpzFixedMultiplier = 0;
        vfCandidates = NULL;
        vfCompositeBiTwin = NULL;
        vfCompositeCunningham1 = NULL;
        vfCompositeCunningham2 = NULL;
        vfCompositeLayerCC1 = NULL;
        vfCompositeLayerCC2 = NULL;
        vfExtendedCandidates = NULL;
        vfExtendedCompositeBiTwin = NULL;
        vfExtendedCompositeCunningham1 = NULL;
        vfExtendedCompositeCunningham2 = NULL;
        vCunningham1Multipliers = NULL;
        vCunningham2Multipliers = NULL;
        nCandidatesWords = 0;
        nCandidatesBytes = 0;
        nCandidatesBytesPrev = 0;
        nSieveExtensionsPrev = 0;
        nMultiplierBytesPrev = 0;
        nPrimeSeq = 0;
        nCandidateCount = 0;
        nCandidateMultiplier = 0;
        nCandidateIndex = 0;
        fCandidateIsExtended = false;
        nCandidateActiveExtension = 0;
        nChainLength = 0;
        nSieveLayers = 0;
        nPrimes = 0;
        nL1CacheElements = 0;
        pindexPrev = NULL;
        fIsReady = false;
        fIsDepleted = true;
    }

    ~CSieveOfEratosthenesHp()
    {
        freeArrays();
    }

    void Reset(unsigned int nSieveSize,
               unsigned int nSieveFilterPrimes,
               unsigned int nSieveExtensions,
               unsigned int nL1CacheSize,
               unsigned int nBits,
               mpz_class& mpzHash,
               mpz_class& mpzFixedMultiplier,
               CBlockIndex* pindexPrev)
    {
        this->nSieveSize = nSieveSize;
        this->nSieveFilterPrimes = nSieveFilterPrimes;
        this->nSieveExtensions = nSieveExtensions;
        nL1CacheElements = nL1CacheSize * 8;
        this->nBits = nBits;
        this->mpzHash = mpzHash;
        this->mpzFixedMultiplier = mpzFixedMultiplier;
        this->pindexPrev = pindexPrev;
        nPrimeSeq = 0;
        nCandidateCount = 0;
        nCandidateMultiplier = 0;
        nCandidateIndex = 0;
        fCandidateIsExtended = false;
        nCandidateActiveExtension = 0;
        nCandidatesWords = (nSieveSize + nWordBits - 1) / nWordBits;
        nCandidatesBytes = nCandidatesWords * sizeof(sieve_word_t);
        nChainLength = TargetGetLength(nBits);

        // Override target length if requested
        if (nSieveTargetLength > 0)
            nChainLength = nSieveTargetLength;
        nSieveLayers = nChainLength + nSieveExtensions;

        // Filter only a certain number of prime factors
        // Most composites are still found
        nPrimes = nSieveFilterPrimes;
        const unsigned int nMultiplierBytes = nPrimes * nSieveLayers * sizeof(unsigned int);

        // Allocate arrays if parameters have changed
        if (nCandidatesBytes != nCandidatesBytesPrev || nSieveExtensions != nSieveExtensionsPrev || nMultiplierBytes != nMultiplierBytesPrev)
        {
            nCandidatesBytesPrev = nCandidatesBytes;
            nSieveExtensionsPrev = nSieveExtensions;
            nMultiplierBytesPrev = nMultiplierBytes;
            freeArrays();
            vfCandidates = (sieve_word_t *)malloc(nCandidatesBytes);
            vfCompositeBiTwin = (sieve_word_t *)malloc(nCandidatesBytes);
            vfCompositeCunningham1 = (sieve_word_t *)malloc(nCandidatesBytes);
            vfCompositeCunningham2 = (sieve_word_t *)malloc(nCandidatesBytes);
            vfCompositeLayerCC1 = (sieve_word_t *)malloc(nCandidatesBytes);
            vfCompositeLayerCC2 = (sieve_word_t *)malloc(nCandidatesBytes);
            vfExtendedCandidates = (sieve_word_t *)malloc(nSieveExtensions * nCandidatesBytes);
            vfExtendedCompositeBiTwin = (sieve_word_t *)malloc(nSieveExtensions * nCandidatesBytes);
            vfExtendedCompositeCunningham1 = (sieve_word_t *)malloc(nSieveExtensions * nCandidatesBytes);
            vfExtendedCompositeCunningham2 = (sieve_word_t *)malloc(nSieveExtensions * nCandidatesBytes);
            vCunningham1Multipliers = (unsigned int *)malloc(nMultiplierBytes);
            vCunningham2Multipliers = (unsigned int *)malloc(nMultiplierBytes);
        }

        // Initialize arrays
        memset(vfCandidates, 0, nCandidatesBytes);
        memset(vfCompositeBiTwin, 0, nCandidatesBytes);
        memset(vfCompositeCunningham1, 0, nCandidatesBytes);
        memset(vfCompositeCunningham2, 0, nCandidatesBytes);
        memset(vfCompositeLayerCC1, 0, nCandidatesBytes);
        memset(vfCompositeLayerCC2, 0, nCandidatesBytes);
        memset(vfExtendedCandidates, 0, nSieveExtensions * nCandidatesBytes);
        memset(vfExtendedCompositeBiTwin, 0, nSieveExtensions * nCandidatesBytes);
        memset(vfExtendedCompositeCunningham1, 0, nSieveExtensions * nCandidatesBytes);
        memset(vfExtendedCompositeCunningham2, 0, nSieveExtensions * nCandidatesBytes);
        memset(vCunningham1Multipliers, 0xFF, nMultiplierBytes);
        memset(vCunningham2Multipliers, 0xFF, nMultiplierBytes);

        fIsReady = true;
        fIsDepleted = false;
    }

    // Get total number of candidates for power test
    unsigned int GetCandidateCount()
    {
        if (nCandidateCount)
            return nCandidateCount;

        unsigned int nCandidates = 0;
#ifdef __GNUC__
        for (unsigned int i = 0; i < nCandidatesWords; i++)
            nCandidates += __builtin_popcountl(vfCandidates[i]);
        for (unsigned int j = 0; j < nSieveExtensions; j++)
            for (unsigned int i = nCandidatesWords / 2; i < nCandidatesWords; i++)
                nCandidates += __builtin_popcountl(vfExtendedCandidates[j * nCandidatesWords + i]);
#else
        for (unsigned int i = 0; i < nCandidatesWords; i++)
        {
            sieve_word_t lBits = vfCandidates[i];
            for (unsigned int j = 0; j < nWordBits; j++)
            {
                nCandidates += (lBits & 1);
                lBits >>= 1;
            }
        }
        for (unsigned int j = 0; j < nSieveExtensions; j++)
        {
            for (unsigned int i = nCandidatesWords / 2; i < nCandidatesWords; i++)
            {
                sieve_word_t lBits = vfExtendedCandidates[j * nCandidatesWords + i];
                for (unsigned int j = 0; j < nWordBits; j++)
                {
                    nCandidates += (lBits & 1);
                    lBits >>= 1;
                }
            }
        }
#endif
        nCandidateCount = nCandidates;
        return nCandidates;
    }

    // Scan for the next candidate multiplier (variable part)
    // Return values:
    //   True - found next candidate; nVariableMultiplier has the candidate
    //   False - scan complete, no more candidate and reset scan
    bool GetNextCandidateMultiplier(unsigned int& nVariableMultiplier, unsigned int& nCandidateType)
    {
        sieve_word_t *vfActiveCandidates;
        sieve_word_t *vfActiveCompositeTWN;
        sieve_word_t *vfActiveCompositeCC1;

        if (fCandidateIsExtended)
        {
            vfActiveCandidates = vfExtendedCandidates + nCandidateActiveExtension * nCandidatesWords;
            vfActiveCompositeTWN = vfExtendedCompositeBiTwin + nCandidateActiveExtension * nCandidatesWords;
            vfActiveCompositeCC1 = vfExtendedCompositeCunningham1 + nCandidateActiveExtension * nCandidatesWords;
        }
        else
        {
            vfActiveCandidates = vfCandidates;
            vfActiveCompositeTWN = vfCompositeBiTwin;
            vfActiveCompositeCC1 = vfCompositeCunningham1;
        }

        // Acquire the current word from the bitmap
        sieve_word_t lBits = vfActiveCandidates[GetWordNum(nCandidateIndex)];

        while (true)
        {
            nCandidateIndex++;
            if (nCandidateIndex >= nSieveSize)
            {
                // Check if extensions are available
                if (!fCandidateIsExtended && nSieveExtensions > 0)
                {
                    fCandidateIsExtended = true;
                    nCandidateActiveExtension = 0;
                    nCandidateIndex = nSieveSize / 2;
                }
                else if (fCandidateIsExtended && nCandidateActiveExtension + 1 < nSieveExtensions)
                {
                    nCandidateActiveExtension++;
                    nCandidateIndex = nSieveSize / 2;
                }
                else
                {
                    // Out of candidates
                    fCandidateIsExtended = false;
                    nCandidateActiveExtension = 0;
                    nCandidateIndex = 0;
                    nCandidateMultiplier = 0;
                    fIsDepleted = true;
                    return false;
                }

                // Fix the pointers
                if (fCandidateIsExtended)
                {
                    vfActiveCandidates = vfExtendedCandidates + nCandidateActiveExtension * nCandidatesWords;
                    vfActiveCompositeTWN = vfExtendedCompositeBiTwin + nCandidateActiveExtension * nCandidatesWords;
                    vfActiveCompositeCC1 = vfExtendedCompositeCunningham1 + nCandidateActiveExtension * nCandidatesWords;
                }
                else
                {
                    vfActiveCandidates = vfCandidates;
                    vfActiveCompositeTWN = vfCompositeBiTwin;
                    vfActiveCompositeCC1 = vfCompositeCunningham1;
                }
            }

            if (nCandidateIndex % nWordBits == 0)
            {
                // Update the current word
                lBits = vfActiveCandidates[GetWordNum(nCandidateIndex)];

                // Check if any bits are set
                if (lBits == 0)
                {
                    // Skip an entire word
                    nCandidateIndex += nWordBits - 1;
                    continue;
                }
            }

            if (lBits & GetBitMask(nCandidateIndex))
            {
                if (fCandidateIsExtended)
                    nCandidateMultiplier = nCandidateIndex * (2 << nCandidateActiveExtension);
                else
                    nCandidateMultiplier = nCandidateIndex;
                nVariableMultiplier = nCandidateMultiplier;
                if (~vfActiveCompositeTWN[GetWordNum(nCandidateIndex)] & GetBitMask(nCandidateIndex))
                    nCandidateType = PRIME_CHAIN_BI_TWIN;
                else if (~vfActiveCompositeCC1[GetWordNum(nCandidateIndex)] & GetBitMask(nCandidateIndex))
                    nCandidateType = PRIME_CHAIN_CUNNINGHAM1;
                else
                    nCandidateType = PRIME_CHAIN_CUNNINGHAM2;
                return true;
            }
        }
    }

    // Get progress percentage of the sieve
    unsigned int GetProgressPercentage();

    // Weave the sieve for the next prime in table
    // Return values:
    //   True  - weaved another prime; nComposite - number of composites removed
    //   False - sieve already completed
    bool Weave();

    bool IsReady() { return fIsReady; }
    bool IsDepleted() { return fIsDepleted; }
    void Deplete() { fIsDepleted = true; }
};
