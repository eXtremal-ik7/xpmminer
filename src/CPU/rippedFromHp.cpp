#include "rippedFromHp.h"
#include <limits>
#include "system.h"

#define error printf

static void *pindexBest = 0;

static unsigned int int_invert(unsigned int a, unsigned int nPrime)
{
    // Extended Euclidean algorithm to calculate the inverse of a in finite field defined by nPrime
    int rem0 = nPrime, rem1 = a % nPrime, rem2;
    int aux0 = 0, aux1 = 1, aux2;
    int quotient, inverse;

    while (1)
    {
        if (rem1 <= 1)
        {
            inverse = aux1;
            break;
        }

        rem2 = rem0 % rem1;
        quotient = rem0 / rem1;
        aux2 = -quotient * aux1 + aux0;

        if (rem2 <= 1)
        {
            inverse = aux2;
            break;
        }

        rem0 = rem1 % rem2;
        quotient = rem1 / rem2;
        aux0 = -quotient * aux2 + aux1;

        if (rem0 <= 1)
        {
            inverse = aux0;
            break;
        }

        rem1 = rem2 % rem0;
        quotient = rem2 / rem0;
        aux1 = -quotient * aux0 + aux2;
    }

    return (inverse + nPrime) % nPrime;
}

void CSieveOfEratosthenesHp::ProcessMultiplier(sieve_word_t *vfComposites, const unsigned int nMinMultiplier, const unsigned int nMaxMultiplier, const PrimeSource& vPrimes, unsigned int *vMultipliers, unsigned int nLayerSeq)
{
    // Wipe the part of the array first
    if (nMinMultiplier < nMaxMultiplier)
        memset(vfComposites + GetWordNum(nMinMultiplier), 0, (nMaxMultiplier - nMinMultiplier + nWordBits - 1) / nWordBits * sizeof(sieve_word_t));

    for (unsigned int nPrimeSeq = 1; nPrimeSeq < nPrimes; nPrimeSeq++)
    {
        const unsigned int nPrime = vPrimes[nPrimeSeq];
        unsigned int nVariableMultiplier = vMultipliers[nPrimeSeq * nSieveLayers + nLayerSeq];
        if (nVariableMultiplier < nMinMultiplier)
            nVariableMultiplier += (nMinMultiplier - nVariableMultiplier + nPrime - 1) / nPrime * nPrime;
#ifdef USE_ROTATE
        const unsigned int nRotateBits = nPrime % nWordBits;
        sieve_word_t lBitMask = GetBitMask(nVariableMultiplier);
        for (; nVariableMultiplier < nMaxMultiplier; nVariableMultiplier += nPrime)
        {
            vfComposites[GetWordNum(nVariableMultiplier)] |= lBitMask;
            lBitMask = (lBitMask << nRotateBits) | (lBitMask >> (nWordBits - nRotateBits));
        }
        vMultipliers[nPrimeSeq * nSieveLayers + nLayerSeq] = nVariableMultiplier;
#else
        for (; nVariableMultiplier < nMaxMultiplier; nVariableMultiplier += nPrime)
        {
            vfComposites[GetWordNum(nVariableMultiplier)] |= GetBitMask(nVariableMultiplier);
        }
        vMultipliers[nPrimeSeq * nSieveLayers + nLayerSeq] = nVariableMultiplier;
#endif
    }
}

bool CSieveOfEratosthenesHp::Weave()
{
#ifdef DEBUG_MINING      
  timeMark beginPoint = getTimeMark();  
#endif
  
    // Check whether fixed multiplier fits in an unsigned long
    bool fUseLongForFixedMultiplier = mpzFixedMultiplier < ULONG_MAX_VALUE;
    unsigned long nFixedMultiplier;
    mpz_class mpzFixedFactor;
    if (fUseLongForFixedMultiplier)
        nFixedMultiplier = mpzFixedMultiplier.get_ui();
    else
        mpzFixedFactor = mpzHash * mpzFixedMultiplier;

    unsigned int nCombinedEndSeq = 1;
    unsigned int nFixedFactorCombinedMod = 0;

    for (unsigned int nPrimeSeqLocal = 1; nPrimeSeqLocal < nPrimes; nPrimeSeqLocal++)
    {
        if (pindexPrev != pindexBest)
            break;  // new block
        unsigned int nPrime = vPrimes[nPrimeSeqLocal];
        if (nPrimeSeqLocal >= nCombinedEndSeq)
        {
            // Combine multiple primes to produce a big divisor
            unsigned int nPrimeCombined = 1;
            while (nPrimeCombined < UINT_MAX_VALUE / vPrimes[nCombinedEndSeq])
            {
                nPrimeCombined *= vPrimes[nCombinedEndSeq];
                nCombinedEndSeq++;
            }

            if (fUseLongForFixedMultiplier)
            {
                nFixedFactorCombinedMod = mpz_tdiv_ui(mpzHash.get_mpz_t(), nPrimeCombined);
                nFixedFactorCombinedMod = (uint64)nFixedFactorCombinedMod * (nFixedMultiplier % nPrimeCombined) % nPrimeCombined;
            }
            else
                nFixedFactorCombinedMod = mpz_tdiv_ui(mpzFixedFactor.get_mpz_t(), nPrimeCombined);
        }

        unsigned int nFixedFactorMod = nFixedFactorCombinedMod % nPrime;
        if (nFixedFactorMod == 0)
        {
            // Nothing in the sieve is divisible by this prime
            continue;
        }
        // Find the modulo inverse of fixed factor
        unsigned int nFixedInverse = int_invert(nFixedFactorMod, nPrime);
        if (!nFixedInverse)
            return error("CSieveOfEratosthenes::Weave(): int_invert of fixed factor failed for prime #%u=%u", nPrimeSeqLocal, vPrimes[nPrimeSeqLocal]);
        unsigned int nTwoInverse = (nPrime + 1) / 2;
        
        // Check whether 32-bit arithmetic can be used for nFixedInverse
        const bool fUse32BArithmetic = (UINT_MAX_VALUE / nTwoInverse) >= nPrime;

        if (fUse32BArithmetic)
        {
            // Weave the sieve for the prime
            for (unsigned int nChainSeq = 0; nChainSeq < nSieveLayers; nChainSeq++)
            {
                // Find the first number that's divisible by this prime
                vCunningham1Multipliers[nPrimeSeqLocal * nSieveLayers + nChainSeq] = nFixedInverse;
                vCunningham2Multipliers[nPrimeSeqLocal * nSieveLayers + nChainSeq] = nPrime - nFixedInverse;

                // For next number in chain
                nFixedInverse = nFixedInverse * nTwoInverse % nPrime;
            }
        }
        else
        {
            // Weave the sieve for the prime
            for (unsigned int nChainSeq = 0; nChainSeq < nSieveLayers; nChainSeq++)
            {
                // Find the first number that's divisible by this prime
                vCunningham1Multipliers[nPrimeSeqLocal * nSieveLayers + nChainSeq] = nFixedInverse;
                vCunningham2Multipliers[nPrimeSeqLocal * nSieveLayers + nChainSeq] = nPrime - nFixedInverse;

                // For next number in chain
                nFixedInverse = (uint64)nFixedInverse * nTwoInverse % nPrime;
            }
        }
    }
    
#ifdef DEBUG_MINING    
    timeMark prepareEndPoint = getTimeMark();
#endif
    
    // Process the array in chunks that fit the L1 cache
    const unsigned int nArrayRounds = (nSieveSize + nL1CacheElements - 1) / nL1CacheElements;

    // Calculate the number of CC1 and CC2 layers needed for BiTwin candidates
    const unsigned int nBiTwinCC1Layers = (nChainLength + 1) / 2;
    const unsigned int nBiTwinCC2Layers = nChainLength / 2;

    // Only 50% of the array is used in extensions
    const unsigned int nExtensionsMinMultiplier = nSieveSize / 2;
    const unsigned int nExtensionsMinWord = nExtensionsMinMultiplier / nWordBits;

    // Loop over each array one at a time for optimal L1 cache performance
    for (unsigned int j = 0; j < nArrayRounds; j++)
    {
        const unsigned int nMinMultiplier = nL1CacheElements * j;
        const unsigned int nMaxMultiplier = std::min(nL1CacheElements * (j + 1), nSieveSize);
        const unsigned int nExtMinMultiplier = std::max(nMinMultiplier, nExtensionsMinMultiplier);
        const unsigned int nMinWord = nMinMultiplier / nWordBits;
        const unsigned int nMaxWord = (nMaxMultiplier + nWordBits - 1) / nWordBits;
        const unsigned int nExtMinWord = std::max(nMinWord, nExtensionsMinWord);
        if (pindexPrev != pindexBest)
            break;  // new block

        // Loop over the layers
        for (unsigned int nLayerSeq = 0; nLayerSeq < nSieveLayers; nLayerSeq++) {
            if (pindexPrev != pindexBest)
                break;  // new block
            if (nLayerSeq < nChainLength)
            {
                ProcessMultiplier(vfCompositeLayerCC1, nMinMultiplier, nMaxMultiplier, vPrimes, vCunningham1Multipliers, nLayerSeq);
                ProcessMultiplier(vfCompositeLayerCC2, nMinMultiplier, nMaxMultiplier, vPrimes, vCunningham2Multipliers, nLayerSeq);
            }
            else
            {
                // Optimize: First halves of the arrays are not needed in the extensions
                ProcessMultiplier(vfCompositeLayerCC1, nExtMinMultiplier, nMaxMultiplier, vPrimes, vCunningham1Multipliers, nLayerSeq);
                ProcessMultiplier(vfCompositeLayerCC2, nExtMinMultiplier, nMaxMultiplier, vPrimes, vCunningham2Multipliers, nLayerSeq);
            }

            // Apply the layer to the primary sieve arrays
            if (nLayerSeq < nChainLength)
            {
                if (nLayerSeq < nBiTwinCC2Layers)
                {
                    for (unsigned int nWord = nMinWord; nWord < nMaxWord; nWord++)
                    {
                        vfCompositeCunningham1[nWord] |= vfCompositeLayerCC1[nWord];
                        vfCompositeCunningham2[nWord] |= vfCompositeLayerCC2[nWord];
                        vfCompositeBiTwin[nWord] |= vfCompositeLayerCC1[nWord] | vfCompositeLayerCC2[nWord];
                    }
                }
                else if (nLayerSeq < nBiTwinCC1Layers)
                {
                    for (unsigned int nWord = nMinWord; nWord < nMaxWord; nWord++)
                    {
                        vfCompositeCunningham1[nWord] |= vfCompositeLayerCC1[nWord];
                        vfCompositeCunningham2[nWord] |= vfCompositeLayerCC2[nWord];
                        vfCompositeBiTwin[nWord] |= vfCompositeLayerCC1[nWord];
                    }
                }
                else
                {
                    for (unsigned int nWord = nMinWord; nWord < nMaxWord; nWord++)
                    {
                        vfCompositeCunningham1[nWord] |= vfCompositeLayerCC1[nWord];
                        vfCompositeCunningham2[nWord] |= vfCompositeLayerCC2[nWord];
                    }
                }
            }

            // Apply the layer to extensions
            for (unsigned int nExtensionSeq = 0; nExtensionSeq < nSieveExtensions; nExtensionSeq++)
            {
                const unsigned int nLayerOffset = nExtensionSeq + 1;
                if (nLayerSeq >= nLayerOffset && nLayerSeq < nChainLength + nLayerOffset)
                {
                    const unsigned int nLayerExtendedSeq = nLayerSeq - nLayerOffset;
                    sieve_word_t *vfExtCC1 = vfExtendedCompositeCunningham1 + nExtensionSeq * nCandidatesWords;
                    sieve_word_t *vfExtCC2 = vfExtendedCompositeCunningham2 + nExtensionSeq * nCandidatesWords;
                    sieve_word_t *vfExtTWN = vfExtendedCompositeBiTwin + nExtensionSeq * nCandidatesWords;

                    if (nLayerExtendedSeq < nBiTwinCC2Layers)
                    {
                        for (unsigned int nWord = nExtMinWord; nWord < nMaxWord; nWord++)
                        {
                            vfExtCC1[nWord] |= vfCompositeLayerCC1[nWord];
                            vfExtCC2[nWord] |= vfCompositeLayerCC2[nWord];
                            vfExtTWN[nWord] |= vfCompositeLayerCC1[nWord] | vfCompositeLayerCC2[nWord];
                        }
                    }
                    else if (nLayerExtendedSeq < nBiTwinCC1Layers)
                    {
                        for (unsigned int nWord = nExtMinWord; nWord < nMaxWord; nWord++)
                        {
                            vfExtCC1[nWord] |= vfCompositeLayerCC1[nWord];
                            vfExtCC2[nWord] |= vfCompositeLayerCC2[nWord];
                            vfExtTWN[nWord] |= vfCompositeLayerCC1[nWord];
                        }
                    }
                    else
                    {
                        for (unsigned int nWord = nExtMinWord; nWord < nMaxWord; nWord++)
                        {
                            vfExtCC1[nWord] |= vfCompositeLayerCC1[nWord];
                            vfExtCC2[nWord] |= vfCompositeLayerCC2[nWord];
                        }
                    }
                }
            }
        }

        // Combine the bitsets
        // vfCandidates = ~(vfCompositeCunningham1 & vfCompositeCunningham2 & vfCompositeBiTwin)
        for (unsigned int i = nMinWord; i < nMaxWord; i++)
            vfCandidates[i] = ~(vfCompositeCunningham1[i] & vfCompositeCunningham2[i] & vfCompositeBiTwin[i]);

        // Combine the extended bitsets
        for (unsigned int j = 0; j < nSieveExtensions; j++)
            for (unsigned int i = nExtMinWord; i < nMaxWord; i++)
                vfExtendedCandidates[j * nCandidatesWords + i] = ~(
                    vfExtendedCompositeCunningham1[j * nCandidatesWords + i] &
                    vfExtendedCompositeCunningham2[j * nCandidatesWords + i] &
                    vfExtendedCompositeBiTwin[j * nCandidatesWords + i]);
    }

    // The sieve has been partially weaved
    this->nPrimeSeq = nPrimes - 1;

#ifdef DEBUG_MINING    
    timeMark endPoint = getTimeMark();
    fprintf(stderr,
            "HP: prepare %.3lfmsec, sieve %.3lfmsec\n",
            usDiff(beginPoint, prepareEndPoint) / 1000.0,
            usDiff(prepareEndPoint, endPoint) / 1000.0);
#endif
    return false;
}