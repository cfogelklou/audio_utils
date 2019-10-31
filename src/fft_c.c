/******************************************************************************
  Copyright 2014 Chris Fogelklou, Applaud Apps (Applicaudia)

  This source code may under NO circumstances be distributed without the 
  express written permission of the author.

  @author: chris.fogelklou@gmail.com
*******************************************************************************/ 
#include "fft_c.h"
#include "audutils_debug.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define MALLOC malloc
#define FREE free
#define REALLOC realloc

//#define FFT_DBG
#ifdef FFT_DBG
#include <stdio.h>
static FILE * pf = 0;
#endif

//#include "DaTunerApi_old.h"

bool_t fft_NumberOfBitsNeeded ( uint_t nPowerOfTwo , uint_t * pnBitsNeeded);
uint_t fft_ReverseBits ( uint_t nIndex, uint_t nNumBits );
static bool_t fft_fftIfft (   
        FftFloat_t * const pFft, 
        const bool_t bInverseTransform,
        const fft_float_t * const pAdRealIn,
        const fft_float_t * const pAdImagIn,
        fft_float_t * const pAdRealOut,
        fft_float_t * const pAdImagOut,
        const uint_t nNumSamples);
static bool_t  fft_InitTwiddles( FftFloat_t *pFft );

//---------------------------------------------------------------------------
// Initializes or re-initializes the FFT.
FftFloat_t *Fft_Init( FftFloat_t *pFft, uint_t nSize )
{
    FftFloat_t *pRval = 0;
    if (!( FFT_IsPowerOfTwo( nSize ) ) )
    {
        return pRval;
    }
    if (NULL == pFft)
    {
        pFft = (FftFloat_t *)MALLOC( sizeof( FftFloat_t ) );
        memset( pFft, 0, sizeof( FftFloat_t ) );
        pFft->allocatedLocally = TRUE;
    }
    ASSERT( pFft != NULL );
    if (pFft != NULL)
    {
        pRval = pFft;
        
        if (!(FFT_IsPowerOfTwo(pRval->lastFftSize)))
        {
            // Assume struct is uninitialized
            pRval->lastFftSize = 0;
            pRval->pReverseBitsLut = 0;
            TRACE_WARNING(("Warning: Uninitialized fft struct passed in."));
        }

        if (nSize != pRval->lastFftSize)
        {
            uint_t i = 0;    
            pRval->lastFftSize = nSize;
            fft_InitTwiddles( pRval );
            pRval->pReverseBitsLut = (uint_t * )REALLOC( pRval->pReverseBitsLut, nSize * sizeof( uint_t ) );
            ASSERT( pRval->pReverseBitsLut );
            ASSERT_FN( fft_NumberOfBitsNeeded ( nSize, &pRval->bitsNeeded ));

            for ( i = 0; i < nSize; i++ ) 
            {
                pRval->pReverseBitsLut[ i ] = fft_ReverseBits ( i, pRval->bitsNeeded );
            }    
        }
    }
#ifdef FFT_DBG
    if (pf == 0)
    {
        pf = fopen("fft_dbg.csv", "w");
    }
    if (pf != 0)
    {
        fprintf( pf, "tIdx, bIdx, pAdRealOut[tIdx], pAdImagOut[tIdx], pAdRealOut[bIdx], pAdImagOut[bIdx]\n" );
    }
#endif
    return pRval;
}


//---------------------------------------------------------------------------
void Fft_DeInit( FftFloat_t *pFft )
{
    if (pFft != NULL)
    {
        const bool_t allocated = pFft->allocatedLocally;
        if (pFft->pReverseBitsLut != NULL)
        {
            FREE( pFft->pReverseBitsLut);
            pFft->pReverseBitsLut = NULL;
        }
        if (pFft->pFftLut != NULL)
        {
            FREE( pFft->pFftLut);
            pFft->pFftLut = NULL;
        }
        if (pFft->pIfftLut != NULL)
        {
            FREE( pFft->pIfftLut);
            pFft->pIfftLut = NULL;
        }
        memset( pFft, 0, sizeof( FftFloat_t ) );
        if (allocated)
        {
            FREE( pFft );
        }
    }
#ifdef FFT_DBG
    if (pf != 0)
    {
        fclose(pf);
        pf = 0;
    }
#endif
}


//---------------------------------------------------------------------------
bool_t  FFT_IsPowerOfTwo ( const int_t nX )
{
    return ((nX & -nX) == nX);
}


//---------------------------------------------------------------------------
bool_t  fft_NumberOfBitsNeeded ( uint_t nPowerOfTwo , uint_t * pnBitsNeeded )
{
    
    uint_t i = 0;
    bool_t success = FALSE;
    uint_t bitsNeeded = 0;

    if ( nPowerOfTwo < 2 )
    {
        bitsNeeded = 0;
    }

    while ((i < 32) && (success == FALSE))
    {
        if ( nPowerOfTwo & (1 << i) ) 
        {
            bitsNeeded = i;
            success = TRUE;
        }
        ++i;
    }
    if (pnBitsNeeded)
    {
        *pnBitsNeeded = bitsNeeded;
    }
    return success;
}



//---------------------------------------------------------------------------
uint_t  fft_ReverseBits ( uint_t nIndex, uint_t nNumBits )
{
    uint_t i = 0;
    uint_t rev = 0;;

    if (nIndex != 0)
    {
        for ( ; i < nNumBits; i++ )
        {
            rev = (rev << 1) | (nIndex & 1);
            nIndex >>= 1;
        }
    }

    return rev;
}


//---------------------------------------------------------------------------
fft_float_t  FFT_IndexToFrequency ( uint_t nNumSamples, uint_t nIndex, fft_float_t fs )
{
    fft_float_t rval = 0.0;
    if ( nIndex < nNumSamples/2 ) 
    {
        rval = (fft_float_t)nIndex / (fft_float_t)nNumSamples;
    }
    else
    {
        rval = -((fft_float_t)(nNumSamples-nIndex) / (fft_float_t)nNumSamples);
    }
    rval *= fs;
    return rval;
}


//---------------------------------------------------------------------------
bool_t   FFT_FFT (             
        FftFloat_t * const pFft, 
        const fft_float_t * const adRealIn, 
        const fft_float_t * const adImagIn, 
        fft_float_t * const adRealOut, 
        fft_float_t * const adImagOut, 
        uint_t const nSize )

{
    bool_t status = fft_fftIfft( pFft, FALSE, adRealIn, adImagIn, adRealOut, adImagOut, nSize);
#ifdef FFT_DBG
    {
        uint_t i;
        for ( i = 0; i < nSize; i++ ) 
        {
            if (pf != 0)
            {
                fprintf(pf, "xr[ i=%d ], xi[ i=%d ] = %d, %d\n", i, i, (int32_t)adRealOut[i], (int32_t)adImagOut[i] );
            }
        }   
    }
#endif
    return status;
}


//---------------------------------------------------------------------------
bool_t   FFT_IFFT ( 
        FftFloat_t *pFft, 
        const fft_float_t * const adRealIn, 
        const fft_float_t * const adImagIn, 
        fft_float_t * const adRealOut, 
        fft_float_t * const adImagOut, 
        uint_t const nSize )
{
    return fft_fftIfft(pFft, TRUE, adRealIn, adImagIn, adRealOut, adImagOut, nSize);
}


//---------------------------------------------------------------------------
bool_t   FFT_Phase( fft_float_t *pAdRealIn, fft_float_t *pAdImagIn, fft_float_t *adPhase,  uint_t nSize )
{
    uint_t i;
    bool_t status = FALSE;
    ASSERT ((pAdRealIn != NULL) && (adPhase != NULL));
    if (pAdImagIn == NULL) 
    {
        for (i = 0; i < nSize; i++) 
        {
            adPhase[i] = 0.0;
        }
        status = TRUE;
    }
    else
    {
        for (i = 0; i < nSize; i++) 
        {
            adPhase[i] = (pAdRealIn[i] == 0) ? MREQ_PI/2 : atan(pAdImagIn[i]/pAdRealIn[i]);
        }
        status = TRUE;
    }
    return status;
}


//---------------------------------------------------------------------------
bool_t   FFT_Magnitude( fft_float_t *pAdRealIn, fft_float_t *pAdImagIn, fft_float_t *adMagnitude,  uint_t nSize )
{
    uint_t i;
    ASSERT ((pAdRealIn != NULL) && (adMagnitude != NULL));
    for (i = 0; i < nSize; i++) 
    {
        adMagnitude[i] = sqrt(pow(pAdRealIn[i], 2) + pow(pAdImagIn[i], 2));
    }
    return TRUE;
}

//---------------------------------------------------------------------------
bool_t FFT_MagnitudePhase( fft_float_t *pAdRealIn, fft_float_t *pAdImagIn, fft_float_t *pAdMagnitude, fft_float_t *pAdPhase, uint_t nSize )
{
    uint_t i;
    ASSERT ((pAdRealIn != NULL) && (pAdMagnitude != NULL) && (pAdPhase != NULL));
    for (i = 0; i < nSize; i++) 
    {
        pAdMagnitude[i] = sqrt(pow(pAdRealIn[i], 2) + pow(pAdImagIn[i], 2));
        pAdPhase[i] = (pAdRealIn[i] == 0) ? MREQ_PI/2 : atan(pAdImagIn[i]/pAdRealIn[i]);
    }
    return TRUE;
}

//---------------------------------------------------------------------------
bool_t fft_fftIfft (   
        FftFloat_t * const pFft, 
        const bool_t bInverseTransform,
        const fft_float_t * const pAdRealIn,
        const fft_float_t * const pAdImagIn,
        fft_float_t * const xr,
        fft_float_t * const xi,
        const uint_t nNumSamples)
{
    if (pFft != 0)
    {
        uint_t i;
        ASSERT( ( nNumSamples == pFft->lastFftSize ) && (pAdRealIn != NULL) && (xr != NULL) && (xi != NULL) );
        const uint_t* const pReverseBitsLut = pFft->pReverseBitsLut;

        // Reverse ordering of samples so FFT can be done in place.
        if (pAdImagIn == NULL)
        {
            for ( i = 0; i < nNumSamples; i++ ) 
            {
                const uint_t rev = pReverseBitsLut[ i ];
                xr[rev] = pAdRealIn[i];
            }
            memset(xi, 0, sizeof(fft_float_t) * nNumSamples);
        }
        else
        {
            for ( i = 0; i < nNumSamples; i++ ) 
            {
                const uint_t rev = pReverseBitsLut[ i ];
                xr[rev] = pAdRealIn[i];
                xi[rev] = pAdImagIn[i];
            }
        }

        // Declare some local variables and start the FFT.
        {
            uint_t nBlockSize;

            const FftLut_t * const pLut = (!bInverseTransform) ? pFft->pFftLut : pFft->pIfftLut;
            uint_t iter = 0;
            uint_t nBlockEnd = 1;

#ifdef FFT_DBG
            uint_t innerloops = 0;
#endif

            for ( nBlockSize = 2; nBlockSize <= nNumSamples; nBlockSize <<= 1 )
            {
                uint_t i;
                const fft_float_t sm2 = pLut[ iter ].sm2;
                const fft_float_t sm1 = pLut[ iter ].sm1;
                const fft_float_t cm2 = pLut[ iter ].cm2;
                const fft_float_t cm1 = pLut[ iter ].cm1;
                const fft_float_t w = 2 * cm1;
                fft_float_t ar[3], ai[3];
                ++iter;

                for ( i=0; i < nNumSamples; i += nBlockSize )
                {
                    uint_t j, n;

                    ar[2] = cm2;
                    ar[1] = cm1;

                    ai[2] = sm2;
                    ai[1] = sm1;

                    for ( j = i, n = 0; n < nBlockEnd; j++, n++ )
                    {

                        fft_float_t tr, ti;     /* temp real, temp imaginary */
                        const uint_t k = j + nBlockEnd;

                        ar[0] = w * ar[1] - ar[2];
                        ar[2] = ar[1];
                        ar[1] = ar[0];

                        ai[0] = w * ai[1] - ai[2];
                        ai[2] = ai[1];
                        ai[1] = ai[0];


                        tr = ar[ 0 ] * xr[ k ] - ai[ 0 ] * xi[ k ];
                        ti = ar[ 0 ] * xi[ k ] + ai[ 0 ] * xr[ k ];

#ifdef FFT_DBG
                        if (pf != 0)
                        {
                            fprintf(pf,"Innerloops: %d\n", innerloops++);
                            fprintf(pf, "tr = ar * xr[ k=%d ] - ai * xi[ k=%d ] = (%f * %f) - (%f * %f) = %f\n", k, k, ar[ 0 ] , xr[ k ], ai[ 0 ] , xi[ k ], tr);
                            fprintf(pf, "ti = ar * xi[ k=%d ] + ai * xr[ k=%d ] = (%f * %f) + (%f * %f) = %f\n", k, k, ar[ 0 ] , xi[ k ], ai[ 0 ] , xr[ k ], ti);
                        }
#endif

                        xr[ k ] = xr[ j ] - tr;
                        xi[ k ] = xi[ j ] - ti;

#ifdef FFT_DBG
                        if (pf != 0)
                        {
                            fprintf(pf, "xr[ k=%d ] = xr[ j=%d ] - tr = %f - (%f) = %f\n", k, j, xr[ j ], tr, xr[k] );
                            fprintf(pf, "xi[ k=%d ] = xi[ j=%d ] - ti = %f - (%f) = %f\n", k, j, xi[ j ], ti, xi[k] );
                        }
#endif

                        {
#ifdef FFT_DBG
                            const fft_float_t oldxrj = xr[ j ];
                            const fft_float_t oldxij = xi[ j ];
#endif
                            xr[ j ] = xr[ j ] + tr;
                            xi[ j ] = xi[ j ] + ti;

#ifdef FFT_DBG
                            if (pf != 0)
                            {
                                fprintf(pf, "xr[ j=%d ] = xr[ j=%d ] + tr = %f + %f = %f\n", j, j, oldxrj, tr, xr[ j ] );
                                fprintf(pf, "xi[ j=%d ] = xi[ j=%d ] + ti = %f + %f = %f\n", j, j, oldxij, ti, xi[ j ] );
                            }
#endif
                        }
                    }
                }

                nBlockEnd = nBlockSize;
            }
        }

        // Normalize
        if ( bInverseTransform )
        {
            fft_float_t dDenom = 1.0 / (fft_float_t)nNumSamples;

            for ( i=0; i < nNumSamples; i++ )
            {
                xr[i] *= dDenom;
                xi[i] *= dDenom;
            }
        }
    }
    return TRUE;
}


//---------------------------------------------------------------------------
static bool_t  fft_InitTwiddles( FftFloat_t *pFft )
{
    uint_t i;
    uint_t nBlockSize;
    uint_t lutSize = 0;
    uint_t numSamples = pFft->lastFftSize;

    for ( nBlockSize = 2; nBlockSize <= numSamples; nBlockSize <<= 1 )
    {
        ++lutSize;
    }

    pFft->pFftLut  = (FftLut_t *)REALLOC( pFft->pFftLut, sizeof( FftLut_t ) * lutSize );
    pFft->pIfftLut = (FftLut_t *)REALLOC( pFft->pIfftLut, sizeof( FftLut_t )  * lutSize );
    
    // Do once for regular transform, once for inverse
    for ( i = 0; i < 2; i++)
    {
        uint_t j = 0;
        // lookup tables for sm and cm.
        FftLut_t *pLut = (i == 0) ? pFft->pFftLut : pFft->pIfftLut;
        // Inverse, not inverse
        fft_float_t dAngleNumerator = ( i == 0 ) ? -MREQ_TWO_PI : MREQ_TWO_PI;
        for ( nBlockSize = 2; nBlockSize <= numSamples; nBlockSize <<= 1 )
        {
            fft_float_t dDeltaAngle = dAngleNumerator / (fft_float_t)nBlockSize;
            pLut[j].sm2 = sin ( -2 * dDeltaAngle );
            pLut[j].sm1 = sin ( -dDeltaAngle );
            pLut[j].cm2 = cos ( -2 * dDeltaAngle );
            pLut[j].cm1 = cos ( -dDeltaAngle );
            ++j;
        }
    }

    return TRUE;
}

uint_t fft_IntegerLog2(uint_t x)
{
  int pos = 0;
  int64_t n = x;
  n &= 0x00000000ffffffff;
  if (n >= (int64_t)((int64_t)1u) << 32) {
    n >>= 32;
    pos += 32;
  }
  if (n >= 1 << 16) {
    n >>= 16;
    pos += 16;
  }
  if (n >= 1 << 8) {
    n >>= 8;
    pos += 8;
  }
  if (n >= 1 << 4) {
    n >>= 4;
    pos += 4;
  }
  if (n >= 1 << 2) {
    n >>= 2;
    pos += 2;
  }
  if (n >= 1 << 1) {
    pos += 1;
  }
  return ((n == 0) ? (-1) : pos);
}
