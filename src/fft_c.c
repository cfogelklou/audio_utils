/******************************************************************************
  Copyright 2014 Chris Fogelklou, Applaud Apps (Applicaudia)

  This source code may under NO circumstances be distributed without the 
  express written permission of the author.

  @author: chris.fogelklou@gmail.com
*******************************************************************************/ 
#include "fft_c.h"
//#include "audutils_debug.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define MALLOC malloc
#define FREE free
#define REALLOC realloc
#define ASSERT(x)
#define ASSERT_FN(x) if(x)do{;}while(0)
#define TRACE_WARNING(x)


//#define FFT_DBG
#ifdef FFT_DBG
#include <stdio.h>
static FILE * pf = 0;
#endif


#ifndef M_PI
#define M_PI (3.141592653589793238462643383279502884197169399375105820974944592307816406286)
#endif

#ifndef MREQ_TWO_PI
#define MREQ_TWO_PI (2*M_PI)
#endif

//#include "DaTunerApi_old.h"

bool fft_NumberOfBitsNeeded ( unsigned int nPowerOfTwo , unsigned int * pnBitsNeeded);
unsigned int fft_ReverseBits ( unsigned int nIndex, unsigned int nNumBits );
static bool fft_fftIfft (   
        FftFloat_t * const pFft, 
        const bool bInverseTransform,
        const fft_float_t * const pAdRealIn,
        const fft_float_t * const pAdImagIn,
        fft_float_t * const pAdRealOut,
        fft_float_t * const pAdImagOut,
        const unsigned int nNumSamples);
static bool fft_fftIfftf(
  FftFloat_t* const pFft,
  const bool bInverseTransform,
  const float* const pAdRealIn,
  fft_float_t* const pAdRealOut,
  fft_float_t* const pAdImagOut,
  const unsigned int nNumSamples);

static bool  fft_InitTwiddles( FftFloat_t *pFft );

//---------------------------------------------------------------------------
// Initializes or re-initializes the FFT.
FftFloat_t *Fft_Init( FftFloat_t *pFft, unsigned int nSize )
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
        pFft->allocatedLocally = true;
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
            unsigned int i = 0;    
            pRval->lastFftSize = nSize;
            fft_InitTwiddles( pRval );
            pRval->pReverseBitsLut = (unsigned int * )REALLOC( pRval->pReverseBitsLut, nSize * sizeof( unsigned int ) );
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
        const bool allocated = pFft->allocatedLocally;
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
bool  FFT_IsPowerOfTwo ( const int nX )
{
    return ((nX & -nX) == nX);
}


//---------------------------------------------------------------------------
bool  fft_NumberOfBitsNeeded(
  unsigned int nPowerOfTwo,
  unsigned int* pnBitsNeeded)
{

  unsigned int i = 0;
  bool success = false;
  unsigned int bitsNeeded = 0;

  if (nPowerOfTwo < 2) {
    bitsNeeded = 0;
  }

  while ((i < 32) && (success == false)){
    if (nPowerOfTwo & (1 << i)) {
      bitsNeeded = i;
      success = true;
    }
    ++i;
  }
  
  if (pnBitsNeeded){
    *pnBitsNeeded = bitsNeeded;
  }
  return success;
}



//---------------------------------------------------------------------------
unsigned int  fft_ReverseBits ( unsigned int nIndex, unsigned int nNumBits )
{
    unsigned int i = 0;
    unsigned int rev = 0;;

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
fft_float_t  FFT_IndexToFrequency ( unsigned int nNumSamples, unsigned int nIndex, fft_float_t fs )
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
bool   FFT_FFT (             
        FftFloat_t * const pFft, 
        const fft_float_t * const adRealIn, 
        const fft_float_t * const adImagIn, 
        fft_float_t * const adRealOut, 
        fft_float_t * const adImagOut, 
        unsigned int const nSize )

{
    bool status = fft_fftIfft( pFft, false, adRealIn, adImagIn, adRealOut, adImagOut, nSize);
#ifdef FFT_DBG
    {
        unsigned int i;
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
bool   FFT_FFTf(
  FftFloat_t* const pFft,
  const float* const adRealIn,
  fft_float_t* const adRealOut,
  fft_float_t* const adImagOut,
  unsigned int const nSize)

{
  bool status = fft_fftIfftf(pFft, false, adRealIn, adRealOut, adImagOut, nSize);
  return status;
}

//---------------------------------------------------------------------------
bool   FFT_IFFT ( 
        FftFloat_t *pFft, 
        const fft_float_t * const adRealIn, 
        const fft_float_t * const adImagIn, 
        fft_float_t * const adRealOut, 
        fft_float_t * const adImagOut, 
        unsigned int const nSize )
{
    return fft_fftIfft(pFft, true, adRealIn, adImagIn, adRealOut, adImagOut, nSize);
}


//---------------------------------------------------------------------------
bool   FFT_Phase( fft_float_t *pAdRealIn, fft_float_t *pAdImagIn, fft_float_t *adPhase,  unsigned int nSize )
{
    unsigned int i;
    bool status = false;
    ASSERT ((pAdRealIn != NULL) && (adPhase != NULL));
    if (pAdImagIn == NULL) 
    {
        for (i = 0; i < nSize; i++) 
        {
            adPhase[i] = 0.0;
        }
        status = true;
    }
    else
    {
        for (i = 0; i < nSize; i++) 
        {
            adPhase[i] = (pAdRealIn[i] == 0) ? M_PI/2 : atan(pAdImagIn[i]/pAdRealIn[i]);
        }
        status = true;
    }
    return status;
}


//---------------------------------------------------------------------------
bool   FFT_Magnitude( fft_float_t *pAdRealIn, fft_float_t *pAdImagIn, fft_float_t *adMagnitude,  unsigned int nSize )
{
    unsigned int i;
    ASSERT ((pAdRealIn != NULL) && (adMagnitude != NULL));
    for (i = 0; i < nSize; i++) 
    {
        adMagnitude[i] = sqrt(pow(pAdRealIn[i], 2) + pow(pAdImagIn[i], 2));
    }
    return true;
}

//---------------------------------------------------------------------------
bool FFT_MagnitudePhase( fft_float_t *pAdRealIn, fft_float_t *pAdImagIn, fft_float_t *pAdMagnitude, fft_float_t *pAdPhase, unsigned int nSize )
{
    unsigned int i;
    ASSERT ((pAdRealIn != NULL) && (pAdMagnitude != NULL) && (pAdPhase != NULL));
    for (i = 0; i < nSize; i++) 
    {
        pAdMagnitude[i] = sqrt(pow(pAdRealIn[i], 2) + pow(pAdImagIn[i], 2));
        pAdPhase[i] = (pAdRealIn[i] == 0) ? M_PI/2 : atan(pAdImagIn[i]/pAdRealIn[i]);
    }
    return true;
}

static bool fft_afterReverse(
  FftFloat_t* const pFft,
  const bool bInverseTransform,
  fft_float_t* const xr,
  fft_float_t* const xi
) {
  if (pFft != 0) {
    const unsigned int nNumSamples = pFft->lastFftSize;
    // Declare some local variables and start the FFT.
    {
      unsigned int nBlockSize;

      const FftLut_t* const pLut = (!bInverseTransform) ? pFft->pFftLut : pFft->pIfftLut;
      unsigned int iter = 0;
      unsigned int nBlockEnd = 1;

#ifdef FFT_DBG
      unsigned int innerloops = 0;
#endif

      for (nBlockSize = 2; nBlockSize <= nNumSamples; nBlockSize <<= 1) {
        unsigned int i;
        const fft_float_t sm2 = pLut[iter].sm2;
        const fft_float_t sm1 = pLut[iter].sm1;
        const fft_float_t cm2 = pLut[iter].cm2;
        const fft_float_t cm1 = pLut[iter].cm1;
        const fft_float_t w = 2 * cm1;
        fft_float_t ar[3], ai[3];
        ++iter;

        for (i = 0; i < nNumSamples; i += nBlockSize) {
          unsigned int j, n;

          ar[2] = cm2;
          ar[1] = cm1;

          ai[2] = sm2;
          ai[1] = sm1;

          for (j = i, n = 0; n < nBlockEnd; j++, n++) {

            fft_float_t tr, ti;     /* temp real, temp imaginary */
            const unsigned int k = j + nBlockEnd;

            ar[0] = w * ar[1] - ar[2];
            ar[2] = ar[1];
            ar[1] = ar[0];

            ai[0] = w * ai[1] - ai[2];
            ai[2] = ai[1];
            ai[1] = ai[0];


            tr = ar[0] * xr[k] - ai[0] * xi[k];
            ti = ar[0] * xi[k] + ai[0] * xr[k];

#ifdef FFT_DBG
            if (pf != 0)
            {
              fprintf(pf, "Innerloops: %d\n", innerloops++);
              fprintf(pf, "tr = ar * xr[ k=%d ] - ai * xi[ k=%d ] = (%f * %f) - (%f * %f) = %f\n", k, k, ar[0], xr[k], ai[0], xi[k], tr);
              fprintf(pf, "ti = ar * xi[ k=%d ] + ai * xr[ k=%d ] = (%f * %f) + (%f * %f) = %f\n", k, k, ar[0], xi[k], ai[0], xr[k], ti);
            }
#endif

            xr[k] = xr[j] - tr;
            xi[k] = xi[j] - ti;

#ifdef FFT_DBG
            if (pf != 0)
            {
              fprintf(pf, "xr[ k=%d ] = xr[ j=%d ] - tr = %f - (%f) = %f\n", k, j, xr[j], tr, xr[k]);
              fprintf(pf, "xi[ k=%d ] = xi[ j=%d ] - ti = %f - (%f) = %f\n", k, j, xi[j], ti, xi[k]);
            }
#endif

            {
#ifdef FFT_DBG
              const fft_float_t oldxrj = xr[j];
              const fft_float_t oldxij = xi[j];
#endif
              xr[j] = xr[j] + tr;
              xi[j] = xi[j] + ti;

#ifdef FFT_DBG
              if (pf != 0)
              {
                fprintf(pf, "xr[ j=%d ] = xr[ j=%d ] + tr = %f + %f = %f\n", j, j, oldxrj, tr, xr[j]);
                fprintf(pf, "xi[ j=%d ] = xi[ j=%d ] + ti = %f + %f = %f\n", j, j, oldxij, ti, xi[j]);
              }
#endif
            }
          }
        }

        nBlockEnd = nBlockSize;
      }
    }

    // Normalize
    if (bInverseTransform)
    {
      fft_float_t dDenom = 1.0 / (fft_float_t)nNumSamples;

      for (unsigned int i = 0; i < nNumSamples; i++)
      {
        xr[i] *= dDenom;
        xi[i] *= dDenom;
      }
    }
  }
  return true;
}

//---------------------------------------------------------------------------
bool fft_fftIfft(
  FftFloat_t* const pFft,
  const bool bInverseTransform,
  const fft_float_t* const pXr,
  const fft_float_t* const pXi,
  fft_float_t* const xr,
  fft_float_t* const xi,
  const unsigned int nNumSamples)
{
  if (pFft != 0) {
    unsigned int i;
    const fft_float_t* const pAdRealIn = pXr;
    const fft_float_t* const pAdImagIn = pXi;
    ASSERT((nNumSamples == pFft->lastFftSize) 
      && (pAdRealIn != NULL) && (xr != NULL) && (xi != NULL));
    const unsigned int* const pReverseBitsLut = pFft->pReverseBitsLut;

    // Reverse ordering of samples so FFT can be done in place.
    if (pAdImagIn == NULL) {
      for (i = 0; i < nNumSamples; i++) {
        const unsigned int rev = pReverseBitsLut[i];
        xr[rev] = pAdRealIn[i];
      }
      memset(xi, 0, sizeof(fft_float_t) * nNumSamples);
    }
    else {
      for (i = 0; i < nNumSamples; i++) {
        const unsigned int rev = pReverseBitsLut[i];
        xr[rev] = pAdRealIn[i];
        xi[rev] = pAdImagIn[i];
      }
    }
  }

  return fft_afterReverse(pFft, bInverseTransform, xr, xi);
}

//---------------------------------------------------------------------------
bool fft_fftIfftf(
  FftFloat_t* const pFft,
  const bool bInverseTransform,
  const float* const pXr,
  fft_float_t* const xr,
  fft_float_t* const xi,
  const unsigned int nNumSamples)
{
  if (pFft != 0) {
    unsigned int i;
    ASSERT((nNumSamples == pFft->lastFftSize)
      && (pAdRealIn != NULL) && (xr != NULL) && (xi != NULL));
    const unsigned int* const pReverseBitsLut = pFft->pReverseBitsLut;

    // Reverse ordering of samples so FFT can be done in place.
    for (i = 0; i < nNumSamples; i++) {
      const unsigned int rev = pReverseBitsLut[i];
      xr[rev] = pXr[i];
    }
    memset(xi, 0, sizeof(fft_float_t) * nNumSamples);
  }

  return fft_afterReverse(pFft, bInverseTransform, xr, xi);
}

//---------------------------------------------------------------------------
static bool  fft_InitTwiddles( FftFloat_t *pFft )
{
    unsigned int i;
    unsigned int nBlockSize;
    unsigned int lutSize = 0;
    unsigned int numSamples = pFft->lastFftSize;

    for ( nBlockSize = 2; nBlockSize <= numSamples; nBlockSize <<= 1 )
    {
        ++lutSize;
    }

    pFft->pFftLut  = (FftLut_t *)REALLOC( pFft->pFftLut, sizeof( FftLut_t ) * lutSize );
    pFft->pIfftLut = (FftLut_t *)REALLOC( pFft->pIfftLut, sizeof( FftLut_t )  * lutSize );
    
    // Do once for regular transform, once for inverse
    for ( i = 0; i < 2; i++)
    {
        unsigned int j = 0;
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

    return true;
}

unsigned int fft_IntegerLog2(unsigned int x)
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
