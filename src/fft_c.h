/******************************************************************************
  Copyright 2014 Chris Fogelklou, Applaud Apps (Applicaudia)

  This source code may under NO circumstances be distributed without the 
  express written permission of the author.

  @author: chris.fogelklou@gmail.com
*******************************************************************************/ 

#ifndef __fft
#define __fft

#include "HelperMath.hpp"

#define fft_float_t double

typedef struct
{
  fft_float_t sm2;
  fft_float_t sm1;
  fft_float_t cm2;
  fft_float_t cm1;
} FftLut_t;

typedef struct _FftFloat_t
{
  uint_t lastFftSize;
  uint_t bitsNeeded;
  uint_t *pReverseBitsLut;
  FftLut_t *pFftLut;
  FftLut_t *pIfftLut;
  bool_t allocatedLocally;
} FftFloat_t;


#if defined( __cplusplus )
extern "C" {
#endif

  bool_t FFT_FFT ( 
    FftFloat_t * const pFft, 
    const fft_float_t * const adRealIn, 
    const fft_float_t * const adImagIn, 
    fft_float_t * const adRealOut, 
    fft_float_t * const adImagOut, 
    uint_t const nSize);

  bool_t FFT_IFFT (             
    FftFloat_t *pFft, 
    const fft_float_t * const adRealIn, 
    const fft_float_t * const adImagIn, 
    fft_float_t * const adRealOut, 
    fft_float_t * const adImagOut, 
    uint_t const nSize);

  bool_t FFT_Phase( fft_float_t *adRealIn, fft_float_t *adImagIn, fft_float_t *adPhase,  uint_t nSize );
  bool_t FFT_Magnitude( fft_float_t *adRealIn, fft_float_t *adImagIn, fft_float_t *adMagnitude,  uint_t nSize );
  bool_t FFT_MagnitudePhase( fft_float_t *pAdRealIn, fft_float_t *pAdImagIn, fft_float_t *pAdMagnitude, fft_float_t *pAdPhase, uint_t nSize );
  bool_t FFT_IsPowerOfTwo ( const int_t nX );
  fft_float_t  FFT_IndexToFrequency ( uint_t nNumSamples, uint_t nIndex, fft_float_t fs );

  // Initializes or re-initializes the FFT.
  FftFloat_t *Fft_Init( FftFloat_t *pFft, uint_t nSize );

  void Fft_DeInit( FftFloat_t *pFft );

  bool_t fft_NumberOfBitsNeeded ( uint_t nPowerOfTwo , uint_t * pnBitsNeeded);
  uint_t fft_ReverseBits ( uint_t nIndex, uint_t nNumBits );

  uint_t fft_IntegerLog2(uint_t x);

#if defined( __cplusplus )
}
#endif

#endif
