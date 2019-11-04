/******************************************************************************
  Copyright 2014 Chris Fogelklou, Applaud Apps (Applicaudia)

  This source code may under NO circumstances be distributed without the 
  express written permission of the author.

  @author: chris.fogelklou@gmail.com
*******************************************************************************/ 

#ifndef __fft
#define __fft

//#include "helper_math.hpp"
#include <stdint.h>
#ifndef __cplusplus
#include <stdbool.h>
#endif

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
  unsigned int lastFftSize;
  unsigned int bitsNeeded;
  unsigned int *pReverseBitsLut;
  FftLut_t *pFftLut;
  FftLut_t *pIfftLut;
  bool allocatedLocally;
} FftFloat_t;


#if defined( __cplusplus )
extern "C" {
#endif

  bool FFT_FFT ( 
    FftFloat_t * const pFft, 
    const fft_float_t * const adRealIn, 
    const fft_float_t * const adImagIn, 
    fft_float_t * const adRealOut, 
    fft_float_t * const adImagOut, 
    unsigned int const nSize);

  bool   FFT_FFTf(
    FftFloat_t* const pFft,
    const float* const adRealIn,
    fft_float_t* const adRealOut,
    fft_float_t* const adImagOut,
    unsigned int const nSize);

  bool FFT_IFFT (             
    FftFloat_t *pFft, 
    const fft_float_t * const adRealIn, 
    const fft_float_t * const adImagIn, 
    fft_float_t * const adRealOut, 
    fft_float_t * const adImagOut, 
    unsigned int const nSize);

  bool FFT_Phase( fft_float_t *adRealIn, fft_float_t *adImagIn, fft_float_t *adPhase,  unsigned int nSize );
  bool FFT_Magnitude( fft_float_t *adRealIn, fft_float_t *adImagIn, fft_float_t *adMagnitude,  unsigned int nSize );
  bool FFT_MagnitudePhase( fft_float_t *pAdRealIn, fft_float_t *pAdImagIn, fft_float_t *pAdMagnitude, fft_float_t *pAdPhase, unsigned int nSize );
  bool FFT_IsPowerOfTwo ( const int nX );
  fft_float_t  FFT_IndexToFrequency ( unsigned int nNumSamples, unsigned int nIndex, fft_float_t fs );

  // Initializes or re-initializes the FFT.
  FftFloat_t *Fft_Init( FftFloat_t *pFft, unsigned int nSize );

  void Fft_DeInit( FftFloat_t *pFft );

  bool fft_NumberOfBitsNeeded ( unsigned int nPowerOfTwo , unsigned int * pnBitsNeeded);
  unsigned int fft_ReverseBits ( unsigned int nIndex, unsigned int nNumBits );

  unsigned int fft_IntegerLog2(unsigned int x);

#if defined( __cplusplus )
}
#endif

#endif
