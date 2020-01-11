/******************************************************************************
  Copyright 2014 Chris Fogelklou, Applaud Apps (Applicaudia)

  This source code may under NO circumstances be distributed without the
  express written permission of the author.

  @author: chris.fogelklou@gmail.com
*******************************************************************************/
#ifndef FFTCLASS_HPP
#define FFTCLASS_HPP

#include <string.h>
#include <math.h>

#include "audiolib_types.h"
#include "fft_c.h"
#include "audutils_debug.h"
#include "pcm_q.h"

//#define USE_FAST_MATH

#ifdef USE_FAST_MATH
typedef float fftc_t;
#else
typedef double fftc_t;
#endif

typedef fftc_t fftin_t;

#ifndef E_FFTSIZE_MAX
#define E_FFTSIZE_MAX 2048
#endif

namespace Fft {

  static const double SQRT2 = sqrt(2.0);

  bool getPower(double pAdRealIn[], double pAdImagIn[], double adPow[], int nSize);

  bool getAmplitude(double pAdRealIn[], double pAdImagIn[], double adPow[], int nSize);

  double getFftPeakAmpFromPowers(double pAdPowers[], int fftSize);

  double getFftPeakAmpFromAmplitudes(double pAdAmplitudes[], int fftSize);

  double getFftRmsFromPowers(double pAdPowers[], int fftSize);

  double getFftRmsFromAmplitudes(double pAdAmplitudes[], int fftSize);

  int findPeak(
    double adMagnitude[], 
    int magnitudeLen, 
    int nMin, 
    int nMax);

  int  bitsNeeded(
    const unsigned int nPowerOfTwo);
};

//---------------------------------------------------------------------------
// An FFT class
class FftClass {


public:

  //---------------------------------------------------------------------------
  // Initializes or re-initializes the FFT.
  FftClass(const int nFftSize = 512);

  virtual ~FftClass();

  void setLutSize(const int fftSize = 512);

  // ---------------------------------------------------------------------------
  int getSize();

  // ---------------------------------------------------------------------------
  void setFftSize(int nFftSize);

  //---------------------------------------------------------------------------
  bool doFft(
    fftin_t adRealIn[], 
    fftin_t adImagIn[], 
    double adRealOut[],
    double adImagOut[], 
    int nSize);

  //---------------------------------------------------------------------------
  bool doIfft(
    fftin_t adRealIn[], 
    fftin_t adImagIn[], 
    double adRealOut[],
    double adImagOut[], 
    int nSize);

private:
  inline fftc_t getSin(const int tblIdx) {
    const int tblIdxMinusFftSizeOver2 = tblIdx - (mFftSize / 2);
    return (tblIdxMinusFftSizeOver2 > 0) ? -mSinLut[tblIdxMinusFftSizeOver2]
                                         : mSinLut[tblIdx];
  }

  inline fftc_t getCos(const int tblIdx) {
    const int newIdx = (tblIdx - (mFftSize / 4)) & (mFftSize - 1);
    return getSin(newIdx);
  }

  //---------------------------------------------------------------------------

  bool doFftOrIfft(
    bool bInverseTransform, 
    fftin_t pAdRealIn[],
    fftin_t pAdImagIn[], 
    double _xr[], 
    double _xi[],
    int nNumSamples);

public:
  // ---------------------------------------------------------------------------
  void doRealInverseFft(fftin_t pAdRealIn[], double xr[], double xi[]);

  // ---------------------------------------------------------------------------
  template <typename TIn, typename Tout>
  bool doRealForwardFft(TIn pAdRealIn[], Tout xr[], Tout xi[]) {
    ASSERT((pAdRealIn != NULL) && (xr != NULL) && (xi != NULL));

    // Reverse ordering of samples so FFT can be done in place.
    for (int i = 0; i < mFftSize; i++) {
      xr[mReverseBitsLut[i]] = (Tout)pAdRealIn[i];
    }
    memset(xi, 0, sizeof(Tout) * mFftSize);

    // int p = mBitsNeeded;
    int Bp = 1 << (mBitsNeeded - 1);
    int Np = 2;

    // First loop only operates on real inputs because xi is all zeroes.
    int twiddleMul = (mFftSize >> 1);
    {
      int BaseT = 0;
      const int Npp = Np >> 1;
      for (int b = Bp; b > 0; b--) {
        const int BaseB = BaseT + Npp;
        for (int k_ = 0; k_ < Npp; k_++) {
          const int j = BaseT + k_;
          const int k = BaseB + k_;
          const int tblIdx = k_ * twiddleMul;
          const Tout xreal = xr[k] * (Tout)getCos(tblIdx);
          const Tout ximag = xr[k] * (Tout)getSin(tblIdx);
          xr[k] = xr[j] - xreal;
          xi[k] = xi[j] - ximag;
          xr[j] = xr[j] + xreal;
          xi[j] = xi[j] + ximag;
        }
        BaseT += Np;
      }
      Bp = Bp >> 1;
      Np = Np << 1;
      twiddleMul >>= 1;
    }
    // for (int P = 1; P < mBitsNeeded; P++)
    for (int P = (mBitsNeeded - 1); P > 0; P--) {
      int BaseT = 0;
      const int Npp = Np >> 1;
      for (int b = Bp; b > 0; b--) {
        const int BaseB = BaseT + Npp;
        for (int k_ = 0; k_ < Npp; k_++) {
          const int j = BaseT + k_;
          const int k = BaseB + k_;
          Tout xreal = xr[k];
          Tout ximag = xi[k];
          const int tblIdx = k_ * twiddleMul;
          const Tout twr = (Tout)getCos(tblIdx);
          const Tout twi = (Tout)getSin(tblIdx);
          const Tout tmp = xreal * twr - ximag * twi;
          ximag = xreal * twi + ximag * twr;
          xreal = tmp;
          xr[k] = xr[j] - xreal;
          xi[k] = xi[j] - ximag;
          xr[j] = xr[j] + xreal;
          xi[j] = xi[j] + ximag;
        }
        BaseT += Np;
      }
      Bp = Bp >> 1;
      Np = Np << 1;
      twiddleMul >>= 1;
    }

#ifdef USE_FAST_MATH
    for (int i = 0; i < mFftSize; i++) {
      _xr[i] = fxr[i];
      _xi[i] = fxi[i];
    }
#endif
    return true;
  }
private:
  // Allows the input arrays to be moved bit-reversed into pWorkingBuf
  int mReverseBitsLut[E_FFTSIZE_MAX + 4]; // = NULL;

#ifndef USE_FAST_MATH

  // Lookup table for sine values
  double mSinLut[(E_FFTSIZE_MAX / 2) + 4]; // = NULL;

#else

  // Lookup table for sine values
  float mSinLut[(E_FFTSIZE_MAX / 2) + 4]; // = NULL;

  fftc_t fxr[E_FFTSIZE_MAX + 4];
  fftc_t fxi[E_FFTSIZE_MAX + 4];
#endif

  // Stores the current size of the FFT
  int mFftSize;

  // Stores the number of bits needed for fftSize
  int mBitsNeeded;
};

#endif // #define FFTCLASS_HPP
