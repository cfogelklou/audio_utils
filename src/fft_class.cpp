#include "fft_class.hpp"
#include "audutils_defines.h"

#ifndef M_PI
#define M_PI (3.141592653589793238462643383279502884197169399375105820974944592307816406286)
#endif

#ifndef MREQ_TWO_PI
#define MREQ_TWO_PI (2*M_PI)
#endif

namespace Fft {

  bool getPower( double pAdRealIn[], double pAdImagIn[], double adPow[],  int nSize )
  {
    int i;
    ASSERT ((pAdRealIn != NULL) && (adPow != NULL));

    if (NULL == pAdImagIn) {
      for (i = 0; i < nSize; i++)  {
        adPow[i] = pAdRealIn[i] * pAdRealIn[i]; 
      }
    }
    else {
      for (i = 0; i < nSize; i++) {
        adPow[i] = pAdRealIn[i] * pAdRealIn[i] + pAdImagIn[i] * pAdImagIn[i];
      }
    }
    return true;
  }

  bool getAmplitude(double pAdRealIn[], double pAdImagIn[], double adPow[], int nSize)
  {
    int i;
    ASSERT((pAdRealIn != NULL) && (adPow != NULL));

    if (NULL == pAdImagIn) {
      for (i = 0; i < nSize; i++)  {
        adPow[i] = sqrt( pAdRealIn[i] * pAdRealIn[i] );
      }
    }
    else {
      for (i = 0; i < nSize; i++) {
        adPow[i] = sqrt( pAdRealIn[i] * pAdRealIn[i] + pAdImagIn[i] * pAdImagIn[i] );
      }
    }
    return true;
  }

  double getFftPeakAmpFromPowers(double pAdPowers[], int fftSize) {
    double sumOfPowers = 0;
    for (int i = 0; i < (fftSize / 2); i++) {
      sumOfPowers += pAdPowers[i];
    }
    const double amp = 2*sqrt(sumOfPowers) / fftSize;
    return amp;
  }

  double getFftPeakAmpFromAmplitudes(double pAdAmplitudes[], int fftSize) {
    double sumOfPowers = 0;
    for (int i = 0; i < (fftSize / 2); i++) {
      sumOfPowers += pAdAmplitudes[i] * pAdAmplitudes[i];
    }
    const double amp = 2*sqrt(sumOfPowers) / fftSize;
    return amp;
  }

  double getFftRmsFromPowers(double pAdPowers[], int fftSize) {
    return getFftPeakAmpFromPowers( pAdPowers, fftSize )/SQRT2 ;
  }

  double getFftRmsFromAmplitudes(double pAdAmplitudes[], int fftSize) {
    return getFftPeakAmpFromAmplitudes( pAdAmplitudes, fftSize )/SQRT2;
  }

  int findPeak( double adMagnitude[], int magnitudeLen, int nMin, int nMax )
  {
    int peakIdx = 0;
    if (magnitudeLen > 0) {
      double peakPwr = 0;
      nMin = MAX(0, nMin);
      nMax = MIN(magnitudeLen, nMax);
      for (int i = nMin; i < nMax; i++) {
        const double mag = adMagnitude[i];
        if (mag > peakPwr) {
          peakPwr = mag;
          peakIdx = i;
        }
      }
    }
    return peakIdx;
  } 

  //---------------------------------------------------------------------------
  int  bitsNeeded(
    const unsigned int nPowerOfTwo)
  {

    unsigned int i = 0;
    bool success = false;
    unsigned int bitsNeeded = 0;

    if (nPowerOfTwo < 2) {
      bitsNeeded = 0;
    }

    while ((i < 32) && (success == false)) {
      if (nPowerOfTwo & (1 << i)) {
        bitsNeeded = i;
        success = true;
      }
      ++i;
    }

    return bitsNeeded;
  }
}

//---------------------------------------------------------------------------
// Initializes or re-initializes the FFT.
FftClass::FftClass(const int nFftSize)
  : mFftSize(nFftSize)
{
  ASSERT(nFftSize <= E_FFTSIZE_MAX);
  if (FFT_IsPowerOfTwo(nFftSize)) {
    mFftSize = nFftSize;
    setFftSize(nFftSize);
  }
}

//---------------------------------------------------------------------------
FftClass::~FftClass() { 
  TRACE_VERBOSE(("Destructor: FftClass();"));
}

//---------------------------------------------------------------------------
int FftClass::getSize() { 
  return mFftSize; 
}

//---------------------------------------------------------------------------
void FftClass::setLutSize(const int fftSize) {
  const double fDenom = -MREQ_TWO_PI / (double)mFftSize;

  for (int i = 0; i < (mFftSize / 2); i++) {
    mSinLut[i] = (fftc_t)sin(((double)i) * fDenom);
  }

  // Create a bit-reversed table
  uint_t bitsNeeded;
  ASSERT_FN(fft_NumberOfBitsNeeded(fftSize, &bitsNeeded));
  for (int i = 0; i < mFftSize; i++) {
    mReverseBitsLut[i] = fft_ReverseBits(i, bitsNeeded);
  }
}

// ---------------------------------------------------------------------------
void FftClass::setFftSize(int nFftSize) {
  ASSERT(nFftSize <= E_FFTSIZE_MAX);
  if (FFT_IsPowerOfTwo(nFftSize)) {
    uint_t n;
    fft_NumberOfBitsNeeded(nFftSize, &n);
    mBitsNeeded = n;
    if (0 != mBitsNeeded) {
      mFftSize = nFftSize;

      // Get / create a lookup table for this FFT size.
      setLutSize(nFftSize);
    }
  }
}

//---------------------------------------------------------------------------
bool FftClass::doFft(fftin_t adRealIn[], fftin_t adImagIn[], double adRealOut[],
  double adImagOut[], int nSize)

{
  return doFftOrIfft(false, adRealIn, adImagIn, adRealOut, adImagOut, nSize);
}

//---------------------------------------------------------------------------
bool FftClass::doIfft(fftin_t adRealIn[], fftin_t adImagIn[], double adRealOut[],
  double adImagOut[], int nSize) {
  return doFftOrIfft(true, adRealIn, adImagIn, adRealOut, adImagOut, nSize);
}

  //---------------------------------------------------------------------------
  bool FftClass::doFftOrIfft(bool bInverseTransform, fftin_t pAdRealIn[],
    fftin_t pAdImagIn[], double _xr[], double _xi[],
    int nNumSamples) {
    int i;

    ASSERT((nNumSamples == mFftSize) && (pAdRealIn != NULL) && (_xr != NULL) &&
      (_xi != NULL));
    // Reverse ordering of samples so FFT can be done in place.
#ifdef USE_FAST_MATH
    fftc_t* xr = fxr;
    fftc_t* xi = fxi;
#else
    fftc_t* xr = _xr;
    fftc_t* xi = _xi;
#endif
    if (pAdImagIn == NULL) {
      for (i = 0; i < nNumSamples; i++) {
        xr[mReverseBitsLut[i]] = (fftc_t)pAdRealIn[i];
      }
      memset(xi, 0, sizeof(fftc_t) * nNumSamples);

    }
    else {
      for (i = 0; i < nNumSamples; i++) {
        const int rev = mReverseBitsLut[i];
        xr[rev] = (fftc_t)pAdRealIn[i];
        xi[rev] = (fftc_t)pAdImagIn[i];
      }
    }

    // int p = mBitsNeeded;
    int Bp = 1 << (mBitsNeeded - 1);
    int Np = 2;
    int twiddleMul = (mFftSize >> 1);

    for (int P = 0; P < mBitsNeeded; P++) {
      int b;
      int BaseT = 0;
      int Npp = Np >> 1;
      for (b = 0; b < Bp; b++) {
        int k_;
        int BaseB = BaseT + Npp;
        for (k_ = 0; k_ < Npp; k_++) {
          fftc_t tmp;
          const int j = BaseT + k_;
          const int k = BaseB + k_;
          fftc_t xreal = xr[k];
          fftc_t ximag = xi[k];

          const int tblIdx = k_ * twiddleMul;
          const fftc_t twr = getCos(tblIdx);
          const fftc_t twi =
            (bInverseTransform) ? -getSin(tblIdx) : getSin(tblIdx);

          // Calculate real part
          tmp = xreal * twr - ximag * twi;

          // Calculate imaginary part
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

    // Normalize
    if (bInverseTransform) {
      double dDenom = (double)1.0 / (double)nNumSamples;
      for (i = 0; i < nNumSamples; i++) {
        _xr[i] = xr[i] * dDenom;
        _xi[i] = xi[i] * dDenom;
      }
    }
#ifdef USE_FAST_MATH
    else {
      for (i = 0; i < nNumSamples; i++) {
        _xr[i] = fxr[i];
        _xi[i] = fxi[i];
      }
    }
#endif

    return true;
  }

  // ---------------------------------------------------------------------------
  void FftClass::doRealInverseFft(fftin_t pAdRealIn[], double xr[], double xi[]) {
    doRealForwardFft(pAdRealIn, xr, xi);
    const double dDenom = (double)1.0 / (double)mFftSize;

    for (int i = 0; i < mFftSize; i++) {
      xr[i] *= dDenom;
      xi[i] *= dDenom;
    }
  }

