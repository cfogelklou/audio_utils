#ifndef HELPERMATH_HPP
#define HELPERMATH_HPP

/******************************************************************************
  Copyright 2014 Chris Fogelklou, Applaud Apps (Applicaudia)

  This source code may under NO circumstances be distributed without the 
  express written permission of the author.

  @author: chris.fogelklou@gmail.com
*******************************************************************************/ 

#include <math.h>
#include "audiolib_types.h"
#include "audutils_defines.h"

#ifndef MREQ_PI
#define MREQ_PI 3.14159265358979323846264338327
#endif

#ifndef MREQ_TWO_PI
#define MREQ_TWO_PI (2 * MREQ_PI)
#endif

#define MREQ_SQRT2 1.41421356237309504880

#ifdef __cplusplus

template <typename T>
bool HP_getPower(T pAdRealIn[], T pAdImagIn[], T adPow[], int nSize) {
  int i;
  ASSERT((pAdRealIn != NULL) && (adPow != NULL));

  if (NULL == pAdImagIn) {
    for (i = 0; i < nSize; i++) {
      adPow[i] = pAdRealIn[i] * pAdRealIn[i];
    }
  } else {
    for (i = 0; i < nSize; i++) {
      adPow[i] = pAdRealIn[i] * pAdRealIn[i] + pAdImagIn[i] * pAdImagIn[i];
    }
  }
  return true;
}

template <typename T>
bool HP_getAmplitude(T pAdRealIn[], T pAdImagIn[], T adPow[], int nSize) {
  int i;
  ASSERT((pAdRealIn != NULL) && (adPow != NULL));

  if (NULL == pAdImagIn) {
    for (i = 0; i < nSize; i++) {
      adPow[i] = sqrt(pAdRealIn[i] * pAdRealIn[i]);
    }
  } else {
    for (i = 0; i < nSize; i++) {
      adPow[i] =
          sqrt(pAdRealIn[i] * pAdRealIn[i] + pAdImagIn[i] * pAdImagIn[i]);
    }
  }
  return true;
}

template <typename T> T HP_getFftPeakAmpFromPowers(T pAdPowers[], int fftSize) {
  T sumOfPowers = 0;
  for (int i = 0; i < (fftSize / 2); i++) {
    sumOfPowers += pAdPowers[i];
  }
  const T amp = 2 * sqrt(sumOfPowers) / fftSize;
  return amp;
}

template <typename T>
T HP_getFftPeakAmpFromAmplitudes(T pAdAmplitudes[], int fftSize) {
  T sumOfPowers = 0;
  for (int i = 0; i < (fftSize / 2); i++) {
    sumOfPowers += pAdAmplitudes[i] * pAdAmplitudes[i];
  }
  const T amp = 2 * sqrt(sumOfPowers) / fftSize;
  return amp;
}

template <typename T> T HP_getFftRmsFromPowers(T pAdPowers[], int fftSize) {
  return HP_getFftPeakAmpFromPowers(pAdPowers, fftSize) / MREQ_SQRT2;
}

template <typename T>
T HP_getFftRmsFromAmplitudes(T pAdAmplitudes[], int fftSize) {
  return HP_getFftPeakAmpFromAmplitudes(pAdAmplitudes, fftSize) / MREQ_SQRT2;
}

template <typename T>
int HP_findPeak(T adMagnitude[], int magnitudeLen, int nMin, int nMax) {
  int peakIdx = 0;
  if (magnitudeLen > 0) {
    T peakPwr = 0;
    nMin = MAX(0, nMin);
    nMax = MIN(magnitudeLen, nMax);
    for (int i = nMin; i < nMax; i++) {
      const T mag = adMagnitude[i];
      if (mag > peakPwr) {
        peakPwr = mag;
        peakIdx = i;
      }
    }
  }
  return peakIdx;
}


template <typename T>
T HP_abs( const T in ) {
  return (in < 0) ? -in : in;
}

template <typename T>
int HP_findAbsPeak(T adMagnitude[], int magnitudeLen, int nMin, int nMax) {
  int peakIdx = 0;
  if (magnitudeLen > 0) {
    T peakPwr = 0;
    nMin = MAX(0, nMin);
    nMax = MIN(magnitudeLen, nMax);
    for (int i = nMin; i < nMax; i++) {
      const T mag = HP_abs(adMagnitude[i]);
      if (mag > peakPwr) {
        peakPwr = mag;
        peakIdx = i;
      }
    }
  }
  return peakIdx;
}

template <typename T>
T HP_findAbsPeakVal(T adMagnitude[], int magnitudeLen, int nMin, int nMax) {
  return adMagnitude[HP_findAbsPeak(adMagnitude, magnitudeLen, nMin, nMax)];
}

#endif // __cplusplus

// The point is for this to be faster than typical int64_t*int64_t since
// it is known that the upper words are 0.  BUT we need to profile to find out
// if it really is.
static __inline int64_t mul3232(const int32_t sa, const int32_t sb) {
  const uint32_t x = ABS(sa);
  const uint32_t y = ABS(sb);
  const uint32_t a = (x >> 16) & 0xffff;
  const uint32_t b = x & 0xffff;
  const uint32_t c = (y >> 16) & 0xffff;
  const uint32_t d = y & 0xffff;
  const int64_t r = ((uint64_t)(d*b)) + (((uint64_t)(d*a)) << 16) + (((uint64_t)(c*b)) << 16) + (((uint64_t)(c*a)) << 32);
  const int signr = ((sa < 0) && (sb > 0)) || ((sa > 0) && (sb < 0));
  return (signr) ? -r : r;
}

#endif //#ifndef HELPERMATH_HPP
