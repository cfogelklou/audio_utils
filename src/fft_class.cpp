#include "fft_class.hpp"

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
}