#ifndef FFT_WINDOW_H__
#define FFT_WINDOW_H__

/*=============================================================================
//	
//  This software has been released under the terms of the GNU Public
//  license. See http://www.gnu.org/copyleft/gpl.html for details.
//
//  Copyright 2001 Anders Johansson ajh@atri.curtin.edu.au
//
//  Modified into a class by C. Fogelklou
//
//=============================================================================
*/

#include <stdint.h>
#include <math.h>

#ifndef M_PI
#define PI        (3.141592653589793238462643383279502884197169399375105820974944592307816406286)
#define M_PI PI
#endif

// Window types that can be applied before the FFT.
typedef enum teFFTWin {
  FFTW_NONE,
  FFTW_HANNING,
  FFTW_HANN,
  FFTW_HAMMING,
  FFTW_BLACKMAN
} FFTWinT;

#ifdef __cplusplus

template<typename FloatT>
void fftw_boxcar(int n, FloatT w[])
{
  int i;
  /* Calculate window coefficients */
  for (i = 0; i < n; i++)
    w[i] = (FloatT)(1.0);
}


/*
* Triang a.k.a Bartlett
*
*               |    (N-1)|
*           2 * |k - -----|
*               |      2  |
* w = 1.0 - ---------------
*                    N+1
* n window length
* w buffer for the window parameters
*/
template<typename FloatT>
void fftw_triang(int n, FloatT w[])
{
  double k1 = (double)(n & 1);
  double k2 = 1 / ((double)n + k1);
  int    end = (n + 1) >> 1;
  int	   i;

  /* Calculate window coefficients */
  for (i = 0; i < end; i++)
    w[i] = w[n - i - 1] = (FloatT)((2.0*((double)(i + 1)) - (1.0 - k1))*k2);
}


/*
* Hanning
*                   2*pi*k
* w = 0.5 - 0.5*cos(------), where 0 < k <= N
*                    N+1
* n window length
* w buffer for the window parameters
*/
template<typename FloatT>
void fftw_hanning(int n, FloatT w[])
{
  int	   i;
  double k = (double)(2 * M_PI / ((double)(n + 1))); /* 2*pi/(N+1) */

  /* Calculate window coefficients */
  for (i = 0; i < n; i++)
    w[i] = (FloatT)(0.5*(1.0 - cos(k*(double)(i + 1))));
}

/*
* Hann (same as matlab function)
*                   2*pi*k
* w = 0.5*(1 - cos(-------)), where 0 <= k <= N
*                    N
* L window length
* w buffer for the window parameters
*/
template<typename FloatT>
void fftw_hann(int L, FloatT w[])
{
  int    n;
  const int N = (L - 1);
  const double coeff = 2 * M_PI / ((double)(N)); /* 2*pi/(N) */

  /* Calculate window coefficients */
  for (n = 0; n <= N; n++)
    w[n] = (FloatT)(0.5*(1.0 - cos(coeff*(double)(n))));
}

/*
* Hamming
*                        2*pi*k
* w(k) = 0.54 - 0.46*cos(------), where 0 <= k < N
*                         N-1
*
* n window length
* w buffer for the window parameters
*/
template<typename FloatT>
void fftw_hamming(int n, FloatT w[])
{
  int      i;
  double k = 2 * M_PI / ((FloatT)(n - 1)); /* 2*pi/(N-1) */

  /* Calculate window coefficients */
  for (i = 0; i < n; i++)
    w[i] = (FloatT)(0.54 - 0.46*cos(k*i));
}

/*
* Blackman
*                       2*pi*k             4*pi*k
* w(k) = 0.42 - 0.5*cos(------) + 0.08*cos(------), where 0 <= k < N
*                        N-1                 N-1
*
* n window length
* w buffer for the window parameters
*/
template<typename FloatT>
void fftw_blackman(int n, FloatT w[])
{
  int      i;
  double k1 = 2 * M_PI / ((double)(n - 1)); /* 2*pi/(N-1) */
  double k2 = 2 * k1; /* 4*pi/(N-1) */

  /* Calculate window coefficients */
  for (i = 0; i < n; i++)
    w[i] = (FloatT)(0.42 - 0.50*cos(k1*(FloatT)i) + 0.08*cos(k2*(double)i));
}

/*
* Flattop
*                                        2*pi*k                     4*pi*k
* w(k) = 0.2810638602 - 0.5208971735*cos(------) + 0.1980389663*cos(------), where 0 <= k < N
*                                          N-1                        N-1
*
* n window length
* w buffer for the window parameters
*/
template<typename FloatT>
void fftw_flattop(int n, FloatT w[])
{
  int      i;
  double k1 = 2 * M_PI / ((FloatT)(n - 1)); /* 2*pi/(N-1) */
  double k2 = 2 * k1;                   /* 4*pi/(N-1) */

  /* Calculate window coefficients */
  for (i = 0; i < n; i++)
    w[i] = (FloatT)(0.2810638602 - 0.5208971735*cos(k1*(double)i) + 0.1980389663*cos(k2*(double)i));
}

/* Computes the 0th order modified Bessel function of the first kind.
* (Needed to compute Kaiser window)
*
* y = sum( (x/(2*n))^2 )
*      n
*/
#define BIZ_EPSILON (1E-21) /* Max error acceptable */

template<typename FloatT>
FloatT besselizero(FloatT x)
{
  double temp;
  double sum = 1.0;
  double u = 1.0;
  double halfx = x / 2.0;
  int      n = 1;

  do {
    temp = halfx / (FloatT)n;
    u *= temp * temp;
    sum += u;
    n++;
  } while (u >= BIZ_EPSILON * sum);
  return((FloatT)sum);
}

/*
* Kaiser
*
* n window length
* w buffer for the window parameters
* b beta parameter of Kaiser window, Beta >= 1
*
* Beta trades the rejection of the low pass filter against the
* transition width from passband to stop band.  Larger Beta means a
* slower transition and greater stop band rejection.  See Rabiner and
* Gold (Theory and Application of DSP) under Kaiser windows for more
* about Beta.  The following table from Rabiner and Gold gives some
* feel for the effect of Beta:
*
* All ripples in dB, width of transition band = D*N where N = window
* length
*
* BETA    D       PB RIP   SB RIP
* 2.120   1.50  +-0.27      -30
* 3.384   2.23    0.0864    -40
* 4.538   2.93    0.0274    -50
* 5.658   3.62    0.00868   -60
* 6.764   4.32    0.00275   -70
* 7.865   5.0     0.000868  -80
* 8.960   5.7     0.000275  -90
* 10.056  6.4     0.000087  -100
*/
//			void kaiser(int n, FloatT w[], FloatT b)
//			{
//			  FloatT tmp;
//			  FloatT k1  = 1.0/besselizero(b);
//			  int	   k2  = 1 - (n & 1);
//			  int      end = (n + 1) >> 1;
//			  int      i; 
//			  
//			  /* Calculate window coefficients */
//			  for (i=0 ; i<end ; i++){
//			    tmp = (FloatT)(2*i + k2) / ((FloatT)n - 1.0);
//			    w[end-(1&(!k2))+i] = w[end-1-i] = k1 * besselizero(b*Math.sqrt(1.0 - tmp*tmp));
//			  }
//			}


// ////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////////////////////
template<typename FloatT>
void FFTW_SetFftWindow(const FFTWinT winType, FloatT arr[], const int size )
{
  switch (winType) {
  case FFTW_HANNING:
    fftw_hanning(size, arr);
    break;
  case FFTW_HANN:
    fftw_hann(size, arr);
    break;
  case FFTW_HAMMING:
    fftw_hamming(size, arr);
    break;
  case FFTW_BLACKMAN:
    fftw_blackman(size, arr);
    break;
  case FFTW_NONE:
  default:
    fftw_boxcar(size, arr);
    break;
  }
}

#endif // __cplusplus

#endif
