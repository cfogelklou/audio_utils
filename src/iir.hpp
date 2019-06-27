#ifndef BIQUAD_H__
#define BIQUAD_H__

/******************************************************************************
  Copyright 2014 Chris Fogelklou, Applaud Apps (Applicaudia)

  This source code may under NO circumstances be distributed without the 
  express written permission of the author.

  @author: chris.fogelklou@gmail.com
*******************************************************************************/

#ifdef __cplusplus

#include <math.h>
#include <float.h>

#include "audiolib_types.h"
#include "audutils_debug.h"
#include "audutils_defines.h"
#include "pcm_q.h"

extern const float IIR_FLT_MIN;
extern const float IIR_FLT_MAX;

typedef struct tIIR_CoefDesignDbl
{
  // Denominator  - Coefficients to past outputs
  double  a2; 
  double  a1;
  double  a0;
  // Nominator    - Coefficients to input and past input samples
  double  b2;
  double  b1; 
  double  b0;    
} IIR_CoefDesignDbl;

typedef struct tIIR_CoefDesignFlt
{
  // Denominator  - Coefficients to past outputs
  float  a2; 
  float  a1;
  float  a0;
  // Nominator    - Coefficients to input and past input samples
  float  b2;
  float  b1; 
  float  b0;    
} IIR_CoefDesignFlt;

//----------------------------------------------------------------------------

// History for a single filter.
typedef struct tIIR_DF2HistoryDbl
{
  double w2;
  double w1;
} IIR_DF2HistoryDbl;

// History for a single filter.
typedef struct tIIR_DF2HistoryFlt
{
  float w2;
  float w1;
} IIR_DF2HistoryFlt;

typedef struct tIIR_DF1HistoryDbl
{
  double x2;
  double x1;
  double y2;
  double y1;
} IIR_DF1HistoryDbl;

typedef struct tIIR_DF1HistoryFlt
{
  float x2;
  float x1;
  float y2;
  float y1;
} IIR_DF1HistoryFlt;

// Type of filter.
typedef enum eIIR_Type {
  IIR_LPF,
  IIR_HPF,
  IIR_BPF,
  IIR_NOTCH,
  IIR_PEQ,
  IIR_LSH,
  IIR_HSH,
} IIR_Type;
 

// Complex type.
typedef struct tIIR_ComplexDbl {
  double re;
  double im;
} IIR_ComplexDbl;

typedef struct tIIR_ComplexFlt {
  float re;
  float im;
} IIR_ComplexFlt;

typedef struct tIIR_ComplexPolDbl {
  double ma;
  double ph;
} IIR_ComplexPolDbl;

typedef struct tIIR_ComplexPolFlt {
  float ma;
  float ph;
} IIR_ComplexPolFlt;


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
template <typename IIR_ComplexRectT, typename IIR_ComplexPolT> 
void bq_Rect2Pol(IIR_ComplexRectT *pRect, IIR_ComplexPolT *pPol)
{
  IIR_ComplexPolT cResult;
  cResult.ma = MA_Sqrt((pRect->re * pRect->re) + (pRect->im * pRect->im));
  cResult.ph = (pRect->re == 0.0) ? M_PI/2 : MA_ATan(MA_Divide(pRect->im, pRect->re));
  pPol->ma = cResult.ma;
  pPol->ph = cResult.ph;
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
template <typename IIR_ComplexPolT, typename IIR_ComplexRectT > 
void 	bq_Pol2Rect(IIR_ComplexPolT *pPol, IIR_ComplexRectT *pRect)
{
  IIR_ComplexRectT cResult;
  cResult.re = pPol->ma * MA_Cos(pPol->ph);
  cResult.im = pPol->ma * MA_Sin(pPol->ph);
  pRect->re = cResult.re;
  pRect->im = cResult.im;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
template <typename IIR_ComplexPolT> 
void     bq_PolDivide2Pol(IIR_ComplexPolT *pNum, IIR_ComplexPolT *pDenom, IIR_ComplexPolT *pResult)
{

  ASSERT(pResult);
  ASSERT(pNum);
  ASSERT(pDenom);

  // Complex divide numerator by denominator to get result.
  if (pDenom->ma == 0.0) {
    if (pNum->ma == 0.0){
      pResult->ma = 0.0;
    }
    else {
      pResult->ma = (pNum->ma > 0.0) ? IIR_FLT_MAX : IIR_FLT_MIN;
    }
  }
  else
  {
    pResult->ma = MA_Divide(pNum->ma, pDenom->ma);
  }
  pResult->ph = pNum->ph - pDenom->ph;
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
template <typename IIR_ComplexRectT, typename IIR_ComplexPolT> 
static void  bq_Divide2Pol(IIR_ComplexRectT *pNum, IIR_ComplexRectT *pDenom, IIR_ComplexPolT *pResult)
{
  if (pResult) {
    IIR_ComplexPolT cNum;
    IIR_ComplexPolT cDenom;
    IIR_ComplexPolT cResult;

    ASSERT(pNum);
    ASSERT(pDenom);

    // Convert to polar.
    bq_Rect2Pol(pNum, &cNum);
    bq_Rect2Pol(pDenom, &cDenom);	

    // Divide
    bq_PolDivide2Pol(&cNum, &cDenom, &cResult);

    *pResult = cResult;
  }
}




///////////////////////////////////////////////////////////////////////////////
// Converts bandwidth in octaves to Q.
///////////////////////////////////////////////////////////////////////////////
template <typename floatT> 
floatT IIR_OctavesToQ(const floatT fOctaves)
{
  floatT fNum, fDenom;

  if ((fOctaves < 0.1) || (fOctaves > 2.0)) {
    return 1.0;
  }

  fNum = (floatT)MA_Pow2(fOctaves / 2.0);     // sqrt(2^octaves)
  fDenom = (floatT)MA_Pow2(fOctaves) - 1.0;   // 2^octaves - 1
  return (floatT)MA_Divide(fNum, fDenom);                 // sqrt(2^octaves) / (2^octaves - 1)

}


///////////////////////////////////////////////////////////////////////////////
// Converts Q to bandwidth in octaves.
//
//     (    1            {   1   }^2       )
// 2*ln(  -----  + sqrt[ { ----- }   + 1 ] )
//     (  2 * Q          { 2 * Q }         )
// ------------------------------------------
//                    ln(2)
///////////////////////////////////////////////////////////////////////////////
template <typename floatT> 
floatT IIR_QToOctaves(floatT fQ)
{
  double fOneOverTwoQ;
  double fTemp;

  if ((fQ <= 0.0)) return 0.0;

  fOneOverTwoQ = MA_Divide(1.0 , (2.0 * fQ));

  fTemp =  MA_Sqrt( (fOneOverTwoQ * fOneOverTwoQ) + 1.0  );
  fTemp += fOneOverTwoQ;
  fTemp =  MA_Ln(  fTemp  );
  fTemp =  fTemp * M_2_OVER_LN2;

  return (floatT)fTemp;
}

///////////////////////////////////////////////////////////////////////////////
// IIR_InitBW
// Based on http://www.musicdsp.org/files/Audio-EQ-Cookbook.txt
// 
// @param type     : The type of filter
// @param dbGain   : The gain of the filter, in dB.
// @param freq     : Center frequency, in Hz.
// @param srate    : Sampling frequency, in Hz.
// @param bandwidth: Bandwidth in octaves.
///////////////////////////////////////////////////////////////////////////////
template <typename IIR_CoefDesignT, typename floatT>
bool IIR_InitBW(
  IIR_CoefDesignT * const pThis, 
  const IIR_Type type, 
  const floatT dbGain, 
  const floatT freq, 
  const floatT srate, 
  const floatT bandwidth,
  bool  scaleA0TO1 = true)
{
  //double A, omega, sn, cs, alpha, beta;
  double a0, a1, a2, b0, b1, b2;

  /* setup variables */
  const double A = (type >= IIR_PEQ) ? MA_Pow(10.0, dbGain/40.0) : MA_Pow(10.0, dbGain/20.0);
  const double omega = 2.0 * M_PI * freq / srate;
  const double sn = MA_Sin(omega);
  const double cs = MA_Cos(omega);
  const double alpha = sn * MA_Sinh(M_LN2/2.0 * bandwidth * omega /sn);
  const double beta = MA_Sqrt(A + A);

  ASSERT(pThis);

  switch (type) {
  case IIR_LPF:
    b0 = (1.0 - cs)/2.0;
    b1 = 1.0 - cs;
    b2 = (1.0 - cs)/2.0;
    a0 = 1.0 + alpha;
    a1 = -2.0 * cs;
    a2 = 1.0 - alpha;
    break;
  case IIR_HPF:
    b0 = (1 + cs) /2;
    b1 = -(1 + cs);
    b2 = (1 + cs) /2;
    a0 = 1 + alpha;
    a1 = -2 * cs;
    a2 = 1 - alpha;
    break;
  case IIR_BPF:
    b0 = alpha;
    b1 = 0;
    b2 = -alpha;
    a0 = 1 + alpha;
    a1 = -2 * cs;
    a2 = 1 - alpha;
    break;
  case IIR_NOTCH:
    b0 = 1;
    b1 = -2 * cs;
    b2 = 1;
    a0 = 1 + alpha;
    a1 = -2 * cs;
    a2 = 1 - alpha;
    break;
  case IIR_PEQ:
    b0 = 1 + (alpha * A);
    b1 = -2 * cs;
    b2 = 1 - (alpha * A);
    a0 = 1 + (alpha /A);
    a1 = -2 * cs;
    a2 = 1 - (alpha /A);
    break;
  case IIR_LSH:
    b0 = A * ((A + 1) - (A - 1) * cs + beta * sn);
    b1 = 2 * A * ((A - 1) - (A + 1) * cs);
    b2 = A * ((A + 1) - (A - 1) * cs - beta * sn);
    a0 = (A + 1) + (A - 1) * cs + beta * sn;
    a1 = -2 * ((A - 1) + (A + 1) * cs);
    a2 = (A + 1) + (A - 1) * cs - beta * sn;
    break;
  case IIR_HSH:
    b0 = A * ((A + 1) + (A - 1) * cs + beta * sn);
    b1 = -2 * A * ((A - 1) + (A + 1) * cs);
    b2 = A * ((A + 1) + (A - 1) * cs - beta * sn);
    a0 = (A + 1) - (A - 1) * cs + beta * sn;
    a1 = 2 * ((A - 1) - (A + 1) * cs);
    a2 = (A + 1) - (A - 1) * cs - beta * sn;
    break;
  default:
    return false;
  }

  if (scaleA0TO1) {
    // a0 becomes 1 by dividing all coeffs by a0.
    pThis->b0 = (floatT)(b0/a0);
    pThis->b1 = (floatT)(b1/a0);
    pThis->b2 = (floatT)(b2/a0);
    pThis->a0 = 1;
    pThis->a1 = (floatT)(a1/a0);
    pThis->a2 = (floatT)(a2/a0);
  }
  else {
    pThis->b0 = (floatT)b0;
    pThis->b1 = (floatT)b1;
    pThis->b2 = (floatT)b2;
    pThis->a0 = (floatT)a0;
    pThis->a1 = (floatT)a1;
    pThis->a2 = (floatT)a2;
  }
  return true;
}


///////////////////////////////////////////////////////////////////////////////
// IIR_InitQ
// 
// @param type    : The type of filter
// @param dbGain  : The gain of the filter, in dB.
// @param freq    : Center frequency, in Hz.
// @param srate   : Sampling frequency, in Hz.
// @param Q       : Q factor
///////////////////////////////////////////////////////////////////////////////
template <typename IIR_CoefDesignT, typename floatT>
bool IIR_InitQ(
  IIR_CoefDesignT * const pThis, 
  const IIR_Type type, 
  const floatT dbGain, 
  const floatT freq, 
  const floatT srate, 
  const floatT Q,
  bool  scaleA0TO1 = true)
{
  floatT fOctaves = IIR_QToOctaves(Q);
  return IIR_InitBW(pThis, type, dbGain, freq, srate, fOctaves);
}

///////////////////////////////////////////////////////////////////////////////
// IIR_GetMagPhaseResponse
// 
//     Arguments: fOmega = 2*pi*f/fs | f & fs are in Hz.
// 
// 
//     Y(z)       b2z^(-2) + b1z^(-1)  + b0  |
//     ---- is    ---------------------      |a0 == 1
//     X(z)       a2z^(-2) + a1z^(-1)  + a0  |
///////////////////////////////////////////////////////////////////////////////
template <typename IIR_CoefDesignT, typename floatT, typename IIR_ComplexRectT>
void IIR_GetMagPhaseResponse(IIR_CoefDesignT *pThis, floatT fOmega, IIR_ComplexRectT *pResult)
{
  IIR_ComplexRectT cA;
  IIR_ComplexRectT cB;

  // Pre-calculate some values that are used quite a bit.
  double cosW = MA_Cos(fOmega);
  double cos2W = MA_Cos(2*fOmega);
  double sinW = MA_Sin(fOmega);
  double sin2W = MA_Sin(2*fOmega);

  if (!(pResult)) return;

  // Calculate the numerator real and imaginary parts
  cB.re = pThis->b0 + (pThis->b1 * cosW) + (pThis->b2 * cos2W);
  cB.im = -(pThis->b1 * sinW) - (pThis->b2 * sin2W);

  // Calculate the denominator real and imaginary parts.
  cA.re = 1.0 + (pThis->a1 * cosW) + (pThis->a2 * cos2W);
  cA.im = -(pThis->a1 * sinW) - (pThis->a2 * sin2W);

  // Divide the B's by the A's.
  bq_Divide2Pol(&cB, &cA, pResult);

}


///////////////////////////////////////////////////////////////////////////////
// Process a single sample of a biquad. (Direct form 2)
//                   1.0      w(n-0)    b0
//x(n) ------+--------*------------------*---------- +------- y(n)
//           |\                 |                   /|
//           | \              delay                / |
//           |  \               |                 /  |
//           |   \   -a1        | w(n-1)   b1    /   |
//           |    <---*---------+----------*---->    |
//            \                 |                   /
//             \              delay                /
//              \               |                 /
//               \   -a2        | w(n-2)   b2    /
//                <---*---------+----------*---->
///////////////////////////////////////////////////////////////////////////////
template <typename IIR_CoefDesignT, typename IIR_DF2HistoryT, typename FloatT>
FloatT IIR_DF2ProcessSingleSampleFlt(
  IIR_CoefDesignT * const pH, 
  IIR_DF2HistoryT * const pHist, 
  const FloatT x)
{
  // compute result
  const FloatT w1 = pHist->w1;
  const FloatT w0 = (FloatT)(x - pHist->w2 * pH->a2 - w1 * pH->a1);
  const FloatT y  = (FloatT)(pHist->w2 * pH->b2 + w1 * pH->b1 + w0 * pH->b0);
  // shift p0 to p1 and p1 to p2
  pHist->w2 = w1;
  pHist->w1 = w0;
  return y;
}

///////////////////////////////////////////////////////////////////////////////
template <typename IIR_CoefDesignT, typename IIR_DF2HistoryT, typename FloatT>
double IIR_DF2ProcessSingleSampleDbl(
  IIR_CoefDesignT * const pH, 
  IIR_DF2HistoryT * const pHist, 
  const FloatT x)
{
  // compute result
  const double w1 = pHist->w1;
  const double w0 = x - pHist->w2 * pH->a2 - w1 * pH->a1;
  const double y  = (double)(pHist->w2 * pH->b2 + w1 * pH->b1 + w0 * pH->b0);
  // shift p0 to p1 and p1 to p2
  pHist->w2 = w1;
  pHist->w1 = w0;
  return y;
}

///////////////////////////////////////////////////////////////////////////////
// Process multiple samples in a biquad.
///////////////////////////////////////////////////////////////////////////////
template <typename IIR_CoefDesignT, typename IIR_DF2HistoryT, typename FloatT>
void IIR_DF2ProcessMultipleSamples(
  IIR_CoefDesignT * const  pH, 
  IIR_DF2HistoryT * const  pHist, 
  const FloatT *pXAry, 
  FloatT *pYAry, 
  int nSamples)
{
  // compute result
  while (nSamples > 0) {
    double w0     = *pXAry++ - pHist->w1 * pH->a1 - pHist->w2 * pH->a2;
    *pYAry++ = (FloatT)(w0 * pH->b0 + pHist->w1 * pH->b1 + pHist->w2 * pH->b2);
    // shift p0 to p1 and p1 to p2
    pHist->w2 = pHist->w1;
    pHist->w1 = (FloatT)w0;
    nSamples--;
  }
}


///////////////////////////////////////////////////////////////////////////////
// Process multiple samples in a biquad.
///////////////////////////////////////////////////////////////////////////////
template <typename IIR_CoefDesignT, typename IIR_DF2HistoryT, typename FloatT>
void IIR_DF2ProcessMultipleSamplesAndAddToOutput(
  IIR_CoefDesignT * const  pH, 
  IIR_DF2HistoryT * const  pHist, 
  const FloatT *pXAry, 
  FloatT *pYAry, 
  int nSamples)
{
  // compute result
  while (nSamples > 0) {
    double w0     = *pXAry++ - pHist->w1 * pH->a1 - pHist->w2 * pH->a2;
    *pYAry += (FloatT)(w0 * pH->b0 + pHist->w1 * pH->b1 + pHist->w2 * pH->b2);
    // shift p0 to p1 and p1 to p2
    pHist->w2 = pHist->w1;
    pHist->w1 = (FloatT)w0;
    ++pYAry;
    --nSamples;
  }
}


///////////////////////////////////////////////////////////////////////////////
// Process single sample using DF1.
///////////////////////////////////////////////////////////////////////////////
template <typename IIR_CoefDesignT, typename IIR_DF1HistoryT, typename FloatT>
FloatT IIR_DF1ProcessSingleSample(
  IIR_CoefDesignT * const pH, 
  IIR_DF1HistoryT * const pHist, 
  const FloatT x)
{

  double v = pH->b0*x + pH->b1*pHist->x1 + pH->b2*pHist->x2;
  double y = v - pH->a1*pHist->y1 - pH->a2*pHist->y2;
  pHist->y2 = pHist->y1;
  pHist->y1 = y;
  pHist->x2 = pHist->x1;
  pHist->x1 = x;

  return (FloatT)y;
}

#endif

#endif //__TCBQPROC_H__

