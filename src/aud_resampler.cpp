/******************************************************************************
Copyright 2014 Chris Fogelklou, Applaud Apps (Applicaudia)

This source code may under NO circumstances be distributed without the 
express written permission of the author.

@author: chris.fogelklou@gmail.com
*******************************************************************************/

#include <math.h>
#include <string.h>
#include "AudResampler.h"
#include "DaTunerDebug.h"
#include "DaTunerDefs.h"
#include "Iir.hpp"

#define float40_t double
#define AR_INIT_COMPLETE 0x7E5A5555
typedef void (*pAudFnResampProcess)( struct tAudResampler * const pThis, AudResampleIO *pIoAry, int numIos );

extern "C" {

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
static void AudUpSampleProcessSamples( AudResampler * const pThis, AudResampleIO *pIoAry, int numIos )
{
  AudResamplerChannel *pChannel = pThis->pChannels;
  ASSERT( numIos == (int)pThis->numChannels );
  do {
    int nSamples = pIoAry->numXSamples;
    ar_float *pyOut = pIoAry->bufs.flt.pyOut;
    const ar_float *pSamples = pIoAry->bufs.flt.pxIn;
    while( nSamples ) {
      int_t xIdx;
      int_t readIdx = pChannel->writeIdx;

      pChannel->arr.pXFltHistory[ pChannel->writeIdx++ ] = *pSamples;
      pSamples += pIoAry->xHop;
      if(pChannel->writeIdx >= pChannel->xHistorySize) pChannel->writeIdx = 0;

      // Use this one x sample with every set of taps.
      for( xIdx = 0; xIdx < pThis->ratio; xIdx++ )  {
        int_t tapIdx;
        float40_t sum = 0;
        int_t histIdx = readIdx;

        for( tapIdx = xIdx; tapIdx < pThis->numTaps; tapIdx += pThis->ratio )
        {
          sum += (float40_t)pThis->coefAry.pCoeffs[ tapIdx ] * (float40_t)pChannel->arr.pXFltHistory[ histIdx-- ];
          if(histIdx < 0) histIdx = pChannel->xHistorySize-1;
        }

        sum *= (float40_t)pThis->ratio;
        pyOut[ xIdx * pIoAry->yHop ] = (ar_float)sum;
      }

      pyOut += pThis->ratio * pIoAry->yHop;
      nSamples--;
    }
    pChannel++;
    pIoAry++;
  }
  while (--numIos > 0);
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
static void AudDnSampleProcessSamplesFlt( AudResampler * const pThis, AudResampleIO *pIoAry, int numIos )
{
  AudResamplerChannel *pChannel = pThis->pChannels;
  ASSERT( numIos == (int)pThis->numChannels );
  do {
    int nSamples = pIoAry->numXSamples;
    ar_float *pyOut = pIoAry->bufs.flt.pyOut;
    ar_float *pxIn  = pIoAry->bufs.flt.pxIn;
    while( nSamples ) {
      int_t readIdx;
      int_t xIdx = 0;

      // Load up the xHistory array to prepare for our next output sample.
      while ( xIdx < pThis->ratio ) {
        pChannel->arr.pXFltHistory[pChannel->writeIdx++] = pxIn[xIdx*pIoAry->xHop];
        if(pChannel->writeIdx >= pChannel->xHistorySize) pChannel->writeIdx = 0;
        xIdx++;
      }

      readIdx = pChannel->writeIdx - 1;
      if (readIdx < 0) {readIdx = pChannel->xHistorySize - 1;}

      // This has nothing to do with the above if.
      {
        int_t tapIdx = 0;

        // Initialize
        float40_t sum = 0;

        // All taps.  Candidate for loop unrolling?
        while (tapIdx < pThis->numTaps)
        {
          sum += (pThis->coefAry.pCoeffs[ tapIdx ] * pChannel->arr.pXFltHistory[ readIdx-- ] );
          if (readIdx < 0) readIdx = pChannel->xHistorySize - 1;
          tapIdx++;
        }

        *pyOut = (ar_float)sum;
        pyOut += pIoAry->yHop;
        pxIn += pThis->ratio*pIoAry->xHop;
        nSamples -= pThis->ratio;
      }
    }
    pChannel++;
    pIoAry++;
  }
  while (--numIos > 0);
}

#define DO_OPT
#ifndef DO_OPT
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
static void AudDnSampleProcessSamplesI16( AudResampler * const pThis, AudResampleIO *pIoAry, int numIos )
{
  AudResamplerChannel *pChannel = pThis->pChannels;
  ASSERT( numIos == (int)pThis->numChannels );
  do {
    int nSamples = pIoAry->numXSamples;
    ar_float *pyOut = pIoAry->bufs.i16.pyOut;
    int16_t *pxIn  = pIoAry->bufs.i16.pxIn;
    while( nSamples > 0 ) {
      int_t readIdx;
      int_t xIdx = 0;

      // Load up the xHistory array to prepare for our next output sample.
      while ( xIdx < pThis->ratio ) {
        pChannel->arr.pXI32History[pChannel->writeIdx++] = pxIn[xIdx*pIoAry->xHop];
        if(pChannel->writeIdx >= pChannel->xHistorySize) pChannel->writeIdx = 0;
        xIdx++;
      }

      readIdx = pChannel->writeIdx - 1;
      if (readIdx < 0) {readIdx = pChannel->xHistorySize - 1;}

      // This has nothing to do with the above if.
      {
        int_t tapIdx = 0;

        // Initialize
        ar_acc_t sum = 0;

        // All taps.  Candidate for loop unrolling?
        while (tapIdx < pThis->numTaps) {
          sum += ((ar_acc_t)pThis->coefAry.i32.pCoeffs32[ tapIdx ] * pChannel->arr.pXI32History[ readIdx-- ] );
          if (readIdx < 0) readIdx = pChannel->xHistorySize - 1;
          tapIdx++;
        }

       *pyOut = (ar_float)(sum * pThis->coefAry.i32.postFilterScale * (1.0/32768.0));
        pyOut += pIoAry->yHop;
        pxIn += pThis->ratio*pIoAry->xHop;
        nSamples -= pThis->ratio;
      }
    }
    pChannel++;
    pIoAry++;
  }
  while (--numIos > 0);
}

#else

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
static void AudDnSampleProcessSamplesI16( AudResampler * const pThis, AudResampleIO *pIoAry, int numIos )
{
  const int numTaps = pThis->numTaps;
  AudResamplerChannel *pChannel = pThis->pChannels;
  ASSERT( numIos == (int)pThis->numChannels );
  do {
    int nSamples = pIoAry->numXSamples;
    ar_float *pyOut = pIoAry->bufs.i16.pyOut;
    int16_t *pxIn  = pIoAry->bufs.i16.pxIn;
    while( nSamples > 0 ) {
      int_t readIdx;

      // Load up the xHistory array to prepare for our next output sample.
      {
        int_t xIdx = 0;
        const int lastItem = pThis->ratio * pIoAry->xHop;
        for (xIdx = 0; xIdx < lastItem; xIdx += pIoAry->xHop) {
          pChannel->arr.pXI32History[pChannel->writeIdx++] = pxIn[xIdx];
          if(pChannel->writeIdx >= pChannel->xHistorySize) {
            memcpy(
              &pChannel->arr.pXI32History[0], 
              &pChannel->arr.pXI32History[pChannel->xHistorySize - numTaps],
              sizeof(pChannel->arr.pXI32History[0]) * numTaps);
            pChannel->writeIdx = numTaps;  // Leave at one past the last item.
          }
        }
      }

      // writeIdx points to the next location to be written.
      // readIdx starts at last written location.
      readIdx = pChannel->writeIdx - 1;
      ASSERT(readIdx >= (numTaps-1));// {readIdx = pChannel->xHistorySize - 1;}

      // This has nothing to do with the above if.
      {
        int_t tapIdx = 0;

        // Initialize
        ar_acc_t sum = 0;

        // All taps.  Candidate for loop unrolling?
        while (tapIdx < numTaps) {
          sum += ((ar_acc_t)pThis->coefAry.i32.pCoeffs32[ tapIdx++ ] * pChannel->arr.pXI32History[ readIdx-- ] );
        }

       *pyOut = (ar_float)(sum * pThis->coefAry.i32.postFilterScale * (1.0/32768.0));
        pyOut += pIoAry->yHop;
        pxIn += pThis->ratio*pIoAry->xHop;
        nSamples -= pThis->ratio;
      }
    }
    ASSERT( nSamples == 0 );
    pChannel++;
    pIoAry++;
  }
  while (--numIos > 0);
}
#endif

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
static ar_coef audresample_getPeak(const ar_coef * const pArr, const int numSamps)
{
  int c;
  ar_coef absPeak = 0;
  for (c = 0; c < numSamps; c++) {
    const ar_coef samp = (ABS(pArr[c]));
    if (samp > absPeak) {
      absPeak = samp;
    }
  }
  return absPeak;
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
bool_t AudResampleInit( 
  AudResampler                 *pThis, 
  bool_t                        isUpsampler, 
  const AudResampleConfig      *pCfg, 
  AudResamplerChannel          *pChannelAry, 
  uint_t                        numChannels,
  AudResamplerFixedPointCfg    *pOptionalFixedPointCfg)
{
  bool_t rval = FALSE;
  ASSERT(( (pThis != NULL) && (pChannelAry != NULL) && (pCfg != NULL) && (numChannels > 0)));
  if ( (pThis != NULL) && (pChannelAry != NULL) && (pCfg != NULL) && (numChannels > 0)  )
  {
    uint_t channel;
    rval = TRUE;
    memset( pThis, 0, sizeof( AudResampler ) );
    pThis->numTaps    = pCfg->numTaps;
    pThis->ratio      = pCfg->ratio;
    pThis->pChannels  = pChannelAry;
    pThis->numChannels = numChannels;
    pThis->isFixedPt = (NULL != pOptionalFixedPointCfg) ? TRUE : FALSE;
    {
      const int historySize0 = pThis->pChannels[0].xHistorySize;
      for (channel = 0; channel < numChannels; channel++) {
        AudResamplerChannel *pCh = &pThis->pChannels[channel];
        ASSERT( pCh->xHistorySize >= pThis->numTaps );
        rval &= (pCh->xHistorySize >= pThis->numTaps);
        rval &= (historySize0 == pCh->xHistorySize);
        ASSERT( rval );
#ifndef DO_OPT
        pCh->writeIdx   = 0;
#else
        pCh->writeIdx   = pThis->numTaps-1;
#endif
        if (!pThis->isFixedPt) {
          memset( pCh->arr.pXFltHistory, 0, sizeof( ar_float )* pCh->xHistorySize );
        }
        else {
          memset( pCh->arr.pXI32History, 0, sizeof( int32_t )* pCh->xHistorySize );
        }
      }
    }
    
    pThis->upsampleInnerLoops = pThis->numTaps / pThis->ratio;
    
    ASSERT( (pThis->upsampleInnerLoops * pThis->ratio) == pThis->numTaps );
    
    if ((rval) && (pThis->coefAry.pCoeffs == NULL)) {
      pThis->coefAry.pCoeffs = pCfg->pFirCoeffs;
      pThis->numTaps = pCfg->numTaps;
    }
    if (pThis->isFixedPt) {
      int c;
      ar_coef coeffMul;
      ASSERT( pOptionalFixedPointCfg->pFixedPtCoeffsBuf );
      ASSERT( pCfg->numTaps <= pOptionalFixedPointCfg->maxCoeffs );

      pThis->coefAry.i32.pCoeffs32 = pOptionalFixedPointCfg->pFixedPtCoeffsBuf;
      pThis->numTaps = MIN( pOptionalFixedPointCfg->maxCoeffs, pCfg->numTaps );
      
      // Get peak sample.
      {
        const double absPeak = audresample_getPeak(pCfg->pFirCoeffs, pThis->numTaps);
        ASSERT( 0 != absPeak );

        // coeffMul is calculated so that the biggest coefficient ends up with a 
        // value of (FIXED_POINT_MUL-1) before rounding
        coeffMul = (FIXED_POINT_MUL-1)/absPeak;
      }

      // Figure out what to multiply each result by to get a value from 0..1;
      pThis->coefAry.i32.postFilterScale = (ar_coef)(1.0/coeffMul);
      
      for (c = 0; c < pThis->numTaps; c++) {
        const ar_acc_t tmp = ROUND(coeffMul * pCfg->pFirCoeffs[c]);
        const ar_acc_t abstmp = ABS(tmp);
        ASSERT( abstmp <= FIXED_POINT_MUL);
        pThis->coefAry.i32.pCoeffs32[c] = (ar_coef_fix_t)(tmp);
      }
    }

    pThis->initCompleteFlag = AR_INIT_COMPLETE;
  }
  return rval;
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void AudResampleDownsample( AudResampler * const pThis, AudResampleIO *pIoAry, int numIos )
{
  if ((pThis) && (pThis->initCompleteFlag == AR_INIT_COMPLETE)) {
    //pAudFnResampProcess  AudResampleDoFn;
    if (pThis->isFixedPt) { 
      AudDnSampleProcessSamplesI16( pThis, pIoAry, numIos );
    }
    else {
      AudDnSampleProcessSamplesFlt( pThis, pIoAry, numIos );
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void AudResampleUpsample( AudResampler * const pThis, AudResampleIO *pIoAry, int numIos )
{
  if ((pThis) && (pThis->initCompleteFlag == AR_INIT_COMPLETE)) {
    if (!pThis->isFixedPt) {
      AudUpSampleProcessSamples( pThis, pIoAry, numIos );
    }
  }
}

}
