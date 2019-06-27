#ifndef AUDRESAMPLER_H__
#define AUDRESAMPLER_H__

/******************************************************************************
  Copyright 2014 Chris Fogelklou, Applaud Apps (Applicaudia)

  This source code may under NO circumstances be distributed without the 
  express written permission of the author.

  @author: chris.fogelklou@gmail.com
*******************************************************************************/


#include "audiolib_types.h"


#ifdef __cplusplus
extern "C" {
#endif 

typedef float ar_float;
typedef double ar_coef;

#define AR_32_BIT_FIXED 1
#if (AR_32_BIT_FIXED > 0)
#define FIXED_POINT_MUL ((double)((1 << 30)-1))
typedef int32_t ar_coef_fix_t;
typedef int64_t ar_acc_t;
#else
#define FIXED_POINT_MUL 32767.0
typedef int16_t ar_coef_fix_t;
typedef int32_t ar_acc_t;
#endif


// Standard configuration.
typedef struct tAudResampleConfig
{
  int_t         numTaps;
  uint_t        ratio;
  int_t         bits;
  const ar_coef  *pFirCoeffs;
} AudResampleConfig;

extern const AudResampleConfig DownsampleX1CoeffsCfg;
extern const AudResampleConfig DownsampleX2CoeffsCfg;
extern const AudResampleConfig DownsampleX3CoeffsCfg;
extern const AudResampleConfig DownsampleX4CoeffsCfg;
extern const AudResampleConfig DownsampleX5CoeffsCfg;
extern const AudResampleConfig DownsampleX6CoeffsCfg;
extern const AudResampleConfig DownsampleX7CoeffsCfg;
extern const AudResampleConfig DownsampleX8CoeffsCfg;
extern const AudResampleConfig DownsampleX9CoeffsCfg;
extern const AudResampleConfig DownsampleX10CoeffsCfg;


typedef struct tAudResampleIOFloat {
    // Pointer to the input samples.
    ar_float *pxIn;

    // pointer to the sample buffer being operated on.
    ar_float *pyOut;
} AudResampleIOFloat;

typedef struct tAudResampleIOInt16 {
    // Pointer to the input samples.
    int16_t *pxIn;

    // pointer to the sample buffer being operated on.
    ar_float *pyOut;
} AudResampleIOInt16;


typedef struct tAudResampleIO {

  union {
    struct {
      // Pointer to the input samples.
      ar_float *pxIn;

      // pointer to the sample buffer being operated on.
      ar_float *pyOut;
    } flt;
    struct {
      // Pointer to the input samples.
      int16_t *pxIn;

      // pointer to the sample buffer being operated on.
      ar_float *pyOut;
    } i16;
  } bufs;

  // Number of samples in the X buffer.
  int   numXSamples;

  // Spacing between input channels (for interleaved inputs.)
  int   xHop;

  // Number of samples in the Y buffer.
  int   numYSamples;

  // Spacing between input channels (for interleaved inputs.)
  int   yHop;
} AudResampleIO;

struct tAudResampler;


typedef struct tAudResamplerChannel {

  union {
    // Pointer to buffer to keep old X samples.
    ar_float              *pXFltHistory;
    // Pointer to buffer to keep old X samples.
    int32_t               *pXI32History;
  } arr;

  // Size of the history buffer.
  int_t               xHistorySize;

  // Current history write index (used for reference version.)
  int_t               writeIdx;


} AudResamplerChannel;


typedef struct tAudResampler
{

  // Number of taps in the filter coefficients.
  int_t               numTaps;

  // Ratio of input / output frequency.  Always >= 1.
  int_t               ratio;

  bool_t              isFixedPt;
  union {
    // Pointer to the filter coefficients.
    const ar_coef        *pCoeffs;
    struct {
      ar_coef_fix_t       *pCoeffs32;
      ar_coef              postFilterScale;
    } i32;
  }coefAry;

  // Internal data for each channel (history of incoming samples, pointer into the history, etc.)
  AudResamplerChannel *pChannels;

  // Number of channels in the pChannels array.
  uint_t              numChannels;

  // Number of inner loops to run.
  int                 upsampleInnerLoops;

  uint32_t            initCompleteFlag;

} AudResampler;


typedef struct tAudResamplerFixedPointCfg {
  ar_coef_fix_t             *pFixedPtCoeffsBuf;
  int_t                maxCoeffs;
} AudResamplerFixedPointCfg;

// Initialize a resampler object
bool_t AudResampleInit( 
  AudResampler                 *pThis, 
  bool_t                        isUpsampler, 
  const AudResampleConfig      *pCfg, 
  AudResamplerChannel          *pChannelAry, 
  uint_t                        numChannels,
  AudResamplerFixedPointCfg    *pOptionalFixedPointCfg);

// Calls the function pointer in the resampler object.
void AudResampleDownsample( AudResampler * const pThis, AudResampleIO *pIoAry, int numIos );


void AudResampleUpsample( AudResampler * const pThis, AudResampleIO *pIoAry, int numIos );

#ifdef __cplusplus
}
#endif 

#endif
