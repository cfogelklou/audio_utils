#ifndef BIQUADFIXED_H__
#define BIQUADFIXED_H__

/******************************************************************************
  Copyright 2014 Chris Fogelklou, Applaud Apps (Applicaudia)

  This source code may under NO circumstances be distributed without the
  express written permission of the author.

  @author: chris.fogelklou@gmail.com
*******************************************************************************/

#include <math.h>
#include "audiolib_types.h"
#include "Iir.hpp"
#include "HelperMath.hpp"

class FixedIir {
private:
  const static int FIXIIR_SHIFTS = 12;
  const static int FIXIIR_ROUND = (1 << (FIXIIR_SHIFTS - 1));

public:
  int32_t mA0;
  int32_t mA1;
  int32_t mA2;
  int32_t mB0;
  int32_t mB1;
  int32_t mB2;
  int32_t mP1;
  int32_t mP2;

  FixedIir(FixedIir &iir)
      : mA0(iir.mA0), mA1(iir.mA1), mA2(iir.mA2), mB0(iir.mB0), mB1(iir.mB1),
        mB2(iir.mB2), mP1(0), mP2(0) {}

  FixedIir(IIR_CoefDesign &iir)
      : mA0(1), mA1(0), mA2(0), mB0(1), mB1(0), mB2(0), mP1(0), mP2(0)

  {
    ConvertFromIir(iir);
  }

  FixedIir()
      : mA0(1), mA1(0), mA2(0), mB0(1), mB1(0), mB2(0), mP1(0), mP2(0)

  {}

  void doReset() {
    mP1 = 0;
    mP2 = 0;
  }

  /*---------------------------------------------------------------------*/ /**
    *
    * Copies and scales the coefficients from the input IIR.  Assumes that
    * the input IIR coefficients are scaled to 1.0 (processing does not
    * do amplitude scaling on output so output amplitude may not match input amp
    * litude, but frequency response should match.)
    *
    * @param samp
    *            : The sample to truncate
    * @return: returns truncated sample.
    */
  /*----------------------------------------------------------------------*/
  void ConvertFromIir(IIR_CoefDesign &iir) {
    double maxCoeff = MAX(HP_abs(1.0), HP_abs(iir.a1));
    maxCoeff = MAX(maxCoeff, HP_abs(iir.a2));
    maxCoeff = MAX(maxCoeff, HP_abs(iir.b0));
    maxCoeff = MAX(maxCoeff, HP_abs(iir.b1));
    maxCoeff = MAX(maxCoeff, HP_abs(iir.b2));
    const double scaling = ((double)(1 << FIXIIR_SHIFTS)) / maxCoeff;
    mA0 = (int)ROUND(scaling * 1.0);
    mA1 = (int)ROUND(scaling * iir.a1);
    mA2 = (int)ROUND(scaling * iir.a2);
    mB0 = (int)ROUND(scaling * iir.b0);
    mB1 = (int)ROUND(scaling * iir.b1);
    mB2 = (int)ROUND(scaling * iir.b2);
    mP1 = 0;
    mP2 = 0;
  }

  /*---------------------------------------------------------------------*/ /**
    *
    * Clips a 32 bit int to a 16-bit int16_t
    *
    * @param samp
    *            : The sample to truncate
    * @return: returns truncated sample.
    */
  /*----------------------------------------------------------------------*/
  static int16_t clip16(int samp) {
    int16_t rval;
    if (samp < -32768) {
      rval = (int16_t)-32768;
    } else if (samp > 32767) {
      rval = (int16_t)32767;
    } else {
      rval = (int16_t)samp;
    }
    return rval;
  }

  /*---------------------------------------------------------------------*/ /**
    *
    * Clips a 64 bit int64_t to a 32-bit int
    *
    * @param samp
    *            : The sample to truncate
    * @return: returns truncated sample.
    */
  /*----------------------------------------------------------------------*/
  static int clip32(int64_t samp) {
    int rval;
    if (samp < -2147483648) {
      rval = (int)-2147483648;
    } else if (samp > 2147483647) {
      rval = (int)2147483647;
    } else {
      rval = (int)samp;
    }
    return rval;
  }

  /*---------------------------------------------------------------------*/ /**
    *
    * Process a single sample of a biquad.
    *
    * @param sample
    *            : The next sample
    * @return: returns the next sample.
    */
  /*----------------------------------------------------------------------*/
  int16_t processSingleSample(int16_t sample) {
    int p0;

    // Note t64 might need to be 64 bits if FIXIIR_SHIFTS is increased!
    int64_t t64 = FIXIIR_ROUND + (sample << FIXIIR_SHIFTS) - (mP1 * mA1) - (mP2 * mA2);
    p0 = (int)(t64 >> FIXIIR_SHIFTS);

    t64 = FIXIIR_ROUND + (p0 * mB0) + (mP1 * mB1) + (mP2 * mB2);

    // shift p0 to p1 and p1 to p2
    mP2 = mP1;
    mP1 = p0;

    return clip16((int)t64 >> FIXIIR_SHIFTS);
  }

  /*---------------------------------------------------------------------*/ /**
    *
    * Process a batch of samples in pSamplesIn. Output to pSamplesOut.
    *
    * @param pSamplesIn
    *            : The array of samples.
    * @param nSamples
    *            : The number of samples in nSamples.
    * @param pSamplesOut
    *            : Spits out the output samples.
    */
  /*----------------------------------------------------------------------*/
  void processSamples(int16_t pSamplesIn[], int nSamples,
                      int16_t pSamplesOut[]) {
    int p0;

    for (int i = 0; i < nSamples; i++) {
      const int sample = pSamplesIn[i];

      // Note t64 might need to be 64 bits if FIXIIR_SHIFTS is increased!
      int64_t t64 = FIXIIR_ROUND + (sample << FIXIIR_SHIFTS) - (mP1 * mA1) - (mP2 * mA2);
      p0 = (int)(t64 >> FIXIIR_SHIFTS);

      t64 = FIXIIR_ROUND + (p0 * mB0) + (mP1 * mB1) + (mP2 * mB2);

      // shift p0 to p1 and p1 to p2
      mP2 = mP1;
      mP1 = p0;

      pSamplesOut[i] = clip16((int)t64 >> FIXIIR_SHIFTS);
    }
  }

  /*---------------------------------------------------------------------*/ /**
    *
    * Process a batch of samples in pSamplesIn. Output to pSamplesOut.
    *
    * @param pSamplesIn
    *            : The array of samples.
    * @param nSamples
    *            : The number of samples in nSamples.
    * @param pSamplesOut
    *            : Spits out the output samples.
    */
  /*----------------------------------------------------------------------*/
  void processSamples(double pSamplesIn[], int nSamples, double pSamplesOut[]) {
    int p0;
    for (int i = 0; i < nSamples; i++) {
      const int sample = (int)ROUND(pSamplesIn[i] * 32768.0);

      // Note t64 might need to be 64 bits if FIXIIR_SHIFTS is increased!
      int64_t t64 = FIXIIR_ROUND + (sample << FIXIIR_SHIFTS) - (mP1 * mA1) - (mP2 * mA2);
      p0 = (int)(t64 >> FIXIIR_SHIFTS);

      t64 = FIXIIR_ROUND + (p0 * mB0) + (mP1 * mB1) + (mP2 * mB2);

      // shift p0 to p1 and p1 to p2
      mP2 = mP1;
      mP1 = p0;

      pSamplesOut[i] = clip16((int)t64 >> FIXIIR_SHIFTS);
    }
  }
}
