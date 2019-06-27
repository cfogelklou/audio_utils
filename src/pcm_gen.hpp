/******************************************************************************
  Copyright 2014 Chris Fogelklou, Applaud Apps (Applicaudia)

  This source code may under NO circumstances be distributed without the
  express written permission of the author.

  @author: chris.fogelklou@gmail.com
*******************************************************************************/
#ifndef PCMGEN_HPP
#define PCMGEN_HPP
#include <math.h>
#include <float.h>
#include "audiolib_types.h"
#include "audutils_defines.h"
#include "audutils_debug.h"

// -----------------------------------------------------------------------------------------------
class PcmGen {

private:
  const static uint64_t M_WAVE_TABLE_SIZE_SHIFTS = 11ULL;
  const static uint64_t M_WAVE_TABLE_SIZE = (1ULL << M_WAVE_TABLE_SIZE_SHIFTS);
  const static uint64_t M_WAVE_MUL_SHIFTS = 44ULL;
  const static uint64_t M_WAVE_TABLE_MUL = (1ULL << M_WAVE_MUL_SHIFTS);
  const static uint64_t M_WAVE_TABLE_MAX_IDX =
      (1ULL << (M_WAVE_TABLE_SIZE_SHIFTS + M_WAVE_MUL_SHIFTS));
  const static int M_AMP_CHANGE_TIME_MS = 500;
  const static int M_AMPLITUDE_SHIFTS = 15;
  const static uint32_t M_MAX_AMPLITUDE = (32768U << M_AMPLITUDE_SHIFTS);
  int32_t mPrivateWaveTable[M_WAVE_TABLE_SIZE];
  uint64_t mWaveTablePhase;
  uint64_t mWaveTablePhasePerSample;
  double mFS;
  double mFreq;
  uint32_t mAmplitude;
  uint32_t mTargetAmplitude;
  uint32_t mIncPerSample;

  // ---------------------------------------------------------------------------------------------
  inline int32_t incAndGetValue() {
    mWaveTablePhase += mWaveTablePhasePerSample;
    while (mWaveTablePhase >= M_WAVE_TABLE_MAX_IDX) {
      mWaveTablePhase -= M_WAVE_TABLE_MAX_IDX;
    }
    uint32_t tblIdx = (uint32_t)((mWaveTablePhase >> M_WAVE_MUL_SHIFTS) & (M_WAVE_TABLE_SIZE - 1));
    return mPrivateWaveTable[tblIdx];
  }

public:
  // ---------------------------------------------------------------------------------------------
  PcmGen()
      : mWaveTablePhase(0), mWaveTablePhasePerSample(0), mFS(44100),
        mFreq(1000), mAmplitude(0), mTargetAmplitude(0) {
    ASSERT((M_WAVE_TABLE_SIZE_SHIFTS + M_WAVE_MUL_SHIFTS) < 62);
    ASSERT(M_WAVE_TABLE_MAX_IDX == (((uint64_t)M_WAVE_TABLE_SIZE) << M_WAVE_MUL_SHIFTS));

    setFs(mFS);
  }

  // ---------------------------------------------------------------------------------------------
  bool setUpWaveProperties(float *const harmonicStrenghtsArr,
                           const int numHarmonics) {
    bool rval = false;
    double totalStrength = 0;
    for (int i = 0; i < numHarmonics; i++) {
      ASSERT((harmonicStrenghtsArr[i] >= 0) &&
             (harmonicStrenghtsArr[i] < FLT_MAX));
      totalStrength += harmonicStrenghtsArr[i];
    }
    ASSERT(totalStrength > 0);
    const double scale = 32760.0 / totalStrength;
    // Clear out the wave table.
    for (int i = 0; i < M_WAVE_TABLE_SIZE; i++) {
      mPrivateWaveTable[i] = 0;
    }

    // Add the lookup table for each index.
    for (int i = 0; i < M_WAVE_TABLE_SIZE; i++) {
      double newAmp = 0;
      for (int h = 0; h < numHarmonics; h++) {
        newAmp +=
            (float)(scale * harmonicStrenghtsArr[h] *
                    sin(i * (h + 1) * 2.0 * M_PI / (double)M_WAVE_TABLE_SIZE));
      }
      mPrivateWaveTable[i] = ROUND(newAmp);
    }

    for (int i = 0; i < M_WAVE_TABLE_SIZE; i++) {
      ASSERT(mPrivateWaveTable[i] <= 32767);
      ASSERT(mPrivateWaveTable[i] >= -32767);
    }
    rval = true;
    return rval;
  }

  // ---------------------------------------------------------------------------------------------
  bool setFs(double fs) {
    bool rval = false;
    if (fs > 0) {
      mFS = fs;
      {
        double sampsForVolumeChange =
            (double)M_AMP_CHANGE_TIME_MS * (double)mFS / 1000.0;
        mIncPerSample = ROUND(((double)M_MAX_AMPLITUDE) / sampsForVolumeChange);
        mIncPerSample = MAX(1, mIncPerSample);
      }
      setFreq(mFreq);
      rval = true;
    }
    return rval;
  }

  // ---------------------------------------------------------------------------------------------
  bool setAmplitude(double amplitude) {
    bool rval = false;
    if (mFS > 0) {
      amplitude = MAX(0.0, amplitude);
      uint64_t newTarget = ROUND(amplitude * M_MAX_AMPLITUDE);
      mTargetAmplitude = (uint32_t)(MIN(newTarget, M_MAX_AMPLITUDE));
      rval = true;
    }
    return rval;
  }

  // ---------------------------------------------------------------------------------------------
  bool setFreq(double freq) {
    bool rval = false;
    if (mFS > 0) {
      mFreq = freq;
      const double fInc =
          ((mFreq / mFS) * (double)(M_WAVE_TABLE_MUL * M_WAVE_TABLE_SIZE));
      mWaveTablePhasePerSample = (uint64_t)ROUNDF(fInc);
      ASSERT(0 == (mWaveTablePhasePerSample & 0x8000000000000000ULL));
      rval = true;
    }
    return rval;
  }

  // ---------------------------------------------------------------------------------------------
  bool setUpGenProperties(double freq, double amplitude) {
    bool rval = false;
    if (mFS > 0) {
      rval = setFreq(freq);
      rval &= setAmplitude(amplitude);
    }
    return rval;
  }

  // ---------------------------------------------------------------------------------------------
  bool doGenerate(short *const x, int startingIdx, int nBufLen) {
    if ((mAmplitude == 0) && (mTargetAmplitude == 0)) {
      return false;
    } else {
      int len = nBufLen - startingIdx;
      ASSERT(len > 0);
      len = MAX(0, len);
      int i = 0;
      const int32_t diff = mTargetAmplitude - mAmplitude;
      if (0 != diff) {
        const int absdiff = ABS(diff);
        const int inc = (diff < 0) ? -(int)mIncPerSample : mIncPerSample;
        const int incIterations = 1 + ((absdiff - 1) / mIncPerSample);
        const int loops = MIN(len, incIterations);
        for (i = 0; i < loops; i++) {
        x[i + startingIdx] = (short)(((mAmplitude>>M_AMPLITUDE_SHIFTS) * incAndGetValue()) >> 15);
          mAmplitude += inc;
        }
        if ((loops == 0) || (i == incIterations)) {
          mAmplitude = mTargetAmplitude;
        }
      }
      {
        const uint32_t amplitude16 = (mAmplitude >> M_AMPLITUDE_SHIFTS);
        for (; i < len; i++) {
          x[i + startingIdx] = (short)((amplitude16 * incAndGetValue()) >> 15);
        }
      }
      return true;
    }
  }
};

#endif
