/******************************************************************************
Copyright 2015 Chris Fogelklou, Applaud Apps (Applicaudia)

This source code may under NO circumstances be distributed without the
express written permission of the author.

@author: chris.fogelklou@gmail.com
*******************************************************************************/
#ifndef SINELUT_HPP
#define SINELUT_HPP
#include <string.h>
#include <math.h>
#include <stddef.h>

#include "HelperMath.hpp"
#include "audiolib_types.h"

#define SINELUT_SIZE 4096
#define SINELUT_MASK (SINELUT_SIZE - 1)

template <typename TSineFloat>
class SineLut {

private:

    // Singleton instance.
    static SineLut *mInst;

    // / Lookup table for sine values
    TSineFloat mSinLut[(SINELUT_SIZE / 2) + 4];

    // Initializes or re-initializes the FFT.
    SineLut() {
        const double mul = MREQ_TWO_PI / SINELUT_SIZE;
        for (int i = 0; i < (SINELUT_SIZE / 2); i++) {
            mSinLut[i] = (TSineFloat)sin((double)i * mul);
        }
    }

public:
    // ---------------------------------------------------------------------------

    static SineLut &inst();


    // ---------------------------------------------------------------------------
    // Uses, as an input, a value from 0 .. SINELUT_SIZE
    inline TSineFloat getSin(const int tblIdx) {
        const unsigned int tblIdxMasked = ((unsigned int)tblIdx) & SINELUT_MASK;
        const int tblIdxMinusFftSizeOver2 = tblIdxMasked - (SINELUT_SIZE / 2);
        return (tblIdxMinusFftSizeOver2 >= 0) ? -mSinLut[tblIdxMinusFftSizeOver2] : mSinLut[tblIdxMasked];
    }

    // ---------------------------------------------------------------------------
    // Uses, as an input, a value from 0 .. SINELUT_SIZE
    inline TSineFloat getCos(const int tblIdx) {
        const int newIdx = (tblIdx + (SINELUT_SIZE / 4)) & SINELUT_MASK;
        return getSin(newIdx);
    }

    // ---------------------------------------------------------------------------
    // Uses, as an input, a value from 0 .. 1
    inline TSineFloat getSinF(const TSineFloat fIdx) {
        const TSineFloat lu0 = SINELUT_SIZE * fIdx;
        const int x = (lu0 < 0) ? (int)(lu0 - 0.5) : (int)(lu0 + 0.5); // round
        return getSin(x);
    }

    // ---------------------------------------------------------------------------
    // Uses, as an input, a value from 0 .. 1
    inline TSineFloat getCosF(const TSineFloat fIdx) {
        const TSineFloat lu0 = SINELUT_SIZE * fIdx;
        const int x = (lu0 < 0) ? (int)(lu0 - 0.5) : (int)(lu0 + 0.5); // round
        return getCos(x);
    }

    // ---------------------------------------------------------------------------
    // fIdx 0..1 in ratio of radians.
    inline TSineFloat getSinFInterp(const TSineFloat fIdx) {
        const TSineFloat luIdx = SINELUT_SIZE * fIdx;
        const TSineFloat lu0 = floor(luIdx); // the lowest index
        const TSineFloat diffx = luIdx - lu0; // how far along from the lowest index.
        const int nlu0 = ((int)lu0);
        TSineFloat y0 = getSin(nlu0);
        TSineFloat slope = getSin(nlu0 + 1) - y0;
        return y0 + (slope * diffx);
    }

    // ---------------------------------------------------------------------------
    // fIdx 0..1 in ratio of radians.
    inline TSineFloat getCosFInterp(const TSineFloat fIdx) {
        const TSineFloat luIdx = SINELUT_SIZE * fIdx;
        const TSineFloat lu0 = floor(luIdx); // the lowest index
        const TSineFloat diffx = luIdx - lu0; // how far along from the lowest index.
        const int nlu0 = ((int)lu0);
        TSineFloat y0 = getCos(nlu0);
        TSineFloat slope = getCos(nlu0 + 1) - y0;
        return y0 + (slope * diffx);
    }

    // ---------------------------------------------------------------------------
    // fIdx 0..2PI translates to appropriate sine value, just like sin(x)
    inline TSineFloat getSinFRad(const TSineFloat fIdx) {
        const TSineFloat radians = fIdx / MREQ_TWO_PI;
        return getSinF(radians);
    }

    // ---------------------------------------------------------------------------
    // fIdx 0..2PI translates to appropriate cosine value, just like cos(x)
    inline TSineFloat getCosFRad(const TSineFloat fIdx) {
        const TSineFloat radians = fIdx / MREQ_TWO_PI;
        return getCosF(radians);
    }


    // ---------------------------------------------------------------------------
    // fIdx 0..2PI translates to appropriate sine value, just like sin(x)
    inline TSineFloat getSinFRadInterp(const TSineFloat fIdx) {
        const TSineFloat radians = fIdx / MREQ_TWO_PI;
        return getSinFInterp(radians);
    }

    // ---------------------------------------------------------------------------
    // fIdx 0..2PI translates to appropriate cosine value, just like cos(x)
    inline TSineFloat getCosFRadInterp(const TSineFloat fIdx) {
        const TSineFloat radians = fIdx / MREQ_TWO_PI;
        return getCosFInterp(radians);
    }
};

#endif // SINELUT_HPP
