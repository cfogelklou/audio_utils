#include "sine_lut.hpp"

// Singleton for float Sine lookups
template<>
SineLut<float> *SineLut<float>::mInst = NULL;

template<>
SineLut<float> & SineLut<float>::inst()
{
    if (NULL == mInst) {
        mInst = new SineLut();
    }
    return *mInst;
}

// Singleton for double Sine lookups
template<>
SineLut<double> *SineLut<double>::mInst = NULL;
template<>
SineLut<double> & SineLut<double>::inst()
{
    if (NULL == mInst) {
        mInst = new SineLut();
    }
    return *mInst;
}
