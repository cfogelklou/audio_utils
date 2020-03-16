/******************************************************************************
  Copyright 2014 Chris Fogelklou, Applaud Apps (Applicaudia)

  This source code may under NO circumstances be distributed without the 
  express written permission of the author.

  @author: chris.fogelklou@gmail.com
*******************************************************************************/ 
#ifndef __Defines
    #define __Defines

#define OSAL_DEBUG_MALLOC

#include <math.h>
#ifndef PI
#define PI        (3.141592653589793238462643383279502884197169399375105820974944592307816406286)
#endif
#ifndef TWOPI
#define TWOPI     (2*PI)
#endif
#define HALFPI    (PI/2)
#define PI2TENTH  (TWOPI/10)
#define PISQUARED (PI*PI)
#define BASELN      (2.71828182845904523536028747135266)
#ifndef LN10
#define LN10        (2.30258509299404568401799145468436)
#endif
#ifndef LN2
#define LN2         (0.693147180559945309417232121458)
#endif
#define LN1000      (3.0 * 2.30258509299404568401799145468436)
#define LOG10E      (0.434294481903251827651128918916605)
#define LOG1000E    (LOG10E/3.0)
#define LN10DIV20   (0.115129254649702284200899572734218)
#define LOG10EMUL20 (8.6858896380650365530225783783321)

#ifndef M_PI
#define M_PI      PI
#endif

#define MA_Pow(x,p)     pow((double)x, (double)p)
#define MA_Pow2(power)	MA_Pow(2, power)
#define MA_Sqrt(x)		  sqrt(x)
#define MA_Ln(x)		    log(x)
#define MA_Log(x)       log10(x)
#define MA_Divide(x, y) ((x) / (y))
#define MA_Sin(x)       sin(x)
#define MA_Cos(x)       cos(x)
#define MA_Sinh(x)      sinh(x)
#define MA_Tan(x)       (MA_Sin(x)/MA_Cos(x))
#define MA_ATan(x)      atan(x)
#define MA_Floor(x)     floor(x)
#define MA_Ceil(x)      ceil(x)
#ifndef M_LN2
#define M_LN2           LN2
#endif
#define M_2_OVER_LN2    (2.0/LN2)


#ifndef MIN
#define MIN(x,y) (((x) < (y)) ? (x) : (y))
#endif

#ifndef MAX
#define MAX(x,y) (((x) > (y)) ? (x) : (y))
#endif

#ifndef ARRSZ
#define ARRSZ(arr) (sizeof(arr)/sizeof(arr[0]))
#endif

#ifndef ROUNDF
#define ROUNDF(x) (((x)>=0)?((x)+0.5):((x)-0.5))
#endif

#ifndef ROUND
#if 0
// Allows us to compile in some checks to ensure that round is always used properly.
#include "DaTunerDebug.h"
static int doround(double x) {
  long long round1 = (long long)(ROUNDF(x));
  long long round2 = (((x)>=0)?(int)((x)+0.5):(int)((x)-0.5));
  ASSERT(round1 == round2);
  return (int)round2;
}
#define ROUND(x) doround((double)x)
#else
#define ROUND(x) (((x)>=0)?(int)((x)+0.5):(int)((x)-0.5))
#endif
#endif // #ifndef roundROUND

#ifndef ABS
#define ABS(x) ( ((x) < 0) ? (-(x)) : (x) )
#endif

#ifndef COMPILE_TIME_ASSERT
// Will not compile if cond is false.
#define COMPILE_TIME_ASSERT(cond) do{switch(((cond)?1:0)){case 0:break;case (cond):break;}}while(0)
#endif

#ifndef MAKE_MASK
#define MAKE_MASK(bits) ((1u << (bits)) - 1)
#endif

#ifndef SET_BITS
#define SET_BITS(x, y) (x |= (y))
#endif

#ifndef CLEAR_BITS
#define CLEAR_BITS(x, y) (x &= ~(y))
#endif

#endif
