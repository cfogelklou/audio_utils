/******************************************************************************
  Copyright 2014 Chris Fogelklou, Applaud Apps (Applicaudia)

  This source code may under NO circumstances be distributed without the 
  express written permission of the author.

  @author: chris.fogelklou@gmail.com
*******************************************************************************/ 

#ifndef __DbMonitor
#define __DbMonitor
#include "audiolib_types.h"
#include "PcmQ.h"
#include "DaTunerApi.h"

#ifdef __cplusplus
extern "C" {
#endif

void DbMonitorGetDefaults( DaTunerDbMonParameters *pInitData );

void DbMonitorGetCurrent( DaTunerDbMonParameters *pInitData );

void DbMonitorInit( 
  const float fs, 
  const uint_t blockSize,
  DaTunerDbMonParameters *pInitData );

void DbMonitorClose( void );

void DbMonitorChangeThresholds( DaTunerDbMonParameters * const pInitData );

void DbMonitorChangeManualSensDb( const float manDb );

void DbMonitorChangeAutoSensDb( const float autoDb ); 

  // Processes the incoming data.
  // Returns -1 if the current state is untriggered, otherwise returns
  // the index of the first "triggered" sample within the buffer.
int32_t DbMonitorProcessSamples( 
  const pcm_t *pSamps, 
  const int32_t samples, 
  float * const pdbFSSignal, 
  float * const pdbFSNoiseFloor,
  float * const pdbFSCurrThreshold);


#ifdef __cplusplus
}
#endif

#endif
