/******************************************************************************
Copyright 2014 Chris Fogelklou, Applaud Apps (Applicaudia)

This source code may under NO circumstances be distributed without the 
express written permission of the author.

@author: chris.fogelklou@gmail.com
*******************************************************************************/ 
#include <assert.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <float.h>
#include "audiolib_types.h"
#include "DbMonitor.h"
#include "MemPools.h"
#include "DaTunerDebug.h"
#include "DaTunerApi.h"
#include "Platform.h"
#include "Iir.hpp"

#ifdef DEBUG_FILES
#include <stdio.h>
#endif

#ifdef DEBUG_FILES
#define SAFE_FWRITE( pbuf, size, count, pfile) if(pfile) { fwrite( (pbuf), (size), (count), (pfile) ); }
#endif

#define  NF_DB_UP_PER_SECOND 1
#define NF_DB_DN_PER_SECOND 12

typedef struct {
  float                   dc;
  float                   fs;
  DaTunerDbMonParameters  prm;
  double                  sampPeriodS;
  bool_t                  signalOn;
  float                   nfDeltaPerSampDnLin;
  float                   nfDeltaPerSampUpLin;
  uint_t                  blockSize;
  IIR_CoefDesignFlt       lpf_slow;
  IIR_DF2HistoryFlt       lpf_slow_state;
  IIR_CoefDesignFlt       lpf_fast;
  IIR_DF2HistoryFlt       lpf_fast_state;
  pcm_t                   noise_floor;
  float                   currTrigLeveldB;
  float                   currTrigLevelLin;
#ifdef DEBUG_FILES
  FILE *pFileTriggers;
  FILE *pRawInput;
#endif
} DbMonitor;

static DbMonitor dbmon;

static const DaTunerDbMonParameters defaults = {
  TRUE,
  12,                       // fAutoThresholdDbAboveTrigger
  -48,                      // fManualThresholdDBFS
  3,                        // fHysteresisDb
  2,                        // fastEnvelopeLpfFreq
  0.75,                     // slowEnvelopeLpfFreq
  NF_DB_UP_PER_SECOND,      // fNoiseFloorDbUpPerSecond
  NF_DB_DN_PER_SECOND,      // fNoiseFloorDbDownPerSecond
  TRUE                      // useFftNrTrigger;
};

extern "C" {

// ////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////////////////////
void DbMonitorGetDefaults( DaTunerDbMonParameters *pInitData )
{
  *pInitData = defaults;
}


// ////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////////////////////
void DbMonitorGetCurrent( DaTunerDbMonParameters *pInitData )
{
  *pInitData = dbmon.prm;
}


// ////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////////////////////
void DbMonitorInit( 
  const float fs, 
  const uint_t blockSize,
  DaTunerDbMonParameters *pInitData )
{
  memset( &dbmon, 0, sizeof(dbmon) );

  DbMonitorGetDefaults( &dbmon.prm );

  dbmon.fs = fs;
  dbmon.sampPeriodS = 1.0f/fs;
  dbmon.blockSize = blockSize;

  // Start with a higher noise floor and let it stabilize.
  dbmon.noise_floor = powf(10, -40/20.0);

  DbMonitorChangeThresholds( pInitData );

  // Create a very slow LPF to capture the DC.
  ASSERT_FN( IIR_InitQ( &dbmon.lpf_slow, IIR_LPF, 0.0f, pInitData->fSlowEnvelopeLpfFreq, fs, 1.0f ) );
  ASSERT_FN( IIR_InitQ( &dbmon.lpf_fast, IIR_LPF, 0.0f, pInitData->fFastEnvelopeLpfFreq, fs, 1.0f ) );

#ifdef DEBUG_FILES
  dbmon.pFileTriggers  = fopen( "../outbin/dbmon_triggers.raw", "wb" );
  dbmon.pRawInput  = fopen( "../outbin/dbmon_insamps.raw", "wb" );
#endif
}


// ////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////////////////////
void DbMonitorClose( void )
{
#ifdef DEBUG_FILES
  fclose(dbmon.pFileTriggers);
  fclose(dbmon.pRawInput);
#endif
}


// ////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////////////////////
void DbMonitorChangeThresholds( DaTunerDbMonParameters * const pInitData )
{
  if (pInitData->fManualThresholdDBFS <= 0) {
    dbmon.prm.fManualThresholdDBFS = pInitData->fManualThresholdDBFS;
  }
  if (pInitData->fAutoThresholdDbAboveTrigger > 0) {
    dbmon.prm.fAutoThresholdDbAboveTrigger = pInitData->fAutoThresholdDbAboveTrigger;
  }

  dbmon.prm.useAutoThreshold = pInitData->useAutoThreshold;

  if (pInitData->fHysteresisDb > 0) {
    dbmon.prm.fHysteresisDb = pInitData->fHysteresisDb;
  }
  {
    float db;
    if(pInitData->fNoiseFloorDbUpPerSecond > 0) {
      dbmon.prm.fNoiseFloorDbUpPerSecond = pInitData->fNoiseFloorDbUpPerSecond;
      db   = (float) (dbmon.sampPeriodS * dbmon.prm.fNoiseFloorDbUpPerSecond); 
      dbmon.nfDeltaPerSampUpLin  = powf( 10, db/20 );
    }
    if(pInitData->fNoiseFloorDbDownPerSecond > 0) {
      dbmon.prm.fNoiseFloorDbDownPerSecond = pInitData->fNoiseFloorDbDownPerSecond;
      db   = (float) (0 - dbmon.sampPeriodS * dbmon.prm.fNoiseFloorDbDownPerSecond); 
      dbmon.nfDeltaPerSampDnLin  = powf( 10, db/20 );
    }
  }
  dbmon.currTrigLeveldB = 0;
  dbmon.currTrigLevelLin = powf( 10,  dbmon.currTrigLeveldB/20 );
}


// ////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////////////////////
void DbMonitorChangeManualSensDb( const float manDb ) 
{
  
  DaTunerDbMonParameters prm;
  DbMonitorGetCurrent(&prm);
  if (manDb >= 0) {
    // Trigger level is at points that can be met.
    prm.useAutoThreshold = TRUE;
    TRACE(("manDb should be in range [-96..0]"));
  } else {
    // User wants to use manual threshold.
    prm.useAutoThreshold = FALSE;
    prm.fManualThresholdDBFS = manDb;
  }

  DbMonitorChangeThresholds( &prm );
}


// ////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////////////////////
void DbMonitorChangeAutoSensDb( const float autoDb ) 
{
  
  DaTunerDbMonParameters prm;
  DbMonitorGetCurrent(&prm);
  if (autoDb <= 0) {
    // Trigger level is at points that can be met.
    prm.useAutoThreshold = FALSE;
    TRACE(("autoDb should be in range [0..96]"));
  } else {
    // User wants to use auto threshold.
    prm.useAutoThreshold = TRUE;
    prm.fAutoThresholdDbAboveTrigger = autoDb;
  }

  DbMonitorChangeThresholds( &prm );
}


// ////////////////////////////////////////////////////////////////////////////
// Processes the incoming data.
// Returns -1 if the current state is untriggered, otherwise returns
// the index of the first "triggered" sample within the buffer.
// ////////////////////////////////////////////////////////////////////////////
int32_t DbMonitorProcessSamples( 
  const pcm_t *pSamps, 
  const int32_t samples, 
  float * const pdbFSSignal, 
  float * const pdbFSNoiseFloor,
  float * const pdbFSCurrThreshold)
{
  int32_t i;
  int32_t triggeredSamp = -1;

  // get the running positive and negative peaks.
  const float MIN_LEVEL = powf(10.0f, -98.0f/20.0f);
  pcm_t signalEnv = MIN_LEVEL;
  float dbFSSignal = -100;
  float dbFSNoiseFloor = -100;

  ASSERT( samples > 0);
  for (i = 0; i < samples; i++) {
    const float absSamp = (pSamps[i] < 0) ? -pSamps[i] : pSamps[i];

    // LPF the slow envelope
    const pcm_t signalEnvSlow = IIR_DF2ProcessSingleSampleFlt( &dbmon.lpf_slow, &dbmon.lpf_slow_state, absSamp );

    // LPF the signal envelope
    const pcm_t signalEnvFast  = IIR_DF2ProcessSingleSampleFlt( &dbmon.lpf_fast, &dbmon.lpf_fast_state, absSamp );

    signalEnv      = MAX( signalEnvSlow, signalEnvFast );
    signalEnv      = MAX( MIN_LEVEL, signalEnv );

    if (signalEnvFast > dbmon.noise_floor ) {
      dbmon.noise_floor *= dbmon.nfDeltaPerSampUpLin;
    } else {
      dbmon.noise_floor *= dbmon.nfDeltaPerSampDnLin;
    }

    dbmon.noise_floor = MAX( MIN_LEVEL, dbmon.noise_floor );
    dbmon.noise_floor = MIN( 1.0f, dbmon.noise_floor );

    {
#ifdef DEBUG_FILES
      float trig = 0.0;
#endif
      if (signalEnv >= dbmon.currTrigLevelLin) {
        triggeredSamp = i;
#ifdef DEBUG_FILES
        trig = 1.0;
#endif
      }
#ifdef DEBUG_FILES
      {
        pcm_t triggersAry[] = {pSamps[i], signalEnvSlow, signalEnvFast, signalEnv, dbmon.noise_floor, trig};
        SAFE_FWRITE( triggersAry, sizeof( pcm_t ), ARRSZ(triggersAry), dbmon.pFileTriggers );
        SAFE_FWRITE( &pSamps[i], sizeof( pcm_t ), 1, dbmon.pRawInput );
      }
#endif
    }
  }

  dbFSSignal  = 20 * log10f( signalEnv );
  dbFSNoiseFloor = 20 * log10f( dbmon.noise_floor );
  if (dbmon.prm.useAutoThreshold) {
    dbmon.currTrigLeveldB = 20 * log10f( dbmon.noise_floor ) + dbmon.prm.fAutoThresholdDbAboveTrigger;
    dbmon.currTrigLevelLin = powf( 10,  dbmon.currTrigLeveldB/20 );
  }
  else {
    dbmon.currTrigLeveldB = dbmon.prm.fManualThresholdDBFS;
    dbmon.currTrigLevelLin = powf( 10, dbmon.prm.fManualThresholdDBFS/20 );
  }


  if (pdbFSSignal) {
    *pdbFSSignal = dbFSSignal;
  }
  if (pdbFSNoiseFloor) {
    *pdbFSNoiseFloor = dbFSNoiseFloor;
  }
  if (pdbFSCurrThreshold) {
    *pdbFSCurrThreshold = dbmon.currTrigLeveldB;
  }

  return triggeredSamp;
}

}
