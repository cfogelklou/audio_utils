/******************************************************************************
  Copyright 2014 Chris Fogelklou, Applaud Apps (Applicaudia)

  This source code may under NO circumstances be distributed without the 
  express written permission of the author.

  @author: chris.fogelklou@gmail.com
*******************************************************************************/
#include <string.h>
#include <math.h>
#include <stdint.h>
#include "pcm_q.h"
#include "audiolib_types.h"
#include "audutils_debug.h"
#include "audutils_defines.h"

#define Platform_EnterOSCritical()
#define Platform_ExitOSCritical()

  typedef int32_t msg_hdr_t;
#define MSG_ALIGN (sizeof(msg_hdr_t))
#define DOALIGN(_LEN)     (((_LEN) + MSG_ALIGN - 1) & (uint_t)(0-MSG_ALIGN))

//lint -emacro( 767, PcmQ_LOCKIFDEFINED)
/* [Internal] Conditional lock, lock only if the mutex exists. */
#define PcmQ_LOCKIFDEFINED( mutexId )

//lint -emacro( 767, PcmQ_UNLOCK_IFDEFINED)
/* [Internal] Conditional unlock, unlock only if the mutex exists. */
#define PcmQ_UNLOCK_IFDEFINED( mutexId )

/* [Internal] Does an insertion at the current write pointer and count, but
does not update write pointer or count variables.
*/
static uint_t PcmQUnprotectedInsert( PcmQ_t * const pQ, const pcm_t *pWrBuf, uint_t nLen, uint_t *pNewWrIdx);

/*
**=============================================================================
**  Abstract:
**    Public function to initialize the PcmQ_t structure.
**
**  Parameters:
**
**  Return values:
**
**=============================================================================
*/
bool_t PcmQCreate( PcmQ_t * const pQ, pcm_t *pBuf, uint_t nBufSz, PcmQ_AccessFlags_t flags )
{

  ASSERT(( NULL != pQ ) && ( NULL != pBuf ) );

  memset( pQ, 0, sizeof( PcmQ_t ) );

  pQ->pfBuf        = pBuf;
  pQ->nBufSz      = nBufSz;

  if ( ( (( flags & PcmQ_NoCountMutexForWrite) != 0 ) &&    //lint !e641 !e655 OK: Just uint flags.
    (( flags & PcmQ_NoCountMutexForWrite) != 0 ) ) ||  //lint !e641 !e655 OK: Just uint flags.
    (( flags & PcmQ_NoMutexes ) != 0 ) )             //lint !e641 !e655 OK: Just uint flags.
  {
  }
  else
  {
    if (( flags & PcmQ_NoCountMutexForWrite) == 0) //lint !e641 !e655 OK: Just uint flags.
    {
      pQ->wrCntProt = TRUE;
    }

    // If reads take place from an interrupt context, don't protect the count variable
    // while reading
    if (( flags & PcmQ_NoCountMutexForRead) == 0 ) //lint !e641 !e655 OK: Just uint flags.
    {
      pQ->rdCntProt = TRUE;
    }
  }

  if ( ((flags & PcmQ_ReentrantWrites) || (flags & PcmQ_ReentrantReads)) && //lint !e641 !e655 OK: Just uint flags.
    (!( flags & PcmQ_NoMutexes )) ) //lint !e641 !e655 OK: Just uint flags.
  {

    // If writes should be re-entrant (more than one process writing)
    // then protect the write functions using a mutex
    if ( flags & PcmQ_ReentrantWrites ) //lint !e641 !e655 OK: Just uint flags.
    {
      ASSERT(( flags & PcmQ_NoCountMutexForWrite) == 0 ); //lint !e641 !e655 OK: Just uint flags.
    }

    // If writes should be re-entrant (more than one process reading)
    // then protect the read functions using a mutex
    if ( flags & PcmQ_ReentrantReads ) //lint !e641 !e655 OK: Just uint flags.
    {
      ASSERT(( flags & PcmQ_NoCountMutexForRead) == 0 ); //lint !e641 !e655 OK: Just uint flags.
    }
  }
  return TRUE;

}

/*
**=============================================================================
**  Abstract:
**  Public function to deallocate the things in the
**   Q that were allocated.
**
**  Parameters:
**
**  Return values:
**
**=============================================================================
*/
bool_t PcmQDestroy( PcmQ_t * const pQ )
{
  ASSERT( NULL != pQ );
  return TRUE;
}


/*
**=============================================================================
**  Abstract:
**  Write function that lacks protection from multiple threads trying to write
**  at the same time.  Use this function if only one thread/task will write
**  to the queue at a time.
**  Parameters:
**
**  Return values:
**
**=============================================================================
*/
uint_t PcmQWrite( PcmQ_t * const pQ, const pcm_t *pWrBuf, uint_t nLen )
{
  uint_t WordsWritten   = 0;

  ASSERT( NULL != pQ );

  if ( nLen )
  {
    // Write nothing if there isn't room in the buffer.
    uint_t WordsToWrite = 0;
    uint_t nWrIdx = pQ->nWrIdx;
    const uint_t nBufSz = pQ->nBufSz;
    pcm_t * const pBuf = pQ->pfBuf;

    PcmQ_LOCKIFDEFINED( pQ->wrMutex );

    WordsToWrite = (nLen <= (nBufSz - pQ->nCount)) ? nLen : 0;

    // We can definitely read WordsToWrite bytes.
    while (WordsToWrite > 0)
    {
      // Calculate how many contiguous bytes to the end of the buffer
      uint_t Words = MIN( WordsToWrite, (nBufSz - nWrIdx) );

      // Copy that many bytes.
      memcpy(&pBuf[ nWrIdx ], &pWrBuf[ WordsWritten ], Words*sizeof(pcm_t) );

      // Circular buffering.
      nWrIdx += Words;
      if (nWrIdx >= nBufSz)
      {
        nWrIdx -= nBufSz;
      }

      // Increment the number of bytes written.
      WordsWritten += Words;
      WordsToWrite -= Words;
    }

    pQ->nWrIdx = nWrIdx;

    // Increment the count.  (protect with mutex)
    if (pQ->wrCntProt){ Platform_EnterOSCritical(); }
    pQ->nCount = pQ->nCount + WordsWritten;
    if (pQ->wrCntProt){ Platform_ExitOSCritical(); }
    PcmQ_UNLOCK_IFDEFINED( pQ->wrMutex );
  }
  return WordsWritten;
}

/*
**=============================================================================
**=============================================================================
*/
static uint_t PcmQUnprotectedInsert( PcmQ_t * const pQ, const pcm_t *pWrBuf, uint_t nLen, uint_t *pNewWrIdx )
{
  uint_t WordsWritten   = 0;
  uint_t nWrIdx = pQ->nWrIdx;

  if ( nLen )
  {
    // Write nothing if there isn't room in the buffer.
    uint_t WordsToWrite = nLen;
    const uint_t nBufSz = pQ->nBufSz;
    pcm_t * const pBuf = pQ->pfBuf;

    while (WordsToWrite > 0)
    {
      uint_t Words = MIN( WordsToWrite, (nBufSz - nWrIdx) );

      memcpy(&pBuf[ nWrIdx ], &pWrBuf[ WordsWritten ], Words*sizeof(pcm_t) );

      // Circular buffering.
      nWrIdx += Words;
      if (nWrIdx >= nBufSz)
      {
        nWrIdx -= nBufSz;
      }

      // Increment the number of bytes written.
      WordsWritten += Words;
      WordsToWrite -= Words;
    }
  }
  if (NULL != pNewWrIdx)
  {
    *pNewWrIdx = nWrIdx;
  }
  return WordsWritten;
}


/*
**=============================================================================
**=============================================================================
*/
uint_t PcmQCommitWrite( PcmQ_t * const pQ, uint_t nLen )
{
  uint_t WordsWritten   = 0;

  ASSERT( NULL != pQ );

  if ( nLen )
  {
    // Write nothing if there isn't room in the buffer.
    PcmQ_LOCKIFDEFINED( pQ->wrMutex );

    // Circular buffering.
    pQ->nWrIdx += nLen;
    if (pQ->nWrIdx >= pQ->nBufSz)
    {
      pQ->nWrIdx -= pQ->nBufSz;
    }

    // Increment the number of bytes written.
    WordsWritten += nLen;

    // Increment the count.  (protect with mutex)
    if (pQ->wrCntProt){ Platform_EnterOSCritical(); }
    pQ->nCount = pQ->nCount + nLen;
    if (pQ->wrCntProt){ Platform_ExitOSCritical(); }

    PcmQ_UNLOCK_IFDEFINED( pQ->wrMutex );
  }
  return WordsWritten;    
}



/*
**=============================================================================
**  Abstract:
**  Read function that lacks protection from multiple threads trying to read
**  at the same time.  Use this function if only one thread/task will read
**  from the queue at a time.
**  Parameters:
**
**  Return values:
**
**=============================================================================
*/
uint_t PcmQRead( PcmQ_t * const pQ, pcm_t *pRdBuf, uint_t nLen )
{
  uint_t WordsRead   = 0;
  ASSERT( NULL != pQ );

  if ( nLen )
  {
    // Calculate how many bytes can be read from the RdBuffer.
    uint_t WordsToRead = 0;
    const uint_t nBufSz = pQ->nBufSz;
    uint_t nRdIdx = pQ->nRdIdx;
    const pcm_t * const pBuf = pQ->pfBuf;

    PcmQ_LOCKIFDEFINED( pQ->rdMutex );

    // No count MUTEX needed because count is native integer (single cycle write or read)
    // and can only get larger if a process writes while we are reading.
    WordsToRead = MIN(pQ->nCount, nLen);

    // We can definitely read WordsToRead bytes.
    while (WordsToRead > 0)
    {
      // Calculate how many contiguous bytes to the end of the buffer
      uint_t Words = MIN( WordsToRead, (nBufSz - nRdIdx) );

      // Copy that many bytes.
      memcpy( &pRdBuf[ WordsRead ], &pBuf[ nRdIdx ], Words*sizeof(pcm_t) );

      // Circular buffering.
      nRdIdx += Words;
      if (nRdIdx >= nBufSz)
      {
        nRdIdx -= nBufSz;
      }

      // Increment the number of bytes read.
      WordsRead += Words;
      WordsToRead -= Words;
    }

    pQ->nRdIdx = nRdIdx;

    // Decrement the count.
    //PcmQ_LOCKIFDEFINED( pQ->rdCntMutex );
    if (pQ->rdCntProt){ Platform_EnterOSCritical(); }
    pQ->nCount = pQ->nCount - WordsRead;
    //PcmQ_UNLOCK_IFDEFINED( pQ->rdCntMutex );
    if (pQ->rdCntProt){ Platform_ExitOSCritical(); }

    PcmQ_UNLOCK_IFDEFINED( pQ->rdMutex );

  }
  return WordsRead;
}


/*
**=============================================================================
**  Abstract:
**  Read function that lacks protection from multiple threads trying to read
**  at the same time.  Use this function if only one thread/task will read
**  from the queue at a time.
**  Parameters:
**
**  Return values:
**
**=============================================================================
*/
uint_t PcmQCommitRead( PcmQ_t * const pQ, uint_t nLen )
{
  uint_t WordsRead   = 0;
  ASSERT( NULL != pQ );

  if ( nLen )
  {
    PcmQ_LOCKIFDEFINED( pQ->rdMutex );

    // No count MUTEX needed because count is native integer (single cycle write or read)
    // and can only get larger if a process writes while we are reading.
    nLen = MIN(pQ->nCount, nLen);

    pQ->nRdIdx += nLen;
    if (pQ->nRdIdx >= pQ->nBufSz)
    {
      pQ->nRdIdx -= pQ->nBufSz;
    }

    // Increment the number of bytes read.
    WordsRead += nLen;

    // Decrement the count.
    if (pQ->rdCntProt){ Platform_EnterOSCritical(); }
    pQ->nCount = pQ->nCount - nLen;
    if (pQ->rdCntProt){ Platform_ExitOSCritical(); }

    PcmQ_UNLOCK_IFDEFINED( pQ->rdMutex );

  }
  return WordsRead;
}



/*
**=============================================================================
**  Abstract:
**
**  Parameters:
**
**  Return values:
**
**=============================================================================
*/
uint_t PcmQGetWriteReady( PcmQ_t * const pQ )
{
  uint_t rval = 0;
  ASSERT( NULL != pQ );

  if (pQ->wrCntProt){ Platform_EnterOSCritical(); }

  rval = (pQ->nBufSz - pQ->nCount);

  if (pQ->wrCntProt){ Platform_ExitOSCritical(); }

  return rval;
}


/*
**=============================================================================
**  Abstract:
**
**  Parameters:
**
**  Return values:
**
**=============================================================================
*/
uint_t PcmQGetContiguousWriteReady( PcmQ_t * const pQ )
{

  uint_t bytesReady = 0;
  ASSERT( NULL != pQ );

  if (pQ->wrCntProt){ Platform_EnterOSCritical(); }

  bytesReady = (pQ->nBufSz - pQ->nCount);
  bytesReady = MIN( bytesReady, pQ->nBufSz - pQ->nWrIdx );

  if (pQ->wrCntProt){ Platform_ExitOSCritical(); }
  return bytesReady;
}

/*
**=============================================================================
**  Abstract:
**
**  Parameters:
**
**  Return values:
**
**=============================================================================
*/
uint_t PcmQGetReadReady( PcmQ_t * const pQ )
{
  uint_t bytesReady = 0;
  ASSERT( NULL != pQ );

  if (pQ->rdCntProt){ Platform_EnterOSCritical(); }

  bytesReady = pQ->nCount;

  if (pQ->rdCntProt){ Platform_ExitOSCritical(); }

  return bytesReady;
}


/*
**=============================================================================
**  Abstract:
**
**  Parameters:
**
**  Return values:
**
**=============================================================================
*/
uint_t PcmQGetContiguousReadReady( PcmQ_t * const pQ )
{

  uint_t bytesReady = 0;
  ASSERT( NULL != pQ );

  if (pQ->rdCntProt){ Platform_EnterOSCritical(); }
  bytesReady = MIN( pQ->nCount, pQ->nBufSz - pQ->nRdIdx );
  if (pQ->rdCntProt){ Platform_ExitOSCritical(); }

  return bytesReady;
}


/*
**=============================================================================
**  Abstract:
**
**  Parameters:
**
**  Return values:
**
**=============================================================================
*/
void PcmQFlush( PcmQ_t * const pQ )
{
  // Get the read mutex to only allow a single thread to read from the
  // queue at a time.

  ASSERT( NULL != pQ );

  PcmQ_LOCKIFDEFINED( pQ->rdMutex );
  PcmQ_LOCKIFDEFINED( pQ->wrMutex );

  // delete the single count mutex
  if ((pQ->rdCntProt) || (pQ->wrCntProt))
  {
    Platform_EnterOSCritical();
  }

  pQ->nCount = 0;
  pQ->nRdIdx = pQ->nWrIdx = 0;

  if ((pQ->rdCntProt) || (pQ->wrCntProt))
  {
    Platform_ExitOSCritical();
  }

  PcmQ_UNLOCK_IFDEFINED( pQ->wrMutex );
  PcmQ_UNLOCK_IFDEFINED( pQ->rdMutex );

}

/*
**=============================================================================
*  Abstract:
*
**=============================================================================
*/
uint_t PcmQPeek( PcmQ_t * const pQ, pcm_t *pRdBuf, uint_t nLen )
{
  uint_t WordsRead   = 0;

  ASSERT( NULL != pQ );
  if (nLen )
  {

    uint_t nRdIdx;
    uint_t WordsToRead;

    PcmQ_LOCKIFDEFINED( pQ->rdMutex );

    nRdIdx = pQ->nRdIdx;

    // Calculate how many bytes can be read from the RdBuffer.
    WordsToRead = MIN(pQ->nCount, nLen);

    // We can definitely read WordsToRead bytes.
    while (WordsToRead > 0)
    {
      // Calculate how many contiguous bytes to the end of the buffer
      uint_t Words = MIN( WordsToRead, (pQ->nBufSz - nRdIdx) );

      // Copy that many bytes.
      memcpy( &pRdBuf[ WordsRead ], &pQ->pfBuf[ nRdIdx ], Words*sizeof(pcm_t) );

      // Circular buffering.
      nRdIdx += Words;
      if (nRdIdx >= pQ->nBufSz)
      {
        nRdIdx -= pQ->nBufSz;
      }

      // Increment the number of bytes read.
      WordsRead += Words;
      WordsToRead -= Words;
    }

    PcmQ_UNLOCK_IFDEFINED( pQ->rdMutex );
  }
  return WordsRead;
}

/*
**=============================================================================
*  Abstract:
*
**=============================================================================
*/
void * PcmQGetWritePtr( PcmQ_t * const pQ )
{
  void *pRVal = 0;

  if (pQ->wrCntProt){ Platform_EnterOSCritical(); }
  pRVal = (void *)&pQ->pfBuf[ pQ->nWrIdx ];
  if (pQ->wrCntProt){ Platform_ExitOSCritical(); }

  return pRVal;
}

/*
**=============================================================================
*  Abstract:
*
**=============================================================================
*/
void * PcmQGetReadPtr( PcmQ_t * const pQ )
{
  void *pRVal = 0;

  if (pQ->rdCntProt){ Platform_EnterOSCritical(); }
  pRVal = (void *)&pQ->pfBuf[ pQ->nRdIdx ];
  if (pQ->rdCntProt){ Platform_ExitOSCritical(); }

  return pRVal;
}

/*
**=============================================================================
*  Abstract:
*
**=============================================================================
*/
void PcmQSetRdIdxFromPointer( PcmQ_t * const pQ, void *pRdPtr )
{
  if (pQ->rdCntProt){ Platform_EnterOSCritical(); }
  {
    const pcm_t *const pRd8 = (const pcm_t *)pRdPtr;
    //(Lint):Note 946: Relational or subtract operator applied to pointers [MISRA 2004 Rule 17.3]
    int_t newRdIdx = (int_t)(pRd8 - pQ->pfBuf); //lint !e946

    // Check for within range.
    if ((newRdIdx >= 0) && (newRdIdx <= (int_t)pQ->nBufSz)) //lint !e574 !e737
    {
      int_t newCount;

      // If last read advanced pointer to end of buffer, this is OK, just set to beginning.
      if (newRdIdx == (int_t)pQ->nBufSz) //lint !e737
      {
        newRdIdx = 0;
      }

      // New count is amount write is ahead of read.
      newCount = (int_t)pQ->nWrIdx - newRdIdx;

      // Assume we are being called from consumer, so wr==rd results in zero count
      if (newCount < 0)
      {
        //(Lint):Info 737: Loss of sign in promotion from int to unsigned int
        //(Lint):Info 713: Loss of precision (assignment) (unsigned int to int)
        newCount += pQ->nBufSz; //lint !e713 !e737
      }

      // Set read index and count.
      pQ->nRdIdx = (uint_t)newRdIdx;
      pQ->nCount = (uint_t)newCount;
    }
  }
  if (pQ->rdCntProt){ Platform_ExitOSCritical(); }
}

/*
**=============================================================================
*  Abstract:
*
**=============================================================================
*/
uint_t PcmQUnread( PcmQ_t * const pQ, uint_t nLen )
{
  uint_t bytesUnread = 0;
  ASSERT( NULL != pQ );

  if ( nLen )
  {
    // Calculate how many bytes can be read from the RdBuffer.
    uint_t WordsToUnRead = 0;
    int_t nReadIdx = 0;

    PcmQ_LOCKIFDEFINED( pQ->rdMutex );

    // No count MUTEX needed because count is native integer (single cycle write or read)
    // and can only get larger if a process writes while we are reading.
    WordsToUnRead = MIN((pQ->nBufSz - pQ->nCount), nLen);

    // We can definitely read WordsToRead bytes.
    nReadIdx = (int_t)pQ->nRdIdx - WordsToUnRead; //lint !e713 !e737
    if (nReadIdx < 0)
    {
      nReadIdx += pQ->nBufSz; //lint !e713 !e737
    }
    pQ->nRdIdx = (uint_t)nReadIdx;

    // Decrement the count.
    if (pQ->rdCntProt){ Platform_EnterOSCritical(); }
    pQ->nCount = pQ->nCount + WordsToUnRead;
    if (pQ->rdCntProt){ Platform_ExitOSCritical(); }

    PcmQ_UNLOCK_IFDEFINED( pQ->rdMutex );

    bytesUnread = WordsToUnRead;
  }
  return bytesUnread;
}

/*
**=============================================================================
*  Abstract:
*
**=============================================================================
*/
uint_t PcmQForceWrite( PcmQ_t * const pQ, const pcm_t * const pWrBuf, uint_t nLen )
{
  uint_t words = 0;
  ASSERT( NULL != pQ );

  if ( nLen )
  {
    int_t  diff = 0;
    uint_t writeableWords = 0;
    uint_t newWrIdx = 0;

    // Lock both read and write mutexes
    PcmQ_LOCKIFDEFINED( pQ->rdMutex );
    PcmQ_LOCKIFDEFINED( pQ->wrMutex );

    // Calculate the number of bytes that can be written
    writeableWords = pQ->nBufSz - pQ->nCount;
    diff = nLen - writeableWords; //lint !e713 !e737

    // If more bytes should be written than there is space for, 
    // force the read pointer forward
    if ( diff > 0 )
    {
      pQ->nRdIdx += diff; //lint !e713 !e737
      if (pQ->nRdIdx >= pQ->nBufSz) { pQ->nRdIdx -= pQ->nBufSz; }

      if (pQ->rdCntProt){ Platform_EnterOSCritical();}
      pQ->nCount -= diff; //lint !e713 !e737
      if (pQ->rdCntProt){ Platform_ExitOSCritical();}
    }
    else
    {
      diff = 0;
    }

    // Insert the data in the buffer.
    ASSERT_FN( nLen == PcmQUnprotectedInsert( pQ, pWrBuf, nLen, &newWrIdx ) );
    pQ->nWrIdx = newWrIdx;

    if (pQ->wrCntProt){ Platform_EnterOSCritical();}
    pQ->nCount += nLen;
    if (pQ->wrCntProt){ Platform_ExitOSCritical();}

    words = nLen;
    PcmQ_UNLOCK_IFDEFINED( pQ->wrMutex );
    PcmQ_UNLOCK_IFDEFINED( pQ->rdMutex );
  }
  return words;
}

/*
**=============================================================================
*  Abstract:
*
**=============================================================================
*/
uint_t PcmQPeekRandom( PcmQ_t * const pQ, pcm_t *pRdBuf, uint_t bytesFromRdIdx, uint_t nLen )
{
  uint_t WordsRead   = 0;

  ASSERT( NULL != pQ );
  if ( nLen )
  {

    uint_t nRdIdx;
    uint_t WordsToRead;
    int_t nCount;

    PcmQ_LOCKIFDEFINED( pQ->rdMutex );

    nRdIdx = pQ->nRdIdx + bytesFromRdIdx;
    if (nRdIdx >= pQ->nBufSz)
    {
      nRdIdx -= pQ->nBufSz;
    }
    nCount = (pQ->nCount - bytesFromRdIdx); //lint !e713 !e737
    nCount = (nCount < 0) ? 0 : nCount;

    // Calculate how many bytes can be read from the RdBuffer.
    WordsToRead = MIN( ((uint_t )nCount), nLen);

    // We can definitely read WordsToRead bytes.
    while (WordsToRead > 0)
    {
      // Calculate how many contiguous bytes to the end of the buffer
      uint_t Words = MIN( WordsToRead, (pQ->nBufSz - nRdIdx) );

      // Copy that many bytes.
      memcpy( &pRdBuf[ WordsRead ], &pQ->pfBuf[ nRdIdx ], Words*sizeof(pcm_t) );

      // Circular buffering.
      nRdIdx += Words;
      if (nRdIdx >= pQ->nBufSz)
      {
        nRdIdx -= pQ->nBufSz;
      }

      // Increment the number of bytes read.
      WordsRead += Words;
      WordsToRead -= Words;
    }

    PcmQ_UNLOCK_IFDEFINED( pQ->rdMutex );
  }
  return WordsRead;
}

/** [Declaration] Inserts data somewhere into the buffer */
uint_t PcmQPokeRandom( PcmQ_t * const pQ, pcm_t *pWrBuf, uint_t bytesFromStart, uint_t nLen )
{
  uint_t WordsWritten   = 0;

  ASSERT( NULL != pQ );
  if ( nLen )
  {

    uint_t nWrIdx;
    uint_t WordsToWrite;
    int_t nCount;

    PcmQ_LOCKIFDEFINED( pQ->wrMutex );

    nWrIdx = pQ->nWrIdx + bytesFromStart;
    if (nWrIdx >= pQ->nBufSz)
    {
      nWrIdx -= pQ->nBufSz;
    }

    // Can only write if it will fit within nCount
    nCount = pQ->nCount - bytesFromStart; //lint !e713 !e737
    nCount = (nCount < 0) ? 0 : nCount;

    // Calculate how many bytes can be written to the WrBuffer.
    WordsToWrite = MIN( ((uint_t )nCount), nLen);

    // We can definitely read WordsToRead bytes.
    while (WordsToWrite > 0)
    {
      // Calculate how many contiguous bytes to the end of the buffer
      uint_t Words = MIN( WordsToWrite, (pQ->nBufSz - nWrIdx) );

      // Copy that many bytes.
      memcpy( &pQ->pfBuf[ nWrIdx ], &pWrBuf[ WordsWritten ], Words*sizeof(pcm_t) );

      // Circular buffering.
      nWrIdx += Words;
      if (nWrIdx >= pQ->nBufSz)
      {
        nWrIdx -= pQ->nBufSz;
      }

      // Increment the number of bytes read.
      WordsWritten += Words;
      WordsToWrite -= Words;
    }

    PcmQ_UNLOCK_IFDEFINED( pQ->wrMutex );
  }
  return WordsWritten;
}


  /** [Declaration] Reads the last nLen words from the buffer */
int_t PcmQDoReadFromEnd(PcmQ_t * const pQ, pcm_t *pRdBuf, int_t nLen)
{
  int_t shortsRead   = 0;

  if ( nLen > 0 )
  {
    // Calculate how many shorts can be read from the RdBuffer.
    int_t shortsToRead = 0;
    int_t nRdIdx = pQ->nWrIdx - nLen;
    if (nRdIdx < 0) {
      nRdIdx += pQ->nBufSz;
    }

    shortsToRead = nLen;

    // We can definitely read ShortsToRead shorts.
    while (shortsToRead > 0) {
      // Calculate how many contiguous shorts to the end of the buffer
      int_t Shorts = MIN( shortsToRead, (int_t)((pQ->nBufSz - nRdIdx) ));

      // Copy that many shorts.
      memcpy( &pRdBuf[shortsRead], &pQ->pfBuf[nRdIdx], Shorts * sizeof( pcm_t ) );

      // Circular buffering.
      nRdIdx += Shorts;
      if (nRdIdx >= (int_t)pQ->nBufSz)  {
        nRdIdx -= pQ->nBufSz;
      }

      // Increment the number of shorts read.
      shortsRead += Shorts;
      shortsToRead -= Shorts;
    }
  }
  return shortsRead;
}



int_t PcmQDoReadToDoubleFromEnd(PcmQ_t * const pQ, double *pRdBuf, int_t nLen)
{
  int_t shortsRead   = 0;

  if ( nLen > 0 )
  {
    // Calculate how many shorts can be read from the RdBuffer.
    int_t shortsToRead = 0;
    int_t nRdIdx = pQ->nWrIdx - nLen;
    if (nRdIdx < 0)
    {
      nRdIdx += pQ->nBufSz;
    }

    shortsToRead = nLen;

    // We can definitely read ShortsToRead shorts.
    while (shortsToRead > 0)
    {
      int_t i;

      // Calculate how many contiguous shorts to the end of the buffer
      int_t Shorts = MIN( shortsToRead, (int_t)((pQ->nBufSz - nRdIdx) ));

      // Copy that many shorts.
      //HelperMath.aMemCpyToDbl( pRdBuf, shortsRead, pQ->pBuf, nRdIdx, Shorts );
      for( i = 0; i < Shorts; i++)
      {
        pRdBuf[shortsRead + i] = (double)pQ->pfBuf[i + nRdIdx];
      }

      // Circular buffering.
      nRdIdx += Shorts;
      if (nRdIdx >= (int_t)pQ->nBufSz)
      {
        nRdIdx -= pQ->nBufSz;
      }

      // Increment the number of shorts read.
      shortsRead += Shorts;
      shortsToRead -= Shorts;
    }
  }
  return shortsRead;
}

int_t PcmQDoReadToDoubleCustomCommit(
  PcmQ_t * const pQ, double *pRdBuf, uint_t nLen, int_t amountToCommit )
{
  uint_t nToRead   = 0;
  // Calculate how many shorts can be read from the RdBuffer.
  uint_t nRead = 0;

  if ( pQ->nCount >= nLen )
  {

    uint_t rdIdx = pQ->nRdIdx;
    // No count MUTEX needed because count is native integer (single cycle write or read)
    // and can only get larger if a process writes while we are reading.
    nToRead = MIN(pQ->nCount, nLen);

    // We can definitely read ShortsToRead shorts.
    while (nToRead > 0)
    {
      uint_t i;
      // Calculate how many contiguous shorts to the end of the buffer
      const uint_t nShorts = MIN( nToRead, (pQ->nBufSz - rdIdx) );

      // Copy that many shorts.
      //HelperMath.aMemCpyToDbl( pRdBuf, nRead, pQ->pBuf, rdIdx, nShorts );
      for( i = 0; i < nShorts; i++)
      {
        pRdBuf[nRead + i] = (double)pQ->pfBuf[i + rdIdx];
      }

      // Circular buffering.
      rdIdx += nShorts;
      if (rdIdx >= pQ->nBufSz) {
        rdIdx -= pQ->nBufSz;
      }

      // Increment the number of shorts read.
      nRead += nShorts;
      nToRead -= nShorts;
    }
    pQ->nRdIdx += amountToCommit;
    if (pQ->nRdIdx >= pQ->nBufSz) {
      pQ->nRdIdx -= pQ->nBufSz;
    }
    pQ->nRdCount += amountToCommit;
    pQ->nCount -= amountToCommit;
  }
  return nRead;
}
