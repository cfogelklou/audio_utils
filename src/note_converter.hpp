#ifndef NOTE_CONVERTER_HPP_
#define NOTE_CONVERTER_HPP_

/******************************************************************************
  Copyright 2014 Chris Fogelklou, Applaud Apps (Applicaudia)

  This source code may under NO circumstances be distributed without the 
  express written permission of the author.

  @author: chris.fogelklou@gmail.com
*******************************************************************************/

#include <math.h>
#include "audiolib_types.h"
#include "DaTunerDefs.h"
#include "DaTunerDebug.h"

class NoteConverter {

private:
  double mReferenceFrequency;
  double mFreq;
  int mNearestNoteOrdinal;
  double mFreqError;
  static const double mOneOverLn2;
  double mErrCents;
  int mOctave;

public:
  enum notesT
  {
    A, // 0
    Ash, // 1
    B, // 2
    C, // 3
    Csh, // 4
    D, // 5
    Dsh, // 6
    E, // 7
    F, // 8
    Fsh, // 9
    G, // 10
    Gsh // 11
  } notes_t;

  static double DEFAULT_REF_FREQ;// = 440;

  NoteConverter( )
  {
    setReferenceFrequency( DEFAULT_REF_FREQ );
    mFreq = DEFAULT_REF_FREQ;
    mNearestNoteOrdinal = 0;
    mFreqError = 0.0;
    mErrCents = 0.0;
    mOctave = 0;
    setFrequency( DEFAULT_REF_FREQ );
  }

  NoteConverter( double frequency )
  {
    setReferenceFrequency( DEFAULT_REF_FREQ );
    mFreq = DEFAULT_REF_FREQ;
    mNearestNoteOrdinal = 0;
    mFreqError = 0.0;
    mErrCents = 0.0;
    mOctave = 0;
    setFrequency( frequency );
  }

  NoteConverter( double referenceFreq, double frequency )
  {
    setReferenceFrequency( referenceFreq );
    mFreq = referenceFreq;
    mNearestNoteOrdinal = 0;
    mFreqError = 0.0;
    mErrCents = 0.0;
    mOctave = 0;
    setFrequency( frequency );
  }

  virtual ~NoteConverter()
  {
    TRACE_VERBOSE(("Destructor: NoteConverter();"));
  }

  static void indexToOctaveAndNote( int index, int *pOctave, int *pNote )
  {
    int octave = (int)index / 12;
    int note = (int)index % 12;
    if (note < 0) {
      note += 12;
      octave--;
    }
    if (pOctave) {
      *pOctave = octave;
    }
    if (pNote) {
      *pNote = note;
    }
  }

  static double getErrorCentsFromFreqAndRef( const double frequencyHz, const double refFrequencyHz )
  {
    if ((frequencyHz <= 0.0) || (refFrequencyHz <= 0)) return 0;
    const double noteRatio = frequencyHz / refFrequencyHz;
    if (noteRatio < 0.0) return 0;
    return 1200.0 * log(noteRatio) * mOneOverLn2;
  }

  void setFrequency( double frequencyHz )
  {
    if (frequencyHz > 0.0)
    {
      mFreq = frequencyHz;
      const double n = 12.0 * log(frequencyHz/mReferenceFrequency) * mOneOverLn2;
      const double nearestN = ROUND(n);
      mOctave = (int)nearestN / 12;
      int nNoteIdx = (int)nearestN % 12;
      if (nNoteIdx < 0) {
        nNoteIdx += 12;
        mOctave--;
      }
      mNearestNoteOrdinal = nNoteIdx;

      // Calculate frequency error in Hz and cents.
      {
        const double noteFreq = mReferenceFrequency * pow(2.0, nearestN / 12.0);
        const double noteRatio = mFreq / noteFreq;
        mFreqError = mFreq - noteFreq;

        if (noteRatio < 0.0) {
          // Prevent divide by or log of zero.
          mErrCents = 0.0;
        }
        else {
          const double oldErrCents = 1200.0 * log(noteRatio) * mOneOverLn2;
          mErrCents = oldErrCents;
        }
      }
    }
  }

  void setReferenceFrequency( double freq )
  {
    mReferenceFrequency = freq * pow(2.0, -4);
  }

  double getReferenceFrequency(  )
  {
    return mReferenceFrequency * pow(2.0, 4);   
  }	

  void setNote (notesT note, int octave)
  {
    setFrequency( getPerfectFrequency( note, octave ) );
  }

  void setNote (int note, int octave)
  {
    setFrequency( getPerfectFrequency( note, octave ) );
  }

  //private static const double errcentsToFreqCoeff = (Math.log(5)+ Math.log(2)) * Math.log(2);
  double getFreqFromNoteAndError (int note, int octave, double errCents )
  {
    double rval;
    //const double perfect = getPerfectFrequency( note, octave );
    //rval = perfect * Math.exp(errCents * errcentsToFreqCoeff / (12*100));
    const double noteIdx = octave*12 + note;
    rval = mReferenceFrequency * pow(2.0, (noteIdx*100.0 + errCents)/1200.0);

    return rval;
  } 

    //private static const double errcentsToFreqCoeff = (Math.log(5)+ Math.log(2)) * Math.log(2);
  static double getFreqFromRefAndNoteAndError (double refFreq, int note, int octave, double errCents )
  {
    double rval;
    //const double perfect = getPerfectFrequency( note, octave );
    //rval = perfect * Math.exp(errCents * errcentsToFreqCoeff / (12*100));
    const double noteIdx = octave*12 + note;
    rval = refFreq * pow(2.0, (noteIdx*100.0 + errCents)/1200.0);

    return rval;
  }

  double getPerfectFrequency (int note, int octave)
  {
    const double noteIdx = octave*12 + note;
    const double freq = mReferenceFrequency * pow(2.0, noteIdx/12.0);
    return freq;
  }

  static double getPerfectFrequency ( double refFreq, int note, int octave)
  {
    const double noteIdx = octave*12 + note;
    const double freq = refFreq * pow(2.0, noteIdx/12.0);
    return freq;
  }	

  static double getPerfectFrequency ( double refFreq, int noteIdx )
  {
    const double freq = refFreq * pow(2.0, noteIdx/12.0);
    return freq;
  }	

  double getFrequency( notesT note, int octave  )
  {
    setNote (note, octave);
    return mFreq;
  }

  double getFrequency( )
  {
    return mFreq;
  }

  notesT getNearestNote()
  {
    return (notesT)(mNearestNoteOrdinal);
  }

  int getNearestNoteInt()
  {
    return mNearestNoteOrdinal;
  }	

  int getOctave()
  {
    return mOctave;
  }

  double getFreqError()
  {
    return mFreqError;
  }

  double getErrCents()
  {
    return mErrCents;
  }

  static const char * NoteToString(const int n) {
    const char *rval = "";
    switch (n) {
    case A: rval = "A";  break;
    case Ash: rval = "A#"; break;
    case B: rval = "B"; break;
    case C: rval = "C"; break;
    case Csh: rval = "C#"; break;
    case D: rval = "D"; break;
    case Dsh: rval = "D#"; break;
    case E: rval = "E"; break;
    case F: rval = "F"; break;
    case Fsh: rval = "F#"; break;
    case G: rval = "G"; break;
    case Gsh: rval = "G#"; break;
    default:
      break;
    }
    return rval;
  }

};


#endif
