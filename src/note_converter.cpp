/******************************************************************************
  Copyright 2014 Chris Fogelklou, Applaud Apps (Applicaudia)

  This source code may under NO circumstances be distributed without the 
  express written permission of the author.

  @author: chris.fogelklou@gmail.com
*******************************************************************************/
#include "note_converter.hpp"
#include "audutils_defines.h"

double NoteConverter::DEFAULT_REF_FREQ = 440;
const double NoteConverter::mOneOverLn2 = 1.0 / LN2;
