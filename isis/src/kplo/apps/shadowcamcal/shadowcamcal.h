#ifndef shadowcamcal_h
#define shadowcamcal_h

#include "Cube.h"
#include "shadowcamcal.h"
#include "UserInterface.h"

namespace Isis{
  extern void shadowcamcal(UserInterface &ui);
  extern void shadowcamcal(Cube *iCube, UserInterface &ui);
}

#endif