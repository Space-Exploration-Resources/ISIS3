#include "Isis.h"
#include "Application.h"
#include "shadowcamcal.h"

using namespace Isis;

/**
  *
  * @brief Shadowcamcal
  *
  * Performs radiometric corrections to images acquired by the ShadowCam
  *    Camera aboard the Korea Pathfinder Lunar Orbiter spacecraft.
  *
  * @author 2016-09-16 Victor Silva
  *
  * @internal
  *   @history 2022-09-21 Victor Silva - Calibration application for ShadowCam
  */
void IsisMain (){
  UserInterface &ui = Application::GetUserInterface();
  shadowcamcal(ui);
}
