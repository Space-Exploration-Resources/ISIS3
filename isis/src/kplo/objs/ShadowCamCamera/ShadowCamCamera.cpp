/** This is free and unencumbered software released into the public domain.

The authors of ISIS do not claim copyright on the contents of this file.
For more details about the LICENSE terms and the AUTHORS, you will
find files of those names at the top level of this repository. **/

/* SPDX-License-Identifier: CC0-1.0 */

#include "ShadowCamCamera.h"

#include <iomanip>

#include "CameraFocalPlaneMap.h"
#include "IException.h"
#include "IString.h"
#include "iTime.h"
#include "LineScanCameraDetectorMap.h"
#include "LineScanCameraGroundMap.h"
#include "LineScanCameraSkyMap.h"
#include "ShadowCamDistortionMap.h"
#include "NaifStatus.h"

using namespace std;
namespace Isis {
  /**
   * Constructor for the LRO NAC Camera Model
   *
   * @param lab Pvl Label to create the camera model from
   *
   * @internal
   *   @history 2022-10-12 Victor Silva - original object
   */
  ShadowCamCamera::ShadowCamCamera(Cube &cube) : LineScanCamera(cube) {
    m_spacecraftNameLong = "KOREA PATHFINDER LUNAR ORBITER";
    m_spacecraftNameShort = "KPLO";
    // SHADOWCAM instrument kernel code = -155151
    if (naifIkCode() == -155151) {
      m_instrumentNameLong = "KOREA PATHFINDER LUNAR ORBITER SHADOWCAM";
      m_instrumentNameShort = "KPLO SHADOWCAM";
    }
    
    else {
      QString msg = "File does not appear to be a Korea Pathfinder Lunar Orbiter ShadowCam Image: ";
      msg += QString::number(naifIkCode());
      msg += " is not a supported instrument kernel code for Korea Pathfinder Lunar Orbiter.";
      throw IException(IException::Programmer, msg, _FILEINFO_);
    }
    NaifStatus::CheckErrors();

    // Set up the camera info from ik/iak kernels
    SetFocalLength();
    SetPixelPitch();

    // ToDo... make iak to hold these four values with preroll set to 0
    double constantTimeOffset = 0.0,
           //additionalPreroll = 0.0,
           additiveLineTimeError = 0.0,
           multiplicativeLineTimeError = 0.0;

    QString ikernKey = "INS" + toString(naifIkCode()) + "_CONSTANT_TIME_OFFSET";
    constantTimeOffset = getDouble(ikernKey);
    
    //ikernKey = "INS" + toString(naifIkCode()) + "_ADDITIONAL_PREROLL";
    //additionalPreroll = getDouble(ikernKey);

    ikernKey = "INS" + toString(naifIkCode()) + "_ADDITIVE_LINE_ERROR";
    additiveLineTimeError = getDouble(ikernKey);

    ikernKey = "INS" + toString(naifIkCode()) + "_MULTIPLI_LINE_ERROR";
    multiplicativeLineTimeError = getDouble(ikernKey);

    // Get the start time from labels
    // ToDo...get the preRoll Time
    Pvl &lab = *cube.label();
    PvlGroup &inst = lab.findGroup("Instrument", Pvl::Traverse);

    SpiceDouble etStart;
    etStart = iTime((QString)inst["PrerollTime"]).Et();

    double preroll;
    preroll = (double) inst["PrerollLines"];
    

    // Get other info from labels
    // Shadwowcam linerate info in milliseconds, make in secs
    double lineRate = (double) inst["LineRate"] / 1000.0;
    // this value is usually 0 then set to one in two lines below
    double ss = inst["SampleFirstPixel"];
    ss += 1.0;
    
    // TDI direction offset from IAK
    /*
      INS-155151_TDI_A_Offset 
      INS-155151_TDI_B_Offset 
    */
    double tdiOffset;
    if ((QString)inst["TDIDirection"] == "A") {
      ikernKey = "INS" + toString(naifIkCode()) + "_TDI_A_OFFSET";
      tdiOffset = getDouble(ikernKey);
    }
    else{
      ikernKey = "INS" + toString(naifIkCode()) + "_TDI_B_OFFSET";
      tdiOffset = getDouble(ikernKey);
    }

    lineRate *= 1.0 + multiplicativeLineTimeError;
    lineRate += additiveLineTimeError;
    etStart += preroll * lineRate;
    etStart += constantTimeOffset;
    etStart += tdiOffset * lineRate;

    setTime(etStart);

    // Setup detector map
    // ShadowCam uses line_rate_ms
    LineScanCameraDetectorMap *detectorMap = new LineScanCameraDetectorMap(this, etStart, lineRate);
    detectorMap->SetStartingDetectorSample(ss);

    // Setup focal plane map
    CameraFocalPlaneMap *focalMap = new CameraFocalPlaneMap(this, naifIkCode());

    //  Retrieve boresight location from instrument kernel (IK) (addendum?)
    ikernKey = "INS" + toString(naifIkCode()) + "_BORESIGHT_SAMPLE";
    double sampleBoreSight = getDouble(ikernKey);

    ikernKey = "INS" + toString(naifIkCode()) + "_BORESIGHT_LINE";
    double lineBoreSight = getDouble(ikernKey);

    focalMap->SetDetectorOrigin(sampleBoreSight, lineBoreSight);
    focalMap->SetDetectorOffset(0.0, 0.0);

    // Setup distortion map
    ShadowCamDistortionMap *distMap = new ShadowCamDistortionMap(this);
    distMap->SetDistortion(naifIkCode());

    // Setup the ground and sky map
    new LineScanCameraGroundMap(this);
    new LineScanCameraSkyMap(this);

    LoadCache();
    NaifStatus::CheckErrors();
  }
}

/**
 * This is the function that is called in order to instantiate a
 * KPLO ShadowCam object.
 *
 * @param lab Cube labels
 *
 * @return Isis::Camera* KploShadowCam Camera
 * @internal
 *    @history 2022-10-12 Victor Silva - original object
 */
extern "C" Isis::Camera *ShadowCamCameraPlugin(Isis::Cube &cube) {
  return new Isis::ShadowCamCamera(cube);
}
