/**
 * @file
 *
 *   Unless noted otherwise, the portions of Isis written by the USGS are public
 *   domain. See individual third-party library and package descriptions for 
 *   intellectual property information,user agreements, and related information.
 *
 *   Although Isis has been used by the USGS, no warranty, expressed or implied,
 *   is made by the USGS as to the accuracy and functioning of such software 
 *   and related material nor shall the fact of distribution constitute any such 
 *   warranty, and no responsibility is assumed by the USGS in connection 
 *   therewith.
 *
 *   For additional information, launch
 *   $ISISROOT/doc//documents/Disclaimers/Disclaimers.html in a browser or see 
 *   the Privacy &amp; Disclaimers page on the Isis website,
 *   http://isis.astrogeology.usgs.gov, and the USGS privacy and disclaimers on
 *   http://www.usgs.gov/privacy.html.
 */
#include "LutiCamera.h"

#include <iomanip>

#include "CameraFocalPlaneMap.h"
#include "IException.h"
#include "IString.h"
#include "iTime.h"
#include "LineScanCameraDetectorMap.h"
#include "LineScanCameraGroundMap.h"
#include "LineScanCameraSkyMap.h"
#include "LutiDistortionMap.h"
#include "NaifStatus.h"

using namespace std;
namespace Isis {
  /**
   * Constructor for the Luti Camera Model
   *
   * @param lab Pvl Label to create the camera model from
   *
   * @internal 
   *   @history 2020-02-10 Mijung Kim, Create File
   */
  LutiCamera::LutiCamera(Cube &cube) : LineScanCamera(cube) {
    m_spacecraftNameLong = "Korea Pathfinder Lunar Orbiter";
    m_spacecraftNameShort = "KPLO";
    // LUTI instrument kernel code = -155000
    if (naifIkCode() == -155000 || naifIkCode() == -155100) {
      m_instrumentNameLong = "LUnar Terrain Imager";
      m_instrumentNameShort = "LUTI";
    }
    else if(naifIkCode() == -155101)
    {
      m_instrumentNameLong = "LUnar Terrain Imager";
      m_instrumentNameShort = "LUTI1";
    }
    else if(naifIkCode() == -155102)
    {
      m_instrumentNameLong = "LUnar Terrain Imager";
      m_instrumentNameShort = "LUTI2";
    }
    else {
      QString msg = "File does not appear to be a Lunar Terrain Image: ";
      msg += QString::number(naifIkCode());
      msg += " is not a supported instrument kernel code for Lunar Terrain Imager.";
      throw IException(IException::Programmer, msg, _FILEINFO_);
    }
    NaifStatus::CheckErrors();

    // Set up the camera info from ik/iak kernels
    SetFocalLength();
    SetPixelPitch();

    double constantTimeOffset = 0.0,
           additionalPreroll = 0.0,
           additiveLineTimeError = 0.0,
           multiplicativeLineTimeError = 0.0;

    QString ikernKey = "INS" + toString(naifIkCode()) + "_CONSTANT_TIME_OFFSET";
    constantTimeOffset = getDouble(ikernKey);

    ikernKey = "INS" + toString(naifIkCode()) + "_ADDITIONAL_PREROLL";
    additionalPreroll = getDouble(ikernKey);

    ikernKey = "INS" + toString(naifIkCode()) + "_ADDITIVE_LINE_ERROR";
    additiveLineTimeError = getDouble(ikernKey);

    ikernKey = "INS" + toString(naifIkCode()) + "_MULTIPLI_LINE_ERROR";
    multiplicativeLineTimeError = getDouble(ikernKey);

    // Get the start time from labels
    Pvl &lab = *cube.label();
    PvlGroup &inst = lab.findGroup("Instrument", Pvl::Traverse);
    QString stime = inst["SpacecraftClockPrerollCount"];
    SpiceDouble etStart;

    if(stime != "NULL") {
      etStart = getClockTime(stime).Et();
    }
    else {
      etStart = iTime((QString)inst["PrerollTime"]).Et();
    }
    
    // Get other info from labels
    double csum = inst["SpatialSumming"];
    double lineRate = (double) inst["LineExposureDuration"] / 1000.0;
    double ss = inst["SampleFirstPixel"];
    ss += 1.0;

    lineRate *= 1.0 + multiplicativeLineTimeError;
    lineRate += additiveLineTimeError;
    etStart += additionalPreroll * lineRate;
    etStart += constantTimeOffset;

    setTime(etStart);

    // Setup detector map
    LineScanCameraDetectorMap *detectorMap = new LineScanCameraDetectorMap(this, etStart, lineRate);
    detectorMap->SetDetectorSampleSumming(csum);
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
    LutiDistortionMap *distMap = new LutiDistortionMap(this);
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
 * Luti object.
 *
 * @param lab Cube labels
 *
 * @return Isis::Camera* LutiCamera
 * @internal
 *   @history 2020-02-10 Mijung Kim, Create File
 */
extern "C" Isis::Camera *LutiCameraPlugin(Isis::Cube &cube) {
  return new Isis::LutiCamera(cube);
}
