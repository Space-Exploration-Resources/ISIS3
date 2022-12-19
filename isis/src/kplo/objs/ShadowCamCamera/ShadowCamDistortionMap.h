#ifndef ShadowCamDistortionMap_h
#define ShadowCamDistortionMap_h

/** This is free and unencumbered software released into the public domain.

The authors of ISIS do not claim copyright on the contents of this file.
For more details about the LICENSE terms and the AUTHORS, you will
find files of those names at the top level of this repository. **/

/* SPDX-License-Identifier: CC0-1.0 */

#include <vector>
#include "CameraDistortionMap.h"

namespace Isis {

  /**
   *  Distort/undistort focal plane coordinates
   *
   * Creates a map for adding/removing optical distortions
   * from the focal plane of a camera.
   *
   * @ingroup SpiceInstrumentsAndCameras
   * @ingroup KoreaPathfinderLunarOrbiter
   *
   * @see ShadowCamCamera
   *
   * @author 2022-10-12 Victor Silva
   * @internal
   */
  class ShadowCamDistortionMap : public CameraDistortionMap {
    public:
      ShadowCamDistortionMap(Camera *parent);

      //! Destroys the ShadowCamDistortionMap object.
      virtual ~ShadowCamDistortionMap() {};

      void SetDistortion(const int naifIkCode);

      virtual bool SetFocalPlane(const double dx, const double dy);

      virtual bool SetUndistortedFocalPlane(const double ux, const double uy);

  };
};
#endif
