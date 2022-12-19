#ifndef ShadowCamCamera_h
#define ShadowCamCamera_h

/** This is free and unencumbered software released into the public domain.

The authors of ISIS do not claim copyright on the contents of this file.
For more details about the LICENSE terms and the AUTHORS, you will
find files of those names at the top level of this repository. **/

/* SPDX-License-Identifier: CC0-1.0 */

#include "LineScanCamera.h"

namespace Isis {
  /**
   * @brief KPLO ShadowCam Camera Model
   *
   * This is the camera model for the Korea Pathfinder Lunar Orbiter ShadowCam
   * camera.
   *
   * @ingroup SpiceInstrumentsAndCameras
   * @ingroup KoreaPathfinderLunarOrbiter
   *
   * @author 2022-10-12 Victor Silva
   *
   * @internal
   *   @history 2022-10-12 Victor Silva, Original Object
   */

  class ShadowCamCamera : public LineScanCamera {
    public:
      ShadowCamCamera(Cube &cube);

      //! Destroys the ShadowCamCamera object
      ~ShadowCamCamera() {};

      /**
       * CK frame ID -  - Instrument Code from spacit run on CK
       *
       * @return @b int The appropriate instrument code for the "Camera-matrix"
       *         Kernel Frame ID
       */
      virtual int CkFrameId() const { return (-155151); }

      /**
       * CK Reference ID - J2000
       *
       * @return @b int The appropriate instrument code for the "Camera-matrix"
       *         Kernel Reference ID
       */
      virtual int CkReferenceId() const { return (1); }

      /**
       *  SPK Reference ID - J2000
       *
       * @return @b int The appropriate instrument code for the Spacecraft
       *         Kernel Reference ID
       */
      virtual int SpkReferenceId() const { return (1); }
  };
};
#endif
