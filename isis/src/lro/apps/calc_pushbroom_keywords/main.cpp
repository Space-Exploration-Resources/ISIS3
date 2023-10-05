#include "Isis.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <math.h>

#include "Pvl.h"
#include "Spice.h"
#include "IString.h"
#include "iTime.h"
#include "Constants.h"
#include "TProjection.h"
#include "Process.h"
#include "FileName.h"
#include "Camera.h"
#include "CameraDetectorMap.h"
#include "CameraDistortionMap.h"
#include "CameraFocalPlaneMap.h"
#include "CameraGroundMap.h"
#include "IException.h"
#include "NaifStatus.h"

using namespace std;
using namespace Isis;

void Load(Isis::PvlKeyword &key);
void Unload_kernels();
std::vector<std::string> p_kernels;

/******************************************************************************
*                                                                             *
*_Title	CALC_PUSHBROOM_KEYWORDS calculates SS generic pushbroom keywords      *
*                                                                             *
*_Desc	This program calculates keyword values specific to SOCET SET's        *
*       Generic Pushbroom sensor model and generates a file of the resulting  *
*       keywords and values for input to SS program 'import_pushbroom'        *
*                                                                             *
*       Missions/Cameras Supported by this application:                       *
*             MOC-NA                                                          *
*             MOC-WA                                                          *
*             MRO-HiRISE                                                      *
*             MRO-CTX                                                         *
*             LRO-NAC                                                         *
*                                                                             *
*_Hist	Oct 23 2008 Elpitha H. Kraus, USGS, Flagstaff Original Version        *
*                           (Note, this is a modification of frame2pushbroom) *
*       Jun 29 2009 EHK - Added LRO NAC camera                                *
*       Aug 24 2009 EHK - For LRO NAC, added tests to error and exit program  *
*                         if boresight sample and/or image number of samples  *
*                         been redefined w.r.t instrument kernel file         *
*                         lro_instruments_v08.ti                              *
*       Oct 12 2009 EHK - Changed USGSAstro_GENERIC_PUSHBROOM to              *
*                         GENERIC_PUSHBROOM for SSv5.5+, and added check for  *
*                         target name in all caps                             *
*       Nov 9  2009 EHK - Added lens distortion coefficients for LRO NAC, and *
*                         corrected the sign of the distortion coefficients   *
*                         for MRO-CTX.                                        *
*                         (Note: for Socet Set we need the negative of the    * 
*                          lens distortion coefficients used in ISIS.)        *
*       Nov 10 2009 EHK - calc_pushbroom_keywords was not getting passed the  *
*                         the camera definition checks for summed image, so   *
*                         made the needed modifications up to csum = 2.       *
*       Nov 12 2009 EHK - Switched Socet Set's along_scan_px_size to be in the*
*                         sample direction, and the cross_scan_pixel_size     *
*                         to be in the line direction.  (LRO images are the   *
*                         first we've encountered that applies summing modes  *
*                         in one direction, and before the switch, SS could   *
*                         not compute partials when initializing the sensor.) *
*       Nov 13 2009 EHK - Switching pixel sizes enabled the sensor model, but *
*                         test LRO images were not coming together properly.  *
*                         Applying the summation factor in both along scan    *
*                         and cross scan directions brought the test images   *
*                         in, so for now, keep the pixels square until I      *
*                         understand how the summation modes are applied in   *
*                         LRO.                                                *
*                         Also need to handle image aquistion when flying the *
*                         LRO s/c backwards!!!                                *
*       Jan 27 2010 EHK - For LRO NAC check for backwards aquired images and  *
*                         set a flip flag in a PVL file to signal that the    *
*                         image needs flipping in the lronac4socet.pl script  *
*                         NOTE: code originally in get_mounting_angles that   *
*                         handles images with kappa of ~180 degrees has been  *
*                         moved to this program in order to determine         *
*                         backwards aquired images                            *
*       Apr 06 2010 EHK - Calculate across_scan and along_scan pixel size     *
*                         based on s/c movement rather than nominal pixel     *
*                         sizes.  This is particularly needed for LRO because *
*                         pixels are not necessarily square for summation=1,  *
*                         and are definitely not square for summation=2.      *
*                         Note that this modification will only work under    *
*                         the LRO isis beta release for now.                  *
*       Aug 24 2010 EHK - On 8/22/2010, updated versions of the LRO ik kernel * 
*                         NAC camera model was released in the LRO beta       *
*                         version of ISIS.  Changes pertinent to this program *
*                         involve a sign change to the lens distortion        *
*                         coefficient in the ik kernel.  However, I hardcoded *
*                         the lens distortion coefficents in this program, so *
*                         no code changes were needed here.  I just updated   *
*                         some comments and added this note to the history.   *
*      Jan 25 2016 VHS -  Updated signs of mounting angles omega, phi         *
*                         to reflect when spacecraft is flying backwards.     *
*      Aug 1 2019 CJL  -  Doubled sampling density of ephemera by halving the *
*                         value of dt_ephem.                                  *
*      Jan 02 2020 CJL -  Changed ephemera to 0.25.                           *
*      Aug 27 2021 CHM -  Re-added to run in the 3.10.2_lroc version of ISIS. *
*_End                                                                         *
*                                                                             *
*******************************************************************************/

int get_naif_matrices (Camera* cam, double et, double ck_et,
                       int sccode, int instcode, int targcode,
                       double BJ[3][3], double PJ[3][3]);

void get_mounting_angles (double *CP, double *BJ, double *PJ,
                          double *pos_USR, double *vel_USR,
                          double *mounting_angles);
void IsisMain() {

#define LINELENGTH 500

  ofstream	pushKeyStrm;	         // pushKey output stream

  string        flipPVLFile;             // Output PVL file containing an
                                         // image flip flag.  The file name
                                         // is built by appending "_flip.pvl"
                                         // to the LRO ProductID (i.e. image
                                         // name)
  ofstream	flipFlagStrm;	         // Flip PVL file output stream
  string        flipFlag;                // Y/N flag signaling if image
                                         // needs to be flipped (x-mirrored)

  double        detectorBoresightLine;   // Boresight line in detector space
  double        detectorBoresightSample; // Boresight sample in detector space
  double        boresightSample;         // Boresight sample in image space
                                         // taking summing modes and 
                                         // starting image sample into account
  bool		isMocNA = false;
  bool		isMocWARed = false;
  bool    isHiRise = false;
  bool    isCTX = false;
  bool    isLroNACL = false;
  bool    isLroNACR = false;
  bool    isShadowCam = false;
  bool    isLUTI = false;


  int		csum;
  int		dsum;
  string        instrumentId;
  string        filter;
  string	      target;
  string        target_frame;
  string        productId;
  string        msg;

	// The following are SOCET SET Generic Pushbroom Keywords:

  double	rectification_terms[6];
  double	focal;
  double	atmco[4];
  double	lensco[10];
  double	iocoef_line[10];
  double	iocoef_samp[10];
  double	tri_params[18];
  double	dt_ephem;
  double	t0_ephem;
  int		  num_ephem;
  double	*ephem_pts;
  double	*ephem_rates;
  double	scan_duration;
  double	int_time;
  double	along_scan_px_size;
  double	cross_scan_px_size;
  double	sun_az;
  double	sun_el;
  double	center_gp[3];
  double	sensor_position[3];
  double	mounting_angles[3];

  int		total_lines;
  int		total_samps=0;

  double	et;
  double	ck_et_center;          /* Ephemeris time with the s/c     */
                                       /* time base added to it of the    */
                                       /* center ephemeris point          */
                                       /* This is the et needed for the   */
				       /* NAIF Camera Kernals             */
  double        ck_time_bias;
  double	roll,pitch,yaw;
  double        CP[3][3];
  double	BJ[3][3];
  double	PJ[3][3];

  double	lt;

  double	center_line_coord;     /* Line coordinate at the center   */
                                       /* of the image */
  double	et_center;             /* Ephemeris seconds past J2000 at */
                                       /* the center of the image         */
  double	et_first_ephem;	       /* Ephemeris seconds past J2000 of */
				       /* the first ephemeris point       */
  double	statevector[6];        /* NAIF statevector in body-fixed  */
				       /* coordinates                     */
  double	radius;	               /* Planet Radius                   */
  double	flattening;	       /* Flattening coeficient           */
  double	center_pos[3];	       /* Body Fixed coordinates of the   */
                                       /* center ephemeris point          */
  double	center_vel[3];	       /* Body Fixed velocity coordinates */
                                       /* of the center ephemeris point   */
  double	lat,lon,height;	       /* Geographic coordinates of a pt  */

  int		i,j;

  double 	deg2rad = Isis::PI / 180.0;

  double expected_boresight;
  int expected_samples;

  // Use a regular Process
  Process p;

  // Open the input cube and initialize the camera
  UserInterface &ui = Application::GetUserInterface();
  Cube *input = p.SetInputCube("FROM");
  Camera *cam = input->camera();
  CameraDetectorMap *detectorMap = cam->DetectorMap();
  CameraFocalPlaneMap *focalMap = cam->FocalPlaneMap();
  CameraDistortionMap *distortionMap = cam->DistortionMap();
  CameraGroundMap *groundMap = cam->GroundMap();

//////////////////////////////////////////////////////////////////////////////
// Get user parameters and error check 
//////////////////////////////////////////////////////////////////////////////
  QString from = ui.GetCubeName("FROM");
  QString pushKey = ui.GetFileName("PUSHKEY");

  // Open input cube and get camera model for it
  Cube cube;
  cube.open(from);

  // Make sure this is a lev1 image (ie, not map projected)
  Projection *proj = 0;
  try {
    proj = input->projection();
  }
  catch (IException &e) {
    proj = 0;
    //e.Clear();
  }
  if (proj != 0) {
    msg = "This is a map projected cube ... not a level 1 image";
    throw IException(IException::User,msg,_FILEINFO_);
  }

  // Open output pushbroom keyword file
  pushKeyStrm.open (FileName(pushKey).expanded().toStdString().c_str(),ios::trunc);
  if (pushKeyStrm.bad()) {
    msg = "Unable to open output PUSHKEY file";
    throw IException(IException::User,msg,_FILEINFO_);
  }

//////////////////////////////////////////////////////////////////////////////
// Get required keywords from instrument and band groups
//////////////////////////////////////////////////////////////////////////////

   PvlGroup inst = cube.label()->findGroup("Instrument", Pvl::Traverse);
   instrumentId = inst["InstrumentId"][0].toStdString();

   target = inst["TargetName"][0].toStdString();

   if (strstr(instrumentId.c_str(),"MOC") != NULL) {
     PvlGroup band = cube.label()->findGroup("BandBin", Pvl::Traverse);
     filter = band["FilterName"][0].toStdString();

     if (strcmp(filter.c_str(),"BROAD_BAND")==0)
       isMocNA = true;
     else if (strcmp(filter.c_str(),"RED")==0)
       isMocWARed = true;
     else if (strcmp(filter.c_str(),"BLUE")==0) {
       msg = "MOC WA Blue filter images not supported for Socet Set mapping";
       throw IException(IException::User,msg,_FILEINFO_);
     }
   }
   else if (strstr(instrumentId.c_str(),"HIRISE") != NULL) isHiRise = true;
   else if (strstr(instrumentId.c_str(),"CTX") != NULL)   isCTX = true;
   else if (strstr(instrumentId.c_str(),"NACL") != NULL ) isLroNACL = true;
   else if (strstr(instrumentId.c_str(),"NACR") != NULL ) isLroNACR = true;
   else if (strstr(instrumentId.c_str(), "ShadowCam") != NULL) isShadowCam = true;
   else if (strstr(instrumentId.c_str(), "LUTI") != NULL) isLUTI = true;  // FIXME: Might need to do LUTIA vs LUTIB
   else {
     msg = "Unsupported instrument: " + instrumentId;
     throw IException(IException::User,msg,_FILEINFO_);
   }

///////////////////////////////////////////////////////////////////////////////
// Now load NAIF kernels and get NAIF codes for this image
//  Modified  to load planetary ephemeris SPKs before s/c SPKs since some
//  missions (e.g., MESSENGER) may augment the s/c SPK with new planet
//  ephemerides. (2008-02-27 (KJB))
///////////////////////////////////////////////////////////////////////////////

  PvlGroup kernels = cube.label()->findGroup("Kernels", Pvl::Traverse);
  Load(kernels["TargetPosition"]);

  if (kernels.hasKeyword("SpacecraftPosition"))
    Load(kernels["SpacecraftPosition"]);
  else
    Load(kernels["InstrumentPosition"]);

  if (kernels.hasKeyword("SpacecraftPointing"))
    Load(kernels["SpacecraftPointing"]);
  else
    Load(kernels["InstrumentPointing"]);
  

  if (kernels.hasKeyword("Frame"))
    Load(kernels["Frame"]);

  if (kernels.hasKeyword("Extra"))
    Load(kernels["Extra"]);

  Load(kernels["TargetAttitudeShape"]);
  Load(kernels["Instrument"]);
  Load(kernels["InstrumentAddendum"]);  // Always load after instrument
  Load(kernels["LeapSecond"]);
  Load(kernels["SpacecraftClock"]);

  // Get NAIF ik, spk, sclk, and ck codes
  //
  //    Use ikcode to get parameters from instrument kernel such as focal
  //    length, distortions, focal plane maps, etc
  //
  //    Use spkCode to get spacecraft position from spk file
  //
  //    Use bodyCode to obtain radii and attitude (pole position/omega0)
  //
  //    Use spkbodycode to read body position from spk

  QString tryKey = "NaifIkCode";
  if (kernels.hasKeyword("NaifFrameCode")) tryKey = "NaifFrameCode";
  int ikCode = (int) kernels[tryKey];

  int spkCode = ikCode / 1000;

  SpiceInt bodyCode;
  SpiceBoolean found;
  bodn2c_c (target.c_str(), &bodyCode, &found);

  if (target == "Mars" || target == "MARS")
     target_frame = "IAU_MARS";
  else if (target == "Moon" || target == "MOON")
     target_frame = "IAU_MOON";
  else {
     msg = target + " is not supported";
     throw Isis::IException(Isis::IException::User,msg,_FILEINFO_);
  }
     
//////////////////////////////////////////////////////////////////////////////
// Get Focal Length.
// NOTE:
//   For MOC Wide Angle, cam->focal_length returns the focal length
//      in pixels, so we must convert from pixels to mm using the PIXEL_SIZE
//      of 0.007 mm gotten from $ISIS3DATA/mgs/kernels/ik/moc20.ti.  (The
//      PIXEL_PITCH value gotten from cam->PixelPitch is 1.0 since the
//      focal length used by ISIS in this case is in pixels)
//      For reference: the MOC WA blue filter pixel size needs an adjustment
//      of 1.000452 (see p_scale in MocWideAngleDistortionMap.cpp), so that the
//      final blue filter pixel size = (0.007 / 1.000452)
//
//   For all other cameras, cam->focal_length returns the focal
//      length in mm, as needed by Socet Set
///////////////////////////////////////////////////////////////////////////////

  focal = cam->FocalLength();  // focal length returned in mm

  if (isMocWARed) 
    focal = focal * 0.007;  // pixel to mm conversion

///////////////////////////////////////////////////////////////////////////////
// Get instrument summing modes
///////////////////////////////////////////////////////////////////////////////

  csum = (int) detectorMap->SampleScaleFactor();
  dsum = (int) detectorMap->LineScaleFactor();

  if (isLroNACL || isLroNACR) dsum = csum;

///////////////////////////////////////////////////////////////////////////////
// Calculate location of boresight in image space, these are zero-based values
//
// Note: For MOC NA, the boresight is at the image center
//       For MOC WA, MRO HiRISE MRO CTX, LRO_NACL and LRO_NACR the boresight is
//       not at the detector center, but we are forcing it there in the raw
//       (and padded) image being imported (see total_samps below). Hence the
//       boresight pixel in the orginal image space, will be the center pixel in
//       the raw image.
///////////////////////////////////////////////////////////////////////////////

  // Get line/samp of boresight pixel in detector space (summing == 1) 
  focalMap->SetFocalPlane(0.0,0.0);
  detectorBoresightSample = focalMap->DetectorSample();
  detectorBoresightLine = focalMap->DetectorLine();

  // Convert sample of boresight pixel in detector into image space
  // (summing, etc., is accounted for.)
  detectorMap->SetDetector(detectorBoresightSample,detectorBoresightLine);
  boresightSample = detectorMap->ParentSample();

///////////////////////////////////////////////////////////////////////////////
// Set Atmospheric corretion coefficients to 0 (for Mars)
///////////////////////////////////////////////////////////////////////////////

  for (i=0; i<=3; i++)
    atmco[i] = 0.0;

///////////////////////////////////////////////////////////////////////////////
// Set Lense Distortion coeffiecients 
//    - MOC NA lens distortion coefficients are set to 0
//    - MOC WA and HIRISE noproj'ed images...distortion coefficents should be 0
//    - MRO-CTX and LRO-NAC, use negative of coeficients used in ISIS
///////////////////////////////////////////////////////////////////////////////

  for (i=0; i<=9; i++) 
    lensco[i] = 0.0;

  //MRO-CTX, values gotten from mroctxAddendum003.ti and negated
  if (isCTX) {
    lensco[0] = 0.00686116;    //r
    lensco[2] = -0.0000282406; //r**3
  }

  // LRO NAC, prior to 8/22/2010: values gotten from lro_instruments_v10.ti and negated
  //          as of 8/22/2010: values gotten from lro_lroc_v14.ti as is.
  if (isLroNACL) lensco[2] = 1.81e-5; //r**3
  if (isLroNACR) lensco[2] = 1.83e-5; //r**3

  // Shadowcam, values gotten from the INS-155151_OD_K keyword in kplo_shadowcam_v00.ti and negated as of 1/26/2023 
  if (isShadowCam) {
    lensco[2] = 1.741e-5;
  }

  // LUTI, values gotten from kplo_luti_v08.ti and negated as of 10/5/2023
  if (isLUTI) {
    lensco[2] = 5.46e-5;
  }
// We no longer need WA lens distortion coefficents
// since we have noproj'ed images for MOC WA...code
// has been commented out for reference...
//  if (pars->NarrowAngle ==1)
//    for (i=0; i<=9; i++)
//      lensco[i] = 0.0;
//
//  else if (pars->WideAngleRed == 1) {
//    for (i=0; i<=9; i++)
//      lensco[i] = 0.0;
//    lensco[0] = -1.8434497e-03;   // r
//    lensco[2] = -3.3479688e-03;   // r**3
//    lensco[4] =  1.7049786e-05;   // r**5
//    lensco[6] = -7.0575441e-08;   // r**7
//    lensco[8] =  1.4522835e-10;   // r**9
//  }
//
// WA Blue images have too much distortion for mapping
// but lens distortion coefs are commented here for
// reference
//  else if (pars->WideAngleBlue == 1) {
//    for (i=0; i<=9; i++)
//      lensco[i] = 0.0;
//    lensco[0] =  3.2182150e-04;   r
//    lensco[2] = -3.5153091e-03;   r**3
//    lensco[4] =  1.8074291e-05;   r**5
//    lensco[6] = -7.3424069e-08;   r**7
//    lensco[8] =  1.4571132e-10;   r**9
//  }

///////////////////////////////////////////////////////////////////////////////
// Get total_lines, total_samps and account for summed images
///////////////////////////////////////////////////////////////////////////////

  total_lines = cube.lineCount();

  if (isMocNA)
     total_samps = 2048/csum;
  else if (isMocWARed)
     total_samps = 3564/csum;  // WA Red full detector size is 3456
                               // pixels, with the boresight at 1674.65, but
                               // the raw image being imported is padded
                               // 108/csum pixels to force the boresight
                               // at image center
  else if (isHiRise)
     total_samps = 20000/csum; // HiRise focal plane is 20048 pixels, but
                               // noproj'ed images are 20000 samples wide, with
                               // the boresight at the image center.
  else if (isCTX)
     total_samps = 5114/csum;  // MRO CTX focal plane is 5000 pixels, with
                               // the boresight at 2557.347 (1-based), but the
                               // raw image being imported is padded 114/csum
                               // pixels to force the boresight at the image
                               // center

  else if (isLroNACL) {
     //Make sure camera definition has not changed.

     expected_samples = 5064/csum;
     switch (csum) {
       case 1:
          expected_boresight = 2548.0;
          break;
       case 2:
          expected_boresight = 1274.25;
          break;
       default:
         stringstream sstr_csum;
         sstr_csum << csum;
         msg = "Unexpected summation mode, csum = " + sstr_csum.str();
         throw IException(IException::Programmer,msg,_FILEINFO_);
     }

     if(cube.sampleCount() ==  expected_samples && boresightSample ==  expected_boresight)
        total_samps = 5098/csum;  // Using lro_instruments_v08.ti, the
                                  // LRO NAC LEFT focal plane is 5064 pixels,
                                  // with the boresight at 2548 (0-based), but
                                  // the raw image being imported is padded
                                  // by 34/csum to force the boresight
                                  // at the image center.
     else {
        stringstream sstr_expected_boresight, sstr_expected_samples;
        stringstream sstr_boresightSample, sstr_samples;

        sstr_expected_boresight << expected_boresight;
        sstr_expected_samples << expected_samples;
        sstr_boresightSample << boresightSample;
        sstr_samples << cube.sampleCount();

        msg = "LRO NAC LEFT camera has changed definition!\nBoresight Sample: Expected = " + sstr_expected_boresight.str() + ", Current Value = " + sstr_boresightSample.str() + "\nImage Number of Samples: Expected = " + sstr_expected_samples.str() + ", Current Value = " + sstr_samples.str();
        throw IException(IException::Programmer,msg,_FILEINFO_);
     }
   }

  else if (isLroNACR) {
     //Make sure camera definition has not changed.

     expected_samples = 5064/csum;
     switch (csum) {
       case 1:
          expected_boresight = 2496.0;
          break;
       case 2:
          expected_boresight = 1248.25;
          break;
       default:
         stringstream sstr_csum;
         sstr_csum << csum;
         msg = "Unexpected summation mode, csum = " + sstr_csum.str();
         throw IException(IException::Programmer,msg,_FILEINFO_);
     }

     if(cube.sampleCount() ==  expected_samples && boresightSample ==  expected_boresight)
        total_samps = 5138/csum;  // Using lro_instruments_v08.ti, the
                                  // LRO NAC RIGHT focal plane is 5064 pixels,
                                  // with the boresight at 2568 (0-based), but
                                  // the raw image being imported is padded
                                  // by 74/summing to force the boresight
                                  // at the image center.

     else {
        stringstream sstr_expected_boresight, sstr_expected_samples;
        stringstream sstr_boresightSample, sstr_samples;

        sstr_expected_boresight << expected_boresight;
        sstr_expected_samples << expected_samples;
        sstr_boresightSample << boresightSample;
        sstr_samples << cube.sampleCount();

        msg = "LRO NAC RIGHT camera has changed definition!\nBoresight Sample: Expected = " + sstr_expected_boresight.str() + ", Current Value = " + sstr_boresightSample.str() + "\nImage Number of Samples: Expected = " + sstr_expected_samples.str() + ", Current Value = " + sstr_samples.str();
        throw IException(IException::Programmer,msg,_FILEINFO_);
     }
  }
  else if (isShadowCam) {
     //Make sure camera definition has not changed.

     expected_samples = 3072;
     // INSR-155151_BORESIGHT_SAMPLE
     expected_boresight = 1558.0;

     if(cube.sampleCount() ==  expected_samples && boresightSample ==  expected_boresight)
        total_samps = 3116;       // images padded by 44 pixels to force the boresight at image center 
     else {
        stringstream sstr_expected_boresight, sstr_expected_samples;
        stringstream sstr_boresightSample, sstr_samples;

        sstr_expected_boresight << expected_boresight;
        sstr_expected_samples << expected_samples;
        sstr_boresightSample << boresightSample;
        sstr_samples << cube.sampleCount();

        msg = "ShadowCam has changed definition!";
        if (boresightSample != expected_boresight)
          msg += "\nBoresight Sample: Expected = " + sstr_expected_boresight.str() + ", Current Value = " + sstr_boresightSample.str();
        if (cube.sampleCount() != expected_samples)
          msg += "\nImage Number of Samples: Expected = " + sstr_expected_samples.str() + ", Current Value = " + sstr_samples.str();
        throw IException(IException::Programmer,msg,_FILEINFO_);
     }
  }
  else if (isLUTI) {
    expected_samples = 2048;
    // INS-155101_BORESIGHT_SAMPLE (currently same as LUTIB's INS-155102_BORESIGHT_SAMPLE)
    expected_boresight = 1024;

    if(cube.sampleCount() ==  expected_samples && boresightSample ==  expected_boresight)
      total_samps = 2048;       // images not padded
    else {
      stringstream sstr_expected_boresight, sstr_expected_samples;
      stringstream sstr_boresightSample, sstr_samples;

      sstr_expected_boresight << expected_boresight;
      sstr_expected_samples << expected_samples;
      sstr_boresightSample << boresightSample;
      sstr_samples << cube.sampleCount();

      msg = "LUTI has changed definition!";
      if (boresightSample != expected_boresight)
        msg += "\nBoresight Sample: Expected = " + sstr_expected_boresight.str() + ", Current Value = " + sstr_boresightSample.str();
      if (cube.sampleCount() != expected_samples)
        msg += "\nImage Number of Samples: Expected = " + sstr_expected_samples.str() + ", Current Value = " + sstr_samples.str();
      throw IException(IException::Programmer,msg,_FILEINFO_);
    }
  }

///////////////////////////////////////////////////////////////////////////////
// Get the Interval Time in seconds
///////////////////////////////////////////////////////////////////////////////

  int_time = detectorMap->LineRate();  //LineRate is in seconds

  // For reference, this is the code if calculating interval time
  // via LineExporsureDuration keyword off image labels:
  //
  // if (isMocNA || isMocWARed)
  //   int_time = exposureDuration * (double) dsum / 1000.0;
  // else if (isHiRise)
  //   int_time = exposureDuration * (double) dsum / 1000000.0;

///////////////////////////////////////////////////////////////////////////////
// Given the interval time and total number of lines, calculate the scan
// duration in seconds
///////////////////////////////////////////////////////////////////////////////

  scan_duration = int_time * total_lines;

///////////////////////////////////////////////////////////////////////////////
// Get along and cross scan pixel size for NA and WA sensors.
// NOTE:  
//     1) The MOC WA pixel size is gotten from moc20.ti and is 7 microns
//     2) For others, cam->PixelPitch() returns the pixel pitch (size) in mm.
///////////////////////////////////////////////////////////////////////////////

  if (isMocWARed) {  /* pixel_pitch is in mm */ 
    along_scan_px_size = csum * 0.007;
    cross_scan_px_size = dsum * 0.007;
  }

  else {

    cross_scan_px_size = dsum * cam->PixelPitch();

    // Get the ephermeris time, ground position and undistorted focal plane X
    // coordinate at the center line/samp of image
    cam->SetImage(cube.sampleCount()/2.0,cube.lineCount()/2.0);

    double t_mid = cam->time().Et();

    const double lat_center = cam->UniversalLatitude();
    const double lon_center = cam->UniversalLongitude();
    const double R_center = cam->LocalRadius().meters();

    double uX_center = distortionMap->UndistortedFocalPlaneX();

    // from the ground position at the image center, increment the ephemeris
    // time by the line rate and map the ground position into the sensor in
    // undistorted focal plane coordinates

    cam->setTime(t_mid + int_time);
    double uX, uY;
    groundMap->GetXY(lat_center,lon_center,R_center, &uX, &uY);

    // the along scan pixel size is the difference in focal plane X coordinates
    along_scan_px_size = abs(uX_center - uX);

  }

///////////////////////////////////////////////////////////////////////////////
// Now that we have total_lines, total_samps, along_scan_px_size and 
// cross_scan_px_size, fill the Interior Orientation Coefficient arrays
///////////////////////////////////////////////////////////////////////////////

  for (i=0; i<=9; i++) {
    iocoef_line[i] = 0.0;
    iocoef_samp[i] = 0.0;
  }

  iocoef_line[0] = total_lines / 2.0;
  iocoef_line[1] = 1.0 / along_scan_px_size;

  iocoef_samp[0] = total_samps / 2.0;
  iocoef_samp[2] = 1.0 / cross_scan_px_size;

///////////////////////////////////////////////////////////////////////////////
// Update the Rectification Terms found in the base sensor class
///////////////////////////////////////////////////////////////////////////////

  rectification_terms[0] = total_lines / 2.0;
  rectification_terms[1] = 0.0;
  rectification_terms[2] = 1.0;
  rectification_terms[3] = total_samps / 2.0;
  rectification_terms[4] = 1.0;
  rectification_terms[5] = 0.0;

///////////////////////////////////////////////////////////////////////////////
// Fill the triangulation parameters array
///////////////////////////////////////////////////////////////////////////////

  for (i=0; i<=17; i++)
     tri_params[i] = 0.0;

  tri_params[15] = focal;

///////////////////////////////////////////////////////////////////////////////
// For now, set the Sun Azimuth and Elevation to zero (these keywords are not
// implemented yet in the SOCET SET Generic Pushbroom Sensor model)
///////////////////////////////////////////////////////////////////////////////

  sun_az = 0.0;
  sun_el = 0.0;

///////////////////////////////////////////////////////////////////////////////
// Get the Center Ground Point using the center lat/lon of the image, in radians
///////////////////////////////////////////////////////////////////////////////

  cam->SetImage(boresightSample,total_lines / 2.0);

  // Convert lat to planetographic
  Distance radii[3];
  cam->radii(radii);
  double oglat = Isis::TProjection::ToPlanetographic(cam->UniversalLatitude(), radii[0].meters() ,radii[2].meters());

  // Fill center_gp
  center_gp[0] = oglat * deg2rad;            // ographic lat in radians
  center_gp[1] = cam->UniversalLongitude();  // ISIS-3 lons are +E by default
  center_gp[2] = 0.0;
  //**** NOTE, in import_pushbroom program, center_gp[2] will be set to
  //**** the project's gp_origin_z;

  if (center_gp[1] < -180.0)    // Make sure lon is between +/- 180 degrees
      center_gp[1] = center_gp[1] + 360.0;
  if (center_gp[1] > 180.0)
      center_gp[1] = center_gp[1] - 360.0;

  center_gp[1] = center_gp[1] * deg2rad;  // +East Lon for SOCET in radians

///////////////////////////////////////////////////////////////////////////////
// Now get keyword values that depend on ephemeris data.
//
// First get the ephemeris time at the center image line coordinate.
///////////////////////////////////////////////////////////////////////////////

  center_line_coord = total_lines/2 + 0.5;  // add 0.5 pixel because the
                                            // upperleft corner of pixel
                                            // of the image has line,samp
                                            // coordinate (0.5,0.5)
  cam->SetImage(1.0,center_line_coord);
  et_center = cam->time().Et();

///////////////////////////////////////////////////////////////////////////////
// Calcualate the number of ephemeris points that are needed, based on the
// value of dt_ephem (Delta-Time-Ephemeris).  SOCET SET needs the ephemeris
// points to exceed the image range for interpolation.  For now, attempt a
// padding of 15 ephemeris points on either side of the image.
///////////////////////////////////////////////////////////////////////////////

  if (isMocNA || isHiRise || isCTX || isLroNACL || isLroNACR || isShadowCam || isLUTI) /* Try increment of 1/4 second for NA */
    dt_ephem = 0.25;  /* Make this a user definable increment? */
  else /* Increase increment for WA images to one second */
    dt_ephem = 1.0;

  num_ephem = (int)(scan_duration / dt_ephem) + 30;  // Pad by 15 ephem pts on
						     // each side of the image

  if ((num_ephem%2) == 0)  // if num_ephem is even, make it odd so that the
    num_ephem++;           // number of ephemeris points is equal on either
                           // side of t_center

///////////////////////////////////////////////////////////////////////////////
// Find the ephemeris time for the first ephemeris point, and from that, get
// to_ephem needed by SOCET (to_ephem is relative to et_center)
///////////////////////////////////////////////////////////////////////////////

  et_first_ephem = et_center - (((num_ephem - 1) / 2) * dt_ephem);

  t0_ephem = et_first_ephem - et_center;

///////////////////////////////////////////////////////////////////////////////
// Starting at et_first_ephem, get the spacecraft position and velocity
// at each ephemeris point incremented by dt_ephem
///////////////////////////////////////////////////////////////////////////////

  ephem_pts = new double [num_ephem * 3];
  ephem_rates = new double [num_ephem * 3];

  et = (double)et_first_ephem;
  for (i=0; i<num_ephem; i++) {
    spkez_c (spkCode, et, target_frame.c_str(), "LT+S", bodyCode, statevector, &lt);
    for (j=0; j<=2; j++) {
      ephem_pts[i*3+j] = statevector[j] * 1000.0; /* convert to meters */
      ephem_rates[i*3+j] = statevector[j+3] * 1000.0;
    }
    et = et + dt_ephem;
  }

///////////////////////////////////////////////////////////////////////////////
// Retrieve the position and velocity of the sensor at the image center
// (to be used for the sensor_position and mounting_angles)
///////////////////////////////////////////////////////////////////////////////

  i = (num_ephem-1)/2; 

  for (j=0; j<=2; j++) {
    center_pos[j] = ephem_pts[i*3+j];
    center_vel[j] = ephem_rates[i*3+j];
  }

///////////////////////////////////////////////////////////////////////////////
// Now calculate the sensor_position in lat, lon, height using naif routine
// recgeo_c
///////////////////////////////////////////////////////////////////////////////

  flattening = 1.0 - (radii[2]/radii[0]);
  radius = radii[0].meters();  /* convert to meters */

  recgeo_c (center_pos, radius, flattening, &lon, &lat, &height);

  sensor_position[0] = lat;  /* lat/lon are in radians */
  sensor_position[1] = lon;  /* Lon is +East out of recgeo_c */
  sensor_position[2] = height;  /* height is in meters */

  // Make sure lon is between +/- 180 degrees
  if (sensor_position[1] < -180.0 * deg2rad)
    sensor_position[1] = sensor_position[1] + 360.0 * deg2rad;
  if (sensor_position[1] > 180.0 * deg2rad)
    sensor_position[1] = sensor_position[1] - 360.0 * deg2rad;

///////////////////////////////////////////////////////////////////////////////
// Finally set up matrices for extracting mounting_angles in OPK system
//
// First convert the euler angles obtained from the INSTPARS file to the matrix
// which rotates from the platform coordinate system to the camera coordinate
// system (CP)
///////////////////////////////////////////////////////////////////////////////

  if (isMocNA || isMocWARed) {
    const QString eulerKey = "INS" + QString(ikCode) + "_EULER_ANGLES";
    yaw = cam->getDouble(eulerKey, 0);
    pitch = cam->Spice::getDouble(eulerKey, 1);
    roll = cam->getDouble(eulerKey, 2);

    eul2m_c ((double)yaw, (double)pitch, (double)roll,
              3, 2, 1, CP);
  }
  else if (isHiRise) {

  // SS readout is closer to ISIS lat/lon readout when *not* applying 'CP'
  // So commented code to get CP from pxform_c for now, and 
  // set CP to the identity matrix.  (Randy also thought we don't 
  // need this matrix for HiRISE because we have set it up as an idealized
  // camera...)

     //pxform_c("MRO_SPACECRAFT","MRO_HIRISE_OPTICAL_AXIS",et_center,CP);
     for (i=0; i<=2; i++)
        for (j=0; j<=2; j++)
           CP[i][j] = 0.0;
     CP[0][0]=1.0;
     CP[1][1]=1.0;
     CP[2][2]=1.0;
  }
  else if (isCTX)
     pxform_c("MRO_SPACECRAFT","MRO_CTX",et_center,CP);

  else if (isLroNACL)
     pxform_c("LRO_SC_BUS","LRO_LROCNACL",et_center,CP);

  else if (isLroNACR)
     pxform_c("LRO_SC_BUS","LRO_LROCNACR",et_center,CP);

  else if (isShadowCam)
     pxform_c("KPLO_SPACECRAFT","KPLO_SHC_A",et_center,CP);
  
  else if (isLUTI) {
    PvlGroup archive = cube.label()->findGroup("Archive", Pvl::Traverse);
    bool isLUTIA = archive["ProductId"][0].toStdString().find("LUTIA") != std::string::npos;
    
    if (isLUTIA) {
      pxform_c("KPLO_SPACECRAFT", "KPLO_LUTIA", et_center, CP);
    }
    else {
      pxform_c("KPLO_SPACECRAFT", "KPLO_LUTIB", et_center, CP);
    }
  }
    

///////////////////////////////////////////////////////////////////////////////
// Calculate the center Ephermeris time for the C-kernal and get the
// J2000->platform and J2000->bodyfixed rotation matrices
///////////////////////////////////////////////////////////////////////////////

  if(isCTX || isLroNACL || isLroNACR || isShadowCam || isLUTI)
     ck_time_bias = 0.0;
  else {
     QString ckTimeBiasKey = "INS" + QString(ikCode) + "_CK_TIME_BIAS";
     ck_time_bias = cam->getDouble(ckTimeBiasKey, 0);
     //ck_time_bias = Isis::Spice::getDouble(ckTimeBiasKey, 0);
  }
  ck_et_center = et_center + ck_time_bias;

  if(get_naif_matrices (cam, et_center,ck_et_center,spkCode,ikCode,bodyCode,
                       BJ, PJ)) {
    msg = "Cannot Generate Keywords File: Camera pointing not available in SPICE C-kernel";
    throw IException(IException::Io,msg,_FILEINFO_);
  }

///////////////////////////////////////////////////////////////////////////////
// Now calculated the mounting angles in OPK, radians
///////////////////////////////////////////////////////////////////////////////

  get_mounting_angles ((double *)CP, (double *)BJ, (double *)PJ,
                       center_pos, center_vel, mounting_angles);

///////////////////////////////////////////////////////////////////////////////
// Check if this is a backwards-acquired LRO NAC image, and populate a
// PVL file accordingly.
///////////////////////////////////////////////////////////////////////////////

  if (isLroNACL || isLroNACR || isShadowCam || isLUTI) {
    flipFlag = "N";
    if ( (isLroNACL && abs(mounting_angles[2])/deg2rad < 20) ||
         (isLroNACR && abs(mounting_angles[2])/deg2rad > 160) ||
         (isShadowCam && abs(mounting_angles[2])/deg2rad < 20) ||
         (isLUTI && abs(mounting_angles[2])/deg2rad > 160)  // Unsure yet what the trigger is, best guess is NACR, watch to verify these are flipped
       ) 
      flipFlag = "Y";

    // Generate the ouput pvl file name and open it...
    PvlGroup archive = cube.label()->findGroup("Archive", Pvl::Traverse);
    if (isShadowCam)
      flipPVLFile = archive["PDSIdentifier"][0].toStdString() + "_flip.pvl";
    else
      flipPVLFile = archive["ProductId"][0].toStdString() + "_flip.pvl";
    flipFlagStrm.open (FileName(QString::fromStdString(flipPVLFile)).expanded().toStdString().c_str(),ios::trunc);
    // Populate the pvl file
    flipFlagStrm << "Group = FlipFlag" << endl;
    flipFlagStrm << "   Flip = " << flipFlag.c_str() << endl;
    flipFlagStrm << "End_Group" << endl;
    flipFlagStrm.close ();
  }

///////////////////////////////////////////////////////////////////////////////
// If kappa is ~180 degrees, subtract 180 degrees from it for SOCET SET
// (When kappa = ~180 degrees, the sensor was flying backwards, this is a
// configuration SOCET SET cannot handle directly.)
///////////////////////////////////////////////////////////////////////////////
  if ((mounting_angles[2] >= 170.0*deg2rad) || (mounting_angles[2] <= -170.0*deg2rad)) {
        mounting_angles[0] = (-1)*mounting_angles[0];
        mounting_angles[1] = (-1)*mounting_angles[1];
        if (mounting_angles[2] >= 170.0*deg2rad)
                mounting_angles[2] = mounting_angles[2] - pi_c();
        else
                mounting_angles[2] = mounting_angles[2] + pi_c();
  }
  
///////////////////////////////////////////////////////////////////////////////
// We are done with computing keyword values, so output the Generic Pushbroom
// Keyword file.
///////////////////////////////////////////////////////////////////////////////

  pushKeyStrm.setf(std::ios::scientific);
  pushKeyStrm << "RECTIFICATION_TERMS" << endl;
  pushKeyStrm << "        " << setprecision(14) << rectification_terms[0] << " " << rectification_terms[1] << " " << rectification_terms[2] << endl;
  pushKeyStrm << "        " << setprecision(14) << rectification_terms[3] << " " << rectification_terms[4] << " " << rectification_terms[5] << endl;

  pushKeyStrm << "GROUND_ZERO ";
  pushKeyStrm << setprecision(14) << center_gp[0] << " " << center_gp[1] << " " << center_gp[2] << endl;

  pushKeyStrm << "LOAD_PT ";
  pushKeyStrm << setprecision(14) << center_gp[0] << " " << center_gp[1] << " " << center_gp[2] << endl;

  pushKeyStrm << "COORD_SYSTEM 1" << endl;

  pushKeyStrm << "IMAGE_MOTION 0" << endl;

///////////////////////////////////////////////////////////////////////////////
// Output SOCET SET generic pushbroom sensor model portion of support file
///////////////////////////////////////////////////////////////////////////////

  pushKeyStrm << "SENSOR_TYPE GENERIC_PUSHBROOM" << endl;
  pushKeyStrm << "SENSOR_MODE UNKNOWN" << endl;

  pushKeyStrm << "FOCAL ";
  pushKeyStrm << setprecision(14) << focal << endl;

  pushKeyStrm << "ATMCO";
  for(i=0; i<4; i++) pushKeyStrm << " " << setprecision(14) << atmco[i];
  pushKeyStrm << endl;

  pushKeyStrm << "LENSCO";
  for(i=0; i<10; i++) pushKeyStrm << " " << setprecision(14) << lensco[i];
  pushKeyStrm << endl;

  pushKeyStrm << "IOCOEF_LINE";
  for(i=0; i<10; i++) pushKeyStrm << " " << setprecision(14) << iocoef_line[i];
  pushKeyStrm << endl;

  pushKeyStrm << "IOCOEF_SAMPLE";
  for(i=0; i<10; i++) pushKeyStrm << " " << setprecision(14) << iocoef_samp[i];
  pushKeyStrm << endl;

  pushKeyStrm << "ABERR    0" << endl;
  pushKeyStrm << "ATMREF   0" << endl;
  pushKeyStrm << "PLATFORM   1" << endl;
  pushKeyStrm << "SOURCE_FLAG  1" << endl;
  pushKeyStrm << "SINGLE_EPHEMERIDE  0" << endl;

  pushKeyStrm << "TRI_PARAMETERS" << endl;
  pushKeyStrm << setprecision(14) << tri_params[0];
  for(i=1; i<18; i++) pushKeyStrm << " " << setprecision(14) << tri_params[i];
  pushKeyStrm << endl;

  pushKeyStrm << "\n\nT_CENTER  0.0000e+000" << endl;

  pushKeyStrm << "DT_EPHEM  ";
  pushKeyStrm << setprecision(14) << dt_ephem << endl;

  pushKeyStrm << "T0_EPHEM  ";
  pushKeyStrm << setprecision(14) << t0_ephem << endl;

  pushKeyStrm << "NUMBER_OF_EPHEM   " << num_ephem << endl;

  pushKeyStrm << "EPHEM_PTS" << endl;
  for(i=0; i<num_ephem; i++) {
    pushKeyStrm << setprecision(14) << " " << ephem_pts[i*3];
    pushKeyStrm << setprecision(14) << " " << ephem_pts[i*3+1];
    pushKeyStrm << setprecision(14) << " " << ephem_pts[i*3+2] << endl;
  }

  pushKeyStrm  << "\n\nEPHEM_RATES" << endl;
  for(i=0; i<num_ephem; i++) {
    pushKeyStrm << setprecision(14) << " " << ephem_rates[i*3];
    pushKeyStrm << setprecision(14) << " " << ephem_rates[i*3+1];
    pushKeyStrm << setprecision(14) << " " << ephem_rates[i*3+2] << endl;
  }

  pushKeyStrm << "\n\nSCAN_DURATION " << setprecision(14) << scan_duration << endl;
  pushKeyStrm << "INT_TIME " << setprecision(14) << int_time << endl;
  pushKeyStrm << "ALONG_SCAN_PIXEL_SIZE  " << setprecision(14) << along_scan_px_size << endl;
  pushKeyStrm << "CROSS_SCAN_PIXEL_SIZE  " << setprecision(14) << cross_scan_px_size << endl;
  pushKeyStrm << "SUN_AZIMUTH  " << setprecision(14) << sun_az << endl;
  pushKeyStrm << "SUN_ELEVATION  " << setprecision(14) << sun_el << endl;

  pushKeyStrm << "\n\nCENTER_GP";
  for(i=0; i<3; i++) 
    pushKeyStrm << " " << setprecision(14) << center_gp[i];
  pushKeyStrm << endl;

  pushKeyStrm << "SENSOR_POSITION";
  for(i=0; i<3; i++) 
    pushKeyStrm << " " << setprecision(14) << sensor_position[i];
  pushKeyStrm << endl;

  pushKeyStrm << "MOUNTING_ANGLES";
  for(i=0; i<3; i++) 
    pushKeyStrm << " " << setprecision(14) << mounting_angles[i];
  pushKeyStrm << endl;

  pushKeyStrm << "TOTAL_LINES " << total_lines << endl;
  pushKeyStrm << "\n\n\nTOTAL_SAMPLES " << total_samps << endl;
  pushKeyStrm << "\n\n\n" << endl;

///////////////////////////////////////////////////////////////////////////////
// All done so cleanup and exit
///////////////////////////////////////////////////////////////////////////////

  Unload_kernels ();
  pushKeyStrm.close ();

  p.EndProcess();
}

// Load/furnish NAIF kernel(s)
void Load(Isis::PvlKeyword &key)
{
  NaifStatus::CheckErrors();

  for (int i=0; i<key.size(); i++) {
    if (key[i] == "") continue;
    if ((key[i]).toUpper() == "NULL") continue;
    if ((key[i]).toUpper() == "NADIR") continue;
    if ((key[i]).toUpper() == "TABLE") continue;
    Isis::FileName file(key[i]);
    if (!file.fileExists()) {
      string msg = "Spice file does not exist [" + file.expanded().toStdString() + "]";
      throw Isis::IException(Isis::IException::Io,QString::fromStdString(msg),_FILEINFO_);
    }
    string fileName(file.expanded().toStdString());
    furnsh_c(fileName.c_str());
    p_kernels.push_back(key[i].toStdString());
  }

  NaifStatus::CheckErrors();
}

// Unload the kernels
void Unload_kernels()
{

  NaifStatus::CheckErrors();

  for (unsigned int i=0; i<p_kernels.size(); i++) {
    Isis::FileName file(QString::fromStdString(p_kernels[i]));
    string fileName(file.expanded().toStdString());
    unload_c(fileName.c_str());
  }

  NaifStatus::CheckErrors();
}
