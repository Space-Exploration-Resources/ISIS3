#include <iostream>
#include <iomanip>
#include <string>
#include "Camera.h"
#include "Spice.h"
#include "IException.h"

using namespace std;

int get_naif_matrices (Isis::Camera* cam, double et, double ck_et,
                       int spkCode, int ikCode, int bodyCode, 
		       double BJ[3][3], double PJ[3][3])

/******************************************************************************
*                           ---------------------                             *
*                           | get_naif_matrices |                             *
*                           ---------------------                             *
*                                                                             *
*_Title - get_naif_matrices  Get rotation matrices via Naif kernels           *
*                                                                             *
*_Args  Type            Variable      I/O    Description                      *
*                                                                             *
*       double          et             I     Ephemeris time of desired states *
*                                                                             *
*       double          ck_et          I     ET to use for getting output     *
*                                            from ck                          *
*       int             spkCode        I     Naif code for spacecraft         *
*                                            carrying instrument              *
*       int             ikCode         I     Naif code for viewing instrument *
*                                                                             *
*       int             bodyCode       I     Naif code for target body        *
*                                                                             *
*       double          *BJ            I     Transformation matrix for rotat- *
*                                            ing J2000 to body-fixed          *
*                                            coordinates                      *
*                                                                             *
*       double          *PJ            I     Transformation matrix for rotat- *
*                                            ing J2000 to platform            *
*                                            coordinates                      *
*                                                                             *
*       int             get_naif_matrices                                     *
*                                      O     Return status                    *
*                                            0  okay                          *
*                                           <0  error                         *
*                                                                             *
*_Desc  This routine is a modification of the ISIS levs routine               *
*       lev1u_naif_get_pointing_etal.  This routine calculates the rotational *
*       matrices BJ and PJ needed to rotate vectors in the camera system of   *
*       the instrument to the body-fixed system of the target body.           *
*                                                                             *
*_Hist  Apr 20 2001 Elpitha H. Kraus, USGS, Flagstaff Original Version        *
*       Feb 25 2002 EHK - additional documentation and the calculation of a   *
*                         new value, ck_matched_time,  was added to the       *
*                         levinst structure in ISIS on Jan 14, 2002, so I     *
*                         mirrored the update in this routine.
*                                                                             *
*                         Also, at some time not documented in the history    *
*                         of lev1u_naif_get_pointing_etal, the interpolation  *
*                         logic in calculating PJ was removed from            *
*                         lev1u_naif_get_pointing_etal, so I applied the      *
*                         updates to match the current version of             *
*                         lev1u_naif_get_pointing_etal's calculations of PJ   *
*                                                                             *
*                         Finally, removed (double)*CP from parameter list of *
*                         get_naif_matrices.  CP was not being used in any    *
*                         calculations.                                       *
*       May  7 2007 EHK - Renamed String class to iString in accordance with  *
*                         changes made in the ISIS3.1.11 release              *
*       Feb 12 2008 EHK - The MRO CTX camera does not have a CK_TIME_TOLERANCE*
*                         keyword in the instrument addemdum kernel, so added *
*                         a try-block to handle this case, and set the default*
*                         tolerance to 1.0 sec                                *
*      Aug 27 2021 CHM -  Re-added to run in the 3.10.2_lroc version of ISIS. *
*                                                                             *
*_End                                                                         *
*******************************************************************************/

{
  SpiceDouble sclkdp, sclkoutdp;
  SpiceDouble tol;
  SpiceBoolean found;

  SpiceDouble ck_etout;

/*****************************************************************************
* First get the transformation matrix from J2000 to body-fixed (BJ)
*****************************************************************************/

  tipbod_c ("J2000", (SpiceInt)bodyCode, (SpiceDouble)et, BJ);

/*****************************************************************************
* Now convert the ephemeris time to continuous spacecraft ticks and try to
* get the transformation matrix from J2000 to the platform (PJ).  Set
* tolerance to 1.0 for a Type 1 c kernel and set tolerance to 0.0 for
* all other c kernel types (2 and 3).
*****************************************************************************/

  sce2c_c ((SpiceInt)spkCode, (SpiceDouble)ck_et, &sclkdp);

  const QString key = "INS" + QString::number(ikCode) +  "_CK_TIME_TOLERANCE";
  try {
      tol = cam->getDouble(key);
      //tol = Isis::Spice::getDouble(key);
  }
  catch (Isis::IException &e) {
    printf("setting tol = 1.0\n");
    tol = 1.0;
  }

//??? ask Debbie about ikCode vs ckCode...
//  ckgp_c ((SpiceInt)ikCode, sclkdp, tol, "J2000", PJ, &sclkoutdp, &found);
  int ckCode = spkCode * 1000;
  ckgp_c ((SpiceInt)ckCode, sclkdp, tol, "J2000", PJ, &sclkoutdp, &found);

/*****************************************************************************
* If PJ wasn't found with a tolerance of 0 then issue an error
*****************************************************************************/

  if (!found) {
    return -1;
  }

/*************************************************************************
 * Get the et time of the matched instance used for getting the camera
 * pointing
 *************************************************************************/

  sct2e_c((SpiceInt)spkCode, (SpiceDouble) sclkoutdp, &ck_etout);
  /*spi->inst.ck_matched_time = ck_etout;*/

/*****************************************************************************
* Ok we've got a good PJ transformation matrix so return
*****************************************************************************/

  return 0;
}
