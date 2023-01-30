#include <iostream>
#include <iomanip>
#include <string>
#include <math.h>
#include "Spice.h"

using namespace std;

void get_mounting_angles (double *CP, double *BJ, double *PJ,
                          double *pos_USR, double *vel_USR,
                          double *mounting_angles)
/******************************************************************************
*                           -----------------------                           *
*                           | get_mounting_angles |                           *
*                           -----------------------                           *
*                                                                             *
*_Title - get_mounting_angles  Get mounting angles for SS pushbroom sensor    *
*                                                                             *
*_Args  Type            Variable         I/O    Description                   *
*                                                                             *
*       double          *CP               I     Platform -> Camera rotation   *
*                                               matrix                        *
*                                                                             *
*       double          *BJ               I     J2000 -> body fixed rotation  *
*                                               matrix                        *
*                                                                             *
*       double          *PJ               I     J2000 -> platform rotation    *
*                                               matrix                        *
*                                                                             *
*       double          *pos_USR          I     Sensor postion at center of   *
*                                               image, in body-fixed, xyz,    *
*                                               coordinates, meters           *
*                                                                             *
*       double          *vel_USR          I     Sensor velocity at center of  *
*                                               image, in body-fixed, xyz,    *
*                                               coordinates, meters/sec       *
*                                                                             *
*       double          *mounting_angles  O     Omega, phi, kappa mounting    *
*                                               angles at center of image,    *
*                                               radians                       * 
*                                                                             *
*_Desc  This routine computes mounting angles for the SS Generic Pushbroom    *
*       sensor model.  Currently, the mounting angles incorporate the full    *
*       s/c pointing at the center of the image, along with platform mounting *
*       angles                                                                *
*                                                                             *
*_Hist  Apr 20 2001 Elpitha H. Kraus, USGS, Flagstaff Original Version        *
*       Jan 24 2008 EHK - We traced a bug (resulting in HiRISE images         *
*                         importing ~3-km off target at 45 degrees latitude)  *
*                         to an ographic latitude vs ocentric latitude mix-up.*
*                         Although Socet Set works in ographic latitudes, the *
*                         pushbroom sensor model bases mounting angles on     *
*                         ocentric coordiantes.  We therefore modified the    *
*                         calculation of the 'E' matrix (Body-Fixed -> LSR)   *
*                         on the ocentric lat/lon at pos_USR.  We still have  *
*                         a 40m 'latitide' and 200m 'longititude' residual    *
*                         error, with respect to ISIS3 calculations, to work  *
*                         out.  This may be related to a wrong CP matrix      *
*       Jan 25 2008 EHK - modified to use same Ground Track to                *
*                         Body-Fixed (USR) matrix rotation as the             *
*                         generic pushbroom sensor (named BG in this routine. *
*                         This matrix is equivalent to our previous ground    *
*                         track to body-fixed rotation matrices, E and GL     *
*                         (Body-Fixed -> LSR & LSR -> Ground Track matrices), *
*                         but is a cleaner approach                           *
*       Jan 29 2008 EHK - added more documentation                            *
*       Jan 28 2010 EHK - In order to determine LRO backwards-aquired images  *
*                         removed 180 degree adjustment to kappa from this    *
*                         routine and placed it in calc_pushbroom_keywords.   *
*      Aug 27 2021 CHM -  Re-added to run in the 3.10.2_lroc version of ISIS. *
*                                                                             *
*_End                                                                         *
*******************************************************************************/

{
  double SC[3][3];     // rotation matrix from camera coordinates
                       // to Socet Set plate/focal-plane coordinates (S)
                       // (For Socet, we need +x = along flight direction
                       //                     +y = left of +x
                       //                     +z = up
  double uI[3];        // spacecraft in-track unit vector
  double uC[3];        // spacecraft cross-track unit vector
  double uR[3];        // spacecraft radial unit vector
  double BG[3][3];     // Ground Track (LSR) to Body-Fixed rotation matrix

  double R[3][3];      // Overall rotation matrix from camera to Ground Track
                       // (LSR) coordinates.  The mounting angles are extracted
                       // from this matrix
  double omega, phi, kappa; // Socet Set omega, phi, kappa rotation angles,
                            //  radians

  double temp_mag;
  double temp1[3][3], temp2[3][3];

//*****************************************************************************
// Initilize matrix SC = [PI]
//*****************************************************************************

  eul2m_c (pi_c(), 0.0, 0.0, 1, 2, 3, SC);

//*****************************************************************************
// Build BG: the Ground Trac -to- Body-Fixed rotation matrix....
//
// BG  = (I C R)
//
// where,
//
// I C and R above are the column vectors of BG and,
//
//          pos_USR           pos_USR X vel_USR
//    R = ------------ ,  C = ----------------- , I = C x R.
//        | pos_USR  |       |pos_USR X vel_USR|
//
// where,
//   pos_USR is spacecraft position and vel_USR is spacecraft velocity
//   I,C,R are the in-track, cross-track and Radial spacecraft vectors
//*****************************************************************************

   unorm_c(pos_USR,uR,&temp_mag);

   ucrss_c(pos_USR, vel_USR, uC);

   vcrss_c(uC, uR, uI);

   BG[0][0] = uI[0];
   BG[0][1] = uC[0];
   BG[0][2] = uR[0];
   BG[1][0] = uI[1];
   BG[1][1] = uC[1];
   BG[1][2] = uR[1];
   BG[2][0] = uI[2];
   BG[2][1] = uC[2];
   BG[2][2] = uR[2];


/******************************************************************************
* We now have all the necessary matrices so compute the rotation matrix for
* the mounting angles
******************************************************************************/

  mxm_c (SC,CP,temp1);           // [SC] [CP]
  mxm_c (temp1,PJ,temp2);        // [SC] [CP] [PJ]
  mxmt_c (temp2,BJ,temp1);       // [SC] [CP] [PJ] [BJ]-transpose
  mxm_c (temp1,BG,R);            // [SC] [CP] [PJ] [BJ]-transpose [BG]

/******************************************************************************
* Decompose R to get mounting angles as omega, phi, kappa
******************************************************************************/

  m2eul_c (R, 1, 2, 3, &omega, &phi, &kappa);

  mounting_angles[0] = omega;
  mounting_angles[1] = phi;
  mounting_angles[2] = kappa;

  return;
}
