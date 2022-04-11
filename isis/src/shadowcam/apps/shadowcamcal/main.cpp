/** This application will calibrate EDR images taken by the ShadowCam TDI Camera

After companding and formatting for ISIS's tiling format, the bias pixel average must be calculated per each channel. 
Once the average is calculated, the non-seen columns for per channel must be cropped. The resulting image will be smaller than the original image.

The pixel width will be reduced from 3144 total of all six channels, to 3072 for all six channels. 

The remaining calibration steps can then be performed on the reduced image.**/

#include "Isis.h"
#include "ProcessByBrick.h"
#include "SpecialPixel.h"
#include "Camera.h"
#include "iTime.h"
#include "IException.h"
#include "TextFile.h"
#include "Brick.h"
#include "Table.h"
#include "Statistics.h"
#include <fstream>
#include <vector>

using namespace std;
using namespace Isis;

// Working functions and parameters
void ResetGlobals();
void CopyCubeIntoArray(QString &fileString, vector<double> &data);
void ReadTextDataFile(QString &fileString, vector<double> &data);
void ReadTextDataFile(QString &fileString, vector<vector<double> > &data);
void Calibrate(Buffer &in, Buffer &out);
void CalculateBiasPixelAverage(Buffer &in));
void CropNonSeenPixelColumns(Buffer &in));
void SubtractBiasPixelAverage(Buffer &line);
//void CorrectDark(Buffer &in);
//void CorrectNonlinearity(Buffer &in);
//void CorrectFlatfield(Buffer &in);
//void RadiometricCalibration(Buffer &in);

#define TOTAL_PIXELS 3144
#define PIXELS_PER_CHANNEL 524
#define PRESCAN_PIXELS_PER_CHANNEL 2
#define OVERSCAN_PIXELS_PER_CHANNEL 2
#define BLACK_PIXELS_PER_CHANNEL 8
#define SEEN_PIXELS_PER_CHANNEL (PIXELS_PER_CHANNEL - PRESCAN_PIXELS - OVERSCAN_PIXELS - BLACK_PIXELS) //512
#define NUM_CHANNELS 6
#define TOTAL_SEEN_PIXELS (SEEN_PIXELS_PER_CHANNEL * NUM_CHANNELS) //3072
//#define MAXNONLIN 600
#define SOLAR_RADIUS 695500
#define KM_PER_AU 149597871
//#define MASKED_PIXEL_VALUES 8

// Main moccal routine
void IsisMain() {
  ResetGlobals();
  // We will be processing by line
  ProcessByBrick p;

  // Setup the input and make sure it is a ctx file
  UserInterface &ui = Application::GetUserInterface();

  g_masked = ui.GetBoolean("BIASPIXELSUBTRACT");
  g_dark = ui.GetBoolean("DARKCORRECT");
  g_nonlinear = ui.GetBoolean("NONLINEARITY");
  g_flatfield = ui.GetBoolean("FLATFIELD");
  g_radiometric = ui.GetBoolean("RADIOMETRIC");
  g_iof = (ui.GetString("RADIOMETRICTYPE") == "IOF");

  Isis::Pvl lab(ui.GetCubeName("FROM"));
  Isis::PvlGroup &inst = lab.findGroup("Instrument", Pvl::Traverse);

  p.EndProcess();
}

void ResetGlobals() {
 /* g_exposure = 1.0; // Exposure duration
  g_solarDistance = 1.01; // average distance in [AU]

  g_maskedPixelsLeft.clear();
  g_maskedPixelsRight.clear();

  g_radianceLeft = 1.0;
  g_radianceRight = 1.0;
  g_iofLeft = 1.0;
  g_iofRight = 1.0;

  g_summed = true;
  g_masked = true;
  g_dark = true;
  g_nonlinear = true;
  g_flatfield = true;
  g_radiometric = true;
  g_iof = true;
  g_isLeftNac = true;
  g_maskedLeftOnly = false;
  g_averageDarkLine.clear();
  g_linearOffsetLine.clear();
  g_flatfieldLine.clear();
  g_linearityCoefficients.clear();*/
}

  // Check if it is a NAC image