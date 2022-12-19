/** This is free and unencumbered software released into the public domain.

The authors of ISIS do not claim copyright on the contents of this file.
For more details about the LICENSE terms and the AUTHORS, you will
find files of those names at the top level of this repository. **/

/* SPDX-License-Identifier: CC0-1.0 */

#include "ProcessByLine.h"
#include "SpecialPixel.h"
#include "Message.h"
#include "Camera.h"
#include "iTime.h"
#include "IException.h"
#include "TextFile.h"
#include "Brick.h"
#include "Table.h"
#include "PvlGroup.h"
#include "Statistics.h"
#include "UserInterface.h"
#include "shadowcamcal.h"
#include <fstream>
#include <QTextStream>
#include <QDir>
#include <QRegExp>
#include <QString>
#include <vector>

using namespace std;
namespace Isis {

  void ResetGlobals();

  void BiasPixelSubtract(Buffer &in );
  void CopyCubeIntoVector(QString &fileString, vector<double> &data);
  void ReadTextDataFile(QString &fileString, vector<double> &data);
  void ReadTextDataFile(QString &fileString, vector<vector<double> > &data);
  void Calibrate(Buffer &in, Buffer &out);
  void RemoveMaskedOffset(Buffer &line);
  void CorrectDark(Buffer &in);
  void CorrectNonlinearity(Buffer &in);
  void CorrectFlatfield(Buffer &in);
  void RadiometricCalibration(Buffer &in);
  void GetNearestDarkFile(QString fileString, QString &file);
  void GetNearestDarkFilePair(QString &fileString, QString &file0, QString &file1);
  void GetCalibrationDirectory(QString calibrationType, QString &calibrationDirectory);
  void GetWeightedDarkAverages();
  bool AllowedSpecialPixelType(double pixelValue);

  #define LINE_SIZE 5064
  #define MAXNONLIN 600
  #define SOLAR_RADIUS 695500
  #define KM_PER_AU 149597871
  #define MASKED_PIXEL_VALUES 8


  double g_radianceLeft, g_radianceRight, g_iofLeft, g_iofRight, g_imgTime;
  double g_exposure; // Exposure duration
  double g_solarDistance; // average distance in [AU]

  bool g_biasSubtract, g_summed, g_masked, g_maskedLeftOnly, g_dark, g_nonlinear, g_flatfield, g_radiometric, g_iof, g_isLeftNac;
  bool g_nearestDark, g_nearestDarkPair, g_customDark;
  vector<int> g_maskedPixelsLeft, g_maskedPixelsRight;
  vector<double> g_avgDarkLineCube0, g_avgDarkLineCube1, g_linearOffsetLine, g_flatfieldLine, g_darkTimes, g_weightedDarkTimeAvgs;
  vector<vector<double> > g_linearityCoefficients;
  Buffer *g_darkCube0, *g_darkCube1;

  /**
    * @brief  Calling method of the application
    *
    * Performs radiometric corrections to images acquired by the Narrow Angle
    *    Camera aboard the Lunar Reconnaissance Orbiter spacecraft.
    *
    * @param ui The user interfact to parse the parameters from. 
    */
  void shadowcamcal(UserInterface &ui){
     //Cube *iCube = p.SetInputCube("FROM", OneBand);
    Cube iCube(ui.GetCubeName("FROM"));
    shadowcamcal(&iCube, ui);
  }

  
  void shadowcamcal(Cube *iCube, UserInterface &ui) {
    ResetGlobals();
    // We will be processing by line
    ProcessByLine p;
    g_biasSubtract = ui.GetBoolean("BIASPIXEL");
    g_masked = ui.GetBoolean("MASKED");
    g_dark = ui.GetBoolean("DARK");
    g_nonlinear = ui.GetBoolean("NONLINEARITY");
    g_flatfield = ui.GetBoolean("FLATFIELD");
    g_radiometric = ui.GetBoolean("RADIOMETRIC");
    g_iof = (ui.GetString("RADIOMETRICTYPE") == "IOF");

    Isis::Pvl lab(ui.GetCubeName("FROM"));
    Isis::PvlGroup &inst = lab.findGroup("Instrument", Pvl::Traverse);

    // Check if it is a NAC image
    QString instId = (QString) inst["InstrumentId"];
    instId = instId.toUpper();

    // And check if it has already run through calibration
    if(lab.findObject("IsisCube").hasGroup("Radiometry")) {
      QString msg = "This image has already been calibrated";
      throw IException(IException::User, msg, _FILEINFO_);
    }

    if(lab.findObject("IsisCube").hasGroup("AlphaCube")) {
      QString msg = "This application can not be run on any image that has been geometrically transformed (i.e. scaled, rotated, sheared, or reflected) or cropped.";
      throw IException(IException::User, msg, _FILEINFO_);
    }


    g_exposure = inst["LineExposureDuration"];

    p.SetInputCube(iCube, OneBand);

    // If there is any pixel in the image with a DN > 1000
    //  then the "left" masked pixels are likely wiped out and useless
    if(iCube->statistics()->Maximum() > 1000)
      g_maskedLeftOnly = true;

    QString flatFile, offsetFile, coefficientFile;
    
    vector <QString> darkFiles;

    // This is reserved for file containing custom bias pixel averages per channel. This is only to be used when testing with specific test values. 
    if(g_biasSubtract){
      QString biasPixelAvgsFile = ui.GetString("BIASPIXELFILE");
    }

    if(g_dark) {
      QString darkFileType = ui.GetString("DARKFILETYPE");
    }

    if(g_nonlinear) {
      offsetFile = ui.GetAsString("OFFSETFILE");
    }

    if(g_flatfield) {
      flatFile = ui.GetAsString("FLATFIELDFILE");
    }

    if(g_radiometric) {
      QString radFile = ui.GetAsString("RADIOMETRICFILE");

      
    }
    // Setup the output cube
    Cube * oCube = p.SetOutputCube(ui.GetCubeName("TO"), ui.GetOutputAttribute("TO")); 
    // Start the line-by-line calibration sequence
    p.StartProcess(Calibrate);

    PvlGroup calgrp("Radiometry");
    if(g_masked) {
      PvlKeyword darkColumns("DarkColumns");
      for(unsigned int i = 0; i < g_maskedPixelsLeft.size(); i++)
        darkColumns += toString(g_maskedPixelsLeft[i]);
      for(unsigned int i = 0; i < g_maskedPixelsRight.size(); i++)
        darkColumns += toString(g_maskedPixelsRight[i]);
      calgrp += darkColumns;
    }

    if(g_dark){
      PvlKeyword darks("DarkFiles");
      darks.addValue(darkFiles[0]);
      if(g_nearestDark)
        calgrp += PvlKeyword("DarkFileType", "NearestDarkFile");
      else if (g_nearestDarkPair){
        calgrp += PvlKeyword("DarkFileType", "NearestDarkFilePair");
        darks.addValue(darkFiles[1]);
      }
      else
        calgrp += PvlKeyword("DarkFileType", "CustomDarkFile");

      calgrp += darks;
    }

    if(g_nonlinear) {
      
    }

    if(g_flatfield)
      calgrp += PvlKeyword("FlatFile", flatFile);
    if(g_radiometric) {
      
    }

    oCube->putGroup(calgrp);
    p.EndProcess();
  }

  /**
  * This method resets global variables
  *
  */
  void ResetGlobals() {
    g_biasSubtract = 0.0;
    g_exposure = 1.0; // Exposure duration
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
    g_nearestDarkPair = false;
    g_nearestDark = false;
    g_customDark = false;
    g_avgDarkLineCube0.clear();
    g_avgDarkLineCube1.clear();
    g_linearOffsetLine.clear();
    g_darkTimes.clear();
    g_weightedDarkTimeAvgs.clear();
    g_flatfieldLine.clear();
    g_linearityCoefficients.clear();
    g_imgTime = 0.0;
  }

  /**
  * This method processes buffer by line to calibrate
  *
  * @param in Buffer to hold 1 line of cube data
  * @param out Buffer to hold 1 line of cube data
  *
  */
  void Calibrate(Buffer &in, Buffer &out) {

  }

  void ReadTextDataFile(QString &fileString, vector<double> &data) {
    FileName filename(fileString);
    if(filename.isVersioned())
      filename = filename.highestVersion();
    if(!filename.fileExists()) {
      QString msg = fileString + " does not exist.";
      throw IException(IException::User, msg, _FILEINFO_);
    }
    TextFile file(filename.expanded());
    QString lineString;
    unsigned int line = 0;
    while(file.GetLine(lineString)) {
      data.push_back(toDouble(lineString.split(QRegExp("[ ,;]")).first()));
      line++;
    }
    fileString = filename.expanded();
  }

  /**
  * Read the text data file - overloaded method
  *
  * @param fileString QString
  * @param data multi-dimensional vector of double
  *
  */
  void ReadTextDataFile(QString &fileString, vector<vector<double> > &data) {
    FileName filename(fileString);
    if(filename.isVersioned())
      filename = filename.highestVersion();
    if(!filename.fileExists()) {
      QString msg = fileString + " does not exist.";
      throw IException(IException::User, msg, _FILEINFO_);
    }
    TextFile file(filename.expanded());
    QString lineString;
    while(file.GetLine(lineString)) {
      vector<double> line;
      lineString = lineString.simplified().remove(QRegExp("^[ ,]*")).trimmed();

      QStringList lineTokens = lineString.split(QRegExp("[ ,]"), QString::SkipEmptyParts);
      foreach (QString value, lineTokens) {
        line.push_back(toDouble(value));
      }

      data.push_back(line);
    }

    fileString = filename.expanded();
  }
  // in.size should be 3144 total pizels including 6 channels
  void Cube2scImageStruct(Buffer &in){
    SDC_Image img = new SDC_Image();
    for(int i = 0; i < in.size(); i + ch.++){
      //channels per sensor
      for(int ch = 0; ch < 6; ch++){
        //pixels per channel
        for(int j = i; j < 524 + i; j++){
          img[ch].data.push_back(in[j]);
          if (j == jlj){
            img[ch].bias.push_back(in[j]);
          else
            img[ch].scene.push_back[in[j]];
          }


        }
      }
      
    }
  }

  struct SDC_Channel{
      DataPerChannel(){
        dataColumns = 524; 
        biasColumns = 8; 
        sceneColumns = 512;
      };  //default constructor
      int rows;
      int dataColumns;
      int biasColumns;
      int sceneColumns;
      vector<double> data;
      vector<double> bias;
      vector<double> scene;

      //DataPerChannel() {
      //  rows = r;
      //  dataColumns = 524; 
      //  biasColumns = 8; 
      //  sceneColumns = 512;
      //};
  };

  struct SDC_Image{
    //SDC_Image(){};
    //double lineTime;
    //string TDI_dir;
    //double gain;
    //double offset;
    //int nRows;
    //int channels;
    SDC_Channel channel[6];
    //vector<vector<DataPerChannel(nRows)>> SDC_Channel;
    //SDC_Image(){ SDC_Channel = {6};}
  };

  

 

  /*
     |-LineTime       [1x1] Numeric
%         |-FPAtemp        [1x1] Numeric
%         |-TDI_dir        [1x6] Char
%         |-Gain           [1x6] Numeric
%         |-Offset         [1x6] Numeric
%         |-SDC_Image      [1×6] Numeric
%             |--channel   [1×1] struct
%             |--data      [Mrows X 524] Numeric 
%             |--bias      [Mrows X 8] Numeric 
%             |--scene     [Mrows X 512] Numeric 
%-------------------------------------------------------------------   
% Last updated: 9/6/2022
%--------------------------------------------------------------------

% initialize
bias_pixel_avg = zeros(1,6);

% get the average bias pixel value for each channel i.e. 6 values
% SDC_Image.bias has the bias columns (3 to 10) for each channel
for channel_no = 0:5
bias_pixel_avg(channel_no+1) = mean2(I_struct.SDC_Image(channel_no+1).bias);
end

%Subtract the bias from scene pixels and bias pixels for each channel of the input image 
for channel_no = 0:5
I_bs(channel_no+1).scene = I_struct.SDC_Image(channel_no+1).scene - bias_pixel_avg(channel_no+1);
I_bs(channel_no+1).bias = I_struct.SDC_Image(channel_no+1).bias - bias_pixel_avg(channel_no+1);
end
  */
  void BiasPixelSubtract(Buffer &in) {
    
    double biasPixelAvg[6] = { 0 };
    for (int i = 0; i < sizeof(biasPixelAvg) / sizeof(int); i++){
      biasPixelAvg[i] = mean2(SDC_Image(i).bias);
    }
    SDC_Image img = new SDC_Image();
    //img.channel[1];
    for (int i = 0; i < img.channel.size(); i++){
      cout << img.channel.dataColuns << endl;
      cout << img.channel.sceneColumns << endl;
      cout << img.channel.biasColumns << endl;

    }



  }

  void RemoveMaskedOffset(Buffer &in) {
    
  }


  void CorrectDark(Buffer &in) {
    
  }

 
  void CorrectNonlinearity(Buffer &in) {
    
  }

  void CorrectFlatfield(Buffer &in) {
    
  }

  void RadiometricCalibration(Buffer &in) {
    
  }

  /**
  * This method returns an QString containing the path of an
  * LRO calibration directory
  *
  * @param calibrationType
  * @param calibrationDirectory Path of the calibration directory
  *
  * @internal
  *   @history 2020-01-06 Victor Silva - Added option for base calibration directory
  */
  void GetCalibrationDirectory(QString calibrationType, QString &calibrationDirectory) {
    PvlGroup &dataDir = Preference::Preferences().findGroup("DataDirectory");
    QString missionDir = (QString) dataDir["LRO"];
    if(calibrationType != "")
      calibrationType += "/";

    calibrationDirectory = missionDir + "/calibration/" + calibrationType;
  }

  void GetNearestDarkFile(QString fileString, QString &file) {
    
  }

  void GetNearestDarkFilePair(QString &fileString, QString &file0, QString &file1) {
    
  }

  /**
  * This method copies cube into vector
  * LRO calibration directory
  *
  * @param fileString QString pointer
  * @param data vector of double
  *
  */
  void CopyCubeIntoVector(QString &fileString, vector<double> &data) {
    Cube cube;
    FileName filename(fileString);
    if(filename.isVersioned())
      filename = filename.highestVersion();
    if(!filename.fileExists()) {
      QString msg = fileString + " does not exist.";
      throw IException(IException::User, msg, _FILEINFO_);
    }
    cube.open(filename.expanded());
    Brick brick(cube.sampleCount(), cube.lineCount(), cube.bandCount(), cube.pixelType());
    brick.SetBasePosition(1, 1, 1);
    cube.read(brick);
    data.clear();
    for(int i = 0; i < cube.sampleCount(); i++)
      data.push_back(brick[i]);

    fileString = filename.expanded();

    if(data.empty()){
      QString msg = "Copy from + " + fileString + " into vector failed.";
      throw IException(IException::User, msg, _FILEINFO_);
    }

  }
  /**
  * Allow special pixel types
  *
  * @param pixelValue double
  *
  * @return bool
  *
  */
  bool AllowedSpecialPixelType(double pixelValue) {
    bool result = false;
    result = result || IsHisPixel(pixelValue);
    result = result || IsLisPixel(pixelValue);
    result = result || IsHrsPixel(pixelValue);
    result = result || IsLrsPixel(pixelValue);
    return result;
  }

  /**
  * Get weighted time average for calculating pixel dark
  * average
  *
  * @param w0 double Weighted Time Average for dark file
  * @param w1 double Weighted time Average for dark file
  *
  */
  void GetWeightedDarkAverages() {

  }
}
