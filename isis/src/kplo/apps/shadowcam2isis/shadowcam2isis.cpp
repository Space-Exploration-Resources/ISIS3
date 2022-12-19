
#include "ProcessImport.h"
#include "UserInterface.h"
#include "FileName.h"
#include "IException.h"
#include "IString.h"
#include "Pvl.h"
#include "OriginalLabel.h"
#include "History.h"
#include "LineManager.h"
#include "Application.h"

#include <vector>

#include "shadowcam2isis.h"

#define MAX_INPUT_VALUE 4095

using namespace std;

namespace Isis {

  //ProcessImportPds importPDS;

  static void ResetGlobals();
  static void Import(Buffer &buf);

  // Global variables for processing functions
  Cube *g_ocube;
  QString g_imgExt = "cub";
  std::vector<double> g_xterm {0, 32, 136, 543, 2207};
  std::vector<double> g_bterm {0, 0, 8, 25, 59, 128, 0};
  //bool g_flip = false;

  void shadowcam2isis(UserInterface &ui) {
    ResetGlobals();

    FileName fileName = ui.GetFileName("FROM");
    ProcessImport p;

    if (fileName.removeExtension().addExtension(g_imgExt).fileExists()) {
        p.SetInputFile(fileName.removeExtension().addExtension(g_imgExt).expanded());
    }
    else {
      QString msg = "Cannot find image file for [" + fileName.name() + "]. Confirm that the "
          ".cub file exists and is located in the same directory.";
        throw IException(IException::User, msg, _FILEINFO_);
    }
   
    
    CubeAttributeOutput &outAtt = ui.GetOutputAttribute("TO");

    g_ocube = new Cube();
    g_ocube->setByteOrder(outAtt.byteOrder());
    g_ocube->setFormat(outAtt.fileFormat());
    g_ocube->setMinMax((double) VALID_MIN2, (double) VALID_MAX2);
 
    g_ocube->setPixelType(Isis::Real);
    g_ocube->create(ui.GetCubeName("TO"));
    // Do 8 bit to 12 bit conversion
    // And if NAC-R, flip the frame
    p.StartProcess(Import);
    p.EndProcess();

    
    g_ocube->close();
    delete g_ocube;
  }
  

  void Import(Buffer &in) {
    //LineManager outLines(*g_ocube);
    //outLines.SetLine(1,1);
    //Buffer buf(20, 20, 1, g_ocube->pixelType());
    cout << "Import begins." << endl;
    // Do the decompanding
    for(int pixin = 0; pixin < 20; pixin++) {

      cout << "pixel: ";
      cout << pixin << endl;

    }
  }

  void ResetGlobals() {
    g_ocube = NULL;
    g_xterm.clear();
    g_bterm.clear();
    //g_flip = false;
  }
}