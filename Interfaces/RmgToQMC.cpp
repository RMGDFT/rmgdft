#include <iostream>
#include <string>
#include <fstream>
#include <sys/stat.h>
#include "XmlRep.h"
#include "WriteEshdf.h"

// RMG includes
#include "const.h"
#include "RmgTimer.h"
#include "RmgException.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "transition.h"
#include "Kpoint.h"
#include "InputKey.h"
#include "rmgthreads.h"

/* Main control structure which is declared extern in main.h so any module */
/* may access it.                                                        */
CONTROL ct;

/* PE control structure which is also declared extern in main.h */
PE_CONTROL pct;

std::unordered_map<std::string, InputKey *> ControlMap;

using namespace std;

void convertToVecStrings(int argc, char* argv[], std::vector<std::string>& vec) {
  for (int i = 1; i < argc; i++) {
    vec.push_back(argv[i]);
  }
}

inline bool file_exists (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

int main(int argc, char* argv[]) {

   // Set branch type and save argc and argv in control structure
   ct.rmg_branch = RMG_TO_QMC;
   ct.argc = argc;
   ct.argv = argv;

   // RMG setup functions will use multiple cores is this is not set
   setenv("OMP_NUM_THREADS", "1", true);
   setenv("OMP_RMG_THREADS", "1", true);

   vector<string> vecParams;
   convertToVecStrings(argc, argv, vecParams);
   string fname;
   for (int i = 0; i < vecParams.size(); i++) {
     // in principle check for other flags here
     fname = vecParams[i];
   }
   if (file_exists(fname) == false) {
     cout << "must specify a valid rmg input file name as an argument" << endl;
     exit(1);
   }

   /* Initialize all I/O including MPI group comms */
   /* Also reads control and pseudopotential files*/
   InitIo (argc, argv, ControlMap);
   MPI_Barrier(MPI_COMM_WORLD);

   /*Exit MPI */
   MPI_Finalize ();
   RmgTerminateThreads();
   exit(0);

  ifstream is(fname.c_str());
  xmlNode qboxSampleXml(&is, 0, true);


  string ofname = "qboxEshdf.h5";
  eshdfFile outFile(ofname);
  outFile.writeBoilerPlate("rmg");
  outFile.writeSupercell(qboxSampleXml);
  outFile.writeAtoms(qboxSampleXml);
  outFile.writeElectrons(qboxSampleXml);
}
