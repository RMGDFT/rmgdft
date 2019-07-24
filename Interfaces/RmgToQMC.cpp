#include <iostream>
#include <string>
#include <fstream>
#include <sys/stat.h>
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

// Pointer to Kpoint class arrays for gamma and non-gamma
Kpoint<double> **Kptr_g;
Kpoint<std::complex<double> > **Kptr_c;


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

    /* Initialize some k-point stuff */
    Kptr_g = new Kpoint<double> * [ct.num_kpts_pe];
    Kptr_c = new Kpoint<std::complex<double> > * [ct.num_kpts_pe];

    ct.is_gamma = true;
    for (int kpt = 0; kpt < ct.num_kpts; kpt++) {
        double v1, v2, v3;
        v1 = twoPI * ct.kp[kpt].kpt[0] / Rmg_L.get_xside();
        v2 = twoPI * ct.kp[kpt].kpt[1] / Rmg_L.get_yside();
        v3 = twoPI * ct.kp[kpt].kpt[2] / Rmg_L.get_zside();

        ct.kp[kpt].kvec[0] = v1;
        ct.kp[kpt].kvec[1] = v2;
        ct.kp[kpt].kvec[2] = v3;
        ct.kp[kpt].kmag = v1 * v1 + v2 * v2 + v3 * v3;

        if(ct.kp[kpt].kmag != 0.0) ct.is_gamma = false;
    }

    if(!ct.is_gamma)
    {
        ct.is_use_symmetry = 1;
    }

    for (int kpt = 0; kpt < ct.num_kpts_pe; kpt++)
    {

        int kpt1 = kpt + pct.kstart;
        if(ct.is_gamma) {

            // Gamma point
            Kptr_g[kpt] = new Kpoint<double> (ct.kp[kpt1].kpt, ct.kp[kpt1].kweight, kpt, pct.grid_comm, Rmg_G, Rmg_T, &Rmg_L, ControlMap);

        }
        else {

            // General case
            Kptr_c[kpt] = new Kpoint<std::complex<double>> (ct.kp[kpt1].kpt, ct.kp[kpt1].kweight, kpt, pct.grid_comm, Rmg_G, Rmg_T, &Rmg_L, ControlMap);

        }
        ct.kp[kpt].kidx = kpt;
    }


   string ofname = "rmgEshdf.h5";
   eshdfFile outFile(ofname);
   outFile.writeBoilerPlate("rmg");
   outFile.writeSupercell();
   outFile.writeAtoms();
   outFile.writeElectrons();

   /*Exit MPI */
   MPI_Finalize ();
   RmgTerminateThreads();
   exit(0);
}
