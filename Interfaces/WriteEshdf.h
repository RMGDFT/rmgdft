#ifndef WRITE_ESHDF_H
#define WRITE_ESHDF_H
#include "hdf5.h"
#include <string>
#include <vector>
#include <cmath>

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

void WriteQmcpackRestart(std::string& name);
void WriteQmcpackRestartLocalized(std::string& name);
void WriteForAFQMC(int ns_occ, int Nchol, int Nup, int Ndown, 
        std::vector<double> eigs, std::vector<double> &CholVec);

class fftContainer;

class eshdfFile {
private:
  hid_t file;
  herr_t error;

  int wrapped(int i, int size) const;
  void writeApplication(const std::string& appName, int major, int minor, int sub);
  void writeCreator();
  void writeFormat();
  void writeVersion();

  void readInEigFcn(std::string& wfname, double& eig_value, double& wf_occ, fftContainer& cont);
  void handleSpinGroup(int kidx, int spin_idx, hid_t groupLoc, double& nocc, fftContainer& cont);
  void handleSpinGroup_ON(int spin_idx, hid_t groupLoc, double& nocc);

  eshdfFile(const eshdfFile& f); // no copy constructor
  eshdfFile& operator=(const eshdfFile& f); // operator= not allowed
public:
  eshdfFile(const std::string& hdfFileName);
  ~eshdfFile();

  void writeBoilerPlate(const std::string& appName);
  void writeSupercell(void);
  void writeAtoms(void);
  void writeElectrons(void);
  void writeLocalizedOrbitals();
  void writeCholVec(int, int, int, int, std::vector<double> eigs, std::vector<double> CholVec);
};

#endif
