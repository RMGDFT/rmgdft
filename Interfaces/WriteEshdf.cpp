#include "WriteEshdf.h"
#include "FftContainer.h"
#include "HdfHelpers.h"
#include <sstream>
#include <map>
#include <fcntl.h>
#include <sys/stat.h>


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

using namespace std;
using namespace hdfHelper;

class kpoint {
public:
  double kx;
  double ky;
  double kz;
  kpoint() { kx = 0; ky = 0; kz = 0; }
  kpoint(double x, double y, double z) { kx = x; ky = y; kz = z; }
  kpoint(const kpoint& kp) { kx = kp.kx; ky = kp.ky; kz = kp.kz; }
  kpoint& operator=(const kpoint& kp) { kx = kp.kx; ky = kp.ky; kz = kp.kz; return *this; }
  bool operator==(const kpoint& kp) const {
    if ((std::abs(kx-kp.kx) < 1e-7) && (std::abs(ky-kp.ky) < 1e-7) && (std::abs(kz-kp.kz) < 1e-7)) {
      return true;
    }
    return false;
  }
  bool operator<(const kpoint& kp) const {
    if (abs(kx-kp.kx) > 1e-7) {
      return kx < kp.kx;
    } else if (abs(ky-kp.ky) > 1e-7) {
      return ky < kp.ky;
    } else if (abs(kz-kp.kz) > 1e-7) {
      return kz < kp.kz;
    }
    return false;
  }
};


eshdfFile::eshdfFile(const string& hdfFileName) {
  remove(hdfFileName.c_str());
  file = H5Fcreate(hdfFileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
}

eshdfFile::~eshdfFile() {
  H5Fclose(file);
}


int eshdfFile::wrapped(int i, int size) const {
  if (i < size/2) {
    return i;
  } else {
    return wrapped(i-size, size);
  }
}


// This handles all orbitals for one spin state and kpoint
void eshdfFile::handleSpinGroup(int kidx, int spin_idx, hid_t groupLoc, double& nocc, fftContainer& cont) {

  vector<double> eigvals;
  nocc = 0.0;

  for (int chIdx = 0; chIdx < ct.num_states; chIdx++)
  {
      //cout << "Working on state " << stateCounter << endl;
      stringstream statess;
      statess << "state_" << chIdx;
      hid_t state_group = makeHDFGroup(statess.str(), groupLoc);
      

      std::string wfname(ct.infile);
      wfname = wfname + "_spin" + std::to_string(spin_idx) +
                        "_kpt" + std::to_string(kidx) +
                        "_wf" + std::to_string(chIdx);

      double eig_val, wf_occ;
      readInEigFcn(wfname, eig_val, wf_occ, cont);

      // RMG Porting note - need to know QMCPACK units
      eigvals.push_back(eig_val);
      nocc += wf_occ;

      // write eigfcn to proper place
      hsize_t psig_dims[]={static_cast<hsize_t>(cont.fullSize),2};

      vector<double> temp;
      for (int i = 0; i < cont.fullSize; i++) {
	temp.push_back(cont.kspace[i][0]);
	temp.push_back(cont.kspace[i][1]);
      }
      writeNumsToHDF("psi_g", temp, state_group, 2, psig_dims);
  }
  writeNumsToHDF("number_of_states", ct.num_states, groupLoc);
  writeNumsToHDF("eigenvalues", eigvals, groupLoc);
}

void eshdfFile::readInEigFcn(std::string& wfname, double& eig_value, double& wf_occ, fftContainer& cont) {
  
    vector<double> values;
    OrbitalHeader H;

    int fhand = open(wfname.c_str(), O_RDWR, S_IREAD | S_IWRITE);
    if (fhand < 0) {
        rmg_printf("Can't open restart file %s", wfname.c_str());
        rmg_error_handler(__FILE__, __LINE__, "Terminating.");
    }
    size_t rsize = read (fhand, &H, sizeof(OrbitalHeader));
    if(rsize != sizeof(OrbitalHeader))
        rmg_error_handler (__FILE__,__LINE__,"error reading");

//    if((H.nx != (size_t)sizes_c[0]) || (H.ny != (size_t)sizes_c[1]) || (H.nz != (size_t)sizes_c[2])) {
//        rmg_printf("Grid size mismatch. %d  %d  %d  %lu  %lu  %lu", sizes_c[0], sizes_c[1], sizes_c[2], H.nx, H.ny, H.nz);
//        rmg_error_handler (__FILE__,__LINE__,"Grid size mismatch.");
//    }

    eig_value = H.eig;
    wf_occ = H.occ;

    int factor = 2;
    if(ct.is_gamma) factor = 1;
    double *tbuf = new double[cont.fullSize * factor];
    rsize = read (fhand, tbuf, cont.fullSize * sizeof(double) * factor);
    if(rsize != cont.fullSize * sizeof(double) * factor)
        rmg_error_handler (__FILE__,__LINE__,"problem reading wf file");
    close(fhand);
    for(int idx=0;idx < cont.fullSize;idx++) values.push_back(tbuf[idx]);
    delete [] tbuf;

  if (!ct.is_gamma) {
    int index = 0;
    for (int ix = 0; ix < cont.getNx(); ix++) {
      for (int iy = 0; iy < cont.getNy(); iy++) {
	for (int iz = 0; iz < cont.getNz(); iz++) {
          // RMG porting note -- C or Fortran storage order again?
	  const int qbx = cont.getQboxIndex(ix,iy,iz);
	  cont.rspace[index][0] = values[2*qbx];
	  cont.rspace[index][1] = values[2*qbx+1];
	  index++;
	}
      }
    }
  } else {
    int index = 0;
    for (int ix = 0; ix < cont.getNx(); ix++) {
      for (int iy = 0; iy < cont.getNy(); iy++) {
	for (int iz = 0; iz < cont.getNz(); iz++) {
	  const int qbx = cont.getQboxIndex(ix,iy,iz);
	  cont.rspace[index][0] = values[qbx];
	  cont.rspace[index][1] = 0.0;
	  index++;
	}
      }
    }
  }
  //cout << "in readInEigFcn, before fft, real space L2 norm = " << cont.getL2NormRS() << endl;
  cont.executeFFT();
  const double fixnorm = sqrt(Rmg_L.get_omega()) / static_cast<double>(cont.fullSize);
  cont.fixKsNorm(fixnorm);
  //cout << "in readInEigFcn, after fft, k space L2 norm = " << cont.getL2NormKS() << endl;
}

void eshdfFile::writeApplication(const string& appName, int major, int minor, int sub) {
  vector<int> version{major, minor, sub};

  hid_t h_app = makeHDFGroup("application", file);
  { 
    writeStringToHDF("code", appName, h_app);
    writeNumsToHDF("version", version, h_app);
  }
}


void eshdfFile::writeVersion() {
  vector<int> version{2,1,0};
  writeNumsToHDF("version", version, file);
}

void eshdfFile::writeCreator() {
  string nameStr("creator");
  vector<int> version{0,1,0};
  hid_t h_creator = makeHDFGroup("creator", file);
  {
    writeStringToHDF("program_name", "rmgConverter", h_creator);
    writeNumsToHDF("version", version, h_creator);
  }
}

void eshdfFile::writeFormat() {
  writeStringToHDF("format", "ES-HDF", file);
}

void eshdfFile::writeBoilerPlate(const string& appName) {
  // Fix this up later
  string desString("RMG");
  string versionStr("4.0.0");
  
  const int firstDotIdx = versionStr.find_first_of('.');
  const int secondDotIdx = versionStr.find_last_of('.');
  const int major = stoi(versionStr.substr(0,firstDotIdx));
  const int minor = stoi(versionStr.substr(firstDotIdx+1,secondDotIdx-firstDotIdx-1));
  const int sub = stoi(versionStr.substr(secondDotIdx+1));

  writeApplication(appName, major, minor, sub);
  writeVersion();
  writeCreator();
  writeFormat();
}

void eshdfFile::writeSupercell(void) {

  vector<double> ptvs;
  double temp;
  ptvs.push_back(Rmg_L.get_a0(0));
  ptvs.push_back(Rmg_L.get_a0(1));
  ptvs.push_back(Rmg_L.get_a0(2));
  ptvs.push_back(Rmg_L.get_a1(0));
  ptvs.push_back(Rmg_L.get_a1(1));
  ptvs.push_back(Rmg_L.get_a1(2));
  ptvs.push_back(Rmg_L.get_a2(0));
  ptvs.push_back(Rmg_L.get_a2(1));
  ptvs.push_back(Rmg_L.get_a2(2));

  // write the ptvs to the supercell group of the hdf file
  hid_t supercell_group = makeHDFGroup("supercell", file);
  hsize_t dims[]={3,3};
  writeNumsToHDF("primitive_vectors", ptvs, supercell_group, 2, dims);  
}


void eshdfFile::writeAtoms(void) {

  //make group
  hid_t atoms_group = makeHDFGroup("atoms", file);
  
  //map<string, int> SpeciesNameToInt;  
  //go through each species, extract: 
  //   atomic_number, mass, name, pseudopotential and valence_charge
  //   then write to a file, also set up mapping between species name and number
  for (int speciesNum = 0; speciesNum < ct.num_species; speciesNum++) {

    SPECIES *sp = &Species[speciesNum];
    // RMG porting note - Not used but is this the full name of the element?
    //string spName = species.getAttribute("name");
    //SpeciesNameToInt[spName] = speciesNum;
    //SpeciesNameToInt[spName] = speciesNum;
    string name(sp->atomic_symbol);
    string pseudopotential(sp->pseudo_filename);

    stringstream gname;
    gname << "species_" << speciesNum;
    hid_t species_group = makeHDFGroup(gname.str(), atoms_group);
    writeNumsToHDF("atomic_number", sp->atomic_number, species_group);
    writeNumsToHDF("mass", sp->atomic_mass, species_group);
    writeNumsToHDF("valence_charge", sp->zvalence, species_group);
    writeStringToHDF("name", name, species_group);
    writeStringToHDF("pseudopotential", pseudopotential, species_group);
    speciesNum++;
  }
  writeNumsToHDF("number_of_species", ct.num_species, atoms_group);
  
  // go through atoms and extract their position and type 
  std::vector<int> species_ids;
  std::vector<double> positions;
  // RMG porting note - are positions sequentially stored 3 doubles per atom?
  for (int i = 0; i < ct.num_ions; i++) {
      ION *iptr = &Atoms[i];
      species_ids.push_back(iptr->species);

//    RMG porting note -- this adds a triplet to the vector
//      atNode.getChild("position").getValue(positions);

      // RMG porting note - need to clarify if these are cell relative or absolute coordinates
      positions.push_back(iptr->crds[0]);
      positions.push_back(iptr->crds[1]);
      positions.push_back(iptr->crds[2]);
  }
  hsize_t dims[]={ct.num_ions,3};  
  writeNumsToHDF("positions", positions, atoms_group, 2, dims);
  writeNumsToHDF("species_ids", species_ids, atoms_group);
  writeNumsToHDF("number_of_atoms", ct.num_ions, atoms_group);
}

void eshdfFile::writeElectrons(void) {
  int nspin, nel;

  //wfnNode.getAttribute("nspin", nspin);
  nspin = (int)(ct.spin_flag + 1);

  //wfnNode.getAttribute("nel", nel);
  // RMG porting note - In RMG this is total not per spin state
  nel = (int)ct.nel;

  // RMG porting note, C or Fortran order?
  int nx = Rmg_G->get_NX_GRID(1);
  int ny = Rmg_G->get_NY_GRID(1);
  int nz = Rmg_G->get_NZ_GRID(1);

  fftContainer fftCont(nx,ny,nz);


  hid_t electrons_group = makeHDFGroup("electrons", file);
  writeNumsToHDF("number_of_kpoints", ct.num_kpts, electrons_group);
  writeNumsToHDF("number_of_spins", nspin, electrons_group);

  double avgNup = 0.0;
  double avgNdn = 0.0;
  // go through kpt by kpt and write 
  for (int i = 0; i < ct.num_kpts; i++) {
    stringstream kptElemName;
    kptElemName << "kpoint_" << i;
    hid_t kpt_group = makeHDFGroup(kptElemName.str(), electrons_group);
    vector<double> kvector;
    kvector.push_back(ct.kp[i].kpt[0]);
    kvector.push_back(ct.kp[i].kpt[1]);
    kvector.push_back(ct.kp[i].kpt[2]);
    //writeNumsToHDF("numsym", 1, kpt_group);
    //writeNumsToHDF("symgroup", 1, kpt_group);
    writeNumsToHDF("reduced_k", kvector, kpt_group);

//    writeNumsToHDF("weight", 1.0/static_cast<double>(kpts.size()), kpt_group); 
    // RMG porting note - all kpoints have the same weight is this right?
    writeNumsToHDF("weight", ct.kp[i].kweight, kpt_group); 

    if (i == 0) {
      // figure out the order of the g-vectors
      // write them to an integer array and then put it in the element
      vector<int> gvectors;
      for (int ix = 0; ix < nx; ix++) {
	for (int iy = 0; iy < ny; iy++) {
	  for (int iz = 0; iz < nz; iz++) {
	    gvectors.push_back(wrapped(ix,nx));
	    gvectors.push_back(wrapped(iy,ny));
	    gvectors.push_back(wrapped(iz,nz));
	  }
	}
      }    
      
      hsize_t dims[]={static_cast<hsize_t>(nx*ny*nz),3};
      writeNumsToHDF("gvectors", gvectors, kpt_group, 2, dims);
      writeNumsToHDF("number_of_gvectors", nx*ny*nz, kpt_group);
    }


    // here is where we will read in both the occupations for the 
    // kpoint (species depenent if needed)
    // and also the wavefunction (in real space)
    // fourier transform it and write it into the appropriate hdf element
    double nup = 0;
    double ndn = 0;

    hid_t up_spin_group = makeHDFGroup("spin_0", kpt_group);
    handleSpinGroup(i, 0, up_spin_group, nup, fftCont);

    if (nspin == 2) {
      hid_t dn_spin_group = makeHDFGroup("spin_1", kpt_group);
      handleSpinGroup(i, nspin-1, dn_spin_group, ndn, fftCont);
    } else {
      ndn = nup;
    }

    avgNup += nup / static_cast<double>(ct.num_kpts);
    avgNdn += ndn / static_cast<double>(ct.num_kpts);
  }
  vector<int> nels;
  nels.push_back(static_cast<int>(std::floor(avgNup+0.1)));
  nels.push_back(static_cast<int>(std::floor(avgNdn+0.1)));
  writeNumsToHDF("number_of_electrons", nels, electrons_group);
}

void WriteQmcpackRestart(std::string& name)
{
   string ofname = name + ".h5";
   eshdfFile outFile(ofname);
   outFile.writeBoilerPlate("rmg");
   outFile.writeSupercell();
   outFile.writeAtoms();
   outFile.writeElectrons();
}
