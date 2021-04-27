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
#include "FiniteDiff.h"
#include "prototypes_on.h"


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
    //        rmg_error_handler (__FILE__,__LINE__,"Grid size mismatch.");//    }

    eig_value = H.eig;
    wf_occ = H.occ;

    int factor = 2;
    if(ct.is_gamma) factor = 1;
    double *tbuf = new double[cont.fullSize * factor];
    rsize = read (fhand, tbuf, cont.fullSize * sizeof(double) * factor);
    if(rsize != cont.fullSize * sizeof(double) * factor)
        rmg_error_handler (__FILE__,__LINE__,"problem reading wf file");
    close(fhand);
    for(int idx=0;idx < cont.fullSize*factor;idx++) values.push_back(tbuf[idx]);
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
    hsize_t dims[]={(hsize_t)ct.num_ions,3};  
    writeNumsToHDF("positions", positions, atoms_group, 2, dims);
    writeNumsToHDF("species_ids", species_ids, atoms_group);
    writeNumsToHDF("number_of_atoms", ct.num_ions, atoms_group);
}

void eshdfFile::writeElectrons(void) {
    int nspin;

    //wfnNode.getAttribute("nspin", nspin);
    nspin = ct.nspin;

    //wfnNode.getAttribute("nel", nel);
    // RMG porting note - In RMG this is total not per spin state

    // RMG porting note, C or Fortran order?
    int nx = Rmg_G->get_NX_GRID(1);
    int ny = Rmg_G->get_NY_GRID(1);
    int nz = Rmg_G->get_NZ_GRID(1);

    fftContainer fftCont(nx,ny,nz);


    hid_t electrons_group = makeHDFGroup("electrons", file);

    handleRho(electrons_group);
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

void WriteForAFQMC(int ns_occ, int Nchol, int Nup, int Ndown, 
        std::vector<double> eigs, std::vector<double> &CholVec, std::vector<double> &Hcore)
{
    string ofname = "afqmc_rmg.h5";
    eshdfFile outFile(ofname);
    outFile.writeCholVec(ns_occ, Nchol, Nup, Ndown, eigs, CholVec, Hcore);

}

void eshdfFile::writeCholVec(int ns_occ, int Nchol, int Nup, int Ndown, 
        std::vector<double> eigs, std::vector<double> CholVec, std::vector<double> Hcore)
{
    hid_t hamiltonian_group = makeHDFGroup("Hamiltonian", file);
    hid_t chol_group = makeHDFGroup("DenseFactorized", hamiltonian_group);
    hsize_t chv_dims[]={static_cast<hsize_t>(ns_occ * ns_occ), static_cast<hsize_t>(Nchol)};
    writeNumsToHDF("L", CholVec, chol_group, 2, chv_dims);
    std::vector<int> dims;
    dims.resize(8, 0);
    dims[3] = ns_occ;
    dims[4] = Nup;
    dims[5] = Ndown;
    dims[7] = Nchol;
    writeNumsToHDF("dims", dims, hamiltonian_group);

    std::vector<double> Energies;

    hsize_t h_dims[]={static_cast<hsize_t>(ns_occ),static_cast<hsize_t>(ns_occ)};
    writeNumsToHDF("hcore", Hcore, hamiltonian_group, 2, h_dims);
    Energies.push_back(ct.II);
    Energies.push_back(0.0);
    writeNumsToHDF("Energies", Energies, hamiltonian_group);



}

void WriteQmcpackRestartLocalized(std::string& name)
{
    string ofname = "RMG_localized.h5";
    eshdfFile outFile(ofname);
    outFile.writeBoilerPlate("rmg_on");
    outFile.writeSupercell();
    outFile.writeAtoms();
    outFile.writeLocalizedOrbitals();
}

void eshdfFile::writeLocalizedOrbitals() {
    //#include "init_var.h"
    int nspin;

    //wfnNode.getAttribute("nspin", nspin);
    nspin = ct.nspin;

    //wfnNode.getAttribute("nel", nel);
    // RMG porting note - In RMG this is total not per spin state

    // RMG porting note, C or Fortran order?

    hid_t electrons_group = makeHDFGroup("basisset", file);
    writeNumsToHDF("NbElements", -1, electrons_group);
    writeNumsToHDF("NbCenters", ct.num_ions, electrons_group);
    writeNumsToHDF("number_of_spins", nspin, electrons_group);
    writeNumsToHDF("number_of_orbitals", ct.num_states, electrons_group);
    std::vector<int> num_orb_centers;
    num_orb_centers.resize(ct.num_ions);
    for(int ion = 0; ion < (int)Atoms.size(); ion++)
    {
        num_orb_centers[ion] = Atoms[ion].num_orbitals;
    }
    writeNumsToHDF("NumOrbCenters", num_orb_centers, electrons_group);
    hid_t kpt_group = makeHDFGroup("KPTS_0", file);
    for(int ispin = 0; ispin < nspin; ispin++)
    {
        vector<double> Cij;
        Cij.resize(ct.num_states * ct.num_states);
        std::string wfname(ct.outfile);
        wfname = wfname + "_spin" + std::to_string(ispin) + "_Cij";

        int fhand = open(wfname.c_str(), O_RDWR);
        if (fhand < 0)
        {
            printf("\n  unable to open file %s", wfname.c_str() );
            exit(0);
        }

        read(fhand, Cij.data(),ct.num_states * ct.num_states * sizeof(double));
        close(fhand);
        std::string Cij_name = "eigenset_" + std::to_string(ispin);
        hsize_t Cij_dims[]={static_cast<hsize_t>(ct.num_states),static_cast<hsize_t>(ct.num_states)};
        writeNumsToHDF(Cij_name, Cij, kpt_group, 2, Cij_dims);
    }



    double nup = 0;
    double ndn = 0;

    hid_t up_spin_group = makeHDFGroup("spin_0", electrons_group);
    handleSpinGroup_ON(0, up_spin_group, nup);

    if (nspin == 2) {
        hid_t dn_spin_group = makeHDFGroup("spin_1", electrons_group);
        handleSpinGroup_ON(nspin-1, dn_spin_group, ndn);
    } else {
        ndn = nup;
    }

    vector<int> nels;
    nels.push_back(static_cast<int>(std::floor(nup+0.1)));
    nels.push_back(static_cast<int>(std::floor(ndn+0.1)));
    writeNumsToHDF("number_of_electrons", nels, electrons_group);
}

void eshdfFile::handleSpinGroup_ON(int spin_idx, hid_t groupLoc, double& nocc) {

    //#include "init_var.h"

    vector<double> eigvals;
    nocc = 0.0;

    std::vector<int> num_orb_centers;
    num_orb_centers.resize(ct.num_ions);
    for(int ion = 0; ion < (int)Atoms.size(); ion++)
    {
        num_orb_centers[ion] = Atoms[ion].num_orbitals;
    }


    vector<double> phi, phi_x, phi_y, phi_z, phi_L;
    phi.resize(ct.states[0].size);
    phi_x.resize(ct.states[0].size);
    phi_y.resize(ct.states[0].size);
    phi_z.resize(ct.states[0].size);
    phi_L.resize(ct.states[0].size);

    FiniteDiff FD(&Rmg_L);
    double hxgrid = Rmg_G->get_hxgrid(1);
    double hygrid = Rmg_G->get_hygrid(1);
    double hzgrid = Rmg_G->get_hzgrid(1);

    int order = 8;
    int item = (ct.max_orbit_nx+order) *(ct.max_orbit_ny+order) *(ct.max_orbit_nz+order);
    double *orbital_border = new double[item];

    std::vector<double> hgrid(3);
    std::vector<int> grid_start(3), grid_dim(3);
    hgrid[0] = Rmg_G->get_hxgrid(1) * Rmg_L.get_xside() ;
    hgrid[1] = Rmg_G->get_hygrid(1) * Rmg_L.get_yside() ;
    hgrid[2] = Rmg_G->get_hzgrid(1) * Rmg_L.get_zside() ;
    int chIdx = 0;

    for(int ic = 0; ic < (int) num_orb_centers.size(); ic++)
    {
        std::string center_group = "orbital_group_" + std::to_string(ic);
        hid_t orbital_group = makeHDFGroup(center_group.c_str(), groupLoc);

        writeNumsToHDF("hgrid", hgrid, orbital_group);

        grid_start[0] = ct.states[chIdx].ixmin;
        grid_start[1] = ct.states[chIdx].iymin;
        grid_start[2] = ct.states[chIdx].izmin;
        grid_dim[0] = ct.states[chIdx].orbit_nx;
        grid_dim[1] = ct.states[chIdx].orbit_ny;
        grid_dim[2] = ct.states[chIdx].orbit_nz;
        writeNumsToHDF("grid_start", grid_start, orbital_group);
        writeNumsToHDF("grid_dim", grid_dim, orbital_group);
        writeNumsToHDF("NumOrbThiscenter", num_orb_centers[ic], orbital_group);

        for (int st = 0; st <num_orb_centers[ic]; st++)
        {
            //cout << "Working on state " << stateCounter << endl;

            std::string wfname(ct.outfile);
            wfname = wfname + "_spin" + std::to_string(spin_idx) +
                ".orbit_" + std::to_string(chIdx+st);

            int fhand = open(wfname.c_str(), O_RDWR);
            if (fhand < 0)
            {
                printf("\n  unable to open file %s", wfname.c_str() );
                fflush(NULL);
                exit(0);
            }

            read(fhand, phi.data(), ct.states[chIdx].size * sizeof(double));
            close(fhand);

            int dimx = ct.states[chIdx].orbit_nx;
            int dimy = ct.states[chIdx].orbit_ny;
            int dimz = ct.states[chIdx].orbit_nz;

            int ix, iy, iz;
            int incx, incy, incxs, incys;
            incys = dimz;
            incxs = dimy * dimz;
            incx = (dimy + order) * (dimz + order);
            incy = dimz + order;

            for (ix = 0; ix < dimx + order; ix++)
                for (iy = 0; iy < dimy + order; iy++)
                    for (iz = 0; iz < dimz + order; iz++)
                        orbital_border[ix * incx + iy * incy + iz] = 0.0;

            /* Load up original values from pg  */

            for (ix = 0; ix < dimx; ix++)
                for (iy = 0; iy < dimy; iy++)
                    for (iz = 0; iz < dimz; iz++)
                        orbital_border[(ix + order/2) * incx + (iy + order/2) * incy + iz + order/2] = phi[ix * incxs + iy * incys + iz];



            FD.app8_del2 (orbital_border, phi_L.data(), dimx, dimy, dimz, hxgrid, hygrid, hzgrid);
            FD.app_gradient_eighth (orbital_border, phi_x.data(), phi_y.data(), phi_z.data(),
                    dimx, dimy, dimz, hxgrid, hygrid, hzgrid);

            std::string orb = "orbital_" + std::to_string(st);
            std::string orb_x = "orbital_x_" + std::to_string(st);
            std::string orb_y = "orbital_y_" + std::to_string(st);
            std::string orb_z = "orbital_z_" + std::to_string(st);
            std::string orb_L = "orbital_L_" + std::to_string(st);
            writeNumsToHDF(orb , phi, orbital_group);
            writeNumsToHDF(orb_x , phi_x, orbital_group);
            writeNumsToHDF(orb_y , phi_y, orbital_group);
            writeNumsToHDF(orb_z , phi_z, orbital_group);
            writeNumsToHDF(orb_L , phi_L, orbital_group);
        }
        chIdx += num_orb_centers[ic];
    }
    delete [] orbital_border;
    //  writeNumsToHDF("eigenvalues", eigvals, groupLoc);
}

// This handles charge density
void eshdfFile::handleRho(hid_t groupLoc) {

    hid_t density_group = makeHDFGroup("density", groupLoc);

    // RMG porting note, C or Fortran order?
    int grid_ratio = Rmg_G->default_FG_RATIO;
    int nx = Rmg_G->get_NX_GRID(grid_ratio);
    int ny = Rmg_G->get_NY_GRID(grid_ratio);
    int nz = Rmg_G->get_NZ_GRID(grid_ratio);


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
    writeNumsToHDF("gvectors", gvectors, density_group, 2, dims);
    writeNumsToHDF("number_of_gvectors", nx*ny*nz, density_group);
    vector<int> mesh;
    mesh.push_back(nx);
    mesh.push_back(ny);
    mesh.push_back(nz);
    writeNumsToHDF("mesh", mesh, density_group);

    fftContainer rho_fftCont(nx,ny,nz);

    std::string rhofname(ct.infile);
    rhofname = rhofname + ".rho";

    int fhand = open(rhofname.c_str(), O_RDWR, S_IREAD | S_IWRITE);
    if (fhand < 0) {
        rmg_printf("Can't open restart file %s", rhofname.c_str());
        rmg_error_handler(__FILE__, __LINE__, "Terminating.");
    }

    int factor = ct.nspin;
    double *rho_buf = new double[rho_fftCont.fullSize * factor];
    size_t rsize = read (fhand, rho_buf, rho_fftCont.fullSize * sizeof(double) * factor);
    if(rsize != rho_fftCont.fullSize * sizeof(double) * factor)
        rmg_error_handler (__FILE__,__LINE__,"problem reading rho file");
    close(fhand);

    // write rho to proper place
    hsize_t rhog_dims[]={static_cast<hsize_t>(rho_fftCont.fullSize),2};

    for(int ispin = 0; ispin < ct.nspin; ispin++) {

        std::string spin_str = "spin_" + std::to_string(ispin);

        hid_t spin_density_group = makeHDFGroup(spin_str.c_str(), density_group);
        double *values = &rho_buf[ispin * rho_fftCont.fullSize];
        int index = 0;
        for (int ix = 0; ix < rho_fftCont.getNx(); ix++) {
            for (int iy = 0; iy < rho_fftCont.getNy(); iy++) {
                for (int iz = 0; iz < rho_fftCont.getNz(); iz++) {
                    const int qbx = rho_fftCont.getQboxIndex(ix,iy,iz);
                    rho_fftCont.rspace[index][0] = values[qbx];
                    rho_fftCont.rspace[index][1] = 0.0;
                    index++;
                }
            }
        }

        //cout << " before fft in rho, real space L2 norm = " << rho_fftCont.getL2NormRS() << endl;
        rho_fftCont.executeFFT();
        const double fixnorm = Rmg_L.get_omega() / static_cast<double>(rho_fftCont.fullSize);
        rho_fftCont.fixKsNorm(fixnorm);
        //cout << " fixnorm = " << fixnorm << endl;
        //cout << " after fft in rho, k space L2 norm = " << rho_fftCont.getL2NormKS() << endl;


        vector<double> temp;
        for (int i = 0; i < rho_fftCont.fullSize; i++) {
            temp.push_back(rho_fftCont.kspace[i][0]);
            temp.push_back(rho_fftCont.kspace[i][1]);
        }
        writeNumsToHDF("density_g", temp, spin_density_group, 2, rhog_dims);
    }
}

