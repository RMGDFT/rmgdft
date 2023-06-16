
#include <exception>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
#include <boost/math/special_functions/erf.hpp>


#include "const.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "MapElements.h"
#include "RmgException.h"
#include "transition.h"
#include "InternalPseudo.h"
#include "Atomic.h"

#include "rmg_mangling.h"
#define sratom       RMG_FC_MODULE(haesratom,sratom,mod_HAESRATOM,SRATOM)

extern "C" void sratom(int *na, int *la, double *ea, double *fa,
                       double *rpk, int *nc, int *ncv, int *it,
                       double *rhoc, double *rho, double *rr, double *vi,
                       double *zz, int *mmax, int *iexc, double *etot, int *ierr, double *uu);

double bparm[12] = {
-3.6442293856e-01,
-1.9653418982e-01,
-1.3433604753e-01,
-1.0200558466e-01,
-8.2208091118e-02,
-6.8842555167e-02,
-5.9213652850e-02,
-5.1947028250e-02,
-4.6268559218e-02,
-4.1708913494e-02,
-3.7967227308e-02,
-3.4841573775e-02,
};

double lrc;

// For Z=1
double VofR(double r, int a_in)
{
  double V = -0.5;
  double b = bparm[a_in-1];
  double a = (double)a_in;
  double hprime = -boost::math::erf (a*r) -
                  2.0 * (a*a*b + a / sqrt(PI)) * r * exp(-a*a*r*r);
  double hdprime = (-2.0*a*a*b - 4.0*a/sqrt(PI) + (4.0*a*a*a*a*b + 4.0*a*a*a/sqrt(PI))*r*r) *
                   exp(-a*a*r*r);
  V += hprime/r + hprime*hprime/2.0 + hdprime/2.0;
  return V;
}

double VofZ(int Z, double r, int a)
{
  double Zr = r*(double)Z;
  return VofR(Zr, a) * (double)(Z*Z);
}
extern "C" double vofz_(double *Z, double *r)
{
  double r1 = r[0];
  double z1 = Z[0];
  double Zr = r1 * z1;
  double V0 = VofR(z1*LOGGRID_START, ct.all_electron_parm) * (z1 * z1);
  double V1 = VofR(Zr, ct.all_electron_parm) * (z1 * z1);
  double rc = fabs(2.0 * z1 / sqrt(PI) / V0);

  V1 = V1 + z1 * boost::math::erf (r1 / rc) / r1;
  return V1;
}


// Generates an all electron potential using the formula developed by F. Gygi
// All-Electron Plane-Wave Electronic Structure Calculations
// Journal of Chemical Theory and Computation 2023 19 (4), 1300-1309
// DOI: 10.1021/acs.jctc.2c01191
void LoadAllElectronPseudo(SPECIES *sp)
{
    Atomic A;
    std::stringstream ss; 
    double_2d_array ddd0;  // Used to read in the PP_DIJ
    double_2d_array qqq;   // Used to read in the norms of the augmentation functions (PP_Q)

    sp->max_l = 0;

    // Atomic symbol, mass, number and zvalence and mesh size
    // Atomic symbol is set in SP by higher level routines
    // Maybe check symbols here
    sp->atomic_mass = GetAtomicMass(sp->atomic_symbol);
    sp->atomic_number = GetAtomicNumber(sp->atomic_symbol);
    sp->zvalence = sp->atomic_number;  // not really of course but treated that way internally
    if(sp->zvalence > ct.max_zvalence) ct.max_zvalence = sp->zvalence;

    sp->generated = std::string("Generated internally using procedure of F.Gygi\nJournal of Chemical Theory and Computation 2023 19 (4), 1300-1309");
    sp->author = std::string("");

    // Store UPF functional string for later processing
    sp->functional = std::string("PBE");

    // Use our internal radial mesh from Atomic
    sp->rg_points = MAX_LOGGRID;
    double *t_r = A.GetRgrid();

    sp->gtype = LOG_GRID;
    sp->r = new double[sp->rg_points];
    for(int i = 0;i < sp->rg_points;i++) {
        sp->r[i] = t_r[i];
    }

    sp->kkbeta = sp->rg_points;
   
    // Get the type of pseudopotential
    sp->is_norm_conserving = true;
    sp->is_spinorb = false;

    // No core corrections!
    sp->nlccflag = false;

    // Determine log mesh parameters directly from the mesh
    sp->aa = (sp->r[0] * sp->r[0]) / (sp->r[1] - 2 * sp->r[0]);
    sp->bb = log (sp->r[1] / sp->r[0] - 1);

    // Generate RAB
    sp->rab = new double[sp->rg_points];
    sp->rab[0] = sp->r[0];
    for(int i = 1;i < sp->rg_points;i++) sp->rab[i] = (sp->r[i] - sp->r[i-1]);

    // All electron potential
    sp->vloc0 = new double[sp->rg_points];

    // Get into our internal units
    for(int ix = 0;ix < sp->rg_points;ix++)
    {
        sp->vloc0[ix] = VofZ(sp->atomic_number, sp->r[ix], ct.all_electron_parm);
    }

    // Get the l-value for the local potential if present
    sp->local = 0;

    // Atomic charge density
    sp->atomic_rho = new double[sp->rg_points]();

    // Number of atomic orbitals
    sp->num_atomic_waves = GetNumberOrbitalsL(sp->atomic_symbol);
    sp->num_atomic_waves_m = GetNumberOrbitalsM(sp->atomic_symbol);
    std::vector<int> pqn;
    sp->atomic_wave.resize(sp->num_atomic_waves);
    SetupAllElectonOrbitals(sp->atomic_symbol,
                            pqn,
                            sp->atomic_wave_l,
                            sp->atomic_wave_j,
                            sp->atomic_wave_oc,
                            sp->atomic_wave_energy,
                            sp->aradius,
                            sp->atomic_wave_label);

    for(int iwf = 0;iwf < sp->num_atomic_waves;iwf++)
    {
        sp->atomic_wave[iwf] = new double[sp->rg_points]();
    }
    ct.max_orbitals = std::max(ct.max_orbitals, sp->num_atomic_waves_m);
    sp->num_orbitals = sp->num_atomic_waves_m;

    int ierr;
    int iexc = 4;  // PBE
    int nc = 0;
    int ncv = sp->num_atomic_waves;
    int its;
    double aenergy;
    double *rpk = new double[sp->num_atomic_waves];
    double *vi = new double[sp->rg_points];
    double *rhoc = new double[sp->rg_points];
    double *uu = new double[sp->rg_points * sp->num_atomic_waves];
    sratom(pqn.data(), sp->atomic_wave_l.data(), sp->atomic_wave_energy.data(),
           sp->atomic_wave_oc.data(), rpk, &nc, &ncv, &its,
                       rhoc, sp->atomic_rho, sp->r, vi,
                       &sp->zvalence, &sp->rg_points, &iexc, &aenergy, &ierr, uu);
printf("IIII  %d  %d\n", its, ierr);
    // Normalize rho for 3D grid
    for(int idx = 0;idx < sp->rg_points;idx++)
    {
        sp->atomic_rho[idx] /= 4.0*PI;
    }

    for(int i = 0;i < sp->num_atomic_waves;i++)
    {
        printf("BBBB  %f  %f\n", sp->atomic_wave_energy[i], aenergy);
        for(int idx = 0;idx < sp->rg_points;idx++)
        {
            double ut = uu[i*sp->rg_points+idx] / sp->r[idx];
            //printf("RRRR  %14.8e  %14.8e\n", sp->r[idx], ut);
            sp->atomic_wave[i][idx] = ut;
        }
    }

    delete [] uu;
    delete [] rhoc;
    delete [] vi;
    delete [] rpk;

    // Number of projectors
    sp->nbeta = 0;
    sp->is_ddd_diagonal = true;
    sp->nqf=0;
    sp->nlc=0;


    // Set the maximum number of non-local projecters needed
    ct.max_nl = 0;

    // Optional stuff next
    sp->description = std::string("Regularized norm conserving all electron.");

    // Stuff not present in the UPF format that RMG requires. 
    // We need to find a consistent way of automatically setting these.
    sp->rc = fabs(2.0 * sp->zvalence / sqrt(PI) / sp->vloc0[0]);
    sp->rc = std::max(sp->rc, 0.01);
    sp->rc = std::min(sp->rc, 1.5);
    sp->lradius = 8.5;
    sp->gwidth = 8.0;
    sp->rwidth = 15.0; 
    sp->agwidth = 10.0;
    sp->arwidth = 25.0;

    // Leftover initializations
    sp->mill_radius = 9.0;

    // Finally adjust sp->rg_points to skip the high end
    int iend = sp->rg_points - 1;
    for(int i = 0;i < sp->rg_points;i++) {
        iend = i;
        if(sp->r[i] > 12.0) break;
    }

    sp->rg_points = iend + 1;
}

