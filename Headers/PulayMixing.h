#ifndef RMG_Pulay_H
#define RMG_Pulay_H 1

#include<functional>
#include<complex>

class PulayMixing {

private:

    double mix_first, beta;
    MPI_Comm comm;
    size_t Nsize;
    int pulay_order, refresh_steps;
    int step;
    int max_order = 10;
    double *hist;
    double *res_hist;
    std::complex<double> *res_histG;
    double *A_mat;
    std::vector<double*> hist_ptr;
    std::vector<double*> res_hist_ptr;
    std::vector<std::complex<double>*> res_histG_ptr;
    std::function<void(double*, int)> Precond;
    bool need_precond;
    int nstates;

    bool Gspace = false;
    bool drho_pre = false;
    double q0, q1;
    // q0 for Kerker mixing  1.5 A^-1
//    double ktf = 1.5/0.529177;  // in unit of au^-1
    double ktf = 0.529177;  // in unit of au^-1
    std::complex<double> *c_fm = NULL;

public:

    PulayMixing(size_t Nsize, int pulay_order, int refresh_steps, double mix_init, double beta, MPI_Comm comm);
    ~PulayMixing(void);
    void Mixing(double *xm, double *fm);
    void Mixing_rhoG(double *xm, double *fm);
    void SetPrecond(std::function<void(double*, int)> precon);
    void SetNstates(int nstates){ this->nstates = nstates;}
    void SetGspace(bool drho_pre, bool Gspace, double q0);
    void Refresh();

};
#endif

