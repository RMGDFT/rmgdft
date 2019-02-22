#ifndef RMG_Pulay_H
#define RMG_Pulay_H 1

#include<functional>
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
    double *A_mat;
    std::vector<double*> hist_ptr;
    std::vector<double*> res_hist_ptr;
    std::function<void(double*)> Precond;
    bool need_precond;

public:

    PulayMixing(size_t Nsize, int pulay_order, int refresh_steps, double mix_init, double beta, MPI_Comm comm);
    ~PulayMixing(void);
    void Mixing(double *xm, double *fm);
    void SetPrecond(std::function<void(double*)> precon);
    void Refresh();

};
#endif

