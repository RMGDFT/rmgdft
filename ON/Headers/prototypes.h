void get_vh (double * rho, double * rhoc, double * vh_eig, int min_sweeps, int max_sweeps, int maxlevel, double rms_target);
void mgrid_solv (REAL *v_mat, REAL *f_mat, REAL *work,
                 int dimx, int dimy, int dimz, REAL gridhx, REAL gridhy,
                 REAL gridhz, int level, int *nb_ids, int max_levels,
                 int *pre_cyc, int *post_cyc, int mu_cyc, REAL step, REAL k);
REAL app_cil_orbital (REAL * a, REAL * b, int dimx, int dimy, int dimz,
              REAL gridhx, REAL gridhy, REAL gridhz);
REAL app_cil_orbital6 (REAL * a, REAL * b, int dimx, int dimy, int dimz,
              REAL gridhx, REAL gridhy, REAL gridhz);
void mg_prolong_6 (double *full, double *half, int dimx, int dimy, int dimz, int half_dimx,
                       int half_dimy, int half_dimz, int grid_ratio, int order);

