#ifndef RMG_vhartree_h
#define RMG_vhartree_h


double CPP_get_vh (BaseGrid *G, Lattice *L, TradeImages *T, double * rho, double *vhartree,
                 int min_sweeps, int max_sweeps, int maxlevel,
                 int global_presweeps, int global_postsweeps, int mucycles,
                 double rms_target, double global_step, double coarse_step, int boundaryflag, int density, bool print_status);

#endif
