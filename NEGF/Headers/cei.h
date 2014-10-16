
struct complex_energy_integral
{
    double EB;
    double KT;
    double GAMMA;
    double DELTA2;
    double DELTA;
    int ncircle;
    int nmax_gq1;
    int nmax_gq2;
    int num_probe;
	int *probe_in_block;
    int num_subsystem;
	int *subsystem_idx;
    int num_probe_window;
    int *probe_window_start;
    int *probe_window_end;
    int num_dos_window;
    int *dos_window_start;
    int *dos_window_end;

    int Npulaysave;
    int Npulayrefresh;
    double pulaymix;
    int probe_noneq;
    int energy_point_insert;
};
typedef struct  complex_energy_integral complex_energy_integral;

extern complex_energy_integral cei;
