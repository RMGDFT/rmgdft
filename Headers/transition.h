#ifndef RMG_transition_h
#define RMG_transition_h
extern "C"
{
rmg_double_t my_crtc (void);
MPI_Comm transition_get_grid_comm(void);
void thread_barrier_wait(void);
int transition_get_gridpe(void);
}
#endif

