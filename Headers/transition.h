#ifndef RMG_transition_h
#define RMG_transition_h

#if __cplusplus
#include "BaseGrid.h"
#include "Lattice.h"
#include "TradeImages.h"

extern BaseGrid *Rmg_G;
extern TradeImages *Rmg_T;
extern Lattice Rmg_L;

extern "C"
{
double my_crtc (void);
MPI_Comm transition_get_grid_comm(void);
void thread_barrier_wait(void);
int transition_get_gridpe(void);
}
#endif
#endif

