
#include "TradeImages.h"
#include "common_prototypes.h"
using namespace std;

#include "transition.h"

// C interfaces for transitional usage
extern "C" void init_TradeImages(void)
{
    Rmg_T = new TradeImages(Rmg_G);
}
extern "C" void set_MPI_comm(MPI_Comm comm)
{
    Rmg_T->set_MPI_comm(comm);
}
extern "C"  void trade_imagesx (rmg_double_t * f, rmg_double_t * w, int dimx, int dimy, int dimz, int images, int type)
{
    Rmg_T->trade_imagesx<double>(f, w, dimx, dimy, dimz, images, type);
}
extern "C" void trade_images (rmg_double_t *mat, int dimx, int dimy, int dimz, int type)
{
    Rmg_T->trade_images<double>(mat, dimx, dimy, dimz, type);
}

