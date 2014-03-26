
#include "TradeImages.h"
#include "common_prototypes.h"
using namespace std;


// C interfaces for transitional usage
extern "C" void init_TradeImages(void)
{
    TradeImages *T;
    T = new TradeImages();
}
extern "C" void set_MPI_comm(MPI_Comm comm)
{
    TradeImages T;
    T.set_MPI_comm(comm);
}
extern "C"  void trade_imagesx (rmg_double_t * f, rmg_double_t * w, int dimx, int dimy, int dimz, int images, int type)
{
    TradeImages T;
    T.trade_imagesx<double>(f, w, dimx, dimy, dimz, images, type);
}
extern "C" void trade_images (rmg_double_t *mat, int dimx, int dimy, int dimz, int type)
{
    TradeImages T;
    T.trade_images<double>(mat, dimx, dimy, dimz, type);
}

