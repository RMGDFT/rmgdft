#ifndef RMG_TradeImages_H
#define RMG_TradeImages_H 1

#include <mpi.h>
#include "BaseGrid.h"
#include "rmg_error.h"


/* Type of image trading */
#define FULL_TRADE 1
#define CENTRAL_TRADE 2

/* Maximum number of images for finite difference routines */
#define MAX_TRADE_IMAGES 5

#if __cplusplus

#include "BaseThread.h"

class TradeImages {


private:

    BaseGrid *G;

    /// Synchronous/asynchronous mode. 0=asnychronous (default) 1=synchronous
    int mode;

    // rank of this node in comm
    int gridpe;

    /// Rank of target node based on offsets from current node. Used by asynchronous comm routines.
    int target_node[3][3][3];

    int max_alloc;

    /// Buffers allocated via MPI_Allocmem
    double *swbuf1x;
    double *swbuf2x;
    double *frdx1, *frdx2, *frdy1, *frdy2, *frdz1, *frdz2;
    double *frdx1n, *frdx2n, *frdy1n, *frdy2n, *frdz1n, *frdz2n;
    double *yzpsms_r, *yzpsps_r, *yzmsms_r, *yzmsps_r;
    double *yzpsms_s, *yzpsps_s, *yzmsms_s, *yzmsps_s;
    double *xzpsms_r, *xzpsps_r, *xzmsms_r, *xzmsps_r;
    double *xzpsms_s, *xzpsps_s, *xzmsms_s, *xzmsps_s;
    double *xypsms_r, *xypsps_r, *xymsms_r, *xymsps_r;
    double *xypsms_s, *xypsps_s, *xymsms_s, *xymsps_s;
    double *m0_s, *m0_r;

    MPI_Request sreqs[26];
    MPI_Request rreqs[26];

    void init_trade_imagesx_async(void);
    template <typename RmgType> void RMG_MPI_trade(RmgType *buf, int count, int type, int pe_x_offset, int pe_y_offset, int pe_z_offset, MPI_Comm comm, int tag, MPI_Request *req);
    template <typename RmgType> void trade_imagesx_async (RmgType * f, RmgType * w, int dimx, int dimy, int dimz, int images);
    template <typename RmgType> void trade_imagesx_central_async (RmgType * f, RmgType * w, int dimx, int dimy, int dimz, int images);
    template <typename RmgType> void trade_images1_central_async (RmgType * f, int dimx, int dimy, int dimz);
    template <typename RmgType> void trade_images1_async (RmgType * f, int dimx, int dimy, int dimz);


public:
    /// MPI communicator to use
    MPI_Comm comm;

    TradeImages(BaseGrid *BG);
    ~TradeImages(void);
    void set_synchronous_mode(void);
    void set_asynchronous_mode(void);
    void set_MPI_comm(MPI_Comm comm);
    MPI_Comm get_MPI_comm(void);
    void set_gridpe(int gridpe);
    template <typename RmgType> void trade_imagesx (RmgType *f, RmgType *w, int dimx, int dimy, int dimz, int images, int type);
    template <typename RmgType> void trade_images (RmgType * mat, int dimx, int dimy, int dimz, int type);




};

#endif
#endif
