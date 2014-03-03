#ifndef RMG_TradeImages_H
#define RMG_TradeImages_H 1

#include <mpi.h>
#include "BaseGrid.h"
#include "rmg_error.h"


/* Type of image trading */
#define FULL_TRADE 1
#define CENTRAL_TRADE 2

#if __cplusplus

#include "BaseThread.h"

class TradeImages : public RmgError {


private:

    /// Synchronous/asynchronous mode. 0=asnychronous (default) 1=synchronous
    static int mode;

    /// Rank of target node based on offsets from current node. Used by asynchronous comm routines.
    static int target_node[3][3][3];

    static int max_alloc;

    /// Buffers allocated via MPI_Allocmem
    static rmg_double_t *swbuf1x;
    static rmg_double_t *swbuf2x;
    static rmg_double_t *frdx1, *frdx2, *frdy1, *frdy2, *frdz1, *frdz2;
    static rmg_double_t *frdx1n, *frdx2n, *frdy1n, *frdy2n, *frdz1n, *frdz2n;
    static rmg_double_t *yzpsms_r, *yzpsps_r, *yzmsms_r, *yzmsps_r;
    static rmg_double_t *yzpsms_s, *yzpsps_s, *yzmsms_s, *yzmsps_s;
    static rmg_double_t *xzpsms_r, *xzpsps_r, *xzmsms_r, *xzmsps_r;
    static rmg_double_t *xzpsms_s, *xzpsps_s, *xzmsms_s, *xzmsps_s;
    static rmg_double_t *xypsms_r, *xypsps_r, *xymsms_r, *xymsps_r;
    static rmg_double_t *xypsms_s, *xypsps_s, *xymsms_s, *xymsps_s;
    static rmg_double_t *m0_s, *m0_r;

    static MPI_Request sreqs[26];
    static MPI_Request rreqs[26];

    void init_trade_imagesx_async(void);
    template <typename RmgType> void RMG_MPI_trade(RmgType *buf, int count, int type, int pe_x_offset, int pe_y_offset, int pe_z_offset, MPI_Comm comm, int tag, MPI_Request *req);
    template <typename RmgType> void trade_imagesx_async (RmgType * f, RmgType * w, int dimx, int dimy, int dimz, int images);
    template <typename RmgType> void trade_imagesx_central_async (RmgType * f, RmgType * w, int dimx, int dimy, int dimz, int images);
    template <typename RmgType> void trade_images1_central_async (RmgType * f, int dimx, int dimy, int dimz);
    template <typename RmgType> void trade_images1_async (RmgType * f, int dimx, int dimy, int dimz);


public:
    TradeImages(void);
    void set_synchronous_mode(void);
    void set_asynchronous_mode(void);
    template <typename RmgType> void trade_imagesx (RmgType *f, RmgType *w, int dimx, int dimy, int dimz, int images, int type);
    template <typename RmgType> void trade_images (RmgType * mat, int dimx, int dimy, int dimz, int type);




};

#endif
#endif
