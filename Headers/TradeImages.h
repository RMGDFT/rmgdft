#ifndef RMG_TradeImages_H
#define RMG_TradeImages_H 1

#include <mpi.h>
#include "const.h"
#include "BaseGrid.h"
#include "BaseThread.h"
#include "rmg_error.h"


class TradeImages {



private:
    static rmg_double_t *swbuf1x;
    static rmg_double_t *swbuf2x;
    static int max_alloc;

public:
    TradeImages(void);
    template <typename RmgType> void CPP_trade_imagesx (RmgType *f, RmgType *w, int dimx, int dimy, int dimz, int images, int type);


};

#endif
