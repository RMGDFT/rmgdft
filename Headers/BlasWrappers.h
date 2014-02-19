#ifndef RMG_BlasWrappers_h
#define RMG_BlasWrappers_h

template <typename RmgType>
void QMD_copy (int n, RmgType * x, int incx, RmgType * y, int incy)
{
      int i, ix=0, iy=0;
      if((incx == 1) && (incy == 1)) {
         for(i = 0;i < n;i++) {
             y[i] = x[i];
         }
         return;
      }

      for(i = 0;i < n;i++) {
          y[iy] = x[ix];
          ix += incx;
          iy += incy;
      }
}

#endif

