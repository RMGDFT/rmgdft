
#include <cmath>
#include <complex>
#include <iostream>
#include "Lattice.h"
#include "LaplacianCoeff.h"

template double FiniteDiffLap(double * __restrict__, double * __restrict__, int , int , int , LaplacianCoeff *);
template double FiniteDiffLap(float * __restrict__, float * __restrict__, int , int , int , LaplacianCoeff *);
template double FiniteDiffLap(std::complex<double> * __restrict__, std::complex<double> * __restrict__, int , int , int , LaplacianCoeff *);
template double FiniteDiffLap(std::complex<float> * __restrict__, std::complex<float> * __restrict__, int , int , int , LaplacianCoeff *);


template <typename T>
double FiniteDiffLap(T * __restrict__ a, T * __restrict__ b, int dimx, int dimy, int dimz, LaplacianCoeff *LC)
{

    int Lorder = LC->GetOrder();
    int dim[3];
    LC->GetDim(dim);

    //recalculate relative index
    if(dim[0] != dimx || dim[1] != dimy || dim[2] != dimz) 
    {
        std::cout << "dimxyz wanted : "<< dimx << "  " <<dimy << "  " <<dimz<< std::endl;
        std::cout << "dimxyz shouldbe: "<< dim[0] << "  " <<dim[1] << "  " <<dim[2]<< std::endl;

        rmg_error_handler(__FILE__, __LINE__, "Generalized finite diffenerence dims not match ");
    }


    for(int idx = 0; idx < dim[0] * dim[1] * dim[2]; idx++) b[idx] = 0.0;

    // This code may seem unduly complicated but for the complex case it is almost twice as fast
    // as the simple templated version. I think the issue is poor vectorization of complex data
    // types but they may be treated as a double length array of floats or doubles instead which
    // does vectorize well.

    if(typeid(T) == typeid(std::complex<double>))
    {
        double *aad = (double *)a;
        double *bbd = (double *)b;
        for(auto coeff:LC->coeff_and_index)
        {

            for (int ix = 0; ix < dim[0]; ix++)
            {
                int ixs = (ix+Lorder/2) * (dim[1] + Lorder) * (dim[2] + Lorder);
                int ixc = ix * dim[1] * dim[2];

                for (int iy = 0; iy < dim[1]; iy++)
                {
                    int iys = (iy+Lorder/2) * (dim[2] + Lorder) + ixs + Lorder/2;
                    int iyc = iy * dim[2] + ixc;

                    for(auto index: coeff.relative_index)
                        for (int iz = 0; iz < dim[2]; iz++)
                        {
                            int idxs = iys + iz + index;
                            int idxc = iyc + iz;
                            bbd[2*idxc] += aad[2*idxs] * (double)coeff.coeff;
                            bbd[2*idxc+1] += aad[2*idxs+1] * (double)coeff.coeff;
                        }
                }
            }
        }

    }
    else if(typeid(T) == typeid(std::complex<float>))
    {
        float *aad = (float *)a;
        float *bbd = (float *)b;
        for(auto coeff:LC->coeff_and_index)
        {

            for (int ix = 0; ix < dim[0]; ix++)
            {
                int ixs = (ix+Lorder/2) * (dim[1] + Lorder) * (dim[2] + Lorder);
                int ixc = ix * dim[1] * dim[2];

                for (int iy = 0; iy < dim[1]; iy++)
                {
                    int iys = (iy+Lorder/2) * (dim[2] + Lorder) + ixs + Lorder/2;
                    int iyc = iy * dim[2] + ixc;

                    for(auto index: coeff.relative_index)
                        for (int iz = 0; iz < dim[2]; iz++)
                        {
                            int idxs = iys + iz + index;
                            int idxc = iyc + iz;
                            bbd[2*idxc] += aad[2*idxs] * (float)coeff.coeff;
                            bbd[2*idxc+1] += aad[2*idxs+1] * (float)coeff.coeff;
                        }
                }
            }
        }

    }
    else
    {
        for(auto coeff:LC->coeff_and_index)
        {

            for (int ix = 0; ix < dim[0]; ix++)
            {
                int ixs = (ix+Lorder/2) * (dim[1] + Lorder) * (dim[2] + Lorder);
                int ixc = ix * dim[1] * dim[2];

                for (int iy = 0; iy < dim[1]; iy++)
                {
                    int iys = (iy+Lorder/2) * (dim[2] + Lorder) + ixs + Lorder/2;
                    int iyc = iy * dim[2] + ixc;

                    for(auto index: coeff.relative_index)
                        for (int iz = 0; iz < dim[2]; iz++)
                        {
                            int idxs = iys + iz + index;
                           int idxc = iyc + iz;
                            b[idxc]+= a[idxs] * (T)coeff.coeff;
                        }
                }
            }
        }
    }


    /* Return the diagonal component of the operator */
    return LC->coeff_and_index[0].coeff;
}
