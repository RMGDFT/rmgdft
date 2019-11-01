
#include <cmath>
#include <complex>
#include <iostream>
#include "Lattice.h"
#include "LaplacianCoeff.h"

template void FiniteDiffGrad(double * __restrict__, double * __restrict__, double * __restrict__, double * __restrict__, int , int , int , LaplacianCoeff *);
template void FiniteDiffGrad(float * __restrict__, float * __restrict__, float * __restrict__, float * __restrict__, int , int , int , LaplacianCoeff *);
template void FiniteDiffGrad(std::complex<double> * __restrict__, std::complex<double> * __restrict__, std::complex<double> * __restrict__, std::complex<double> * __restrict__, int , int , int , LaplacianCoeff *);
template void FiniteDiffGrad(std::complex<float> * __restrict__, std::complex<float> * __restrict__, std::complex<float> * __restrict__, std::complex<float> * __restrict__, int , int , int , LaplacianCoeff *);


template <typename T>
void FiniteDiffGrad(T * __restrict__ a, T * __restrict__ gx, T * __restrict__ gy, T * __restrict__ gz, int dimx, int dimy, int dimz, LaplacianCoeff *LC)
{

    int Lorder = LC->GetOrder();
    int dim[3];
    LC->GetDim(dim);

    // recalculate relative index
    if(dim[0] != dimx || dim[1] != dimy || dim[2] != dimz) 
    {
        std::cout << "dimxyz wanted : "<< dimx << "  " <<dimy << "  " <<dimz<< std::endl;
        std::cout << "dimxyz shouldbe: "<< dim[0] << "  " <<dim[1] << "  " <<dim[2]<< std::endl;
        rmg_error_handler(__FILE__, __LINE__, "Generalized finite diffenerence dims not match ");
    }


    // This code may seem unduly complicated but for the complex case it is almost twice as fast
    // as the simple templated version. I think the issue is poor vectorization of complex data
    // types but they may be treated as a double length array of floats or doubles instead which
    // does vectorize well.
    if(typeid(T) == typeid(std::complex<double>))
    {
        double *aad = (double *)a;
        double *gxd = (double *)gx;
        for(int idx = 0; idx < 2 * dim[0] * dim[1] * dim[2]; idx++) gxd[idx] = 0.0;
        for(auto coeff:LC->gx_coeff_and_index)
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
                            gxd[2*idxc]+= aad[2*idxs] * (double)coeff.coeff;
                            gxd[2*idxc+1]+= aad[2*idxs+1] * (double)coeff.coeff;
                        }
                }
            }
        }
    }
    else if(typeid(T) == typeid(std::complex<float>))
    {
        float *aad = (float *)a;
        float *gxd = (float *)gx;
        for(int idx = 0; idx < 2 * dim[0] * dim[1] * dim[2]; idx++) gxd[idx] = 0.0;
        for(auto coeff:LC->gx_coeff_and_index)
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
                            gxd[2*idxc]+= aad[2*idxs] * (float)coeff.coeff;
                            gxd[2*idxc+1]+= aad[2*idxs+1] * (float)coeff.coeff;
                        }
                }
            }
        }
    }
    else
    {
        for(int idx = 0; idx < dim[0] * dim[1] * dim[2]; idx++) gx[idx] = 0.0;
        for(auto coeff:LC->gx_coeff_and_index)
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
                            gx[idxc]+= a[idxs] * (T)coeff.coeff;
                        }
                }
            }
        }
    }

    if(typeid(T) == typeid(std::complex<double>))
    {
        double *aad = (double *)a;
        double *gyd = (double *)gy;
        for(int idx = 0; idx < 2 * dim[0] * dim[1] * dim[2]; idx++) gyd[idx] = 0.0;
        for(auto coeff:LC->gy_coeff_and_index)
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
                            gyd[2*idxc]+= aad[2*idxs] * (double)coeff.coeff;
                            gyd[2*idxc+1]+= aad[2*idxs+1] * (double)coeff.coeff;
                        }
                }
            }
        }
    }
    else if(typeid(T) == typeid(std::complex<float>))
    {
        float *aad = (float *)a;
        float *gyd = (float *)gy;
        for(int idx = 0; idx < 2 * dim[0] * dim[1] * dim[2]; idx++) gyd[idx] = 0.0;
        for(auto coeff:LC->gy_coeff_and_index)
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
                            gyd[2*idxc]+= aad[2*idxs] * (float)coeff.coeff;
                            gyd[2*idxc+1]+= aad[2*idxs+1] * (float)coeff.coeff;
                        }
                }
            }
        }
    }
    else
    {
        for(int idx = 0; idx < dim[0] * dim[1] * dim[2]; idx++) gy[idx] = 0.0;
        for(auto coeff:LC->gy_coeff_and_index)
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
                            gy[idxc]+= a[idxs] * (T)coeff.coeff;
                        }
                }
            }
        }
    }

    if(typeid(T) == typeid(std::complex<double>))
    {
        double *aad = (double *)a;
        double *gzd = (double *)gz;
        for(int idx = 0; idx < 2 * dim[0] * dim[1] * dim[2]; idx++) gzd[idx] = 0.0;
        for(auto coeff:LC->gz_coeff_and_index)
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
                            gzd[2*idxc]+= aad[2*idxs] * (double)coeff.coeff;
                            gzd[2*idxc+1]+= aad[2*idxs+1] * (double)coeff.coeff;
                        }
                }
            }
        }
    }
    else if(typeid(T) == typeid(std::complex<float>))
    {
        float *aad = (float *)a;
        float *gzd = (float *)gz;
        for(int idx = 0; idx < 2 * dim[0] * dim[1] * dim[2]; idx++) gzd[idx] = 0.0;
        for(auto coeff:LC->gz_coeff_and_index)
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
                            gzd[2*idxc]+= aad[2*idxs] * (float)coeff.coeff;
                            gzd[2*idxc+1]+= aad[2*idxs+1] * (float)coeff.coeff;
                        }
                }
            }
        }
    }
    else
    {
        for(int idx = 0; idx < dim[0] * dim[1] * dim[2]; idx++) gz[idx] = 0.0;
        for(auto coeff:LC->gz_coeff_and_index)
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
                            gz[idxc]+= a[idxs] * (T)coeff.coeff;
                        }
                }
            }
        }
    }

}
