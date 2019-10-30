
#include <cmath>
#include <complex>
#include <iostream>
#include "Lattice.h"
#include "LaplacianCoeff.h"

double FiniteDiffLap(double * __restrict__ a, double * __restrict__ b, int dimx, int dimy, int dimz, LaplacianCoeff *LC)
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
             //      daxpy(&dim[2], &coeff.coeff, &a[iys+index], &ione, &b[iyc], &ione);
                    for (int iz = 0; iz < dim[2]; iz++)
                    {
                        int idxs = iys + iz + index;
                       int idxc = iyc + iz;
                        b[idxc]+= a[idxs] * coeff.coeff;
                    }
            }
        }
    }

    /* Return the diagonal component of the operator */
    return LC->coeff_and_index[0].coeff;
}


double FiniteDiffLap(float * __restrict__ a, float * __restrict__ b, int dimx, int dimy, int dimz, LaplacianCoeff *LC)
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
             //      daxpy(&dim[2], &coeff.coeff, &a[iys+index], &ione, &b[iyc], &ione);
                    for (int iz = 0; iz < dim[2]; iz++)
                    {
                        int idxs = iys + iz + index;
                       int idxc = iyc + iz;
                        b[idxc]+= a[idxs] * (float)coeff.coeff;
                    }
            }
        }
    }

    /* Return the diagonal component of the operator */
    return LC->coeff_and_index[0].coeff;
}


double FiniteDiffLap(std::complex<double> * __restrict__ a, std::complex<double> * __restrict__ b, int dimx, int dimy, int dimz, LaplacianCoeff *LC)
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
             //      daxpy(&dim[2], &coeff.coeff, &a[iys+index], &ione, &b[iyc], &ione);
                    for (int iz = 0; iz < dim[2]; iz++)
                    {
                        int idxs = iys + iz + index;
                       int idxc = iyc + iz;
                        b[idxc]+= a[idxs] * coeff.coeff;
                    }
            }
        }
    }

    /* Return the diagonal component of the operator */
    return LC->coeff_and_index[0].coeff;
}


double FiniteDiffLap(std::complex<float> * __restrict__ a, std::complex<float> * __restrict__ b, int dimx, int dimy, int dimz, LaplacianCoeff *LC)
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
             //      daxpy(&dim[2], &coeff.coeff, &a[iys+index], &ione, &b[iyc], &ione);
                    for (int iz = 0; iz < dim[2]; iz++)
                    {
                        int idxs = iys + iz + index;
                       int idxc = iyc + iz;
                        b[idxc]+= a[idxs] * (float)coeff.coeff;
                    }
            }
        }
    }

    /* Return the diagonal component of the operator */
    return LC->coeff_and_index[0].coeff;
}

