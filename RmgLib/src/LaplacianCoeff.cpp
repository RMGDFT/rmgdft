#include <vector>
#include <algorithm>
#include <math.h>
#include <iostream>
#include <csignal>
#include <float.h>
#include "LaplacianCoeff.h"
#include "blas.h"

struct {
    bool operator()(GridPoint a, GridPoint b) const
    {   
        return a.dist < b.dist;
    }   
} customLess_dist;
struct {
    bool operator()(GridPoint a, GridPoint b) const
    {   
        return std::abs(a.index[0]) > std::abs(b.index[0]) ;
    }   
} customLess_x;
struct {
    bool operator()(GridPoint a, GridPoint b) const
    {   
        return std::abs(a.index[1]) > std::abs(b.index[1]);
    }   
} customLess_y;
struct {
    bool operator()(GridPoint a, GridPoint b) const
    {   
        return std::abs(a.index[2]) > std::abs(b.index[2]);
    }   
} customLess_z;

struct {
    bool operator()(GridPoint a, GridPoint b) const
    {   
        return std::abs(a.ijk) < std::abs(b.ijk);
    }   
} customLess_ijk;

LaplacianCoeff::LaplacianCoeff(double a[3][3], int Ngrid[3], int Lorder, int dim[3]){

    // Figure out the rescaling to improve numerical stability
    scale1 = sqrt(a[0][0]*a[0][0] + a[0][1]*a[0][1] + a[0][2]*a[0][2]) / Ngrid[0];
    scale1 = 1.0 / scale1;
    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
            this->a[i][j] = scale1*a[i][j];
    for(int i = 0; i < 3; i++){
        this->Ngrid[i] = Ngrid[i];
        this->dim[i] = dim[i];
        this->Lorder = Lorder;
    }
}
void LaplacianCoeff::CalculateCoeff()
{
    this->CalculateCoeff(this->a, this->Ngrid, this->Lorder, this->dim);
}

void LaplacianCoeff::CalculateCoeff(double a[3][3], int Ngrid[3], int Lorder, int dim[3])
{

    for(int i=0;i<13;i++)
    {
        // 0=x,1=y,2=z,3=xy,4=xz,5=yz,6=nxy,7=nxz,8=nyz
        this->include_axis[i] = false;
        this->plane_centers[i] = 0.0;
        this->plane_dists[i] = 10000000.0;
        for(int j=0;j<12;j++)
        {
            this->axis_lc[i][j] = 0.0;
            this->axis_gc_x[i][j] = 0.0;
            this->axis_gc_y[i][j] = 0.0;
            this->axis_gc_z[i][j] = 0.0;
        }
    }

    std::vector<GridPoint>  points, points1;

    std::vector<GridPoint> der_list;  // only i,j,k are needed.
    int dimension;
    // for laplacian, 2n order is also 2n+1 order since the odd partianl derivatives cancelled automatically.
    Lorder = Lorder /2 * 2;
    if( (this->ibrav == ORTHORHOMBIC_PRIMITIVE || 
        this->ibrav == CUBIC_PRIMITIVE  ||
        this->ibrav == TETRAGONAL_PRIMITIVE) && !this->offdiag)
    {
        dimension = 1;
        points.clear();
        points1.clear();
        GetDerList(der_list, Lorder, dimension, 0);
        GetPointList1D(points1, a[0][0], Ngrid[0], Lorder, 0);
        this->BuildSolveLinearEq(points1, der_list, dimension);
        points.insert(std::end(points), std::begin(points1), std::end(points1));

        points1.clear();
        GetDerList(der_list, Lorder, dimension, 1);
        GetPointList1D(points1, a[1][1], Ngrid[1], Lorder, 1);
        this->BuildSolveLinearEq(points1, der_list, dimension);
        points.insert(std::end(points), std::begin(points1), std::end(points1));

        points1.clear();
        GetDerList(der_list, Lorder, dimension, 2);
        GetPointList1D(points1, a[2][2], Ngrid[2], Lorder, 2);
        this->BuildSolveLinearEq(points1, der_list, dimension);
        points.insert(std::end(points), std::begin(points1), std::end(points1));

        std::stable_sort(points.begin(), points.end(), customLess_dist);
    }
    else if ((this->ibrav == HEXAGONAL) || (this->ibrav == HEXAGONAL2))
    {
        double x[20], w1[20], w2[20];
        for(int i=0;i<20;i++) x[i] = (double)i;
        FiniteDiff::gen_weights(Lorder+1, 2, (double)Lorder/2, x, w1);
        FiniteDiff::gen_weights(Lorder+1, 1, (double)Lorder/2, x, w2);
        double id=1.0;
        if(this->ibrav == HEXAGONAL) id = -1.0;
        dimension = 2;
        double a2d[2][2];
        a2d[0][0] = a[0][0];
        a2d[0][1] = a[0][1];
        a2d[1][0] = a[1][0];
        a2d[1][1] = a[1][1];
        double h = a2d[0][0] / (double)Ngrid[0];
        double h2 = h*h;
        double hf = sqrt(3.0/4.0);
        double hf1 = sqrt(3.0);

        GetPointList2D(points1, a2d, Ngrid, Lorder);
        for(auto &a:points1)
        {
            if(a.index[0] != 0 && a.index[1]==0)
            {
                a.coeff = 2.0*w1[Lorder/2+a.index[0]]/3.0/h2;
                a.coeff_gx = (double)a.index[0]*w2[Lorder/2+a.index[0]]/3.0/h;
                a.coeff_gy = 0.0;
                a.coeff_gz = 0.0;
            }
            else if(a.index[0] == 0 && a.index[1]!=0)
            {
                a.coeff = 2.0*w1[Lorder/2+a.index[1]]/3.0/h2;
                a.coeff_gx = hf*id*(double)a.index[1]*w2[Lorder/2+a.index[1]]/3.0/h/hf1;
                a.coeff_gy = hf*(double)a.index[1]*w2[Lorder/2+a.index[1]]/3.0/h;
                a.coeff_gz = 0.0;
            }
            else if(a.index[0] == a.index[1])
            {
                a.coeff = 2.0*w1[Lorder/2+a.index[0]]/3.0/h2;
                a.coeff_gx = hf*(double)a.index[1]*w2[Lorder/2+a.index[1]]/3.0/h/hf1;
                a.coeff_gy = hf*(double)a.index[1]*w2[Lorder/2+a.index[1]]/3.0/h;
                a.coeff_gz = 0.0;
            }
            else if(a.index[0] == -a.index[1])
            {
                a.coeff = 2.0*w1[Lorder/2+a.index[0]]/3.0/h2;
                a.coeff_gx = -hf*id*(double)a.index[1]*w2[Lorder/2+a.index[1]]/3.0/h/hf1;
                a.coeff_gy = hf*(double)a.index[1]*w2[Lorder/2+a.index[1]]/3.0/h;
                a.coeff_gz = 0.0;
            }
        }

        for(auto &a:points1)
        {
            a.coeff *= scale1*scale1;
            a.coeff_gx *= scale1;
            a.coeff_gy *= scale1;
            a.coeff_gz *= scale1;
        }

        points.insert(std::end(points), std::begin(points1), std::end(points1));

        points1.clear();
        dimension = 1;
        GetDerList(der_list, Lorder, dimension, 2);
        GetPointList1D(points1, a[2][2], Ngrid[2], Lorder, 2);
        this->BuildSolveLinearEq(points1, der_list, dimension);
        points.insert(std::end(points), std::begin(points1), std::end(points1));

        std::stable_sort(points.begin(), points.end(), customLess_dist);

    }
    else if((this->ibrav == MONOCLINIC_PRIMITIVE) && !this->offdiag)
    {
        dimension = 2;
        GetDerList(der_list, Lorder, dimension, 0);
        double a2d[2][2];
        a2d[0][0] = a[0][0];
        a2d[0][1] = a[0][1];
        a2d[1][0] = a[1][0];
        a2d[1][1] = a[1][1];

        GetPointList2D(points1, a2d, Ngrid, Lorder);
        this->BuildSolveLinearEq(points1, der_list, dimension);
        points.insert(std::end(points), std::begin(points1), std::end(points1));

        points1.clear();
        dimension = 1;
        GetDerList(der_list, Lorder, dimension, 2);
        GetPointList1D(points1, a[2][2], Ngrid[2], Lorder, 2);
        this->BuildSolveLinearEq(points1, der_list, dimension);
        points.insert(std::end(points), std::begin(points1), std::end(points1));

        std::stable_sort(points.begin(), points.end(), customLess_dist);
    }
    else if( this->ibrav == CUBIC_FC && !this->offdiag )
    {
        
        dimension = 3;
        der_list.clear();
        points.clear();
        points1.clear();
        GetDerListFCC(der_list, Lorder);
        GetPointListFCC(points1, a, Ngrid, Lorder);
        this->BuildSolveLinearEq(points1, der_list, dimension);
        points.insert(std::end(points), std::begin(points1), std::end(points1));
        std::stable_sort(points.begin(), points.end(), customLess_dist);

    }
    else if( abs(this->ibrav) == CUBIC_BC && !this->offdiag )
    {
        double x[20], w1[20], w2[20];
        for(int i=0;i<20;i++) x[i] = (double)i;
        FiniteDiff::gen_weights(Lorder+1, 2, (double)Lorder/2, x, w1);
        FiniteDiff::gen_weights(Lorder+1, 1, (double)Lorder/2, x, w2);
        // First derivative requires this
        if(this->ibrav > 0)
            for(int i=0;i <= Lorder;i++) w2[i] = -w2[i] * (double)std::abs(i-Lorder/2);
        else
            for(int i=0;i <= Lorder;i++) w2[i] = w2[i] * (double)std::abs(i-Lorder/2);

        points.clear();
        points1.clear();
        double h = a[0][0] / (double)Ngrid[0];
        double h2 = h*h;
        double fac = 1.0 / 8.0 / h;
        GetPointListBCC(points1, a, Ngrid, Lorder);
        for(auto &a:points1)
        {
            if(a.index[0] == a.index[1] && a.index[0] == a.index[2])
            {
              a.coeff = w1[Lorder/2+a.index[0]]/4.0/h2;
              int i1 = a.index[0];
              a.coeff_gx = fac*w2[Lorder/2+i1];
              a.coeff_gy = fac*w2[Lorder/2+i1];
              a.coeff_gz = fac*w2[Lorder/2+i1];
            }
            else if(a.index[0] != 0 && a.index[1] == 0 && a.index[2] == 0)
            {
//A
              a.coeff = w1[Lorder/2+a.index[0]]/4.0/h2;
              int i1 = a.index[0];
              if(this->ibrav > 0) {
                  a.coeff_gx = fac*w2[Lorder/2+i1];
                  a.coeff_gy = fac*w2[Lorder/2+i1];
                  a.coeff_gz = -fac*w2[Lorder/2+i1];
              }
              else {
                  a.coeff_gx = -fac*w2[Lorder/2+i1];
                  a.coeff_gy = fac*w2[Lorder/2+i1];
                  a.coeff_gz = fac*w2[Lorder/2+i1];
              }
            
            } 
            else if(a.index[0] == 0 && a.index[1] != 0 && a.index[2] == 0)
            {
//C
              a.coeff = w1[Lorder/2+a.index[1]]/4.0/h2;
              int i1 = a.index[1];
              if(this->ibrav > 0) {
                  a.coeff_gx = -fac*w2[Lorder/2+i1];
                  a.coeff_gy = fac*w2[Lorder/2+i1];
                  a.coeff_gz = fac*w2[Lorder/2+i1];
              }
              else {
                  a.coeff_gx = fac*w2[Lorder/2+i1];
                  a.coeff_gy = -fac*w2[Lorder/2+i1];
                  a.coeff_gz = fac*w2[Lorder/2+i1];
              }
            }
            else if(a.index[0] == 0 && a.index[1] ==0 && a.index[2] != 0)
            {
//B
              a.coeff = w1[Lorder/2+a.index[2]]/4.0/h2;
              int i1 = a.index[2];
              if(this->ibrav > 0) {
                  a.coeff_gx = fac*w2[Lorder/2+i1];
                  a.coeff_gy = -fac*w2[Lorder/2+i1];
                  a.coeff_gz = fac*w2[Lorder/2+i1];
              }
              else {
                  a.coeff_gx = fac*w2[Lorder/2+i1];
                  a.coeff_gy = fac*w2[Lorder/2+i1];
                  a.coeff_gz = -fac*w2[Lorder/2+i1];
              }
            }
        }

        for(auto &a:points1)
        {
            a.coeff *= scale1*scale1;
            a.coeff_gx *= scale1;
            a.coeff_gy *= scale1;
            a.coeff_gz *= scale1;
        }

        points.insert(std::end(points), std::begin(points1), std::end(points1));
        std::stable_sort(points.begin(), points.end(), customLess_dist);
    }
    else       // TRICLINIC
    {
        dimension = 3;
        der_list.clear();
        points.clear();
        points1.clear();
        GetDerList(der_list, Lorder, dimension, 0);
        GetPointList3D(points1, a, Ngrid, Lorder);

        this->BuildSolveLinearEq(points1, der_list, dimension);
        points.insert(std::end(points), std::begin(points1), std::end(points1));
        std::stable_sort(points.begin(), points.end(), customLess_dist);
    }

    // Put coeffs into axis arrays for efficient application.
    for(auto a:points)
    {   
        if((a.index[0] != 0) && (a.index[1] == 0) && (a.index[2] == 0))
        {
            this->plane_dists[0] = std::min(this->plane_dists[0], std::abs(a.dist));
            this->axis_lc[0][a.index[0]+Lorder/2] = a.coeff;
            this->plane_centers[0] -= a.coeff;
            this->axis_gc_x[0][a.index[0]+Lorder/2] = a.coeff_gx;
            this->axis_gc_y[0][a.index[0]+Lorder/2] = a.coeff_gy;
            this->axis_gc_z[0][a.index[0]+Lorder/2] = a.coeff_gz;
            this->include_axis[0] = true;
        }
        else if((a.index[0] == 0) && (a.index[1] != 0) && (a.index[2] == 0))
        {
            this->plane_dists[1] = std::min(this->plane_dists[1], std::abs(a.dist));
            this->axis_lc[1][a.index[1]+Lorder/2] = a.coeff;
            this->plane_centers[1] -= a.coeff;
            this->axis_gc_x[1][a.index[1]+Lorder/2] = a.coeff_gx;
            this->axis_gc_y[1][a.index[1]+Lorder/2] = a.coeff_gy;
            this->axis_gc_z[1][a.index[1]+Lorder/2] = a.coeff_gz;
            this->include_axis[1] = true;
        }
        else if((a.index[0] == 0) && (a.index[1] == 0) && (a.index[2] != 0))
        {
            this->plane_dists[2] = std::min(this->plane_dists[2], std::abs(a.dist));
            this->axis_lc[2][a.index[2]+Lorder/2] = a.coeff;
            this->plane_centers[2] -= a.coeff;
            this->axis_gc_x[2][a.index[2]+Lorder/2] = a.coeff_gx;
            this->axis_gc_y[2][a.index[2]+Lorder/2] = a.coeff_gy;
            this->axis_gc_z[2][a.index[2]+Lorder/2] = a.coeff_gz;
            this->include_axis[2] = true;
        }
        else if(a.index[0] == a.index[1] && a.index[2] == 0)
        {
            this->plane_dists[3] = std::min(this->plane_dists[3], std::abs(a.dist));
            this->axis_lc[3][a.index[1]+Lorder/2] = a.coeff;
            this->plane_centers[3] -= a.coeff;
            this->axis_gc_x[3][a.index[1]+Lorder/2] = a.coeff_gx;
            this->axis_gc_y[3][a.index[1]+Lorder/2] = a.coeff_gy;
            this->axis_gc_z[3][a.index[1]+Lorder/2] = a.coeff_gz;
            this->include_axis[3] = true;
        }
        else if(a.index[0] == a.index[2] && a.index[1] == 0)
        {
            this->plane_dists[4] = std::min(this->plane_dists[4], std::abs(a.dist));
            this->axis_lc[4][a.index[2]+Lorder/2] = a.coeff;
            this->plane_centers[4] -= a.coeff;
            this->axis_gc_x[4][a.index[2]+Lorder/2] = a.coeff_gx;
            this->axis_gc_y[4][a.index[2]+Lorder/2] = a.coeff_gy;
            this->axis_gc_z[4][a.index[2]+Lorder/2] = a.coeff_gz;
            this->include_axis[4] = true;
        }
        else if(a.index[1] == a.index[2] && a.index[0] == 0)
        {
            this->plane_dists[5] = std::min(this->plane_dists[5], std::abs(a.dist));
            this->axis_lc[5][a.index[2]+Lorder/2] = a.coeff;
            this->plane_centers[5] -= a.coeff;
            this->axis_gc_x[5][a.index[2]+Lorder/2] = a.coeff_gx;
            this->axis_gc_y[5][a.index[2]+Lorder/2] = a.coeff_gy;
            this->axis_gc_z[5][a.index[2]+Lorder/2] = a.coeff_gz;
            this->include_axis[5] = true;
        }
        else if(a.index[0] == -a.index[1] && a.index[2] == 0)
        {
            this->plane_dists[6] = std::min(this->plane_dists[6], std::abs(a.dist));
            this->axis_lc[6][a.index[1]+Lorder/2] = a.coeff;
            this->plane_centers[6] -= a.coeff;
            this->axis_gc_x[6][a.index[1]+Lorder/2] = a.coeff_gx;
            this->axis_gc_y[6][a.index[1]+Lorder/2] = a.coeff_gy;
            this->axis_gc_z[6][a.index[1]+Lorder/2] = a.coeff_gz;
            this->include_axis[6] = true;
        }
        else if(a.index[0] == -a.index[2] && a.index[1] == 0)
        {
            this->plane_dists[7] = std::min(this->plane_dists[7], std::abs(a.dist));
            this->axis_lc[7][a.index[2]+Lorder/2] = a.coeff;
            this->plane_centers[7] -= a.coeff;
            this->axis_gc_x[7][a.index[2]+Lorder/2] = a.coeff_gx;
            this->axis_gc_y[7][a.index[2]+Lorder/2] = a.coeff_gy;
            this->axis_gc_z[7][a.index[2]+Lorder/2] = a.coeff_gz;
            this->include_axis[7] = true;
        }
        else if(a.index[1] == -a.index[2] && a.index[0] == 0)
        {
            this->plane_dists[8] = std::min(this->plane_dists[8], std::abs(a.dist));
            this->axis_lc[8][a.index[2]+Lorder/2] = a.coeff;
            this->plane_centers[8] -= a.coeff;
            this->axis_gc_x[8][a.index[2]+Lorder/2] = a.coeff_gx;
            this->axis_gc_y[8][a.index[2]+Lorder/2] = a.coeff_gy;
            this->axis_gc_z[8][a.index[2]+Lorder/2] = a.coeff_gz;
            this->include_axis[8] = true;
        }

        else if(a.index[0] == a.index[1] && a.index[0] == a.index[2])
        {
            // -1,-1,-1 to 1,1,1
            this->plane_dists[9] = std::min(this->plane_dists[9], std::abs(a.dist));
            this->axis_lc[9][a.index[2]+Lorder/2] = a.coeff;
            this->plane_centers[9] -= a.coeff;
            this->axis_gc_x[9][a.index[2]+Lorder/2] = a.coeff_gx;
            this->axis_gc_y[9][a.index[2]+Lorder/2] = a.coeff_gy;
            this->axis_gc_z[9][a.index[2]+Lorder/2] = a.coeff_gz;
            this->include_axis[9] = true;
        }
        else if(a.index[0] == a.index[1] && a.index[1] == -a.index[2])
        {
            // -1,-1,1 to 1,1,-1
            this->plane_dists[10] = std::min(this->plane_dists[10], std::abs(a.dist));
            this->axis_lc[10][a.index[2]+Lorder/2] = a.coeff;
            this->plane_centers[10] -= a.coeff;
            this->axis_gc_x[10][a.index[2]+Lorder/2] = a.coeff_gx;
            this->axis_gc_y[10][a.index[2]+Lorder/2] = a.coeff_gy;
            this->axis_gc_z[10][a.index[2]+Lorder/2] = a.coeff_gz;
            this->include_axis[10] = true;
        }
        else if(a.index[0] == -a.index[1] && a.index[0] == a.index[2])
        {
            // -1,1,-1 to 1,-1,1
            this->plane_dists[11] = std::min(this->plane_dists[11], std::abs(a.dist));
            this->axis_lc[11][a.index[2]+Lorder/2] = a.coeff;
            this->plane_centers[11] -= a.coeff;
            this->axis_gc_x[11][a.index[2]+Lorder/2] = a.coeff_gx;
            this->axis_gc_y[11][a.index[2]+Lorder/2] = a.coeff_gy;
            this->axis_gc_z[11][a.index[2]+Lorder/2] = a.coeff_gz;
            this->include_axis[11] = true;
        }
        else if(a.index[0] == -a.index[1] && a.index[1] == a.index[2])
        {
            // -1,1,1 to 1,-1,-1
            this->plane_dists[12] = std::min(this->plane_dists[12], std::abs(a.dist));
            this->axis_lc[12][a.index[0]+Lorder/2] = a.coeff;
            this->plane_centers[12] -= a.coeff;
            this->axis_gc_x[12][a.index[0]+Lorder/2] = a.coeff_gx;
            this->axis_gc_y[12][a.index[0]+Lorder/2] = a.coeff_gy;
            this->axis_gc_z[12][a.index[0]+Lorder/2] = a.coeff_gz;
            this->include_axis[12] = true;
        }

    }
}

void LaplacianCoeff::GetPointList1D(std::vector<GridPoint>& points, double a, int Ngrid, int Lorder, int direction){
    GridPoint point;
    double dx;
    for(int i = 1; i <= Lorder/2; i++){
        dx = i* a/Ngrid;
        point.dist = dx;
        point.index[0] = 0;
        point.index[1] = 0;
        point.index[2] = 0;
        point.index[direction] = i;

        point.delta[0] = 0.0;
        point.delta[1] = 0.0;
        point.delta[2] = 0.0;
        point.delta[direction] = dx;
        point.weight_factor = 1.0/std::pow(dx, this->weight_power);
        point.coeff = 0.0;
        points.push_back(point);

        point.delta[direction] = -dx;
        point.index[direction] = -i;
        points.push_back(point);

    }
}

void LaplacianCoeff::GetPointList2D(std::vector<GridPoint>& points, double a[2][2], int Ngrid[2], int Lorder){
    GridPoint point;
    double dx, dy, dist;    

    for(int i = -Lorder/2; i <= Lorder/2; i++){
        for(int j = -Lorder/2; j <= Lorder/2; j++){
            if(i == 0 && j == 0) continue;

            dx = i*a[0][0]/Ngrid[0] + j*a[1][0]/Ngrid[1];
            dy = i*a[0][1]/Ngrid[0] + j*a[1][1]/Ngrid[1];
            if((i!=0) && (j!=0))
            {
                if((this->ibrav == HEXAGONAL) && (i != j)) continue;
                if((this->ibrav == HEXAGONAL2) && (i != -j)) continue;
                if((this->ibrav == MONOCLINIC_PRIMITIVE) && (i*i != j*j)) continue;
                if((this->ibrav == -MONOCLINIC_PRIMITIVE) && (i*i != j*j)) continue;
            }
            dist = sqrt(dx * dx  + dy * dy);
            point.dist = dist;
            point.delta[0] = dx;
            point.delta[1] = dy;
            point.delta[2] = 0.0;
            point.index[0] = i;
            point.index[1] = j;
            point.index[2] = 0;
            point.ijk = std::abs(i) + std::abs(j) ;
            point.weight_factor = 1.0/std::pow(dist, this->weight_power);
            point.coeff = 0.0;
            points.push_back(point);
        }
    }
    std::stable_sort(points.begin(), points.end(), customLess_dist);

    std::vector<GridPoint> points1, points2;
    for(auto p:points){
       if((p.index[0] == 0 || p.index[1] == 0 || p.index[0]+p.index[1] == 0) && p.dist -1.0e-4 < Lorder/2 * a[0][0]/Ngrid[0])  points1.push_back(p);
        else points2.push_back(p);
   }

    points.clear();
    points.insert(points.end(), points1.begin(), points1.end());
    points.insert(points.end(), points2.begin(), points2.end());

   
}
void LaplacianCoeff::ijk_to_point(int i, int j, int k, GridPoint &point, double a[3][3], int Ngrid[3])
{
    double dx = i*a[0][0]/Ngrid[0] + j*a[1][0]/Ngrid[1] + k*a[2][0]/Ngrid[2];
    double dy = i*a[0][1]/Ngrid[0] + j*a[1][1]/Ngrid[1] + k*a[2][1]/Ngrid[2];
    double dz = i*a[0][2]/Ngrid[0] + j*a[1][2]/Ngrid[1] + k*a[2][2]/Ngrid[2];

    double dist = sqrt(dx * dx  + dy * dy + dz * dz);
    point.dist = dist;
    point.delta[0] = dx;
    point.delta[1] = dy;
    point.delta[2] = dz;
    point.index[0] = i;
    point.index[1] = j;
    point.index[2] = k;
    point.ijk = std::abs(i) + std::abs(j) + std::abs(k);
    point.weight_factor = 1.0/std::pow(dist, this->weight_power);
    point.coeff = 0.0;
}

void LaplacianCoeff::GetPointList3D(std::vector<GridPoint>& points, double a[3][3], int Ngrid[3], int Lorder){
    GridPoint point;
    points.clear();
    if(this->ibrav != TRICLINIC_PRIMITIVE)
    {
        printf("\n WARNING this is only work for Triclinic \n");
        std::raise(SIGTERM);
    }
    for(int ia = -Lorder/2; ia <= Lorder/2; ia++){
        if(ia == 0 ) continue;
        int i, j, k;
        //axis 0: a1
        i = ia;
        j = 0;
        k = 0;
        ijk_to_point(i,j,k,point, a, Ngrid);
        points.push_back(point);

        //axis 1: a2
        i = 0;
        j = ia;
        k = 0;
        ijk_to_point(i,j,k,point, a, Ngrid);
        points.push_back(point);

        //axis 2: a3
        i = 0;
        j = 0;
        k = ia;
        ijk_to_point(i,j,k,point, a, Ngrid);
        points.push_back(point);

        //axis 3: a1 + a2
        i = ia;
        j = ia;
        k = 0;
        ijk_to_point(i,j,k,point, a, Ngrid);
        points.push_back(point);

        //axis 4: a1 + a3
        i = ia;
        j = 0;
        k = ia;
        ijk_to_point(i,j,k,point, a, Ngrid);
        points.push_back(point);

        //axis 5: a2 + a3
        i = 0;
        j = ia;
        k = ia;
        ijk_to_point(i,j,k,point, a, Ngrid);
        points.push_back(point);

        //axis 6: a1 - a2
        i = ia;
        j = -ia;
        k = 0;
        ijk_to_point(i,j,k,point, a, Ngrid);
        points.push_back(point);

        //axis 7: a1 - a3
        i = ia;
        j = 0;
        k = -ia;
        ijk_to_point(i,j,k,point, a, Ngrid);
        points.push_back(point);

        //axis 8: a2 - a3
        i = 0;
        j = ia;
        k = -ia;
        ijk_to_point(i,j,k,point, a, Ngrid);

        //axis 9: a1 + a2 + a3
        i = ia;
        j = ia;
        k = ia;
        ijk_to_point(i,j,k,point, a, Ngrid);

        //axis 10: a1 + a2 - a3
        i = ia;
        j = ia;
        k = -ia;
        ijk_to_point(i,j,k,point, a, Ngrid);

        //axis 11: a1 - a2 + a3
        i = ia;
        j = -ia;
        k = ia;
        ijk_to_point(i,j,k,point, a, Ngrid);

        //axis 12: -a1 + a2 + a3
        i = -ia;
        j = ia;
        k = ia;
        ijk_to_point(i,j,k,point, a, Ngrid);
        points.push_back(point);

    }

    std::stable_sort(points.begin(), points.end(), customLess_z);
    std::stable_sort(points.begin(), points.end(), customLess_y);
    std::stable_sort(points.begin(), points.end(), customLess_x);
    std::stable_sort(points.begin(), points.end(), customLess_ijk);
    std::stable_sort(points.begin(), points.end(), customLess_dist);



}

void LaplacianCoeff::BuildSolveLinearEq(std::vector<GridPoint>& points, std::vector<GridPoint>& der_list, int dimension)
{

    int num_derivative = der_list.size();

    bool iprint = false;
    int info, index;
    std::vector<int> num_points;
    for(int i = 2*num_derivative; i < (int)points.size()-1; i++)
    {
        if( (std::abs(points[i].dist - points[i+1].dist) > 1.0e-3 ))
        {
            num_points.push_back(i+1);
        }
    }

    num_points.push_back((int)points.size());


    if(iprint)
    {
        for(auto a:num_points) std::cout << "numpoint" << "  " << a<<std::endl;

        std::cout << "num_derivative " << num_derivative << std::endl;
        index = 0;
        for(auto a:points)
        {   
            // std::cout << index <<"  "<<a.dist << "    "<<a.index[0]<<"  " <<a.index[1]<<"  " <<a.index[2] << "  "<<std::endl;
            std::cout << index <<"  "<<a.dist << "    "<<a.delta[0]<<"  " <<a.delta[1]<<"  " <<a.delta[2] << "  "<<std::endl;
            index++;
        }
        for(auto a:der_list)
        {   
            std::cout << index <<"  "<<a.dist << "    "<<a.index[0]<<"  " <<a.index[1]<<"  " <<a.index[2] << "  "<<std::endl;
            index++;
        }
    }



    double *Am = new double[num_derivative * num_derivative];
    double *Bm = new double[num_derivative * num_derivative];
    int *ipvt = new int[num_derivative];
    for(int i = 0; i < num_derivative * num_derivative; i++) 
    {
        Am[i] = 0.0;
        Bm[i] = 0.0;
    }

    double dx, dy, dz, delta_r, delta_c;

    int point_start = 0, point_end = num_points[0];
    bool Uii = false;
    index = 0;
    bool der_list_reduced = false;
    while(!Uii)
    {
        //std::cout << point_start << " point " << point_end<< std::endl;
        for(int ip = point_start; ip < point_end; ip++)
        {
            dx = points[ip].delta[0];
            dy = points[ip].delta[1];
            dz = points[ip].delta[2];

            for(int irow = 0; irow < num_derivative; irow++)
            {
                delta_r = points[ip].weight_factor;

                for(int i = 1; i <= der_list[irow].index[0]; i++) delta_r *= dx;
                for(int i = 1; i <= der_list[irow].index[1]; i++) delta_r *= dy;
                for(int i = 1; i <= der_list[irow].index[2]; i++) delta_r *= dz;

                for(int icol = 0; icol < num_derivative; icol++)
                {
                    delta_c = points[ip].weight_factor;

                    for(int i = 1; i <= der_list[icol].index[0]; i++) delta_c *= dx;
                    for(int i = 1; i <= der_list[icol].index[1]; i++) delta_c *= dy;
                    for(int i = 1; i <= der_list[icol].index[2]; i++) delta_c *= dz;

                    Am[irow * num_derivative + icol] += delta_r * delta_c;
                }
            } 

        }

        double maxd=0.0, mind=DBL_MAX;
        for(int i=0;i < num_derivative;i++)
        {
            maxd = std::abs(std::max(maxd, Am[i + i*num_derivative]));
            mind = std::abs(std::min(maxd, Am[i + i*num_derivative]));
        }
        if(mind == 0)
        {
            printf("Matrix is singular\n");
        }

        for(int i = 0; i < num_derivative * num_derivative; i++) Bm[i] = Am[i];
        dgetrf(&num_derivative, &num_derivative, Am, &num_derivative, ipvt,  &info );

        //printf("\n point_end %d", point_end);
        //for(int i = 0; i < num_derivative; i++)
        //{   printf("\n");
        //    for(int j = 0; j < num_derivative; j++)
        //        printf("%5.2f", Bm[i*num_derivative + j]);
        //}

        Uii = true;
        if(info != 0)
        {
            Uii = false;
            index++;
            if(index >= (int)num_points.size()) 
            {
                printf ("error in LaplacianCoeff.cpp dgetrf INFO = %d \n", info);
                if(info <0 || der_list_reduced)
                {
                    printf ("not enough grid points in LaplacianCoeff.cpp\n");
                    fflush (NULL);
                    std::raise(SIGTERM);
                }
                else
                {
                    printf("\nreduce number of der_list \n");

                    reduce_der_list(der_list);
                    num_derivative =der_list.size();
                    for(int i = 0; i < num_derivative * num_derivative; i++) Am[i] = 0.0;
                    //for(int i = 0; i < num_derivative; i++)
                    //{
                    //    printf("\n der_list %d %d %d %d", i, der_list[i].index[0], der_list[i].index[1], der_list[i].index[2] );
                    //}
                    point_start =0; 
                    point_end = num_points[0];
                    der_list_reduced = true;
                }
            }
            else
            {
                point_start = point_end;
                point_end = num_points[index];
                for(int i = 0; i < num_derivative * num_derivative; i++) Am[i] = Bm[i];
            }
        }

    }


    int lwork = num_derivative * num_derivative;
    dgetri(&num_derivative, Am, &num_derivative, ipvt, Bm, &lwork, &info );

    if (info != 0)
    {
        printf ("error in zgesv in LaplacianCoeff.cpp with INFO = %d \n", info);
        fflush (NULL);
        std::raise(SIGTERM);
    }

    //     for(int j = 0; j < num_derivative; j++)
    //     {   printf("\n %d ", j);
    //         for(int i = 0; i < num_derivative; i++)
    //             printf(" %f",  Am[i*num_derivative + j]);
    //     }

    for(int ip = 0; ip < point_end; ip++)
    {
        dx = points[ip].delta[0];
        dy = points[ip].delta[1];
        dz = points[ip].delta[2];
        points[ip].coeff = 0.0;
        points[ip].coeff_gx = 0.0;
        points[ip].coeff_gy = 0.0;
        points[ip].coeff_gz = 0.0;

        for(int irow = 0; irow < num_derivative; irow++)
        {
            delta_r = std::pow(points[ip].weight_factor, 2);

            for(int i = 1; i <= der_list[irow].index[0]; i++) delta_r *= dx;
            for(int i = 1; i <= der_list[irow].index[1]; i++) delta_r *= dy;
            for(int i = 1; i <= der_list[irow].index[2]; i++) delta_r *= dz;

            double tem;
            //  the Lapalaciation operators in this case is the sum of the first 3 derivatives 
            // factor of 2 is for 1/2!, expand coefficient in Taylor expansion
            tem = 0.0;
            for(int id = 0; id < num_derivative; id++)
            {
                if(der_list[id].index[0] == 2 && der_list[id].index[1] == 0 && der_list[id].index[2] == 0)
                    tem += Am[irow * num_derivative +id] *2;
                if(der_list[id].index[0] == 0 && der_list[id].index[1] == 2 && der_list[id].index[2] == 0)
                    tem += Am[irow * num_derivative +id] *2;
                if(der_list[id].index[0] == 0 && der_list[id].index[1] == 0 && der_list[id].index[2] == 2)
                    tem += Am[irow * num_derivative +id] *2;
            }

            points[ip].coeff += tem * delta_r;

            tem = 0.0;
            for(int id = 0; id < num_derivative; id++)
            {
                if(der_list[id].index[0] == 1 && der_list[id].index[1] == 0 && der_list[id].index[2] == 0)
                    tem += Am[irow * num_derivative +id];
            }

            points[ip].coeff_gx += tem * delta_r;

            tem = 0.0;
            for(int id = 0; id < num_derivative; id++)
            {
                if(der_list[id].index[0] == 0 && der_list[id].index[1] == 1 && der_list[id].index[2] == 0)
                    tem += Am[irow * num_derivative +id];
            }

            points[ip].coeff_gy += tem * delta_r;

            tem = 0.0;
            for(int id = 0; id < num_derivative; id++)
            {
                if(der_list[id].index[0] == 0 && der_list[id].index[1] == 0 && der_list[id].index[2] == 1)
                    tem += Am[irow * num_derivative +id];
            }

            points[ip].coeff_gz += tem * delta_r;
            //if(der_list[id].index[0] == 1 && der_list[id].index[1] == 0 && der_list[id].index[2] == 0)
            //    tem += Am[irow * num_derivative +id] ;
        }
    }

    points.resize(point_end);

    if(iprint)
    {
        index  =0;
        for(auto a:points)
        {   
            printf("COEFF  %2d    %14.8f    %d  %d  %d   %14.8f\n",index, a.dist, 
                    a.index[0], a.index[1], a.index[2], a.coeff);
            std::cout << index <<"  "<<a.dist << "    "<<a.index[0]<<"  " <<a.index[1]<<"  " <<a.index[2] << "  "<<a.coeff<<std::endl;
            index++;
        }
    }

    for(auto &a:points)
    {
        a.dist *= scale1;
        a.delta[0] *= scale1;
        a.delta[1] *= scale1;
        a.delta[2] *= scale1;

        a.coeff *= scale1*scale1;
        a.coeff_gx *= scale1;
        a.coeff_gy *= scale1;
        a.coeff_gz *= scale1;
    }

    delete []Am;
    delete []Bm;
    delete []ipvt;

}

void LaplacianCoeff::GetDerList(std::vector<GridPoint>& der_list, int Lorder, int dimension, int direction){

    GridPoint point;
    der_list.clear();
    if( dimension == 3){
        for(int k = 0; k <= Lorder; k+=1){
            for(int j = 0; j <= Lorder; j+=1){
                for(int i = 0; i <= Lorder; i+=1){
                    if(i+j+k <= Lorder && i+j+k >0)
                    {
                        point.index[0] = i;
                        point.index[1] = j;
                        point.index[2] = k;
                        point.dist = i+j+k;
                        point.ijk = 1;
                        for(int ix = 1; ix <= i; ix++) point.ijk *=ix;
                        for(int ix = 1; ix <= j; ix++) point.ijk *=ix;
                        for(int ix = 1; ix <= k; ix++) point.ijk *=ix;
                        der_list.push_back(point);
                    }
                }
            }
        }
    }

    if( dimension == 2){
        for(int j = 0; j <= Lorder; j++){
            for(int i = 0; i <= Lorder; i++){
                // if(this->ibrav == HEXAGONAL && i !=0 && j !=0) continue;
                if(i+j <= Lorder && i+j >0  )
                {
                    point.index[0] = i;
                    point.index[1] = j;
                    point.index[2] = 0;
                    point.dist = i+j;
                    point.ijk = 1;
                    for(int ix = 1; ix <= i; ix++) point.ijk *=ix;
                    for(int ix = 1; ix <= j; ix++) point.ijk *=ix;
                    der_list.push_back(point);
                }
            }
        }
    }

    if( dimension == 1){
        for(int i = 1; i <= Lorder; i+=1){
            {
                point.index[0] = 0;
                point.index[1] = 0;
                point.index[2] = 0;
                point.index[direction] = i;
                point.dist = i;
                point.ijk = 1;
                for(int ix = 1; ix <= i; ix++) point.ijk *=ix;
                der_list.push_back(point);
            }
        }
    }

    std::stable_sort(der_list.begin(), der_list.end(), customLess_ijk);
    std::stable_sort(der_list.begin(), der_list.end(), customLess_dist);

}

void LaplacianCoeff::GetDerListFCC(std::vector<GridPoint>& der_list, int Lorder)
{
    der_list.clear();

    std::vector<GridPoint> der_list1;
    GetDerList(der_list1, Lorder, 1, 0);
    der_list.insert(der_list.end(), der_list1.begin(), der_list1.end());
    GetDerList(der_list1, Lorder, 1, 1);
    der_list.insert(der_list.end(), der_list1.begin(), der_list1.end());
    GetDerList(der_list1, Lorder, 1, 2);
    der_list.insert(der_list.end(), der_list1.begin(), der_list1.end());

    std::stable_sort(der_list.begin(), der_list.end(), customLess_ijk);
    std::stable_sort(der_list.begin(), der_list.end(), customLess_dist);

}

void LaplacianCoeff::GetDerListBCC(std::vector<GridPoint>& der_list, int Lorder)
{
    der_list.clear();

#if 0
    std::vector<GridPoint> der_list1;
    GetDerList(der_list1, Lorder, 3, 0);
    der_list.insert(der_list.end(), der_list1.begin(), der_list1.end());
    GetDerList(der_list1, Lorder, 1, 1);
    der_list.insert(der_list.end(), der_list1.begin(), der_list1.end());
    GetDerList(der_list1, Lorder,  2);
    der_list.insert(der_list.end(), der_list1.begin(), der_list1.end());
#endif
    GetDerList(der_list, Lorder, 3, 0);
    printf("SIZE = %lu\n",der_list.size());
    std::stable_sort(der_list.begin(), der_list.end(), customLess_ijk);
    std::stable_sort(der_list.begin(), der_list.end(), customLess_dist);

}

void LaplacianCoeff::GetPointListFCC(std::vector<GridPoint>& points, double a[3][3], int Ngrid[3], int Lorder){
    GridPoint point;
    double dx, dy,dz, dist;    
    // along lattuce vectirs a1 = (0,5, 0.5, 0), a2 = (0.5, 0.0, 0.5), and a3 = (0.0, 0.5, 0.5)
    for(int direction = 0; direction < 3; direction++){
        for(int i = -Lorder/2; i <= Lorder/2; i++){

            if (i == 0) continue;
            dx = i*a[direction][0]/Ngrid[0];
            dy = i*a[direction][1]/Ngrid[0];
            dz = i*a[direction][2]/Ngrid[0];

            dist = sqrt(dx * dx  + dy * dy + dz * dz);
            point.dist = dist;
            point.delta[0] = dx;
            point.delta[1] = dy;
            point.delta[2] = dz;
            point.index[0] = 0;
            point.index[1] = 0;
            point.index[2] = 0;
            point.index[direction] = i;
            point.ijk = std::abs(i);
            point.weight_factor = 1.0/std::pow(dist, this->weight_power);
            point.coeff = 0.0;
            points.push_back(point);

            int direction1 = (direction + 1) %3;

            dx = i*a[direction][0]/Ngrid[0] - i*a[direction1][0]/Ngrid[0];
            dy = i*a[direction][1]/Ngrid[0] - i*a[direction1][1]/Ngrid[0];
            dz = i*a[direction][2]/Ngrid[0] - i*a[direction1][2]/Ngrid[0];

            dist = sqrt(dx * dx  + dy * dy + dz * dz);
            point.dist = dist;
            point.delta[0] = dx;
            point.delta[1] = dy;
            point.delta[2] = dz;
            point.index[0] = 0;
            point.index[1] = 0;
            point.index[2] = 0;
            point.index[direction] = i;
            point.index[direction1] = -i;
            point.ijk = std::abs(i);
            point.weight_factor = 1.0/std::pow(dist, this->weight_power);
            point.coeff = 0.0;
            points.push_back(point);

        }
    }

    std::stable_sort(points.begin(), points.end(), customLess_dist);

}

void LaplacianCoeff::GetPointListBCC(std::vector<GridPoint>& points, double a[3][3], int Ngrid[3], int Lorder){
    GridPoint point;
    double dx, dy,dz, dist;    
    for(int i = -Lorder/2; i <= Lorder/2; i++){
        for(int j = -Lorder/2; j <= Lorder/2; j++){
            for(int k = -Lorder/2; k <= Lorder/2; k++){
                if(i == 0 && j == 0 && k == 0) continue;
                dx = i*a[0][0]/Ngrid[0] + j*a[1][0]/Ngrid[1] + k*a[2][0]/Ngrid[2];
                dy = i*a[0][1]/Ngrid[0] + j*a[1][1]/Ngrid[1] + k*a[2][1]/Ngrid[2];
                dz = i*a[0][2]/Ngrid[0] + j*a[1][2]/Ngrid[1] + k*a[2][2]/Ngrid[2];
                if(
                        (i==j && j==k) ||
                        (i!=0 && j==0 && k==0) ||
                        (i==0 && j!=0 && k==0) ||
                        (i==0 && j==0 && k!=0))
                {
                    dist = sqrt(dx * dx  + dy * dy + dz * dz);
                    point.dist = dist;
                    point.delta[0] = dx;
                    point.delta[1] = dy;
                    point.delta[2] = dz;
                    point.index[0] = i;
                    point.index[1] = j;
                    point.index[2] = k;
                    point.ijk = std::abs(i) + std::abs(j) + std::abs(k);
                    point.weight_factor = 1.0/std::pow(dist, this->weight_power);
                    point.coeff = 0.0;
                    points.push_back(point);
                }
            }
        }
    }

    std::stable_sort(points.begin(), points.end(), customLess_z);
    std::stable_sort(points.begin(), points.end(), customLess_y);
    std::stable_sort(points.begin(), points.end(), customLess_x);
    std::stable_sort(points.begin(), points.end(), customLess_ijk);
    std::stable_sort(points.begin(), points.end(), customLess_dist);
}

void LaplacianCoeff::reduce_der_list(std::vector<GridPoint>& der_list)
{
    int num_derivative = der_list.size();
    for(int idx = num_derivative -1; idx >=0; idx--)
    {
        int i = der_list[idx].index[0];
        int j = der_list[idx].index[1];
        int k = der_list[idx].index[2];

        if (i * j * k != 0)
        {
            if (i == j && i == k) continue;
            der_list.erase(der_list.begin() + idx);
        }
        else if( i * j != 0 && i != j)
        {
            der_list.erase(der_list.begin() + idx);
        }
        else if( i * k != 0 && i != k)
        {
            der_list.erase(der_list.begin() + idx);
        }
        else if( j * k != 0 && j != k)
        {
            der_list.erase(der_list.begin() + idx);
        }
    }
}

LaplacianCoeff::~LaplacianCoeff(void)
{
}
