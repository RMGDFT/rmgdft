#ifndef RMG_LaplacianCoeff_H
#define RMG_LaplacianCoeff_H 1
#include "Lattice.h"
#include <vector>
#include <complex>
// with input of base vectors a[3][3] and required order of Lapalacian operator, find needed neighbors and their coefficients. 

struct GridPoint{
    double dist;
    double delta[3];
    int index[3];
    int ijk;
    double weight_factor;
    double coeff;
    double coeff_gx;
    double coeff_gy;
    double coeff_gz;
};
typedef GridPoint GridPoint;


class LaplacianCoeff {

private:

    double a[3][3];
    int Ngrid[3];
    int dim[3];
    int weight_power = 6;
    bool offdiag = true;
    int ibrav = ORTHORHOMBIC_PRIMITIVE;
    
    int iprint, world_rank;

public:
    

    int Lorder;
    double gen_hxgrid;    // a[0][0] / Ngrid[0] used to identify scale factor in multigrid routines
    LaplacianCoeff ();
    LaplacianCoeff (double a[3][3], int Ngrid[3], int Lorder, int dim[3]);
    void CalculateCoeff (double a[3][3], int Ngrid[3], int Lorder, int dim[3]);
    void CalculateCoeff ();
    ~LaplacianCoeff(void);

    // For laplacian
    // 0=x,1=y,2=z,3=xy,4=xz,5=yz,6=nxy,7=nxz,8=nyz
    double axis_lc[13][12];
    double plane_centers[13];
    double plane_dists[13];
    bool include_axis[13];

    // For gradients
    double axis_gc_x[13][12];
    double axis_gc_y[13][12];
    double axis_gc_z[13][12];

    // For rescaling to improve numercal stability
    double scale1;

    void SetLattice(double a[3][3])
    {
        for(int i = 0; i < 3; i++)
            for(int j = 0; j < 3; j++)
                this->a[i][j] = a[i][j];
    }
    void SetWeightPower(int wp){ this->weight_power = wp;}

    void SetOffdiag(bool flag){
         this->offdiag = flag;
        }
    void SetBrav(int ibrav){
         this->ibrav = ibrav;
        }
    void SetOrder(int Lorder){this->Lorder = Lorder;}
    void SetNgrid(int Ngrid[3]){
        for(int i = 0; i < 3; i++) this->Ngrid[i] = Ngrid[i];
    }
    void SetDim(int dim[3]){
        for(int i = 0; i < 3; i++) this->dim[i] = dim[i];
    }

    int GetOrder(){return this->Lorder;}
    void GetDim(int *dim){
        for(int i = 0; i < 3; i++) dim[i] = this->dim[i];
    }
    
    void UpdateIndex(int dim[3]);

    void BuildSolveLinearEq(std::vector<GridPoint>& points, std::vector<GridPoint>& der_list, int dimension);
    void GetDerList(std::vector<GridPoint>& der_list, int Lorder, int dimension, int direction);
    void GetPointList3D (std::vector<GridPoint>& points, double a[3][3], int Ngrid[3], int Lorder);
    void GetPointList2D (std::vector<GridPoint>& points, double a[3][3], int Ngrid[3], int Lorder);
    void GetPointList1D (std::vector<GridPoint>& points, double a, int Ngrid, int Lorder, int direction);

    void GetPointListFCC (std::vector<GridPoint>& points, double a[3][3], int Ngrid[3], int Lorder);
    void GetPointListBCO (std::vector<GridPoint>& points, double a[3][3], int Ngrid[3], int Lorder);
    void GetPointListBCC (std::vector<GridPoint>& points, double a[3][3], int Ngrid[3], int Lorder);
    void GetDerListFCC(std::vector<GridPoint>& der_list, int Lorder);
    void GetDerListBCO(std::vector<GridPoint>& der_list, int Lorder);
    void GetDerListBCC(std::vector<GridPoint>& der_list, int Lorder);
    void ijk_to_point(int i, int j, int k, GridPoint &point, double a[3][3], int Ngrid[3]);
    void reduce_der_list(std::vector<GridPoint>& der_list);

};

template <typename T>
double FiniteDiffLap(T * __restrict__ a, T * __restrict__ b, int dimx, int dimy, int dimz, LaplacianCoeff *LC);

template <typename T>
void FiniteDiffGrad(T * __restrict__ a, T * __restrict__ gx, T * __restrict__ gy, T * __restrict__ gz, int dimx, int dimy, int dimz, LaplacianCoeff *LC);

extern LaplacianCoeff *LC;
extern LaplacianCoeff *LC_6;
extern LaplacianCoeff *LC_4;

#endif

