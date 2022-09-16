#ifndef RMG_LaplacianCoeff_H
#define RMG_LaplacianCoeff_H 1
#include "Lattice.h"
#include <vector>
#include <complex>
// with input of base vectors a[3][3] and required order of Lapalacian operator, find needed neighbors and their coefficients. 

typedef struct {
    double coeff;
    std::vector<int> i,j,k,relative_index;
} CoeffList; 

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
    double coeff_gxy;
    double coeff_gxz;
    double coeff_gyz;
};
typedef GridPoint GridPoint;


class LaplacianCoeff {

private:

    double a[3][3];
    int Ngrid[3];
    int Lorder;
    int dim[3];
    int weight_power = 3;
    bool offdiag = true;
    int ibrav = ORTHORHOMBIC_PRIMITIVE;
    

public:
    

    LaplacianCoeff ();
    LaplacianCoeff (double a[3][3], int Ngrid[3], int Lorder, int dim[3]);
    void CalculateCoeff (double a[3][3], int Ngrid[3], int Lorder, int dim[3]);
    void CalculateCoeff ();
    ~LaplacianCoeff(void);

    // Used for monoclinic
    double plane_dist_x, plane_dist_y, plane_dist_z;
    double plane_dist_xy, plane_dist_xz, plane_dist_yz;
    double plane_dist_nxy, plane_dist_nxz, plane_dist_nyz;
    double plane_center_x, plane_center_y, plane_center_z;
    double plane_center_xy, plane_center_xz, plane_center_yz;
    double plane_center_nxy, plane_center_nxz, plane_center_nyz;

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

    std::vector<double> axis_x;
    std::vector<double> axis_y;
    std::vector<double> axis_z;
    std::vector<double> axis_xy;
    std::vector<double> axis_xz;
    std::vector<double> axis_yz;
    std::vector<double> axis_nxy;
    std::vector<double> axis_nxz;
    std::vector<double> axis_nyz;

    std::vector<double> axis_x_gx;
    std::vector<double> axis_x_gy;
    std::vector<double> axis_x_gz;
    std::vector<double> axis_y_gx;
    std::vector<double> axis_y_gy;
    std::vector<double> axis_y_gz;
    std::vector<double> axis_z_gx;
    std::vector<double> axis_z_gy;
    std::vector<double> axis_z_gz;
    std::vector<double> axis_xy_gx;
    std::vector<double> axis_xy_gy;
    std::vector<double> axis_xy_gz;
    std::vector<double> axis_xz_gx;
    std::vector<double> axis_xz_gy;
    std::vector<double> axis_xz_gz;
    std::vector<double> axis_yz_gx;
    std::vector<double> axis_yz_gy;
    std::vector<double> axis_yz_gz;

    std::vector<double> axis_nxy_gx;
    std::vector<double> axis_nxy_gy;
    std::vector<double> axis_nxy_gz;
    std::vector<double> axis_nxz_gx;
    std::vector<double> axis_nxz_gy;
    std::vector<double> axis_nxz_gz;
    std::vector<double> axis_nyz_gx;
    std::vector<double> axis_nyz_gy;
    std::vector<double> axis_nyz_gz;

    std::vector<CoeffList> coeff_and_index;
    std::vector<CoeffList> gx_coeff_and_index;
    std::vector<CoeffList> gy_coeff_and_index;
    std::vector<CoeffList> gz_coeff_and_index;
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

    void GenerateList(const std::vector<GridPoint>& points);
    void BuildSolveLinearEq(std::vector<GridPoint>& points, const std::vector<GridPoint>& der_list, int dimension);
    void GetDerList(std::vector<GridPoint>& der_list, int Lorder, int dimension, int direction);
    void GetPointList3D (std::vector<GridPoint>& points, double a[3][3], int Ngrid[3], int Lorder);
    void GetPointList2D (std::vector<GridPoint>& points, double a[2][2], int Ngrid[2], int Lorder);
    void GetPointList1D (std::vector<GridPoint>& points, double a, int Ngrid, int Lorder, int direction);

    void GetPointListFCC (std::vector<GridPoint>& points, double a[3][3], int Ngrid[3], int Lorder);
    void GetPointListBCC (std::vector<GridPoint>& points, double a[3][3], int Ngrid[3], int Lorder);
    void GetDerListFCC(std::vector<GridPoint>& der_list, int Lorder);
    void GetDerListBCC(std::vector<GridPoint>& der_list, int Lorder);

};

template <typename T>
double FiniteDiffLap(T * __restrict__ a, T * __restrict__ b, int dimx, int dimy, int dimz, LaplacianCoeff *LC);

template <typename T>
void FiniteDiffGrad(T * __restrict__ a, T * __restrict__ gx, T * __restrict__ gy, T * __restrict__ gz, int dimx, int dimy, int dimz, LaplacianCoeff *LC);

extern LaplacianCoeff *LC;
extern LaplacianCoeff *LC_6;
extern LaplacianCoeff *HLC;

#endif

