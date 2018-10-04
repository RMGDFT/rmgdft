#ifndef RMG_LaplacianCoeff_H
#define RMG_LaplacianCoeff_H 1
#include "const.h"
#include <vector>
// with input of base vectors a[3][3] and required order of Lapalacian operator, find needed neighbors and their coefficients. 

typedef struct {
    double coeff;
    std::vector<int> i,j,k,relative_index;
} CoeffList; 

struct GridPoint{
    double dist;
    double dx,dy,dz;
    int i, j, k, eq_num;
    int ijk;
    double weight_factor;
    double coeff;
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

    std::vector<CoeffList> coeff_and_index;
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
    void GetDerList(std::vector<GridPoint>& der_list, int Lorder, int dimension);
    void GetPointList3D (std::vector<GridPoint>& points, double a[3][3], int Ngrid[3], int Lorder);
    void GetPointList2D (std::vector<GridPoint>& points, double a[2][2], int Ngrid[2], int Lorder);
    void GetPointList1D (std::vector<GridPoint>& points, double a, int Ngrid, int Lorder);


};
#endif

