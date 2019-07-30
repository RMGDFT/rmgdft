#include <vector>
#include <algorithm>
#include <math.h>
#include <iostream>
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
        return std::abs(a.dx) > std::abs(b.dx) ;
    }   
} customLess_x;
struct {
    bool operator()(GridPoint a, GridPoint b) const
    {   
        return std::abs(a.dy) > std::abs(b.dy);
    }   
} customLess_y;
struct {
    bool operator()(GridPoint a, GridPoint b) const
    {   
        return std::abs(a.dz) > std::abs(b.dz);
    }   
} customLess_z;

struct {
    bool operator()(GridPoint a, GridPoint b) const
    {   
        return std::abs(a.ijk) < std::abs(b.ijk);
    }   
} customLess_ijk;

LaplacianCoeff::LaplacianCoeff(double a[3][3], int Ngrid[3], int Lorder, int dim[3]){
    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
            this->a[i][j] = a[i][j];
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

    std::vector<GridPoint>  points;

    std::vector<GridPoint> der_list;  // only i,j,k are needed.
    int dimension;
    // for laplacian, 2n order is also 2n+1 order since the odd partianl derivatives cancelled automatically.
    Lorder = Lorder /2 * 2;

    if( (this->ibrav == ORTHORHOMBIC_PRIMITIVE || this->ibrav == CUBIC_PRIMITIVE) && !this->offdiag)
    {
        std::vector<GridPoint>  points1;
        GridPoint point;
        dimension = 1;
        GetDerList(der_list, Lorder, dimension);
        GetPointList1D(points1, a[0][0], Ngrid[0], Lorder);
        this->BuildSolveLinearEq(points1, der_list, dimension);

        for(auto p:points1)
        {
            point.i = p.i;
            point.j = 0;
            point.k = 0;
            point.coeff = p.coeff;
            point.dist = p.dist;
            points.push_back(point);
        }


        points1.clear();
        GetPointList1D(points1, a[1][1], Ngrid[1], Lorder);
        this->BuildSolveLinearEq(points1, der_list, dimension);

        for(auto p:points1)
        {
            point.i = 0;
            point.j = p.i;
            point.k = 0;
            point.coeff = p.coeff;
            point.dist = p.dist;
            points.push_back(point);
        }


        points1.clear();
        GetPointList1D(points1, a[2][2], Ngrid[2], Lorder);
        this->BuildSolveLinearEq(points1, der_list, dimension);

        for(auto p:points1)
        {
            point.i = 0;
            point.j = 0;
            point.k = p.i;
            point.coeff = p.coeff;
            point.dist = p.dist;
            points.push_back(point);
        }

        std::stable_sort(points.begin(), points.end(), customLess_dist);

    }

    else if ( (this->ibrav == HEXAGONAL) && !this->offdiag)
    {
        std::vector<GridPoint>  points1;
        GridPoint point;
        dimension = 2;
        GetDerList(der_list, Lorder, dimension);
        double a2d[2][2];
        a2d[0][0] = a[0][0];
        a2d[0][1] = a[0][1];
        a2d[1][0] = a[1][0];
        a2d[1][1] = a[1][1];

        GetPointList2D(points1, a2d, Ngrid, Lorder);
        this->BuildSolveLinearEq(points1, der_list, dimension);

        for(auto p:points1)
        {
            point.i = p.i;
            point.j = p.j;
            point.k = 0;
            point.coeff = p.coeff;
            point.dist = p.dist;
            points.push_back(point);
        }


        points1.clear();
        dimension = 1;
        GetDerList(der_list, Lorder, dimension);
        GetPointList1D(points1, a[2][2], Ngrid[2], Lorder);
        this->BuildSolveLinearEq(points1, der_list, dimension);

        for(auto p:points1)
        {
            point.i = 0;
            point.j = 0;
            point.k = p.i;
            point.coeff = p.coeff;
            point.dist = p.dist;
            points.push_back(point);
        }

        std::stable_sort(points.begin(), points.end(), customLess_dist);
    }
    else
    {
        dimension = 3;
        GetDerList(der_list, Lorder, dimension);

        GetPointList3D(points, a, Ngrid, Lorder);

        this->BuildSolveLinearEq(points, der_list, dimension);
    }

    this->GenerateList(points);
}

void LaplacianCoeff::GetPointList1D(std::vector<GridPoint>& points, double a, int Ngrid, int Lorder){
    GridPoint point;
    double dx, dy,dz, dist;    
    for(int i = 1; i <= Lorder/2; i++){
        dx = i* a/Ngrid;
        point.dist = dx;
        point.i = i;
        point.j = 0;
        point.k = 0;

        point.dx = dx;
        point.dy = 0.0;
        point.dz = 0.0;
        point.weight_factor = 1.0/std::pow(dx, this->weight_power);
        point.coeff = 0.0;
        points.push_back(point);

        point.dx = -dx;
        point.i = -i;
        points.push_back(point);

    }
}

void LaplacianCoeff::GetPointList2D(std::vector<GridPoint>& points, double a[2][2], int Ngrid[2], int Lorder){
    GridPoint point;
    double dx, dy, dist;    
    for(int i = -Lorder; i <= Lorder; i++){
        for(int j = -Lorder; j <= Lorder; j++){
            if(i == 0 && j == 0) continue;

            dx = i*a[0][0]/Ngrid[0] + j*a[1][0]/Ngrid[1];
            dy = i*a[0][1]/Ngrid[0] + j*a[1][1]/Ngrid[1];

            dist = sqrt(dx * dx  + dy * dy);
            point.dist = dist;
            point.dx = dx;
            point.dy = dy;
            point.dz = 0.0;
            point.i = i;
            point.j = j;
            point.k = 0;
            point.ijk = std::abs(i) + std::abs(j) ;
            point.weight_factor = 1.0/std::pow(dist, this->weight_power);
            point.coeff = 0.0;
            points.push_back(point);
        }
    }
    std::stable_sort(points.begin(), points.end(), customLess_dist);

    std::vector<GridPoint> points1, points2;
    for(auto p:points){
        if((p.i == 0 || p.j == 0 || p.i+p.j == 0) && p.dist -1.0e-4 < Lorder/2 * a[0][0]/Ngrid[0])  points1.push_back(p);
        else points2.push_back(p);
   }

    points.clear();
    points.insert(points.end(), points1.begin(), points1.end());
    points.insert(points.end(), points2.begin(), points2.end());

   
}


void LaplacianCoeff::GetPointList3D(std::vector<GridPoint>& points, double a[3][3], int Ngrid[3], int Lorder){
    GridPoint point;
    double dx, dy,dz, dist;    
    for(int i = -Lorder; i <= Lorder; i++){
        for(int j = -Lorder; j <= Lorder; j++){
            for(int k = -Lorder; k <= Lorder; k++){
                if(i == 0 && j == 0 && k == 0) continue;
                if( (this->ibrav == ORTHORHOMBIC_PRIMITIVE || this->ibrav == CUBIC_PRIMITIVE) && 
                        Lorder > 6 && (i +j == 0 || i+k == 0 || j+k  == 0)) continue;

                dx = i*a[0][0]/Ngrid[0] + j*a[1][0]/Ngrid[1] + k*a[2][0]/Ngrid[2];
                dy = i*a[0][1]/Ngrid[0] + j*a[1][1]/Ngrid[1] + k*a[2][1]/Ngrid[2];
                dz = i*a[0][2]/Ngrid[0] + j*a[1][2]/Ngrid[1] + k*a[2][2]/Ngrid[2];

                dist = sqrt(dx * dx  + dy * dy + dz * dz);
                point.dist = dist;
                point.dx = dx;
                point.dy = dy;
                point.dz = dz;
                point.i = i;
                point.j = j;
                point.k = k;
                point.ijk = std::abs(i) + std::abs(j) + std::abs(k);
                point.weight_factor = 1.0/std::pow(dist, this->weight_power);
                point.coeff = 0.0;
                points.push_back(point);
            }
        }
    }
    std::stable_sort(points.begin(), points.end(), customLess_z);
    std::stable_sort(points.begin(), points.end(), customLess_y);
    std::stable_sort(points.begin(), points.end(), customLess_x);
    std::stable_sort(points.begin(), points.end(), customLess_ijk);
    std::stable_sort(points.begin(), points.end(), customLess_dist);

    if( (this->ibrav == ORTHORHOMBIC_PRIMITIVE || this->ibrav == CUBIC_PRIMITIVE) && Lorder > 6)
        for(int i = Lorder/2; i > 0; i--)
        {
            for(int j = 0; j < 3; j++){

                dx = i*a[j][0]/Ngrid[j];
                dy = i*a[j][1]/Ngrid[j];
                dz = i*a[j][2]/Ngrid[j];

                dist = sqrt(dx * dx  + dy * dy + dz * dz);
                point.dist = dist;
                point.dx = dx;
                point.dy = dy;
                point.dz = dz;
                if(j == 0){
                    point.i = i;
                    point.j = 0;
                    point.k = 0;
                }
                if(j == 1){
                    point.i = 0;
                    point.j = i;
                    point.k = 0;
                }
                if(j == 2){
                    point.i = 0;
                    point.j = 0;
                    point.k = i;
                }
                point.ijk = std::abs(i);
                point.weight_factor = 1.0/std::pow(dist, this->weight_power);
                point.coeff = 0.0;
                points.insert(points.begin(), point);

                dx = -i*a[j][0]/Ngrid[j];
                dy = -i*a[j][1]/Ngrid[j];
                dz = -i*a[j][2]/Ngrid[j];

                dist = sqrt(dx * dx  + dy * dy + dz * dz);
                point.dist = dist;
                point.dx = dx;
                point.dy = dy;
                point.dz = dz;
                if(j == 0){
                    point.i = -i;
                    point.j = 0;
                    point.k = 0;
                }
                if(j == 1){
                    point.i = 0;
                    point.j = -i;
                    point.k = 0;
                }
                if(j == 2){
                    point.i = 0;
                    point.j = 0;
                    point.k = -i;
                }
                point.ijk = std::abs(i);
                point.weight_factor = 1.0/std::pow(dist, this->weight_power);
                point.coeff = 0.0;
                points.insert(points.begin(), point);
            }
        }

    if( (this->ibrav == ORTHORHOMBIC_PRIMITIVE || this->ibrav == CUBIC_PRIMITIVE) && Lorder == 6)
    {
        points.insert(points.begin()+79, points[97]);
        points.insert(points.begin()+79, points[97]);
        points.insert(points.begin()+79, points[97]);
        points.insert(points.begin()+79, points[97]);
        points.insert(points.begin()+79, points[97]);
        points.insert(points.begin()+79, points[97]);
    }


}

void LaplacianCoeff::BuildSolveLinearEq(std::vector<GridPoint>& points, const std::vector<GridPoint>& der_list, int dimension){

    int num_derivative = der_list.size();
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


//      for(auto a:num_points) std::cout << "numpoint" << "  " << a<<std::endl;
//
//            index = 0;
//           for(auto a:points)
//          {   
//             std::cout << index <<"  "<<a.dist << "    "<<a.i<<"  " <<a.j<<"  " <<a.k << "  "<<std::endl;
//            index++;
//       }



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
    while(!Uii)
    {
        for(int ip = point_start; ip < point_end; ip++)
        {
            dx = points[ip].dx;
            dy = points[ip].dy;
            dz = points[ip].dz;

            for(int irow = 0; irow < num_derivative; irow++)
            {
                delta_r = points[ip].weight_factor;

                for(int i = 1; i <= der_list[irow].i; i++) delta_r *= dx;
                for(int i = 1; i <= der_list[irow].j; i++) delta_r *= dy;
                for(int i = 1; i <= der_list[irow].k; i++) delta_r *= dz;


                for(int icol = 0; icol < num_derivative; icol++)
                {
                    delta_c = points[ip].weight_factor;

                    for(int i = 1; i <= der_list[icol].i; i++) delta_c *= dx;
                    for(int i = 1; i <= der_list[icol].j; i++) delta_c *= dy;
                    for(int i = 1; i <= der_list[icol].k; i++) delta_c *= dz;

                    Am[irow * num_derivative + icol] += delta_r * delta_c;
                }
            } 

        }

        for(int i = 0; i < num_derivative * num_derivative; i++) Bm[i] = Am[i];
        dgetrf(&num_derivative, &num_derivative, Am, &num_derivative, ipvt,  &info );

        // printf("\n point_end %d", point_end);
        //  for(int i = 0; i < num_derivative; i++)
        //  {   printf("\n");
        //  for(int j = 0; j < num_derivative; j++)
        //      printf("%d %e", i, Bm[i*num_derivative + j]);
        //   }

        Uii = true;
        for(int i = 0; i < num_derivative; i++)
        {
            if(std::abs(Am[i * num_derivative + i] ) < 1.0e-10)
            {
                index++;
                point_start = point_end;
                point_end = num_points[index];
                Uii = false;
                for(int i = 0; i < num_derivative * num_derivative; i++) Am[i] = Bm[i];
                break;
            }
        }

    }


    int lwork = num_derivative * num_derivative;
    dgetri(&num_derivative, Am, &num_derivative, ipvt, Bm, &lwork, &info );

    if (info != 0)
    {
        printf ("error in zgesv in LaplacianCoeff.cppwith INFO = %d \n", info);
        fflush (NULL);
        exit (0);
    }

    //    for(int j = 0; j < num_derivative; j++)
    //    {   printf("\n %d ", j);
    //        for(int i = 0; i < 3; i++)
    //            printf(" %f",  Am[i*num_derivative + j]);
    //    }

    double coeff0 = 0.0;
    for(int ip = 0; ip < point_end; ip++)
    {
        dx = points[ip].dx;
        dy = points[ip].dy;
        dz = points[ip].dz;

        for(int irow = 0; irow < num_derivative; irow++)
        {
            delta_r = std::pow(points[ip].weight_factor, 2);

            for(int i = 1; i <= der_list[irow].i; i++) delta_r *= dx;
            for(int i = 1; i <= der_list[irow].j; i++) delta_r *= dy;
            for(int i = 1; i <= der_list[irow].k; i++) delta_r *= dz;

            double tem;
            //  the Lapalaciation operators in this case is the sum of the first 3 derivatives 
            // factor of 2 is for 1/2!, expand coefficient in Taylor expansion
            tem = 0.0;
            for(int id = 0; id < dimension; id++)
            {
                tem += Am[irow * num_derivative +id] *2;
            }

            points[ip].coeff += tem * delta_r;
            coeff0 -= tem * delta_r;
        }
    }

    points.resize(point_end);

//        index  =0;
//        for(auto a:points)
//        {   
//            std::cout << index <<"  "<<a.dist << "    "<<a.i<<"  " <<a.j<<"  " <<a.k << "  "<<a.coeff<<std::endl;
//            index++;
//         }


    delete []Am;
    delete []Bm;
    delete []ipvt;

}
void LaplacianCoeff::GenerateList(const std::vector<GridPoint>&  points)
{
    this->coeff_and_index.clear();
    CoeffList coeff_list;
    int  index;
    coeff_list.coeff = points[0].coeff;
    coeff_list.i.push_back(points[0].i);
    coeff_list.j.push_back(points[0].j);
    coeff_list.k.push_back(points[0].k);
    index = points[0].i * (dim[1]+Lorder) * (dim[2]+Lorder) + points[0].j * (dim[2]+Lorder) + points[0].k;
    coeff_list.relative_index.push_back(index);

    for(int ip = 1; ip < (int)points.size(); ip++)
    {
        if(std::abs(points[ip].coeff) < 1.0e-8) continue;

        if(std::abs(points[ip].coeff - points[ip-1].coeff) > 1.0e-8)
        {
            this->coeff_and_index.push_back(coeff_list);

            //start a new coeff and its index(s) 
            coeff_list.relative_index.clear();
            coeff_list.i.clear();
            coeff_list.j.clear();
            coeff_list.k.clear();
            coeff_list.coeff = points[ip].coeff;
        }

        index = points[ip].i * (dim[1]+Lorder) * (dim[2]+Lorder) + points[ip].j * (dim[2]+Lorder) + points[ip].k;
        coeff_list.relative_index.push_back(index);
        coeff_list.i.push_back(points[ip].i);
        coeff_list.j.push_back(points[ip].j);
        coeff_list.k.push_back(points[ip].k);
    }

    this->coeff_and_index.push_back(coeff_list);

    double coeff0 = 0;
    index = 0;
    for(auto a:this->coeff_and_index)
    {
        coeff0 += a.coeff * a.relative_index.size();
        index += a.relative_index.size();
        //    std::cout << a.coeff << "  "<< a.relative_index.size() <<std::endl;
        //    for(auto b:a.relative_index) std::cout <<"            "<< b<<std::endl;
    }


    // add the original (0,0,0) point
    coeff_list.i.clear();
    coeff_list.j.clear();
    coeff_list.k.clear();
    coeff_list.relative_index.clear();

    coeff_list.coeff = - coeff0;
    coeff_list.i.push_back(0);
    coeff_list.j.push_back(0);
    coeff_list.k.push_back(0);
    coeff_list.relative_index.push_back(0);
    this->coeff_and_index.insert(this->coeff_and_index.begin(), coeff_list);

    for(auto a:this->coeff_and_index)
    {
        //std::cout << a.coeff << std::endl;
        //for (int ip = 0; ip < (int)a.i.size(); ip++)
        //    std::cout<<"      " << a.i[ip] <<"  "<<a.j[ip]<<"  "<<a.k[ip]<<std::endl;
    } 
    //    std::cout << "total points= " <<this->coeff_and_index.size() <<std::endl;
    //    std::cout << "sum of coeffs = " <<coeff0 <<std::endl;
}

void LaplacianCoeff::UpdateIndex(int dim[3])
{

    //printf("\n dim %d %d %d", dim[0], dim[1], dim[2]);
    for(int ip = 0; ip < (int)this->coeff_and_index.size(); ip++)
    {

        for(int idx = 0; idx < (int)this->coeff_and_index[ip].i.size(); idx++)
        {

            this->coeff_and_index[ip].relative_index[idx] =
                this->coeff_and_index[ip].i[idx] * (dim[1]+this->Lorder) * (dim[2]+this->Lorder)
                + this->coeff_and_index[ip].j[idx] * (dim[2]+this->Lorder) + this->coeff_and_index[ip].k[idx];
        }

    }

    this->SetDim(dim);
    //  for(auto a:this->coeff_and_index)
    //  {
    //  printf("\n coff %e ", a.coeff);
    //  for(auto b:a.relative_index) printf("\n    %d ", b);
    //  }

}

void LaplacianCoeff::GetDerList(std::vector<GridPoint>& der_list, int Lorder, int dimension){

    GridPoint point;
    der_list.clear();
    if( dimension == 3){
        for(int k = 0; k <= Lorder; k+=2){
            for(int j = 0; j <= Lorder; j+=2){
                for(int i = 0; i <= Lorder; i+=2){
                    if(i+j+k <= Lorder && i+j+k >0)
                    {
                        point.i = i;
                        point.j = j;
                        point.k = k;
                        point.dist = i+j+k;
                        point.ijk = 1;
                        for(int ix = 2; ix <= i; ix++) point.ijk *=ix;
                        for(int ix = 2; ix <= j; ix++) point.ijk *=ix;
                        for(int ix = 2; ix <= k; ix++) point.ijk *=ix;
                        der_list.push_back(point);
                    }
                }
            }
        }
    }

    if( dimension == 2){
        for(int j = 0; j <= Lorder; j+=2){
            for(int i = 0; i <= Lorder; i+=2){
                if(this->ibrav == HEXAGONAL && i !=0 && j !=0) continue;
                if(i+j <= Lorder && i+j >0  )
                {
                    point.i = i;
                    point.j = j;
                    point.k = 0;
                    point.dist = i+j;
                    point.ijk = 1;
                    for(int ix = 2; ix <= i; ix++) point.ijk *=ix;
                    for(int ix = 2; ix <= j; ix++) point.ijk *=ix;
                    der_list.push_back(point);
                }
            }
        }
    }

    if( dimension == 1){
        for(int i = 2; i <= Lorder; i+=2){
            {
                point.i = i;
                point.j = 0;
                point.k = 0;
                point.dist = i;
                point.ijk = 1;
                for(int ix = 2; ix <= i; ix++) point.ijk *=ix;
                der_list.push_back(point);
            }
        }
    }

    std::stable_sort(der_list.begin(), der_list.end(), customLess_ijk);
    std::stable_sort(der_list.begin(), der_list.end(), customLess_dist);

}
