#include <vector>
#include <algorithm>
#include <math.h>
#include <iostream>
#include "LaplacianCoeff.h"
#include "blas.h"

struct GridPoint{
    double dist;
    double dx,dy,dz;
    int i, j, k, eq_num;
    int ijk;
    double weight_factor;
    double coeff;
    int relative_address;
};
typedef GridPoint GridPoint;

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
    GridPoint point;
    int index = 0;

    std::vector<GridPoint> der_list;  // only i,j,k are needed.
    // for laplacian, 2n order is also 2n+1 order since the odd partianl derivatives cancelled automatically.
    Lorder = Lorder /2 * 2;

    // determine how many neghbors are needed for this order.

    std::vector<int> num_neighbor;
    int num_derivative = 0;
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
                    num_derivative++;
                }
            }
        }
    }

    std::stable_sort(der_list.begin(), der_list.end(), customLess_ijk);
    std::stable_sort(der_list.begin(), der_list.end(), customLess_dist);

//    printf("\n bbb %d\n", num_derivative);
//   index = 0;
//    for(auto a:der_list)
//    {
//        std::cout <<index<< "    "<<a.i<<"  " <<a.j<<"  " <<a.k << "  "<<a.ijk<<std::endl;
//        index++;
//    }

    double dx, dy,dz, dist;    
    for(int i = -Lorder/2; i <= Lorder/2; i++){
        for(int j = -Lorder/2; j <= Lorder/2; j++){
            for(int k = -Lorder/2; k <= Lorder/2; k++){
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

                dx = i*a[j][0]/Ngrid[0];
                dy = i*a[j][1]/Ngrid[0];
                dz = i*a[j][2]/Ngrid[0];

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

                dx = -i*a[j][0]/Ngrid[0];
                dy = -i*a[j][1]/Ngrid[0];
                dz = -i*a[j][2]/Ngrid[0];

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

    if(!this->offdiag) 
    {
        if( (this->ibrav == ORTHORHOMBIC_PRIMITIVE || this->ibrav == CUBIC_PRIMITIVE) )
            std::stable_sort(points.begin(), points.end(), customLess_ijk);

        if(this->ibrav == HEXAGONAL)
        {

            std::vector<GridPoint>  points1, points2;

            for(auto point:points)
            {
                if((point.i == 0 || point.j == 0 || point.i == -point.j) && point.k == 0)
                    points1.push_back(point);
                else if(point.k == 0 || (point.i == 0 && point.j == 0) )
                    points1.push_back(point);
                else
                    points2.push_back(point);
            }

            points.clear();
            for(auto point:points1) points.push_back(point);
            for(auto point:points2) points.push_back(point);
        }
    }

    int info;

    dx = points[0].dx;
    dy = points[0].dy;
    dz = points[0].dz;
    points[0].eq_num = 0;
    int eq_num = 0 ;
    //  all sane |dx|, |dy| , |dz| give a same equation 
    //    points.erase(points.begin()+56,points.begin()+80);
    for(int i = 1; i < (int)points.size(); i++)
    {
        if(fabs(fabs(dx) - fabs(points[i].dx)) < 1.0e-3 
                && fabs(fabs(dy) - fabs(points[i].dy)) <1.0e-3 
                && fabs(fabs(dz) - fabs(points[i].dz)) < 1.0e-3)
        {
            points[i].eq_num = eq_num;
        }
        else
        {
            dx = points[i].dx;
            dy = points[i].dy;
            dz = points[i].dz;
            eq_num++;
            points[i].eq_num = eq_num;
        }

    }

    std::vector<int> num_points;
    for(int i = num_derivative; i < (int)points.size()-1; i++)
    {
        if(
                // (points[i].eq_num+1 >= num_derivative)  &&
                //(points[i].ijk != points[i+1].ijk) )
            (std::abs(points[i].dist - points[i+1].dist) > 1.0e-3 ))
            {
                num_points.push_back(i+1);
            }
    }


    // for(auto a:num_points) std::cout << "numpoint" << "  " << a<<std::endl;

    //    index  =0;
    //    for(auto a:points)
    //    {   
    //        std::cout << index <<"  "<<a.dist << "    "<<a.i<<"  " <<a.j<<"  " <<a.k << "  "<<a.eq_num<<std::endl;
    //        index++;
    //    }



    double *Am = new double[num_derivative * num_derivative];
    double *Bm = new double[num_derivative * num_derivative];
    int *ipvt = new int[num_derivative];
    for(int i = 0; i < num_derivative * num_derivative; i++) 
    {
        Am[i] = 0.0;
        Bm[i] = 0.0;
    }

    double delta_r, delta_c;

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

        //for(int i = 0; i < num_derivative; i++)
        //{   printf("\n");
        //    printf("%d %e", i, Am[i*num_derivative + i]);
        // }
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
            tem  = Am[irow * num_derivative +0] *2;
            tem += Am[irow * num_derivative +1] *2;
            tem += Am[irow * num_derivative +2] *2;
            points[ip].coeff += tem * delta_r;
            coeff0 -= tem * delta_r;
        }
    }

    points.resize(point_end);

    //index  =0;
    //for(auto a:points)
    //{   
    //    std::cout << index <<"  "<<a.dist << "    "<<a.i<<"  " <<a.j<<"  " <<a.k << "  "<<a.coeff<<std::endl;
    //    index++;
   // }


    delete []Am;
    delete []Bm;
    delete []ipvt;

    this->coeff_and_index.clear();
    CoeffList coeff_list;
    coeff_list.coeff = coeff0;
    coeff_list.relative_index.push_back(0);
    this->coeff_and_index.push_back(coeff_list);

    coeff_list.relative_index.clear();
    coeff_list.coeff = points[0].coeff;

    index = points[0].i * dim[1] * dim[2] + points[0].j * dim[2] + points[0].k;
    coeff_list.relative_index.push_back(index);

    for(int ip = 1; ip < point_end; ip++)
    {
        if(std::abs(points[ip].coeff) < 1.0e-8) continue;

        if(std::abs(points[ip].coeff - points[ip-1].coeff) < 1.0e-8)
        {
            index = points[ip].i * dim[1] * dim[2] + points[ip].j * dim[2] + points[ip].k;
            coeff_list.relative_index.push_back(index);
        }
        else
        {

            std::sort(coeff_list.relative_index.begin(), coeff_list.relative_index.end());
            this->coeff_and_index.push_back(coeff_list);

            //start a new coeff and its index(s) 
            coeff_list.relative_index.clear();
            coeff_list.coeff = points[ip].coeff;
            index = points[ip].i * (dim[1]+Lorder) * (dim[2]+Lorder) + points[ip].j * (dim[2]+Lorder) + points[ip].k;
            coeff_list.relative_index.push_back(index);
        }
    }

    this->coeff_and_index.push_back(coeff_list);

    coeff0 = 0;
    index = 0;
    for(auto a:this->coeff_and_index)
    {
        coeff0 += a.coeff * a.relative_index.size();
        index += a.relative_index.size();
        // std::cout << a.coeff << "  "<< a.relative_index.size() <<std::endl;
        // for(auto b:a.relative_index) std::cout <<"            "<< b<<std::endl;
    }

    std::cout << "total points= " <<index <<std::endl;
    std::cout << "sum of coeffs = " <<coeff0 <<std::endl;
}

