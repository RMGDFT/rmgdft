/*
 *
 * Copyright (c) 2014, Emil Briggs
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.

 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
*/

#ifndef RMG_Lattice_H
#define RMG_Lattice_H 1

#include "rmg_error.h"


/* Crystal lattice types */
#define No_Lattice                      0
#define CUBIC_PRIMITIVE                 1
#define CUBIC_FC                        2
#define CUBIC_BC                        3
#define HEXAGONAL                       4
#define TRIGONAL_PRIMITIVE              5
#define TETRAGONAL_PRIMITIVE            6
#define TETRAGONAL_BC                   7
#define ORTHORHOMBIC_PRIMITIVE          8
#define ORTHORHOMBIC_BASE_CENTRED       9
#define ORTHORHOMBIC_BC                 10
#define ORTHORHOMBIC_FC                 11
#define MONOCLINIC_PRIMITIVE            12
#define MONOCLINIC_BASE_CENTRED         13
#define TRICLINIC_PRIMITIVE             14

// Declared outside class since we may want to use it in other places
void TransformVector(double *vec, double *a0, double *a1, double *a2,
                                  double *b0, double *b1, double *b2);

// HEXAGONAL is for 60 degree case HEXAGONAL2 is for 120
#define HEXAGONAL2                      15
class Lattice {

private:

    // Grid bravais lattice type 
    int ibrav;

    // lengths of the sides of the supercell
    double xside;
    double yside;
    double zside;

public:

    // lattice vectors
    double a0[3];
    double a1[3];
    double a2[3];
    double at[9];  // For passing to Fortran routines

    // original input lattice vectors
    double a0i[3];
    double a1i[3];
    double a2i[3];

    // reciprocal lattice vectors
    double b0[3];
    double b1[3];
    double b2[3];
    double bg[9];  // For passing to Fortran routines

    // cell force dE/da = -omega sum (a^-1) sigma
    // sigma: the stress tensor
    double stress_tensor[9];
    double cell_force[9]; 
    double cell_velocity[9]; 

    // cell dimensions
    double celldm[6];
    double a, b, c, cosab, cosac, cosbc;

    // Total cell volume
    double omega;

    void latgen (double * celldm, double * OMEGAI, double *a0, double *a1, double *a2, bool flag);

    void move_cell(double dt, int *cell_movable);
    void cross_product (double * a, double * b, double * c);
    double dot_product (double * a, double * b);
    void to_crystal (double *crystal, double *cartesian);
    void to_crystal_vector (double *crystal, double *cartesian);
    void to_crystal_half (double *crystal, double *cartesian);
    void to_cartesian (double *crystal, double *cartesian);
    void to_cartesian_input (double *crystal, double *cartesian);
    void recips (void);
    int get_ibrav_type(void);
    void set_ibrav_type(int newtype);
    double metric (double * crystal);
    double get_omega(void);
    double get_celldm(int which);
    double get_a0(int which);
    double get_a1(int which);
    double get_a2(int which);
    double get_b0(int which);
    double get_b1(int which);
    double get_b2(int which);

    double get_xside(void);
    double get_yside(void);
    double get_zside(void);
    void lat2abc(double *a0, double *a1, double *a2);
    void lat2celldm (int ibrav, double alat, double *a1, double *a2, double *a3);
    void abc2celldm(void);
    int lat2ibrav (double *a0, double *a1, double *a2);
    void rotate_vectors(double *a0, double *a1, double *a2);
    void save_vectors(double *a0, double *a1, double *a2);
    void remake_cell(int ibrav, double alat, double *a1, double *a2, double *a3, double *nalat);
};

#endif
