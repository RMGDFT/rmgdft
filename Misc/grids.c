/*

  Data and functions related to grid layout and dimensions.

*/

#include "common_prototypes.h"

/* Grid sizes on each PE */
static int PX0_GRID;
static int PY0_GRID;
static int PZ0_GRID;

/* Grid offsets on each PE */
static int PX_OFFSET;
static int PY_OFFSET;
static int PZ_OFFSET;

/* Basis size on each PE */
static int P0_BASIS;

/* Fine grid sizes on each PE */
static int FPX0_GRID;
static int FPY0_GRID;
static int FPZ0_GRID;

/* Fine Grid offsets on each PE */
static int FPX_OFFSET;
static int FPY_OFFSET;
static int FPZ_OFFSET;

/* Fine grid basis size on each PE */
static int FP0_BASIS;

/* Grid bravais lattice type */
static int ibrav;

static int neighbor_first=0;
static int neighbors[6];

void set_grids(void)
{

    int rem, ii, jj, kk;

    // Compute grid sizes for each node.
    PX0_GRID = NX_GRID / PE_X;
    rem = NX_GRID % PE_X;
    if(rem && (ii < rem)) PX0_GRID++;

    PY0_GRID = NY_GRID / PE_Y;
    rem = NY_GRID % PE_Y;
    if(rem && (jj < rem)) PY0_GRID++;

    PZ0_GRID = NZ_GRID / PE_Z;
    rem = NZ_GRID % PE_Z;
    if(rem && (kk < rem)) PZ0_GRID++;

    find_node_sizes(pct.gridpe, NX_GRID, NY_GRID, NZ_GRID, &PX0_GRID, &PY0_GRID, &PZ0_GRID);
    find_node_sizes(pct.gridpe, FNX_GRID, FNY_GRID, FNZ_GRID, &FPX0_GRID, &FPY0_GRID, &FPZ0_GRID);

    P0_BASIS = PX0_GRID * PY0_GRID * PZ0_GRID;
    FP0_BASIS = FPX0_GRID * FPY0_GRID * FPZ0_GRID;

    // Now compute the global grid offset of the first point of the coarse and fine node grids
    find_node_offsets(pct.gridpe, NX_GRID, NY_GRID, NZ_GRID,
                      &PX_OFFSET, &PY_OFFSET, &PZ_OFFSET);

    find_node_offsets(pct.gridpe, FNX_GRID, FNY_GRID, FNZ_GRID,
                      &FPX_OFFSET, &FPY_OFFSET, &FPZ_OFFSET);

}
int get_PX0_GRID(void)
{
    return PX0_GRID;
}
int get_PY0_GRID(void)
{
    return PY0_GRID;
}
int get_PZ0_GRID(void)
{
    return PZ0_GRID;
}

int get_FPX_OFFSET(void)
{
    return FPX_OFFSET;
}
int get_FPY_OFFSET(void)
{
    return FPY_OFFSET;
}
int get_FPZ_OFFSET(void)
{
    return FPZ_OFFSET;
}
int get_P0_BASIS(void)
{
    return P0_BASIS;
}
int get_FP0_BASIS(void)
{
    return FP0_BASIS;
}
int get_FPX0_GRID(void)
{
    return FPX0_GRID;
}
int get_FPY0_GRID(void)
{
    return FPY0_GRID;
}
int get_FPZ0_GRID(void)
{
    return FPZ0_GRID;
}
int get_ibrav_type(void)
{
    return ibrav;
}
void set_ibrav_type(int value)
{
    ibrav = value;
}

void set_neighbors(int *list)
{
    int idx;
    
    for(idx = 0;idx < 6;idx++) 
        neighbors[idx] = list[idx];

    neighbor_first = 1;

}

// Returns a pointer to the neighbors structure which contains the rank
// of neighboring processors in three-dimensional space.
int *get_neighbors(void)
{
    if(!neighbor_first)
        error_handler("Neighbor list not initialized. Please call set_neighbors first");

    return neighbors;
}
