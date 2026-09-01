/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef RMG_pack_H
#define RMG_pack_H 1

// loop counters for doing a pack/unpack

struct pack_plan_3d {
  int nfast;                 // # of elements in fast index
  int nmid;                  // # of elements in mid index
  int nslow;                 // # of elements in slow index
  int nstride_line;          // stride between successive mid indices
  int nstride_plane;         // stride between successive slow indices
  int nqty;                  // # of values/element
};

template <typename PACK_DATA> void pack_3d(PACK_DATA *data, PACK_DATA *buf, struct pack_plan_3d *plan);
template <typename PACK_DATA> void unpack_3d(PACK_DATA *buf, PACK_DATA *data, struct pack_plan_3d *plan);
template <typename PACK_DATA> void unpack_3d_permute1_1(PACK_DATA *buf, PACK_DATA *data, struct pack_plan_3d *plan);
template <typename PACK_DATA> void unpack_3d_permute1_2(PACK_DATA *buf, PACK_DATA *data, struct pack_plan_3d *plan);
template <typename PACK_DATA> void unpack_3d_permute1_n(PACK_DATA *buf, PACK_DATA *data, struct pack_plan_3d *plan);
template <typename PACK_DATA> void unpack_3d_permute2_1(PACK_DATA *buf, PACK_DATA *data, struct pack_plan_3d *plan);
template <typename PACK_DATA> void unpack_3d_permute2_2(PACK_DATA *buf, PACK_DATA *data, struct pack_plan_3d *plan);
template <typename PACK_DATA> void unpack_3d_permute2_n(PACK_DATA *buf, PACK_DATA *data, struct pack_plan_3d *plan);
template <typename PACK_DATA> void pack_3d(PACK_DATA *data, PACK_DATA *buf, struct pack_plan_3d *plan);
template <typename PACK_DATA> void unpack_3d(PACK_DATA *buf, PACK_DATA *data, struct pack_plan_3d *plan);
template <typename PACK_DATA> void unpack_3d_permute1_1(PACK_DATA *buf, PACK_DATA *data, struct pack_plan_3d *plan);
template <typename PACK_DATA> void unpack_3d_permute1_2(PACK_DATA *buf, PACK_DATA *data, struct pack_plan_3d *plan);
template <typename PACK_DATA> void unpack_3d_permute1_n(PACK_DATA *buf, PACK_DATA *data, struct pack_plan_3d *plan);
template <typename PACK_DATA> void unpack_3d_permute2_1(PACK_DATA *buf, PACK_DATA *data, struct pack_plan_3d *plan);
template <typename PACK_DATA> void unpack_3d_permute2_2(PACK_DATA *buf, PACK_DATA *data, struct pack_plan_3d *plan);
template <typename PACK_DATA> void unpack_3d_permute2_n(PACK_DATA *buf, PACK_DATA *data, struct pack_plan_3d *plan);
template <typename PACK_DATA> void pack_3d(PACK_DATA *data, PACK_DATA *buf, struct pack_plan_3d *plan);
template <typename PACK_DATA> void unpack_3d(PACK_DATA *buf, PACK_DATA *data, struct pack_plan_3d *plan);
template <typename PACK_DATA> void unpack_3d_permute1_1(PACK_DATA *buf, PACK_DATA *data, struct pack_plan_3d *plan);
template <typename PACK_DATA> void unpack_3d_permute1_2(PACK_DATA *buf, PACK_DATA *data, struct pack_plan_3d *plan);
template <typename PACK_DATA> void unpack_3d_permute1_n(PACK_DATA *buf, PACK_DATA *data, struct pack_plan_3d *plan);
template <typename PACK_DATA> void unpack_3d_permute2_1(PACK_DATA *buf, PACK_DATA *data, struct pack_plan_3d *plan);
template <typename PACK_DATA> void unpack_3d_permute2_2(PACK_DATA *buf, PACK_DATA *data, struct pack_plan_3d *plan);
template <typename PACK_DATA> void unpack_3d_permute2_n(PACK_DATA *buf, PACK_DATA *data, struct pack_plan_3d *plan);


#if !defined(PACK_POINTER) && !defined(PACK_MEMCPY)
#define PACK_ARRAY
#endif


#endif
