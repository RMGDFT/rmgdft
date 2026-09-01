/*  This file is part of cuTranspose.

    cuTranspose is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    cuTranspose is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with cuTranspose.  If not, see <http://www.gnu.org/licenses/>.

    Copyright 2016 Ibai Gurrutxaga, Javier Muguerza, Jose L. Jodra.
*/

#ifndef KERNELS_012_H_
#define KERNELS_012_H_

/********************************************
 * Public function prototypes               *
 ********************************************/

__global__
void dev_copy( float*       out,
               const float* in,
               int           np0,
               int           np1,
               int           elements_per_thread );
#endif /* KERNELS_012_H_ */
/*  This file is part of cuTranspose.

    cuTranspose is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    cuTranspose is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with cuTranspose.  If not, see <http://www.gnu.org/licenses/>.

    Copyright 2016 Ibai Gurrutxaga, Javier Muguerza, Jose L. Jodra.
*/

#ifndef KERNELS_021_H_
#define KERNELS_021_H_

/********************************************
 * Public function prototypes               *
 ********************************************/

__global__
void dev_transpose_021_ept1( float*       out,
                             const float* in,
                             int           np0,
                             int           np1,
                             int           np2 );

__global__
void dev_transpose_021_ept2( float*       out,
                             const float* in,
                             int           np0,
                             int           np1,
                             int           np2 );

__global__
void dev_transpose_021_ept4( float*       out,
                             const float* in,
                             int           np0,
                             int           np1,
                             int           np2 );

__global__
void dev_transpose_021_in_place( float* data,
                                 int     np0,
                                 int     np1 );
#endif /* KERNELS_021_H_ */
/*  This file is part of cuTranspose.

    cuTranspose is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    cuTranspose is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with cuTranspose.  If not, see <http://www.gnu.org/licenses/>.

    Copyright 2016 Ibai Gurrutxaga, Javier Muguerza, Jose L. Jodra.
*/

#ifndef KERNELS_102_H_
#define KERNELS_102_H_

/********************************************
 * Public function prototypes               *
 ********************************************/

__global__
 void dev_transpose_102_ept1( float*       out,
                              const float* in,
                              int           np0,
                              int           np1,
                              int           np2 );

__global__
void dev_transpose_102_ept2( float*       out,
                             const float* in,
                             int           np0,
                             int           np1,
                             int           np2 );

__global__
void dev_transpose_102_ept4( float*       out,
                             const float* in,
                             int           np0,
                             int           np1,
                             int           np2 );

__global__
void dev_transpose_102_in_place( float* data,
                                 int     np0,
                                 int     np2 );
#endif /* KERNELS_102_H_ */
/*  This file is part of cuTranspose.

    cuTranspose is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    cuTranspose is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with cuTranspose.  If not, see <http://www.gnu.org/licenses/>.

    Copyright 2016 Ibai Gurrutxaga, Javier Muguerza, Jose L. Jodra.
*/

#ifndef KERNELS_120_H_
#define KERNELS_120_H_

/********************************************
 * Public function prototypes               *
 ********************************************/

__global__
void dev_transpose_120_ept1( float*       out,
                             const float* in,
                             int           np0,
                             int           np1,
                             int           np2 );

__global__
void dev_transpose_120_ept2( float*       out,
                             const float* in,
                             int           np0,
                             int           np1,
                             int           np2 );

__global__
void dev_transpose_120_ept4( float*       out,
                             const float* in,
                             int           np0,
                             int           np1,
                             int           np2 );

__global__
void dev_transpose_120_in_place( float* data,
                                 int     np0 );
#endif /* KERNELS_120_H_ */
/*  This file is part of cuTranspose.

    cuTranspose is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    cuTranspose is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with cuTranspose.  If not, see <http://www.gnu.org/licenses/>.

    Copyright 2016 Ibai Gurrutxaga, Javier Muguerza, Jose L. Jodra.
*/

#ifndef KERNELS_201_H_
#define KERNELS_201_H_

/********************************************
 * Public function prototypes               *
 ********************************************/

__global__
void dev_transpose_201_ept1( float*       out,
                             const float* in,
                             int           np0,
                             int           np1,
                             int           np2 );

__global__
void dev_transpose_201_ept2( float*       out,
                              const float* in,
                              int           np0,
                              int           np1,
                              int           np2 );

__global__
void dev_transpose_201_ept4( float*       out,
                             const float* in,
                             int           np0,
                             int           np1,
                             int           np2 );

__global__
void dev_transpose_201_in_place( float* data,
                                 int     np0 );
#endif /* KERNELS_201_H_ */
/*  This file is part of cuTranspose.

    cuTranspose is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    cuTranspose is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with cuTranspose.  If not, see <http://www.gnu.org/licenses/>.

    Copyright 2016 Ibai Gurrutxaga, Javier Muguerza, Jose L. Jodra.
*/

#ifndef KERNELS_210_H_
#define KERNELS_210_H_

/********************************************
 * Public function prototypes               *
 ********************************************/

__global__
void dev_transpose_210_ept1( float*       out,
                             const float* in,
                             int           np0,
                             int           np1,
                             int           np2 );

__global__
void dev_transpose_210_ept2( float*       out,
                             const float* in,
                             int           np0,
                             int           np1,
                             int           np2 );

__global__
void dev_transpose_210_ept4( float*       out,
                             const float* in,
                             int           np0,
                             int           np1,
                             int           np2 );

__global__
void dev_transpose_210_in_place( float* data,
                                 int     np0,
                                 int     np1 );

#endif /* KERNELS_210_H_ */
