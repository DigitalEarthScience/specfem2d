/*
!========================================================================
!
!                   S P E C F E M 2 D  Version 7 . 0
!                   --------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
!
! This software is a computer program whose purpose is to solve
! the two-dimensional viscoelastic anisotropic or poroelastic wave equation
! using a spectral-element method (SEM).
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
! The full text of the license is available in file "LICENSE".
!
!========================================================================
*/


__global__ void compute_acoustic_seismogram_kernel(int nrec_local,
                                                   realw* pressure,
                                                   int* d_ibool,
                                                   realw* hxir, realw* hgammar,
                                                   realw* seismograms,
                                                   int* ispec_selected_rec_loc,
                                                   int it,
                                                   int NSTEP){
  int irec_local = blockIdx.x + blockIdx.y*gridDim.x;
  int tx = threadIdx.x;

  int J = (tx/NGLLX);
  int I = (tx-J*NGLLX);

  __shared__ realw sh_dxd[NGLL2_PADDED];

  if (irec_local < nrec_local) {

    int ispec = ispec_selected_rec_loc[irec_local]-1;

    sh_dxd[tx] = 0;

    if (tx < NGLL2) {

      int iglob = d_ibool[tx+NGLL2_PADDED*ispec]-1;

      realw hlagrange = hxir[irec_local + nrec_local*I]*hgammar[irec_local + nrec_local*J];
      sh_dxd[tx] = hlagrange*pressure[iglob];
    }
    __syncthreads();

    for (unsigned int s=1; s<NGLL2_PADDED ; s *= 2) {
      if (tx % (2*s) == 0) sh_dxd[tx] += sh_dxd[tx + s];
      __syncthreads();
    }

    // Signe moins car pression = -potential_dot_dot
   if (tx == 0) {seismograms[irec_local*NSTEP + it ]                   = -sh_dxd[0];}
   if (tx == 1) {seismograms[irec_local*NSTEP + it + nrec_local*NSTEP] = 0;}
  }
}

