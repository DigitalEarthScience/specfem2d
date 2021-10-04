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

// this file has been automatically generated by script utils/create_specfem2D_gpu_cuda_method_stubs.pl

#include <stdio.h>
#include <stdlib.h>

#include "config.h"

typedef float realw;



//
// src/gpu/assemble_MPI_scalar_cuda.cu
//

void FC_FUNC_(transfer_boun_pot_from_device,
              TRANSFER_BOUN_POT_FROM_DEVICE)(long* Mesh_pointer,
                                             realw* h_send_potential_dot_dot_buffer,
                                             const int* FORWARD_OR_ADJOINT){}

void FC_FUNC_(transfer_asmbl_pot_to_device,
              TRANSFER_ASMBL_POT_TO_DEVICE)(long* Mesh_pointer,
                                            realw* h_buffer_recv_scalar_gpu,
                                            const int* FORWARD_OR_ADJOINT) {}


//
// src/gpu/assemble_MPI_vector_cuda.cu
//

void FC_FUNC_(transfer_boun_accel_from_device,
              TRANSFER_BOUN_ACCEL_FROM_DEVICE)(long* Mesh_pointer,
                                               realw* h_send_accel_buffer,
                                               const int* FORWARD_OR_ADJOINT){}

void FC_FUNC_(transfer_boundary_from_device_a,
              TRANSFER_BOUNDARY_FROM_DEVICE_A)(long* Mesh_pointer) {}

void FC_FUNC_(prepare_boundary_on_device,
              PREPARE_BOUNDARY_ON_DEVICE)(long* Mesh_pointer) {}

void FC_FUNC_(transfer_boundary_to_device_a,
              TRANSFER_BOUNDARY_TO_DEVICE_A)(long* Mesh_pointer,
                                             realw* buffer_recv_vector_gpu,
                                             const int* max_nibool_interfaces_ext_mesh) {}

void FC_FUNC_(transfer_asmbl_accel_to_device,
              TRANSFER_ASMBL_ACCEL_TO_DEVICE)(long* Mesh_pointer,
                                              realw* buffer_recv_vector_gpu,
                                              const int* max_nibool_interfaces_ext_mesh,
                                              const int* nibool_interfaces_ext_mesh,
                                              const int* ibool_interfaces_ext_mesh,
                                              const int* FORWARD_OR_ADJOINT) {}

void FC_FUNC_(sync_copy_from_device,
              SYNC_copy_FROM_DEVICE)(long* Mesh_pointer,
                                     int* iphase,
                                     realw* send_buffer) {}


//
// src/gpu/check_fields_cuda.cu
//

void FC_FUNC_(output_free_device_memory,
              OUTPUT_FREE_DEVICE_MEMORY)(int* myrank_f) {}

void FC_FUNC_(get_free_device_memory,
              get_FREE_DEVICE_MEMORY)(realw* free, realw* used, realw* total) {}

void FC_FUNC_(get_norm_acoustic_from_device,
              GET_NORM_ACOUSTIC_FROM_DEVICE)(realw* norm,long* Mesh_pointer,const int* FORWARD_OR_ADJOINT) {}

void FC_FUNC_(get_norm_elastic_from_device,
              GET_NORM_ELASTIC_FROM_DEVICE)(realw* norm,long* Mesh_pointer,const int* FORWARD_OR_ADJOINT) {}

void FC_FUNC_(get_max_accel,
              GET_MAX_ACCEL)(int* itf,int* sizef,long* Mesh_pointer) {}

void FC_FUNC_(check_max_norm_displ_gpu,
              CHECK_MAX_NORM_DISPL_GPU)(int* size, realw* displ,long* Mesh_pointer,int* announceID) {}

void FC_FUNC_(check_max_norm_vector,
              CHECK_MAX_NORM_VECTOR)(int* size, realw* vector1, int* announceID) {}

void FC_FUNC_(check_max_norm_displ,
              CHECK_MAX_NORM_DISPL)(int* size, realw* displ, int* announceID) {}

void FC_FUNC_(check_max_norm_b_displ_gpu,
              CHECK_MAX_NORM_B_DISPL_GPU)(int* size, realw* b_displ,long* Mesh_pointer,int* announceID) {}

void FC_FUNC_(check_max_norm_b_accel_gpu,
              CHECK_MAX_NORM_B_ACCEL_GPU)(int* size, realw* b_accel,long* Mesh_pointer,int* announceID) {}

void FC_FUNC_(check_max_norm_b_veloc_gpu,
              CHECK_MAX_NORM_B_VELOC_GPU)(int* size, realw* b_veloc,long* Mesh_pointer,int* announceID) {}

void FC_FUNC_(check_max_norm_b_displ,
              CHECK_MAX_NORM_B_DISPL)(int* size, realw* b_displ,int* announceID) {}

void FC_FUNC_(check_max_norm_b_accel,
              CHECK_MAX_NORM_B_ACCEL)(int* size, realw* b_accel,int* announceID) {}

void FC_FUNC_(check_error_vectors,
              CHECK_ERROR_VECTORS)(int* sizef, realw* vector1,realw* vector2) {}


//
// src/gpu/compute_add_sources_viscoacoustic_cuda.cu
//

void FC_FUNC_(compute_add_sources_ac_cuda,
              COMPUTE_ADD_SOURCES_AC_CUDA)(long* Mesh_pointer,
                                           int* iphasef,
                                           int * itf) {}

void FC_FUNC_(compute_add_sources_ac_s3_cuda,
              COMPUTE_ADD_SOURCES_AC_s3_CUDA)(long* Mesh_pointer,
                                              int* iphasef,
                                              int* itf) {}

void FC_FUNC_(compute_add_moving_sources_ac_cuda,
              COMPUTE_ADD_MOVING_SOURCES_AC_CUDA)(long* Mesh_pointer,
                                                  int* iphase_f,
                                                  int* nsources_local_moving,
                                                  int* itf,
                                                  int* NSTEP_f,
                                                  int* nsources_f) {}

void FC_FUNC_(add_sources_ac_sim_2_or_3_cuda,
              ADD_SOURCES_AC_SIM_2_OR_3_CUDA)(long* Mesh_pointer,
                                               int* iphasef,
                                               int* itf,
                                               int* nadj_rec_local,
                                               int* NSTEP) {}


//
// src/gpu/compute_add_sources_viscoelastic_cuda.cu
//

void FC_FUNC_(compute_add_sources_el_cuda,
              COMPUTE_ADD_SOURCES_EL_CUDA)(long* Mesh_pointer,
                                           int* iphasef,
                                           int* itf) {}

void FC_FUNC_(compute_add_sources_el_s3_cuda,
              COMPUTE_ADD_SOURCES_EL_S3_CUDA)(long* Mesh_pointer,
                                              int* iphasef,
                                              int* itf) {}

void FC_FUNC_(add_sources_el_sim_type_2_or_3,
              ADD_SOURCES_EL_SIM_TYPE_2_OR_3)(long* Mesh_pointer,
                                               int* iphasef,
                                               int* itf,
                                               int* nadj_rec_local,
                                               int* NSTEP) {}


//
// src/gpu/compute_coupling_cuda.cu
//

void FC_FUNC_(compute_coupling_ac_el_cuda,
              COMPUTE_COUPLING_AC_EL_CUDA)(long* Mesh_pointer,
                                           int* iphasef,
                                           int* num_coupling_ac_el_facesf,
                                           int* FORWARD_OR_ADJOINT) {}

void FC_FUNC_(compute_coupling_el_ac_cuda,
              COMPUTE_COUPLING_EL_AC_CUDA)(long* Mesh_pointer,
                                           int* iphasef,
                                           int* num_coupling_ac_el_facesf,
                                           int* FORWARD_OR_ADJOINT) {}


//
// src/gpu/compute_forces_acoustic_cuda.cu
//

void FC_FUNC_(compute_forces_acoustic_cuda,
              COMPUTE_FORCES_ACOUSTIC_CUDA)(long* Mesh_pointer,
                                            int* iphase,
                                            int* nspec_outer_acoustic,
                                            int* nspec_inner_acoustic,
                                            int* ATTENUATION_VISCOACOUSTIC,
                                            int* compute_wavefield_1,
                                            int* compute_wavefield_2) {}


//
// src/gpu/compute_forces_viscoelastic_cuda.cu
//

void FC_FUNC_(compute_forces_viscoelastic_cuda,
              COMPUTE_FORCES_VISCOELASTIC_CUDA)(long* Mesh_pointer,
                                                int* iphase,
                                                realw* deltat,
                                                int* nspec_outer_elastic,
                                                int* nspec_inner_elastic,
                                                int* ANISOTROPY,
                                                int* ATTENUATION_VISCOELASTIC,
                                                int* compute_wavefield_1,
                                                int* compute_wavefield_2) {}


//
// src/gpu/compute_kernels_cuda.cu
//

void FC_FUNC_(compute_kernels_elastic_cuda,
              COMPUTE_KERNELS_ELASTIC_CUDA)(long* Mesh_pointer,realw * deltat) {}

void FC_FUNC_(compute_kernels_acoustic_cuda,
              COMPUTE_KERNELS_ACOUSTIC_CUDA)(long* Mesh_pointer,realw * deltat) {}

void FC_FUNC_(compute_kernels_hess_cuda,
              COMPUTE_KERNELS_HESS_CUDA)(long* Mesh_pointer,
                                         int* ELASTIC_SIMULATION,
                                         int* ACOUSTIC_SIMULATION) {}


//
// src/gpu/compute_stacey_acoustic_cuda.cu
//

void FC_FUNC_(compute_stacey_acoustic_cuda,
              COMPUTE_STACEY_ACOUSTIC_CUDA)(long* Mesh_pointer,
                                            int* iphasef,
                                            realw* h_b_absorb_potential_left,
                                            realw* h_b_absorb_potential_right,
                                            realw* h_b_absorb_potential_top,
                                            realw* h_b_absorb_potential_bottom,
                                            int* compute_wavefield_1,
                                            int* compute_wavefield_2,
                                            int* UNDO_ATTENUATION_AND_OR_PML) {}


//
// src/gpu/compute_stacey_viscoelastic_cuda.cu
//

void FC_FUNC_(compute_stacey_viscoelastic_cuda,
              COMPUTE_STACEY_VISCOELASTIC_CUDA)(long* Mesh_pointer,
                                                int* iphasef,
                                                realw* h_b_absorb_elastic_left,
                                                realw* h_b_absorb_elastic_right,
                                                realw* h_b_absorb_elastic_top,
                                                realw* h_b_absorb_elastic_bottom,
                                                int* compute_wavefield_1,
                                                int* compute_wavefield_2,
                                                int *UNDO_ATTENUATION_AND_OR_PML) {}


//
// src/gpu/enforce_acoustic_free_surface_cuda.cu
//

void FC_FUNC_(acoustic_enforce_free_surf_cuda,
              ACOUSTIC_ENFORCE_FREE_SURF_CUDA)(long* Mesh_pointer,int* compute_wavefield_1,int* compute_wavefield_2) {}


//
// src/gpu/helper_functions.cu
//

void FC_FUNC_(pause_for_debug,PAUSE_FOR_DEBUG)() {}


//
// src/gpu/initialize_cuda.cu
//

void FC_FUNC_(initialize_cuda_device,
              INITIALIZE_CUDA_DEVICE)(int* myrank_f,int* ncuda_devices) {
 fprintf(stderr,"ERROR: GPU_MODE enabled without GPU/CUDA Support. To enable GPU support, reconfigure with --with-cuda flag.\n");
 exit(1);
}

void FC_FUNC_(initialize_cuda_aware_mpi,
              INITIALIZE_CUDA_AWARE_MPI)() {}


//
// src/gpu/pml_compute_cuda.cu
//

void FC_FUNC_(pml_boundary_acoustic_cuda,
              PML_BOUNDARY_ACOUSTIC_CUDA)(long* Mesh_pointer,int* compute_wavefield_1,int* compute_wavefield_2) {}


//
// src/gpu/prepare_mesh_constants_cuda.cu
//

void FC_FUNC_(prepare_constants_device,
              PREPARE_CONSTANTS_DEVICE)(long* Mesh_pointer,
                                        int* h_NGLLX, int* NSPEC_AB, int* NGLOB_AB,
                                        realw* h_xix, realw* h_xiz,
                                        realw* h_gammax, realw* h_gammaz,
                                        realw* h_kappav, realw* h_muv,
                                        int* h_ibool,
                                        int* num_interfaces_ext_mesh, int* max_nibool_interfaces_ext_mesh,
                                        int* h_nibool_interfaces_ext_mesh, int* h_ibool_interfaces_ext_mesh,
                                        realw* h_hprime_xx, realw* h_hprimewgll_xx,
                                        realw* h_wxgll,
                                        int* STACEY_BOUNDARY_CONDITIONS,
                                        int* PML_BOUNDARY_CONDITIONS,
                                        int* h_ispec_is_inner,
                                        int* nsources_local_f,
                                        realw* h_sourcearrays, realw * h_source_time_function,
                                        int* NSTEP,
                                        int* h_ispec_selected_source,
                                        int* h_ispec_selected_rec_loc,
                                        int* nrec_local,
                                        realw * h_cosrot,realw * h_sinrot,
                                        int* SIMULATION_TYPE,
                                        int* P_SV,
                                        int* nspec_acoustic,int* nspec_elastic,
                                        int* ispec_is_acoustic, int* ispec_is_elastic,
                                        int* h_myrank,
                                        int* SAVE_FORWARD,
                                        realw* h_xir_store, realw* h_gammar_store,
                                        int* h_NSIGTYPE, int* h_seismotypeVec,
                                        int* nlength_seismogram) {}

void FC_FUNC_(prepare_fields_acoustic_device,
              PREPARE_FIELDS_ACOUSTIC_DEVICE)(long* Mesh_pointer,
                                              realw* rmass_acoustic, realw* rhostore, realw* kappastore,
                                              int* num_phase_ispec_acoustic, int* phase_ispec_inner_acoustic,
                                              int* num_free_surface_faces,
                                              int* free_surface_ispec,
                                              int* free_surface_ijk,
                                              int* ELASTIC_SIMULATION,
                                              int* num_coupling_ac_el_faces,
                                              int* coupling_ac_el_ispec,
                                              int* coupling_ac_el_ijk,
                                              realw* coupling_ac_el_normal,
                                              realw* coupling_ac_el_jacobian2Dw,
                                              int * h_ninterface_acoustic,int * h_inum_interfaces_acoustic,
                                              int* ATTENUATION_VISCOACOUSTIC,
                                              realw* h_A_newmark,realw* h_B_newmark,
                                              int* NO_BACKWARD_RECONSTRUCTION,realw* h_no_backward_acoustic_buffer) {}

void FC_FUNC_(prepare_fields_acoustic_adj_dev,
              PREPARE_FIELDS_ACOUSTIC_ADJ_DEV)(long* Mesh_pointer,
                                               int* APPROXIMATE_HESS_KL,
                                               int* ATTENUATION_VISCOACOUSTIC,
                                               int* NO_BACKWARD_RECONSTRUCTION) {}

void FC_FUNC_(prepare_fields_elastic_device,
              PREPARE_FIELDS_ELASTIC_DEVICE)(long* Mesh_pointer,
                                             realw* rmassx, realw* rmassz,
                                             int* num_phase_ispec_elastic,
                                             int* phase_ispec_inner_elastic,
                                             int* ispec_is_anisotropic,
                                             int* ANISOTROPY,
                                             realw *c11store,realw *c12store,realw *c13store,
                                             realw *c15store,
                                             realw *c23store,
                                             realw *c25store,realw *c33store,
                                             realw *c35store,
                                             realw *c55store,
                                             int* h_ninterface_elastic,int * h_inum_interfaces_elastic,
                                             int* ATTENUATION_VISCOELASTIC,
                                             realw* h_A_newmark_mu,realw* h_B_newmark_mu,
                                             realw* h_A_newmark_kappa,realw* h_B_newmark_kappa) {}

void FC_FUNC_(prepare_fields_elastic_adj_dev,
              PREPARE_FIELDS_ELASTIC_ADJ_DEV)(long* Mesh_pointer,
                                              int* size_f,
                                              int* APPROXIMATE_HESS_KL,
                                              int* ATTENUATION_VISCOELASTIC,
                                              int* NO_BACKWARD_RECONSTRUCTION){}

void FC_FUNC_(prepare_sim2_or_3_const_device,
              PREPARE_SIM2_OR_3_CONST_DEVICE)(long* Mesh_pointer,
                                              int* nadj_rec_local,
                                              realw* h_source_adjoint,
                                              int* NSTEP) {}

void FC_FUNC_(prepare_pml_device,
              PREPARE_PML_DEVICE)(long* Mesh_pointer,
                                  int* NSPEC_PML,
                                  int* NSPEC_PML_X,
                                  int* NSPEC_PML_Z,
                                  int* NSPEC_PML_XZ,
                                  int* h_spec_to_pml,
                                  realw* h_abs_normalized,
                                  realw* ALPHA_MAX_PML,
                                  realw* d0_max,
                                  realw* deltat,
                                  realw* h_alphax_store,
                                  realw* h_alphaz_store,
                                  realw* h_betax_store,
                                  realw* h_betaz_store,
                                  int *PML_nglob_abs_acoustic_f,
                                  int *h_PML_abs_points_acoustic){}

void FC_FUNC_(prepare_stacey_device,
              PREPARE_STACEY_DEVICE)(long* Mesh_pointer,
                                     int* ACOUSTIC_SIMULATION,
                                     int* ELASTIC_SIMULATION,
                                     realw* rho_vp, realw* rho_vs,
                                     int* h_nspec_bottom,
                                     int* h_nspec_left,
                                     int* h_nspec_right,
                                     int* h_nspec_top,
                                     int* h_abs_boundary_ispec, int* h_abs_boundary_ij,
                                     realw* h_abs_boundary_normal,
                                     realw* h_abs_boundary_jacobian1Dw,
                                     int* h_num_abs_boundary_faces,
                                     int* h_edge_abs,
                                     int* h_ib_bottom,
                                     int* h_ib_left,
                                     int* h_ib_right,
                                     int* h_ib_top){}

void FC_FUNC_(prepare_moving_sources_cuda,
              PREPARE_MOVING_SOURCES_CUDA)(long* Mesh_pointer,
                                           int* h_nsources_local_f_moving,
                                           int* NSOURCES,
                                           realw* h_sourcearrays_moving,
                                           int* h_ispec_selected_source_moving,
                                           int* NSTEP,
                                           realw* h_source_time_function_moving) {}

//void FC_FUNC_(recompute_source_position_cuda,
//              RECOMPUTE_SOURCE_POSITION_CUDA)(long* Mesh_pointer,
//                                        int* nsources_local_f,
//                                        realw* h_sourcearrays,
//                                        int* h_ispec_selected_source) {}

void FC_FUNC_(prepare_cleanup_device,
              PREPARE_CLEANUP_DEVICE)(long* Mesh_pointer,
                                      int* ACOUSTIC_SIMULATION,
                                      int* ELASTIC_SIMULATION,
                                      int* ANISOTROPY,
                                      int* APPROXIMATE_HESS_KL,
                                      int* ATTENUATION_VISCOACOUSTIC,
                                      int* ATTENUATION_VISCOELASTIC,
                                      int* NO_BACKWARD_RECONSTRUCTION,
                                      realw * h_no_backward_acoustic_buffer) {}


//
// src/gpu/smooth_cuda.cu
//

void FC_FUNC_(prepare_arrays_gpu,
              PREPARE_arrays_GPU)(long * Container,
                                  realw * xstore_me,
                                  realw * zstore_me,
                                  realw * sigma_h2_inv,
                                  realw * sigma_v2_inv,
                                  realw * h_criterion,
                                  realw * v_criterion,
                                  int * nspec_me,
                                  int * nker,
                                  realw * wgll_sq){}

void FC_FUNC_(compute_smooth,
              COMPUTE_SMOOTH)(long * smooth_pointer,
                              realw * jacobian,
                              realw * xstore_other,
                              realw * zstore_other,
                              realw * data_other,
                              const int * nspec_other){}

void FC_FUNC_(get_smooth,
              GET_SMOOTH)(long * smooth_pointer,realw * data_smooth){}


//
// src/gpu/transfer_fields_cuda.cu
//

void FC_FUNC_(transfer_fields_el_to_device,
              TRANSFER_FIELDS_EL_TO_DEVICE)(int* size, realw* displ, realw* veloc, realw* accel,long* Mesh_pointer) {}

void FC_FUNC_(transfer_fields_el_from_device,
              TRANSFER_FIELDS_EL_FROM_DEVICE)(int* size, realw* displ, realw* veloc, realw* accel,long* Mesh_pointer) {}

void FC_FUNC_(transfer_b_fields_to_device,
              TRANSFER_B_FIELDS_TO_DEVICE)(int* size, realw* b_displ, realw* b_veloc, realw* b_accel,
                                           long* Mesh_pointer) {}

void FC_FUNC_(transfer_b_fields_from_device,
              TRANSFER_B_FIELDS_FROM_DEVICE)(int* size, realw* b_displ, realw* b_veloc, realw* b_accel,long* Mesh_pointer) {}

void FC_FUNC_(transfer_accel_to_device,
              TRNASFER_ACCEL_TO_DEVICE)(int* size, realw* accel,long* Mesh_pointer) {}

void FC_FUNC_(transfer_accel_from_device,
              TRANSFER_ACCEL_FROM_DEVICE)(int* size, realw* accel,long* Mesh_pointer) {}

void FC_FUNC_(transfer_b_accel_from_device,
              TRNASFER_B_ACCEL_FROM_DEVICE)(int* size, realw* b_accel,long* Mesh_pointer) {}

void FC_FUNC_(transfer_b_displ_from_device,
              TRANSFER_B_DISPL_FROM_DEVICE)(int* size, realw* displ,long* Mesh_pointer) {}

void FC_FUNC_(transfer_displ_from_device,
              TRANSFER_DISPL_FROM_DEVICE)(int* size, realw* displ,long* Mesh_pointer) {}

void FC_FUNC_(transfer_kernels_el_to_host,
              TRANSFER_KERNELS_EL_TO_HOST)(long* Mesh_pointer,
                                            realw* h_rho_kl,
                                            realw* h_mu_kl,
                                            realw* h_kappa_kl,
                                            int* NSPEC_AB) {}

void FC_FUNC_(transfer_fields_ac_to_device,
              TRANSFER_FIELDS_AC_TO_DEVICE)(int* size,
                                            realw* potential_acoustic,
                                            realw* potential_dot_acoustic,
                                            realw* potential_dot_dot_acoustic,
                                            long* Mesh_pointer) {}

void FC_FUNC_(transfer_b_fields_ac_to_device,
              TRANSFER_B_FIELDS_AC_TO_DEVICE)(int* size,
                                              realw* b_potential_acoustic,
                                              realw* b_potential_dot_acoustic,
                                              realw* b_potential_dot_dot_acoustic,
                                              long* Mesh_pointer) {}

void FC_FUNC_(transfer_fields_ac_from_device,
              TRANSFER_FIELDS_AC_FROM_DEVICE)(int* size,
                                              realw* potential_acoustic,
                                              realw* potential_dot_acoustic,
                                              realw* potential_dot_dot_acoustic,
                                              long* Mesh_pointer) {}

void FC_FUNC_(transfer_b_fields_ac_from_device,
              TRANSFER_B_FIELDS_AC_FROM_DEVICE)(int* size,
                                                realw* b_potential_acoustic,
                                                realw* b_potential_dot_acoustic,
                                                realw* b_potential_dot_dot_acoustic,
                                                long* Mesh_pointer) {}

void FC_FUNC_(transfer_b_potential_ac_from_device,
              TRANSFER_B_POTENTIAL_AC_FROM_DEVICE)(int* size,
                                                realw* b_potential_acoustic,
                                                long* Mesh_pointer) {}

void FC_FUNC_(transfer_b_potential_ac_to_device,
              TRANSFER_B_POTENTIAL_AC_TO_DEVICE)(int* size,
                                                 realw* b_potential_acoustic,
                                                 long* Mesh_pointer) {}

void FC_FUNC_(transfer_dot_dot_from_device,
              TRNASFER_DOT_DOT_FROM_DEVICE)(int* size, realw* potential_dot_dot_acoustic,long* Mesh_pointer) {}

void FC_FUNC_(transfer_b_dot_dot_from_device,
              TRNASFER_B_DOT_DOT_FROM_DEVICE)(int* size, realw* b_potential_dot_dot_acoustic,long* Mesh_pointer) {}

void FC_FUNC_(transfer_kernels_ac_to_host,
              TRANSFER_KERNELS_AC_TO_HOST)(long* Mesh_pointer,realw* h_rho_ac_kl,realw* h_kappa_ac_kl,int* NSPEC_AB) {}

void FC_FUNC_(transfer_kernels_hess_el_tohost,
              TRANSFER_KERNELS_HESS_EL_TOHOST)(long* Mesh_pointer,realw* h_hess_kl,int* NSPEC_AB) {}

void FC_FUNC_(transfer_kernels_hess_ac_tohost,
              TRANSFER_KERNELS_HESS_AC_TOHOST)(long* Mesh_pointer,realw* h_hess_ac_kl,int* NSPEC_AB) {}

void FC_FUNC_(transfer_viscoacoustic_b_var_to_device,
              TRANSFER_VISCOACOUSTIC_b_VAR_TO_DEVICE)(int* size,
                                                      realw* b_e1_acous_sf,
                                                      realw* b_sum_forces_old,
                                                      long* Mesh_pointer) {}

void FC_FUNC_(transfer_viscoacoustic_var_from_device,
              TRANSFER_VISCOACOUSTIC_VAR_FROM_DEVICE)(int* size,
                                                      realw* e1_acous_sf,
                                                      realw* sum_forces_old,
                                                      long* Mesh_pointer) {}

void FC_FUNC_(transfer_viscoelastic_b_var_to_device,
              TRANSFER_VISCOELASTIC_b_VAR_TO_DEVICE)(int* size,
                                                     realw* b_e1,
                                                     realw* b_e11,
                                                     realw* b_e13,
                                                     realw* b_dux_dxl_old,
                                                     realw* b_duz_dzl_old,
                                                     realw* b_dux_dzl_plus_duz_dxl_old,
                                                     long* Mesh_pointer) {}

void FC_FUNC_(transfer_viscoelastic_var_from_device,
              TRANSFER_VISCOELASTIC_VAR_FROM_DEVICE)(int* size,
                                                     realw* e1,
                                                     realw* e11,
                                                     realw* e13,
                                                     realw* dux_dxl_old,
                                                     realw* duz_dzl_old,
                                                     realw* dux_dzl_plus_duz_dxl_old,
                                                     long* Mesh_pointer) {}

void FC_FUNC_(transfer_async_pot_ac_from_device,
              TRANSFER_ASYNC_POT_AC_FROM_DEVICE)(realw* pot_buffer,long* Mesh_pointer) {}

void FC_FUNC_(transfer_async_pot_ac_to_device,
              TRANSFER_ASYNC_POT_AC_TO_DEVICE)(realw* pot_buffer,
                                               long* Mesh_pointer) {}

void FC_FUNC_(transfer_compute_kernel_answers_from_device,
              TRANSFER_COMPUTE_KERNEL_ANSWERS_FROM_DEVICE)(long* Mesh_pointer,
                                                           realw* rho_kl,int* size_rho,
                                                           realw* mu_kl, int* size_mu,
                                                           realw* kappa_kl, int* size_kappa) {}

void FC_FUNC_(transfer_compute_kernel_fields_from_device,
              TRANSFER_COMPUTE_KERNEL_FIELDS_FROM_DEVICE)(long* Mesh_pointer,
                                                          realw* accel, int* size_accel,
                                                          realw* b_displ, int* size_b_displ,
                                                          realw* epsilondev_xx,realw* epsilondev_yy,realw* epsilondev_xy,
                                                          realw* epsilondev_xz,realw* epsilondev_yz,
                                                          int* size_epsilondev,
                                                          realw* b_epsilondev_xx,realw* b_epsilondev_yy,realw* b_epsilondev_xy,
                                                          realw* b_epsilondev_xz,realw* b_epsilondev_yz,
                                                          int* size_b_epsilondev,
                                                          realw* rho_kl,int* size_rho,
                                                          realw* mu_kl, int* size_mu,
                                                          realw* kappa_kl, int* size_kappa,
                                                          realw* epsilon_trace_over_3,
                                                          realw* b_epsilon_trace_over_3,
                                                          int* size_epsilon_trace_over_3) {}


//
// src/gpu/update_displacement_cuda.cu
//

void FC_FUNC_(update_displacement_cuda,
              UPDATE_DISPLACMENT_CUDA)(long* Mesh_pointer,
                                          realw* deltat_F,
                                          realw* deltatsqover2_F,
                                          realw* deltatover2_F,
                                          realw* b_deltat_F,
                                          realw* b_deltatsqover2_F,
                                          realw* b_deltatover2_F) {}

void FC_FUNC_(update_displacement_ac_cuda,
              UPDATE_DISPLACEMENT_AC_CUDA)(long* Mesh_pointer,
                                           realw* deltat_F,
                                           realw* deltatsqover2_F,
                                           realw* deltatover2_F,
                                           realw* b_deltat_F,
                                           realw* b_deltatsqover2_F,
                                           realw* b_deltatover2_F,
                                           int* compute_b_wavefield,
                                           int* UNDO_ATTENUATION_AND_OR_PML) {}

void FC_FUNC_(kernel_3_a_cuda,
              KERNEL_3_A_CUDA)(long* Mesh_pointer,
                               realw* deltatover2_F,
                               realw* b_deltatover2_F,
                               int* compute_wavefield_1,
                               int* compute_wavefield_2) {}

void FC_FUNC_(kernel_3_b_cuda,
              KERNEL_3_B_CUDA)(long* Mesh_pointer,
                               realw* deltatover2_F,
                               realw* b_deltatover2_F) {}

void FC_FUNC_(kernel_3_acoustic_cuda,
              KERNEL_3_ACOUSTIC_CUDA)(long* Mesh_pointer,
                                      realw* deltatover2,
                                      realw* b_deltatover2,
                                      int* compute_wavefield_1,
                                      int* compute_wavefield_2) {}


//
// src/gpu/write_seismograms_cuda.cu
//

void FC_FUNC_(compute_seismograms_cuda,
              COMPUTE_SEISMOGRAMS_CUDA)(long* Mesh_pointer_f,
                                        int* i_sigf,
                                        double* sisux, double* sisuz,
                                        int* seismo_currentf,
                                        int* nlength_seismogramf,
                                        int* ELASTIC_SIMULATION,
                                        int* ACOUSTIC_SIMULATION,
                                        int* USE_TRICK_FOR_BETTER_PRESSURE,
                                        int* ATTENUATION_VISCOELASTIC,
                                        int* itf,
                                        int* it_endf) {}

