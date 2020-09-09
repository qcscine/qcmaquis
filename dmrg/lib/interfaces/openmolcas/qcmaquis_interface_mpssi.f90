!!  dmrg-interface-utils: interface to the Maquis DMRG program for various
!!                        quantum-chemistry program packages.
!!  Copyright 2013-2020 Leon Freitag, Erik Hedegaard, Sebastian Keller,
!!                      Stefan Knecht, Yingjin Ma, Christopher Stein
!!                      and Markus Reiher
!!                      Laboratory for Physical Chemistry, ETH Zurich
!!  dmrg-interface-utils is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU Lesser General Public License as published by
!!  the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  dmrg-interface-utils is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!!  GNU Lesser General Public License for more details.
!!
!!  You should have received a copy of the GNU Lesser General Public License
!!  along with dmrg-interface-utils. If not, see <http://www.gnu.org/licenses/>.

module qcmaquis_interface_mpssi
! Module for Fortran-C interoperability with the new (non-Python) QCMaquis interface
! MPSSI part

  use iso_c_binding
  use qcmaquis_interface_cfg
  implicit none
#ifdef _MOLCAS_MPP_
#include "mafdecls.fh"
#endif

  contains
      ! Initialise the MPSSI interface
      ! project_names: array of project names for all jobiphs
      ! states: array containing indices of states to be considered for each jobiph
      ! nstates: number of states considered in each jobiph
      ! nprojects: number of jobiphs
      subroutine qcmaquis_mpssi_init(project_names, states, nstates, n_projects)

      interface
        subroutine qcmaquis_mpssi_init_c(project_names, states, nstates, n_projects) bind(C, name='qcmaquis_mpssi_init')
          import c_int, c_double, c_char, c_ptr, c_loc
          integer(c_int), value :: n_projects
          integer(c_int), dimension(*) :: states, nstates
          type(c_ptr),dimension(*) :: project_names
        end subroutine qcmaquis_mpssi_init_c
      end interface

        integer :: n_projects
        integer, dimension(*) :: states
        integer, dimension(n_projects) :: nstates
        character(len=*),intent(in), dimension(*) :: project_names

        !
        integer :: i, sum_states

        integer(c_int), dimension(:), allocatable :: states_
        integer(c_int), dimension(n_projects) :: nstates_


        character(len=255),dimension(n_projects), target :: project_names_str
        type(c_ptr),dimension(n_projects) :: project_names_ptr



        ! same for the states and nstates arrays
        sum_states = sum(nstates)
        allocate(states_(sum_states))

        do i = 1, sum_states
          states_(i) = states(i)-1 ! States in QCMaquis start indexing from 0 but in MOLCAS from 1
        end do


        do i = 1, n_projects
          nstates_(i) = int(nstates(i), c_int)
        end do

        ! convert suffixes from a fortran array of strings to C
        do i=1, n_projects
          project_names_str(i) = trim(project_names(i))//c_null_char
          project_names_ptr(i) = c_loc(project_names_str(i))
        end do

        call qcmaquis_mpssi_init_c(project_names_ptr, &
                                   states_, &
                                   nstates_,&
                                   int(n_projects, c_int))

        if(allocated(states_)) deallocate(states_)

      end subroutine qcmaquis_mpssi_init

      ! Overlap function
      ! if su2u1 flag is set, SU2U1 checkpoints will be used whenever possible
      ! i.e. for the same multiplicity

      ! Also convert state numbers from 1-based in OpenMOLCAS to 0-based in QCMaquis
      function qcmaquis_mpssi_overlap(bra_name, bra_state, &
                              ket_name, ket_state, su2u1) result(res)

        interface
          function qcmaquis_mpssi_overlap_c(bra_name, bra_state, &
                              ket_name, ket_state, su2u1)          &
                              bind(C, name='qcmaquis_mpssi_overlap') result(res)
            import c_int, c_bool, c_double, c_char
            real(c_double) :: res
            integer(c_int), value, intent(in) :: bra_state, ket_state
            character(len=1,kind=c_char), dimension(*), intent(in) :: bra_name, ket_name
            logical(c_bool), value :: su2u1
          end function qcmaquis_mpssi_overlap_c
        end interface

        integer, intent(in) :: bra_state, ket_state
        logical, intent(in) :: su2u1
        character(len=*), intent(in) :: bra_name, ket_name
        real(c_double) :: res

        res = qcmaquis_mpssi_overlap_c(trim(bra_name)//c_null_char,  &
                                        int(bra_state-1, c_int), &
                                        trim(ket_name)//c_null_char,  &
                                        int(ket_state-1, c_int), &
                                        logical(su2u1, c_bool))
      end function qcmaquis_mpssi_overlap

      ! Get aa and bb components of 1-TDM, from which we can construct 1-TDM and 1-spin TDM
      ! bra_name and ket_name are prefixes for bra and ket, respectively
      ! bra_state and ket_state are indexes, starting from 1 (will be converted to indexes starting from 0 here!)
      ! size: L*L with L= number of active orbitals
      ! output: indices -- array of integers of size 2*L*L = 2*size, with indices
      ! values: array of reals of size 2*L*L = size with TDM elements
      ! TDM should be reconstructed from indices and values as follows:
      ! do k=1,size,2
      !   i=indices(k)
      !   j=indices(k+1)
      !   TDMaa(i+1,j+1) = values(k)
      !   TDMbb(i+1,j+1) = values(k+1)
      ! end do
      ! TODO: check if above is true!!!
      ! if bra_name == ket_name and bra_state == ket_state then only the triangular elements are returned
      ! and you are supposed to take care of it properly!
      subroutine qcmaquis_mpssi_get_onetdm_spin(bra_name, bra_state, ket_name, ket_state, tdmaa, tdmbb, size)
        interface
          subroutine qcmaquis_mpssi_get_onetdm_spin_C(bra_name, bra_state, ket_name, ket_state, tdmaa, tdmbb, size) &
            bind(C, name='qcmaquis_mpssi_get_onetdm_spin')
            import c_int, c_double, c_char
            character(len=1,kind=c_char), dimension(*), intent(in) :: bra_name, ket_name
            integer(c_int), value, intent(in) :: bra_state, ket_state
            real(c_double), dimension(*), intent(out) :: tdmaa, tdmbb
            integer(c_int), value :: size
          end subroutine qcmaquis_mpssi_get_onetdm_spin_C

        end interface

        character(len=*), intent(in) :: bra_name, ket_name
        integer, intent(in) :: bra_state, ket_state
        integer, intent(in) :: size
        real*8, dimension(*), intent(out) ::  tdmaa, tdmbb

        call qcmaquis_mpssi_get_onetdm_spin_C(trim(bra_name)//c_null_char,  &
                                        int(bra_state-1, c_int), &
                                        trim(ket_name)//c_null_char,  &
                                        int(ket_state-1, c_int), &
                                        tdmaa, tdmbb, int(size,c_int))

      end subroutine qcmaquis_mpssi_get_onetdm_spin

      !! Perform SU2U1->2U1 QCMaquis checkpoint transformation
      !! pname: project prefix
      !! state: state index
      subroutine qcmaquis_mpssi_transform(pname, state, Ms)
        interface
          subroutine qcmaquis_mpssi_transform_C(pname, state, Ms) bind(C, name='qcmaquis_mpssi_transform')
            import c_int, c_char
            character(len=1,kind=c_char), dimension(*), intent(in) :: pname
            integer(c_int), value :: state
            integer(c_int), value :: Ms
          end subroutine
        end interface

        character(len=*), intent(in) :: pname
        integer, intent(in) :: state

        integer, intent(in), optional :: Ms
        integer(c_int) :: Ms_
#ifdef _MOLCAS_MPP_
    if(dmrg_host_program_settings%myrank == 0)then
#endif
        if (present(Ms)) then
          Ms_ = Ms
        else
          Ms_ = 0
        endif

        call qcmaquis_mpssi_transform_C(trim(pname)//c_null_char, int(state-1, c_int), Ms_)
#ifdef _MOLCAS_MPP_
    end if
#endif
      end subroutine

      !! QCMaquis MPSSI counterrotation
      !! pname: project prefix
      !! state: state index
      !! t: rotation matrix (from eqs 43-45 in the paper), in square form (it's non-symmetric)
      !! t_size: size of t
      subroutine qcmaquis_mpssi_rotate(pname, state, t, t_size, scale_inactive, Ms)

        interface
          subroutine qcmaquis_mpssi_rotate_C(pname, state, t, t_size, scale_inactive, Ms) bind(C,name='qcmaquis_mpssi_rotate')
            import c_char, c_int, c_double
            character(len=1,kind=c_char), dimension(*), intent(in) :: pname
            integer(c_int), value :: state
            real(c_double), dimension(*), intent(in) :: t
            integer(c_int), value :: t_size
            real(c_double), value :: scale_inactive
            integer(c_int), value :: Ms

          end subroutine qcmaquis_mpssi_rotate_C
        end interface

      character(len=*), intent(in) :: pname
      integer, intent(in) :: state, t_size, Ms
      real*8, dimension(*) :: t
      real*8 :: scale_inactive

#ifdef _MOLCAS_MPP_
    if(dmrg_host_program_settings%myrank == 0)then
#endif
      call qcmaquis_mpssi_rotate_C(trim(pname)//c_null_char, int(state-1, c_int), t, int(t_size, c_int), &
                                   dble(scale_inactive), int(Ms, c_int))
#ifdef _MOLCAS_MPP_
    end if
#endif
      end subroutine qcmaquis_mpssi_rotate
  end module qcmaquis_interface_mpssi
