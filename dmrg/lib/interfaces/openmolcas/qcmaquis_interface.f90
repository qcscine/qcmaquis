!!  dmrg-interface-utils: interface to the Maquis DMRG program for various
!!                        quantum-chemistry program packages.
!!  Copyright 2019-2020 Leon Freitag
!!
!!  Contains parts of the old interface:
!!        (C) 2013-2019 Leon Freitag, Erik Hedegaard, Sebastian Keller,
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

module qcmaquis_interface
! Module for Fortran-C interoperability with the new (non-Python) QCMaquis interface

  use iso_c_binding

  use qcmaquis_interface_cfg
  use qcmaquis_interface_utility_routines, only: str, fiedlerorder_length, file_name_generator
  implicit none
#ifdef _MOLCAS_MPP_
#include "mafdecls.fh"
#endif

  ! enum for higher RDM QCMaquis templates for NEVPT2
  enum, bind(C)
    enumerator :: TEMPLATE_4RDM, TEMPLATE_TRANSITION_3RDM
  end enum

  interface qcmaquis_interface_get_1rdm
    module procedure qcmaquis_interface_get_1rdm_compat, qcmaquis_interface_get_1rdm_full
  end interface

  interface qcmaquis_interface_get_2rdm
    module procedure qcmaquis_interface_get_2rdm_compat, qcmaquis_interface_get_2rdm_full
  end interface
  interface

    subroutine qcmaquis_interface_update_integrals_C(integral_indices, integral_values, size) &
    bind(C,name='qcmaquis_interface_update_integrals')
        import c_int, c_double
        integer(c_int), dimension(*) :: integral_indices
        real(c_double), dimension(*) :: integral_values
        integer(c_int), value :: size
    end subroutine

    subroutine qcmaquis_interface_set_state(state) bind(C)
      import c_int
      integer(c_int), value :: state
    end subroutine

    subroutine qcmaquis_interface_set_nsweeps(nsweeps) bind(C)
      import c_int
      integer(c_int), value :: nsweeps
    end subroutine

    subroutine qcmaquis_interface_reset() bind(C)
    end subroutine

    subroutine qcmaquis_interface_optimize() bind(C)
    end subroutine

    function qcmaquis_interface_get_energy() result(res) bind(C)
      import c_double
      real(c_double) :: res
    end function

    subroutine qcmaquis_interface_get_1rdm_C(indices, values, size) bind(C, name='qcmaquis_interface_get_1rdm')
      import c_int, c_double
      integer(c_int), dimension(*) :: indices
      real(c_double), dimension(*) :: values
      integer(c_int), value :: size
    end subroutine

    subroutine qcmaquis_interface_get_spdm_C(indices, values, size) bind(C, name='qcmaquis_interface_get_spdm')
      import c_int, c_double
      integer(c_int), dimension(*) :: indices
      real(c_double), dimension(*) :: values
      integer(c_int), value :: size
    end subroutine

    subroutine qcmaquis_interface_get_2rdm_C(indices, values, size) bind(C,  name='qcmaquis_interface_get_2rdm')
      import c_int, c_double
      integer(c_int), dimension(*) :: indices
      real(c_double), dimension(*) :: values
      integer(c_int), value :: size
    end subroutine

  end interface
  contains

    ! Initialises QCMaquis interface with values from OpenMOLCAS
    ! parameter description see below
    subroutine qcmaquis_interface_init(nel, L, spin, irrep, nsym, L_per_sym, conv_thresh, &
                                        m, nsweeps, project_name,   &
                                        meas_2rdm, &
                                        ! Parameters from old initialize_dmrg:
                                        nroots,    &        ! Number of roots             DEFAULT:       1
                                        lroot,    &        ! Max   roots                 DEFAULT:       1
                                        iroot,     &        ! Target root                 DEFALUT:       1
                                        thre,      &        ! convergence threshold
                                        weight,    &        ! Root Weights
#ifdef _MOLCAS_MPP_
                                        nprocs,    &        ! number of MPI processes
                                        myrank,    &        ! MPI rank
#endif
                                        initial_occ,&       ! Initial HF occupation (optional)
                                        sweep_m    &
                                        )
      interface
        subroutine qcmaquis_interface_init_c(nel, L, spin, irrep, site_types, conv_thresh, &
                                        m, nsweeps, sweep_m, nsweepm, project_name,   &
                                        meas_2rdm) &
                                        bind(C, name='qcmaquis_interface_preinit')
            import c_char, c_int, c_double, c_bool, c_ptr
            integer(c_int), value :: nel, L, spin, irrep, m, nsweeps, nsweepm
            integer(c_int), dimension(*) :: site_types
            type(c_ptr), value :: sweep_m
            character(len=1,kind=c_char), dimension(*), intent(in) :: project_name
            real(c_double), value :: conv_thresh
            logical(c_bool), value :: meas_2rdm

        end subroutine
      end interface

        integer :: nel, & ! number of active electrons
                     L, & ! number of active orbitals (in all symmetries)
                  spin, & ! spin
                 irrep, & ! spatial symmetry
                  nsym, & ! total number of irreps
                     m, & ! bond dimension (if sweep_m is not used)
                 nsweeps  ! number of sweeps
        integer(c_int), dimension(:), optional, intent(in), target :: sweep_m ! list of m if several m values are desired

        ! Parameters from old initialize_dmrg
        integer, intent(in) :: nroots, lroot
        integer , intent(in), dimension(lroot) :: iroot
        integer, optional, intent(in) :: initial_occ(L,nroots)
        real*8 , intent(in)                   :: thre
        real*8 , intent(in), dimension(lroot):: weight

        character(len=*), intent(in) :: project_name
                 real*8 :: conv_thresh ! convergence threshold
        logical, intent(in), value :: meas_2rdm ! should we measure 2-RDM?

        integer, dimension(nsym) :: L_per_sym ! number of orbitals in each symmetry
#ifdef _MOLCAS_MPP_
        integer, intent(in) :: nprocs
        integer, intent(in) :: myrank
#endif


        type(c_ptr) :: sweep_m_
        integer(c_int) :: nsweepm = 0

        integer :: err

        integer(c_int), dimension(L) :: site_types

        integer :: i, j, n

        integer :: ierr
        ! Copy-paste from old initialize_dmrg

            !> Threshold for energy
        E_threshold = thre



        dmrg_symmetry%nirrep = nsym

        allocate(dmrg_state%iroot(lroot), stat=ierr); if( ierr /= 0 )stop ' Error in allocation: iroot(:)'
        allocate(dmrg_state%weight(lroot), stat=ierr); if( ierr /= 0 )stop ' Error in allocation: weight(:)'

        dmrg_state    = type_state             (                                &
                                                irrep,                    &
                                                nel,                  &
                                                spin,                 &
                                                nroots,                  &
                                                lroot,                 &
                                                iroot(1:lroot), &
                                                weight(1:lroot) &
                                              )

        allocate(dmrg_orbital_space%initial_occ(L,lroot), stat=ierr); if( ierr /= 0 ) &
        stop ' Error in allocation: initial_occ(:,:)'
        dmrg_orbital_space = type_orbital_space(L_per_sym, initial_occ)

        allocate(dmrg_energy%dmrg_state_specific(lroot), stat=ierr); if( ierr /= 0 ) &
        stop ' Error in allocation: dmrg_state_specific(:)'
        dmrg_energy%dmrg_state_specific = 0.0d0

        allocate(dmrg_energy%num_sweeps(lroot), stat=ierr); if( ierr /= 0 ) &
        stop ' Error in allocation: num_sweeps(:)'
        dmrg_energy%num_sweeps     = 0

        allocate(dmrg_energy%max_truncW(lroot), stat=ierr); if( ierr /= 0 ) &
        stop ' Error in allocation: max_truncW(:)'
        dmrg_energy%max_truncW     = 0

        allocate(dmrg_file%qcmaquis_checkpoint_file(lroot), stat=ierr); if( ierr /= 0 ) &
        stop ' Error in allocation: qcmaquis_checkpoint_file(:)'
        dmrg_file%qcmaquis_checkpoint_file = ''


        !> initialize parallel settings from the host program
#ifdef _MOLCAS_MPP_
        dmrg_host_program_settings%nprocs = nprocs
        dmrg_host_program_settings%myrank = myrank
        if(dmrg_host_program_settings%nprocs > 1) dmrg_host_program_settings%runs_parallel = .true.
#endif

        ! process sweep_bond_dimensions
        if(present(sweep_m)) then
            nsweepm = size(sweep_m)
            sweep_m_ = c_loc(sweep_m(1))
        else
            sweep_m_ = c_null_ptr
        end if

        ! process site_types
!         if(nsym.eq.1) then
            ! site_types can be nullptr if we have no symmetry
!             site_types = c_null_ptr
!         else

        site_types = 0
        n = 1
        do i=1,nsym
          do j=1,L_per_sym(i)
            site_types(n) = i-1
            n = n+1
            if (n.gt.L+1) stop "number of orbitals in symmetries exceeds total number of orbitals"
          end do
        end do
!         end if

        ! save project name (only the short version!!! Warning, this might be different from project_name

        call getenv('Project', qcmaquis_param%project_name)
        call getenv('CurrDir', qcmaquis_param%currdir)

        ! save wavefunction parameters
        qcmaquis_param%nactel = nel
        qcmaquis_param%ms2 = spin
        qcmaquis_param%L = L

#ifdef _MOLCAS_MPP_
        if(dmrg_host_program_settings%myrank == 0)then
#endif
        ! convert all ints neatly before passing to C
        call qcmaquis_interface_init_c(int(nel, c_int),    &
                                       int(L, c_int),      &
                                       int(spin, c_int),   &
                                       int(irrep, c_int),  &
                                       site_types,         &
                                       conv_thresh,        &
                                       int(m, c_int),      &
                                       int(nsweeps, c_int),&
                                       sweep_m_,           &
                                       nsweepm,            &
                                       trim(project_name)//c_null_char,  &
                                       logical(meas_2rdm, c_bool))
#ifdef _MOLCAS_MPP_
    endif
#endif

    end subroutine qcmaquis_interface_init


    ! Init interface by reading parameters from QCMaquis checkpoint file
    subroutine qcmaquis_interface_init_checkpoint(checkpoint)
        interface
          subroutine qcmaquis_interface_init_checkpoint_c(checkpoint) &
              bind(C, name='qcmaquis_interface_preinit_checkpoint')

              import c_char
              character(len=1,kind=c_char), dimension(*), intent(in) :: checkpoint
          end subroutine
        end interface
        character(len=*), intent(in) :: checkpoint
#ifdef _MOLCAS_MPP_
        if(dmrg_host_program_settings%myrank == 0)then
#endif
          call qcmaquis_interface_init_checkpoint_c(trim(checkpoint)//c_null_char)
#ifdef _MOLCAS_MPP_
        endif
#endif
    end subroutine

    ! Set an arbitrary QCMaquis parameter, required for parsing QCMaquis input from MOLCAS
    subroutine qcmaquis_interface_set_param(key,value_)

        interface
            subroutine qcmaquis_interface_set_param_c(key,value_) &
                bind(C, name='qcmaquis_interface_set_param')
            import c_char
            character(len=1,kind=c_char), dimension(*), intent(in) :: key, value_
            end subroutine
        end interface
        character(len=*), intent(in) :: key, value_
#ifdef _MOLCAS_MPP_
        if(dmrg_host_program_settings%myrank == 0)then
#endif
        call qcmaquis_interface_set_param_c(trim(key)//c_null_char, trim(value_)//c_null_char)
#ifdef _MOLCAS_MPP_
        endif
#endif
    end subroutine qcmaquis_interface_set_param

    ! Remove an arbitrary QCMaquis parameter
      subroutine qcmaquis_interface_remove_param(key)

        interface
            subroutine qcmaquis_interface_remove_param_c(key) &
                bind(C, name='qcmaquis_interface_remove_param')
            import c_char
            character(len=1,kind=c_char), dimension(*), intent(in) :: key
            end subroutine
        end interface
        character(len=*), intent(in) :: key
#ifdef _MOLCAS_MPP_
        if(dmrg_host_program_settings%myrank == 0)then
#endif
          call qcmaquis_interface_remove_param_c(trim(key)//c_null_char)
#ifdef _MOLCAS_MPP_
        endif
#endif
        end subroutine qcmaquis_interface_remove_param


    ! Obtain sweep statistics of the last sweep
    ! wrapper with type conversion
    ! returns:
    ! nsweeps: number of sweeps
    ! m: employed m value
    ! for twosite optimization:
    ! truncated_weight, truncated_fraction, smallest_ev:
    ! truncated weight, truncated fraction and the smallest truncated eigenvalue
    ! (otherwise the last values are zero)
    subroutine qcmaquis_interface_get_iteration_results(nsweeps, m, truncated_weight, &
                                            truncated_fraction, smallest_ev)
      interface
        subroutine qcmaquis_interface_get_iteration_results_C(nsweeps, m, truncated_weight, &
                                            truncated_fraction, smallest_ev) &
                   bind(C,name='qcmaquis_interface_get_iteration_results')
        import c_int, c_size_t, c_double
        integer(c_int) :: nsweeps
        integer(c_size_t) :: m
        real(c_double) :: truncated_weight, truncated_fraction, smallest_ev
        end subroutine
      end interface
        integer,intent(inout) :: nsweeps, m
        real*8,intent(inout) :: truncated_weight, truncated_fraction, smallest_ev
        integer(c_int) :: nsweeps_ = 0, m_ = 0
        real(c_double) :: truncated_weight_ = 0.0d0, truncated_fraction_ = 0.0d0, smallest_ev_ = 0.0d0

        ! call C interface
        call qcmaquis_interface_get_iteration_results_C(nsweeps_, int(m_, c_size_t), truncated_weight_, &
                                            truncated_fraction_, smallest_ev_)
        ! convert types
        nsweeps = int(nsweeps_)
        m = int(m_)
        truncated_weight = dble(truncated_weight_)
        truncated_fraction = dble(truncated_fraction_)
        smallest_ev = dble(smallest_ev_)
    end subroutine

    ! Prepares and updates integrals for the QCMaquis interface
    ! Input as for the FCI dumper
    ! oneint -- one-electron integrals
    ! twoint -- two-electron integrals
    ! corenergy -- core energy

    ! Now it's mostly a copy-paste from the old dmrg_task_fcidump()
    ! TODO: move to the new FCI dumper!
    subroutine qcmaquis_interface_update_integrals(oneint,twoint,corenergy)
      real*8 , intent(in), dimension(*) :: oneint
      real*8 , intent(in), dimension(*) :: twoint
      real*8 , intent(in)               :: corenergy

      ! Arrays that get passed to QCMaquis
      integer(c_int), dimension(:), allocatable :: indices
      real*8, dimension(:), allocatable :: values

      ! calculate the size of the arrays
      integer :: arr_size = 0
      integer :: max_index2

      integer :: isym, ksym, lsym, jsym, ijsym, klsym
      integer :: i, j, k, l, ij, kl, ijkl
      integer :: ndummy
      integer :: norbtot, noccend, noccendi, noccendj
      integer :: offset, offseti, offsetj, offsetk, offsetl

      ! Offset for values and indices array
      integer(c_int) :: offset_integrals

      real*8 , parameter                :: threshold = 1.0d-16
#ifdef _MOLCAS_MPP_
        if(dmrg_host_program_settings%myrank == 0)then
#endif
      offset_integrals = 1
      !     calculate total number of active orbitals
      norbtot = 0
      do ksym = 1, dmrg_symmetry%nirrep
        norbtot = norbtot + dmrg_orbital_space%nash(ksym)
      end do

      ! calculate the maximum size of the integral array, can be less b/c of symmetry
      max_index2 = norbtot*(norbtot+1)/2

      arr_size = 1 & ! core energy
                   + max_index2 & ! one-electron integrals
                   + max_index2*(max_index2+1)/2 ! two-electron integrals

      allocate(indices(4*arr_size), stat=err)
      if (err.ne.0) stop "Error in memory allocation for integral indices"
      indices=0
      allocate(values(arr_size), stat=err)
      if (err.ne.0) stop "Error in memory allocation for integrals"
      values=0.0d0
      ! 2e integrals
      if(dmrg_state%nactel > 1)then
          !> integrals are sorted in (KL|IJ) order
          offset = 0

!         KL IJ
!         define orbital offset for K
          offsetk = 0

!         sum over all irreducible representations
          do ksym = 1, dmrg_symmetry%nirrep

            if(dmrg_orbital_space%nash(ksym) == 0) cycle

            do k = 1, dmrg_orbital_space%nash(ksym)

              offsetl = 0 !define orbital offset for L

              do lsym = 1, ksym !restrict summation to prevent double counting

                if(dmrg_orbital_space%nash(lsym) == 0)cycle

!               set upper summation bound for orbital index (prevent double counting):
!               if not the same irrep l goes from 1 to number of orbitals
                if(ksym == lsym)then
                  noccend = k
                else
                  noccend = dmrg_orbital_space%nash(lsym)
                end if

                do l = 1, noccend

!                 orbital offset for I
                  offseti = 0

!                 restrict summation to prevent double counting for both irrep ISYM and orbital indices i
                  do isym = 1, ksym

                    if(dmrg_orbital_space%nash(isym) == 0)cycle

                    if(isym == ksym)then
                      noccendi = k
                    else
                      noccendi = dmrg_orbital_space%nash(isym)
                    end if

                    do i = 1, noccendi

!                     set orbital offset J
                      offsetj = 0

!                     double counting issue: irrep of J must be smaller or equal to irrep of I
                      do jsym = 1, isym
!                       fetch integrals which are nonzero by symmetry
!                       two cases have to be distinguished: IJ|KL  and IK|JL
                        if(dmrg_orbital_space%nash(jsym) == 0)cycle

                        if(isym == jsym .and. ksym == lsym)then ! first case
                          ijsym   = dmrg_symmetry%multiplication_table(isym,jsym)
                          klsym   = dmrg_symmetry%multiplication_table(ksym,lsym)
                        else                                    ! second case
                          ijsym   = dmrg_symmetry%multiplication_table(isym,ksym)
                          klsym   = dmrg_symmetry%multiplication_table(jsym,lsym)
                        end if

                        if(dmrg_symmetry%multiplication_table(ijsym,klsym) == 1)then
                          offset = (offsetk+k)*(offsetk+k-1)/2*(norbtot*(norbtot+1)/2)+ &
                                   (offsetl+l-1)*(norbtot*(norbtot+1)/2)+1

!                         prevent double counting of symmetry redundant indices: set upper summation index
!                         if IJKL in same irrep, restrict j to at most i
                          if(jsym == isym .and. ksym == lsym)then !.and.ISYM == KSYM
                            noccendj = i
!                           if LJ in irrep1 and IK in irrep2
                          else if(lsym == jsym .and. ksym == isym .and. lsym /= isym)then
!                           second restriction to prevent double counting, J<=L in KL IJ
                            if(k == i)then
                              noccendj = l
                            else ! otherwise all J are needed
                              noccendj = dmrg_orbital_space%nash(jsym)
                            end if
                          else
                            noccendj = dmrg_orbital_space%nash(jsym)
                          end if
                          offset = offset + (offseti+i)*(offseti+i-1)/2+offsetj
                          do j=1,noccendj,1
!                           check for redundant integrals
                            if(JSYM == ISYM .and. KSYM == lsym .and. ISYM == KSYM)then
                              if(k == i .and. l == i .and. j < i) cycle
                              if(k == i .and. j < l) cycle
                            end if
                            ij   = max(i+offseti,j+offsetj)*(max(i+offseti,j+offsetj)-1)/2+min(i+offseti,j+offsetj)
                            kl   = max(k+offsetk,l+offsetl)*(max(k+offsetk,l+offsetl)-1)/2+min(k+offsetk,l+offsetl)
                            ijkl = max(ij,kl)*(max(ij,kl)-1)/2+min(ij,kl)
                            if(dabs(twoint(ijkl)) < threshold)then

                              cycle
                            else
!                              write(fcidump,form1) twoint(ijkl), i+offseti, j+offsetj, k+offsetk, l+offsetl
                            values(offset_integrals) = twoint(ijkl)
                            indices(4*(offset_integrals-1)+1:4*(offset_integrals-1)+4) = &
                                (/ i+offseti, j+offsetj, k+offsetk, l+offsetl /)
                            offset_integrals = offset_integrals+1
                            end if
                          end do ! do j
                        end if ! if dmrg_symmetry%multiplication_table(ijsym,klsym) == 1

                        offsetj = offsetj + dmrg_orbital_space%nash(jsym) !update orbital offset J
                      end do ! do jsym
                    end do ! do i
                    offseti = offseti + dmrg_orbital_space%nash(isym) !update orbital offset I
                  end do ! do isym
                end do ! do l
                offsetl = offsetl + dmrg_orbital_space%nash(lsym) !update orbital offset L
              end do ! do lsym
            end do ! do k
            offsetk = offsetk + dmrg_orbital_space%nash(ksym) !update orbital offset K
          end do ! do ksym

      end if ! dmrg_state%nactel > 1

      ! one-electron integrals
      offset = 0
!     keep track of dummy indices to be ignored in one-electron
!     integrals because the symmetry ISYM < actual ISYM
      ndummy = 0
!
      do isym = 1, dmrg_symmetry%nirrep

        if(dmrg_orbital_space%nash(isym) == 0) cycle

        offset = offset + ndummy

        do i = 1, dmrg_orbital_space%nash(isym) ! only loop through same irrep

          do j = 1, i

            offset = offset + 1

            if(j == i)then
              if(dabs(oneint(offset)-(corenergy/dble(dmrg_state%nactel))) < threshold)then
              cycle
              else
                  ! subtract scaled iLive energy from diagonal elements
!                 write(fcidump,form1) oneint(offset)-(corenergy/dble(dmrg_state%nactel)), &
!                                                     i+ndummy, j+ndummy,0, 0
                values(offset_integrals) = oneint(offset)-(corenergy/dble(dmrg_state%nactel))
                indices(4*(offset_integrals-1)+1:4*(offset_integrals-1)+4) = (/ i+ndummy, j+ndummy, 0, 0 /)
                offset_integrals=offset_integrals+1
              end if
            else
              if(dabs(oneint(offset)) < threshold)then
                cycle
              else
!                 write(fcidump,form1) oneint(offset), i+ndummy, j+ndummy,0, 0
                values(offset_integrals) = oneint(offset)
                indices(4*(offset_integrals-1)+1:4*(offset_integrals-1)+4) = (/ i+ndummy, j+ndummy, 0, 0 /)
                offset_integrals=offset_integrals+1
              end if
            end if
          end do
!         add orbital offset if irrep changes
          if(i < dmrg_orbital_space%nash(isym)) offset = offset + ndummy
        end do

!       update orbital offset
        ndummy = ndummy + dmrg_orbital_space%nash(isym)
      end do

!     last step: core energy
!       write(fcidump,form3) corenergy , 0,0 ,0,0
      values(offset_integrals) = corenergy
      indices(4*(offset_integrals-1)+1:4*(offset_integrals-1)+4) = (/ 0, 0, 0, 0 /)

      ! done with integral dumping: now pass integrals to QCMaquis
      call qcmaquis_interface_update_integrals_C(indices, values, int(offset_integrals,c_int))
      if (allocated(indices)) deallocate(indices)
      if (allocated(values)) deallocate(values)
#ifdef _MOLCAS_MPP_
    endif
#endif
    end subroutine

    ! prepare CI-DEAS and/or Fiedler order guess
    ! fiedler_order_str: if do_fiedler is set, will return the comma-separated list of orbitals in the Fiedler ordering
    ! Must be initialised to the correct length!
    ! hf_occupations(optional, required for CI-DEAS): HF occupations in QCMaquis format (4-docc,3-up,2-down,1-empty)
    ! size: nact x nstates, flattened as row-major order
    subroutine qcmaquis_interface_run_starting_guess(nstates, do_fiedler, do_cideas, fiedler_order_str, hf_occupations)
      interface
        subroutine qcmaquis_interface_run_starting_guess_C(nstates, project_name, do_fiedler, do_cideas, &
                                                           fiedler_order_str, hf_occupations) &
          bind(C,name='qcmaquis_interface_run_starting_guess')
          import c_int, c_char, c_bool, c_ptr
          integer(c_int), value, intent(in) :: nstates
          logical(c_bool), value, intent(in) :: do_fiedler, do_cideas
          character(len=1,kind=c_char), dimension(*) :: project_name
          character(len=1,kind=c_char), dimension(*) :: fiedler_order_str
          type(c_ptr), value :: hf_occupations
        end subroutine qcmaquis_interface_run_starting_guess_C
      end interface
        integer, intent(in) :: nstates
        logical,intent(in) :: do_fiedler, do_cideas
        character(len=:),allocatable :: fiedler_order_str ! make sure this has the correct length
        character(len=:),allocatable :: fiedler_order_str_c
        character(len=512) :: project_name
        integer :: len_fiedler_str
        integer, intent(in), optional, dimension(*) :: hf_occupations
        integer :: L

        ! HF occupations must be converted into integer(c_int) this only works with a copy
        integer(c_int), dimension(:), allocatable, target :: hf_occupations_
        type(c_ptr) :: hf_occupations_ptr
#ifdef _MOLCAS_MPP_
  if(dmrg_host_program_settings%myrank == 0)then
#endif

        L = sum(dmrg_orbital_space%nash(1:dmrg_symmetry%nirrep))

        if (present(hf_occupations)) then
          allocate(hf_occupations_(L*nstates))
          hf_occupations_(1:L*nstates) = &
            int(hf_occupations(1:L*nstates), c_int)

          hf_occupations_ptr = c_loc(hf_occupations_)
        else
          hf_occupations_ptr = c_null_ptr
        end if

        len_fiedler_str = len(fiedler_order_str)
        project_name = trim(qcmaquis_param%project_name)

        fiedler_order_str_c = fiedler_order_str//c_null_char
        call qcmaquis_interface_run_starting_guess_C(int(nstates, c_int), &
                                                     trim(project_name)//c_null_char, &
                                                     logical(do_fiedler, c_bool), &
                                                     logical(do_cideas, c_bool), &
                                                     fiedler_order_str_c, hf_occupations_ptr)

        fiedler_order_str(:) = fiedler_order_str_c(1:len_fiedler_str)
        if (allocated(hf_occupations_)) deallocate(hf_occupations_)
#ifdef _MOLCAS_MPP_
    endif
    if(dmrg_host_program_settings%nprocs > 1 .and. dmrg_host_program_settings%runs_parallel)then
      if (do_fiedler) then
        len_fiedler_str = len(fiedler_order_str)
        call GA_Brdcst(MT_INT, [len_fiedler_str], storage_size(len_fiedler_str)/8,0)
        call GA_Brdcst(MT_BYTE, fiedler_order_str, len_fiedler_str, 0)
      end if
    end if
#endif
    end subroutine qcmaquis_interface_run_starting_guess

! Run optimization and obtain energies for all states
! replaces old run_dmrg_driver()
! Because we optimize state by state, to avoid saving and loading the MPS, we need to obtain 1- and 2-RDM immediately
! so we will return them here
! Parameters: nstates: number of states to optimize
! d1 and d2 (out): state-specific 1- and 2-RDMs as 2D arrays
! spd(out): (if requested) state-specific spin-1-RDM
!  1st dimension: state-specific RDM flattened into 1D in OpenMOLCAS format
!  2nd dimension: root
! entanglement: (optional, boolean) -- whether entanglement should be calculated and saved into HDF5 file
  subroutine qcmaquis_interface_run_dmrg(nstates, d1, d2, spd, entanglement)
    integer, intent(in) :: nstates
    real*8, intent(inout), optional :: d1(:,:), d2(:,:), spd(:,:)
    logical, intent(in), optional :: entanglement
    integer :: i, ii, nsweeps_prev
    integer :: nsweeps, m
    real*8 :: truncated_weight, truncated_fraction, smallest_ev
    logical :: hf_guess ! Check if we have a HF guess
    character(len=512) :: hf_guess_string

    logical :: entanglement_
#ifdef _MOLCAS_MPP_
    if(dmrg_host_program_settings%myrank == 0)then
#endif
    entanglement_ = .false.
    if (present(entanglement)) entanglement_ = entanglement

    ! check if HF guess is available
    hf_guess_string = ""

    hf_guess = allocated(dmrg_orbital_space%initial_occ)
    if (hf_guess) hf_guess = dmrg_orbital_space%initial_occ(1,1) > 0

    ! Enable HF guess by setting the appropriate QCMaquis parameter
    if (hf_guess) then
        call qcmaquis_interface_set_param("init_type","hf")
    end if

    ! Enable entanglement calculation if requested
    ! also enable spin density calculation
    if (entanglement_) then
        call qcmaquis_interface_set_param("MEASURE[ChemEntropy]", "1");
        call qcmaquis_interface_set_param("MEASURE[1spdm]", "1");
    else
        call qcmaquis_interface_remove_param("MEASURE[ChemEntropy]");
        call qcmaquis_interface_remove_param("MEASURE[1spdm]");
    end if

    do i=1,nstates
      ! Save the QCMaquis checkpoint name in the MOLCAS rasscf/dmrgscf.h5 file
      dmrg_file%qcmaquis_checkpoint_file(i)=trim(qcmaquis_param%project_name)//'.checkpoint_state.'// &
        trim(str(i-1))//".h5"
      ! Set the HF occupation string
      if (hf_guess) then
        hf_guess_string=""
        do ii=1,qcmaquis_param%L
          ! Construct a comma-separated string with all occupations
          hf_guess_string=trim(hf_guess_string)//trim(str(int(dmrg_orbital_space%initial_occ(ii,i),kind=8)))
          if ((ii.ne.qcmaquis_param%L)) then
            hf_guess_string=trim(hf_guess_string)//","
          end if
        end do
        ! Pass it as a hf_occ parameter to QCMaquis library
        call qcmaquis_interface_set_param("hf_occ", trim(hf_guess_string))
      end if

      call qcmaquis_interface_set_state(int(i-1,c_int))

      ! get the real number of sweeps at the beginning of the iteration
      ! if we load from a previous checkpoint, nsweeps will be nonzero
      call qcmaquis_interface_get_iteration_results(nsweeps, m, truncated_weight, &
                                                    truncated_fraction, smallest_ev)
      ! add total number of sweeps to the current number of sweeps
      nsweeps_prev = nsweeps
      call qcmaquis_interface_set_nsweeps(int(qcmaquis_param%num_sweeps+nsweeps_prev,c_int))

      call qcmaquis_interface_optimize()
      dmrg_energy%dmrg_state_specific(i) = qcmaquis_interface_get_energy()

      ! update # of sweeps and truncated weight
      ! nsweeps may not be equal to qcmaquis_param%num_sweeps+nsweeps_prev because
      ! the calculation might have needed less sweeps for convergence
      call qcmaquis_interface_get_iteration_results(nsweeps, m, truncated_weight, &
                                                    truncated_fraction, smallest_ev)

      ! get RDMs
      if (present(d1)) call qcmaquis_interface_get_1rdm_compat(d1(:,i))
      if (present(d2)) call qcmaquis_interface_get_2rdm_compat(d2(:,i))
      if (present(spd)) call qcmaquis_interface_get_spdm_compat(spd(:,i))

      ! save number of sweeps
      dmrg_energy%num_sweeps(i) = nsweeps-nsweeps_prev
      dmrg_energy%max_truncW(i) = truncated_weight

    end do
    ! SA energy
    ! If weights are present, use them, otherwise equal weights
    if (allocated(dmrg_state%weight)) then
      dmrg_energy%dmrg = dot_product(dmrg_energy%dmrg_state_specific,dmrg_state%weight)
    else
      dmrg_energy%dmrg = sum(dmrg_energy%dmrg_state_specific)/nstates
    end if
#ifdef _MOLCAS_MPP_
    endif
    if(dmrg_host_program_settings%nprocs > 1 .and. dmrg_host_program_settings%runs_parallel)then
      call GA_sync()
      call GA_Brdcst(MT_DBL, [dmrg_energy%dmrg], storage_size(dmrg_energy%num_sweeps)/8, 0)
      call GA_Brdcst(MT_DBL, [dmrg_energy%dmrg_state_specific], nstates*storage_size(dmrg_energy%dmrg_state_specific)/8, 0)
      if (present(d1)) call GA_Brdcst(MT_DBL, d1, size(d1)*storage_size(d1(1,1))/8, 0)
      if (present(d2)) call GA_Brdcst(MT_DBL, d2, size(d2)*storage_size(d2(1,1))/8, 0)
      if (present(spd)) call GA_Brdcst(MT_DBL, spd, size(spd)*storage_size(spd(1,1))/8, 0)
      call GA_Brdcst(MT_INT, dmrg_energy%num_sweeps, size(dmrg_energy%num_sweeps)*storage_size(dmrg_energy%num_sweeps)/8, 0)
      call GA_Brdcst(MT_DBL, dmrg_energy%max_truncW, size(dmrg_energy%max_truncW)*storage_size(dmrg_energy%max_truncW)/8, 0)
    endif
#endif
  end subroutine

  ! Get 1-RDM and save it into a 2D array (Used by NEVPT2)
  ! Copy-paste from 1_rdm full
  subroutine qcmaquis_interface_get_1rdm_full(d1)
    real*8, intent(inout) :: d1(:,:)
    integer(c_int) :: sz
    integer(c_int), allocatable :: indices(:)
    real*8, allocatable :: values(:)
    integer :: v,i ! counters for values and indices
    integer :: nact

    nact = qcmaquis_param%L

    sz = nact*(nact+1)/2
    d1 = 0.0d0
    allocate(values(sz))
    values = 0.0d0
    allocate(indices(2*sz))

    ! Workaround for symmetry support
    ! In case of symmetry, QCMaquis actually returns us less values
    ! So instead of allocating the arrays with the correct size
    ! (which should be the right solution)
    ! we initialise the indices with -1 and skip them later
    indices = -1

    ! obtain the rdms from qcmaquis
    call qcmaquis_interface_get_1rdm_C(indices, values, sz)

    ! copy the values into the matrix
    do v=0,sz-1
      i = 2*v
      ! TODO: make sure we don't get out of bounds here, i.e. add some
      ! checks for size(d1) eventually
      if ((indices(i+1).ge.0).and.(indices(i+2).ge.0)) then
        ! fill d1 symmetrically
        d1(int(indices(i+1),8)+1, int(indices(i+2),8)+1) = values(v+1)
        d1(int(indices(i+2),8)+1, int(indices(i+1),8)+1) = values(v+1)
      end if
    end do

    if (allocated(values)) deallocate(values)
    if (allocated(indices)) deallocate(indices)
  end subroutine

  ! Get 1-RDM or spin density in MOLCAS-compatible format
  subroutine qcmaquis_interface_get_generic1rdm_compat(method,d1)

    interface ! function that gets passed here requests either 1-RDM or spin-DM from QCMaquis
      subroutine method(indices, values, size) bind(C)
        import c_int, c_double
        integer(c_int), dimension(*) :: indices
        real(c_double), dimension(*) :: values
        integer(c_int), value :: size
      end subroutine
    end interface
    real*8, intent(inout) :: d1(:)
    integer(c_int) :: sz
    integer(c_int), allocatable :: indices(:)
    real*8, allocatable :: values(:)
    integer :: v,i ! counters for values and indices
    sz = size(d1)
    d1 = 0.0d0
    allocate(values(sz))
    values = 0.0d0
    allocate(indices(2*sz))

    ! Workaround for symmetry support
    ! In case of symmetry, QCMaquis actually returns us less values
    ! So instead of allocating the arrays with the correct size
    ! (which should be the right solution)
    ! we initialise the indices with -1 and skip them later
    indices = -1

    ! obtain the rdms from qcmaquis
    call method(indices, values, sz)

    ! copy the values into the matrix
    do v=0,sz-1
      i = 2*v
      ! TODO: make sure we don't get out of bounds here, i.e. add some
      ! checks for size(d1) eventually
      if ((indices(i+1).ge.0).and.(indices(i+2).ge.0)) then
        d1(tri(indices(i+1)+int(1,kind=4), indices(i+2)+int(1,kind=4))) = values(v+1)
      end if
    end do

    if (allocated(values)) deallocate(values)
    if (allocated(indices)) deallocate(indices)
  end subroutine

  ! Specializations of the above
  ! Get 1-RDM
  subroutine qcmaquis_interface_get_1rdm_compat(d1)
    real*8, intent(inout) :: d1(:)

    call qcmaquis_interface_get_generic1rdm_compat(qcmaquis_interface_get_1rdm_C,d1)
  end subroutine

  ! Get spin-1-RDM
  subroutine qcmaquis_interface_get_spdm_compat(d1)
    real*8, intent(inout) :: d1(:)

    call qcmaquis_interface_get_generic1rdm_compat(qcmaquis_interface_get_spdm_C,d1)
  end subroutine

  ! Get 2-RDM and save it into an 4-dimensional array. (Used by NEVPT2)
  subroutine qcmaquis_interface_get_2rdm_full(d2)
    real*8, intent(inout) :: d2(:,:,:,:)
    integer(c_int) :: sz ! size

    ! indices and values that are obtained from QCMaquis interface
    integer(c_int), allocatable :: indices(:)
    real*8, allocatable :: values(:)
    integer :: nact
    integer :: vv,ii ! counters for values and indices
    integer :: ij,jk,kl,li
    integer :: i,j,k,l
    ! temporary rdms
    real*8, allocatable :: rdm2(:,:,:,:)
    nact = qcmaquis_param%L
    ! calculate the size of 2-RDM in QCMaquis
    ! TODO: this should be of the value n2*(n2+1)/2 with n2=nact*(nact+1)/2
    d2 = 0.0d0
    sz=0
    do i=1,nact
    do j=1,nact
    do k=min(i,j),nact
    do l=k,nact
      sz=sz+1
    end do
    end do
    end do
    end do

    allocate(values(sz))
    values = 0.0d0
    allocate(indices(4*sz))
    ! initialise indices to -1, see in 1RDM code why
    indices = -1
    allocate(rdm2(nact,nact,nact,nact))
    rdm2 = 0.0d0
    d2 = 0.0d0
    ! obtain the rdms from qcmaquis
    call qcmaquis_interface_get_2rdm_C(indices, values, sz)

    ! copy the values into the matrix
    do vv=0,sz-1
      ii = 4*vv
      ! TODO: make sure we don't get out of bounds here, i.e. add some
      ! checks for size(d2) eventually
      ij = indices(ii+1)+1
      jk = indices(ii+2)+1
      kl = indices(ii+3)+1
      li = indices(ii+4)+1
      if ((ij+jk+kl+li).eq.0) cycle ! skip empty indices
      rdm2(ij,jk,kl,li) = values(vv+1)
      rdm2(kl,li,ij,jk) = values(vv+1)
      rdm2(jk,ij,li,kl) = values(vv+1)
      rdm2(li,kl,jk,ij) = values(vv+1)
    end do

    ! now symmetrise the RDM for compatibility with MOLCAS (copy-paste from old interface)
    do i=1,nact
      do j=1,nact
        do k=1,nact
          do l=1,nact
            if(j.eq.k)then
              d2(i,j,k,l)=rdm2(i,j,k,l)
            else
              d2(i,j,k,l)=rdm2(i,j,k,l)+rdm2(i,k,j,l)
            end if
          end do
        end do
      end do
    end do

    if (allocated(values)) deallocate(values)
    if (allocated(indices)) deallocate(indices)
    if (allocated(rdm2)) deallocate(rdm2)
  end subroutine qcmaquis_interface_get_2rdm_full
  ! Get 2-RDM in MOLCAS-compatible format (packed form)
  ! Warning, a lot of copy-paste from the old interface
  subroutine qcmaquis_interface_get_2rdm_compat(d2)
    real*8, intent(inout) :: d2(:)
    integer :: nact
    integer :: ij,kl
    integer :: i,j,k,l,ijkl
    ! temporary rdms
    real*8, allocatable :: rdm2T(:,:,:,:)

    nact = qcmaquis_param%L
    d2 = 0.0d0

    allocate(rdm2T(nact,nact,nact,nact))
    call qcmaquis_interface_get_2rdm_full(rdm2T)
    ! copy symmetrised RDM to the output
    ij=0
    ijkl=1
    do i=1,nact
      do j=1,i
        ij=ij+1
        kl=0
        do k=1,nact
          do l=1,k
            kl=kl+1
            if(ij.ge.kl)then
              d2(ijkl)=rdm2T(i,k,l,j)*0.5d0
              ! We need a factor of 0.5 for 2-RDM here
              ! as used in the previous interface
              ijkl=ijkl+1
            end if
          end do
        end do
      end do
    end do

    if (allocated(rdm2T)) deallocate(rdm2T)

  end subroutine

    ! redirect QCMaquis stdout to file
  subroutine qcmaquis_interface_stdout(filename)
    interface
      subroutine qcmaquis_interface_stdout_C(filename) bind(C, name='qcmaquis_interface_stdout')
        import c_char
        character(len=1,kind=c_char), dimension(*), intent(in) :: filename
      end subroutine
    end interface
    character(len=*), intent(in) :: filename

#ifdef _MOLCAS_MPP_
  if(dmrg_host_program_settings%myrank == 0)then
#endif
    call qcmaquis_interface_stdout_C(trim(filename)//c_null_char)
#ifdef _MOLCAS_MPP_
  endif
#endif
  end subroutine

  function qcmaquis_interface_get_overlap(state) result(res)

    interface
      function qcmaquis_interface_get_overlap_c(filename) result(res)  bind(C, name='qcmaquis_interface_get_overlap')
        import c_double, c_char
        real(c_double) :: res
        character(len=1,kind=c_char), dimension(*), intent(in) :: filename
      end function
    end interface

    real(c_double) :: res
    integer,intent(in) :: state
    integer :: state_real
    character(len=2300) :: filename
    ! TODO: get rid of the old function

#ifdef _MOLCAS_MPP_
  if(dmrg_host_program_settings%myrank == 0)then
#endif
    state_real = state - 1
    call file_name_generator(state_real,"checkpoint_state.",".h5", filename)

    res = qcmaquis_interface_get_overlap_c(trim(filename)//c_null_char)
#ifdef _MOLCAS_MPP_
  endif
  if(dmrg_host_program_settings%nprocs > 1 .and. dmrg_host_program_settings%runs_parallel)then
    call GA_sync()
    Call GA_Brdcst(MT_DBL, [res], storage_size(res)/8, 0)
  endif
#endif
  end function

  integer(kind=4) function tri(i,j)
    integer(kind=4), intent(in) :: i,j
    tri=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
  end function tri

  ! remove the comments from the string starting with comment_string
  subroutine remove_comment(instring, comment_string)
    character(len=*),intent(inout) :: instring
    character(len=*),intent(in) :: comment_string

    integer :: idx
    idx = index(trim(instring), trim(comment_string))
    if (idx.gt.0) then
      if (idx.gt.1) then
          instring = trim(adjustl(instring(1:idx-1)))
      else
          instring = ""
      end if
    end if
  end subroutine

  ! Delete the checkpoint corresponding to state state
  subroutine qcmaquis_interface_delete_chkp(state)
    integer, intent(in) :: state
    integer :: state_real
    character(len=2300) :: filename
    ! TODO: get rid of the old function
#ifdef _MOLCAS_MPP_
    if(dmrg_host_program_settings%myrank == 0)then
#endif
      state_real = state - 1
      call file_name_generator(state_real,"checkpoint_state.",".h5",filename)

      call system("rm -rf "//trim(filename))
#ifdef _MOLCAS_MPP_
    endif
#endif
  end subroutine


  ! NEVPT2
  subroutine qcmaquis_interface_measure_and_save_4rdm(state)
    interface
      subroutine qcmaquis_interface_measure_and_save_4rdm_C(state) &
        bind(C,name='qcmaquis_interface_measure_and_save_4rdm')
        import c_int
        integer(c_int), value :: state
      end subroutine
    end interface

    integer,intent(in) :: state
#ifdef _MOLCAS_MPP_
    if(dmrg_host_program_settings%myrank == 0)then
#endif
      call qcmaquis_interface_measure_and_save_4rdm_C(int(state,c_int))
#ifdef _MOLCAS_MPP_
    endif
#endif
  end subroutine

  subroutine qcmaquis_interface_measure_and_save_3rdm(state)
    interface
      subroutine qcmaquis_interface_measure_and_save_3rdm_C(state) &
        bind(C,name='qcmaquis_interface_measure_and_save_3rdm')
        import c_int
        integer(c_int), value :: state
      end subroutine
    end interface

    integer,intent(in) :: state
#ifdef _MOLCAS_MPP_
    if(dmrg_host_program_settings%myrank == 0)then
#endif
      call qcmaquis_interface_measure_and_save_3rdm_C(int(state,c_int))
#ifdef _MOLCAS_MPP_
    endif
#endif
  end subroutine

  subroutine qcmaquis_interface_measure_and_save_trans3rdm(state, bra_state)
    interface
      subroutine qcmaquis_interface_measure_and_save_trans3rdm_C(state, bra_state) &
        bind(C,name='qcmaquis_interface_measure_and_save_trans3rdm')
        import c_int
        integer(c_int), value :: state, bra_state
      end subroutine
    end interface

    integer,intent(in) :: state, bra_state
#ifdef _MOLCAS_MPP_
    if(dmrg_host_program_settings%myrank == 0)then
#endif
      call qcmaquis_interface_measure_and_save_trans3rdm_C(int(state,c_int), int(bra_state,c_int))
#ifdef _MOLCAS_MPP_
    endif
#endif
  end subroutine

  ! get number of 4-RDM elements measured by QCMaquis
  ! if slice (array with 4 indices) is present, get number of
  ! elements with those indices as the first 4 indices
  function qcmaquis_interface_get_4rdm_elements(slice) result(res)
    interface
      function qcmaquis_interface_get_4rdm_elements_C(L, slice) result(res) &
        bind(C,name='qcmaquis_interface_get_4rdm_elements')
        import c_int, c_ptr
        integer(c_int), value :: L
        type(c_ptr), value :: slice
        integer(c_int) :: res
      end function
    end interface

    integer, dimension(4), optional, intent(in), target :: slice
    integer :: res
    type(c_ptr) :: slice_ptr

#ifdef _MOLCAS_MPP_
    if(dmrg_host_program_settings%myrank == 0)then
#endif
    if (present(slice)) then
        slice_ptr = c_loc(slice)
      else
        slice_ptr = c_null_ptr
    end if

    res = qcmaquis_interface_get_4rdm_elements_C(int(qcmaquis_param%L, c_int), slice_ptr)
#ifdef _MOLCAS_MPP_
    endif
#endif
  end function qcmaquis_interface_get_4rdm_elements

  ! same for 3-RDM elements
    function qcmaquis_interface_get_3rdm_elements(bra_neq_ket,slice) result(res)
    interface
      function qcmaquis_interface_get_3rdm_elements_C(L, bra_neq_ket, slice) result(res) &
        bind(C,name='qcmaquis_interface_get_3rdm_elements')
        import c_int, c_ptr, c_bool
        integer(c_int), value :: L
        logical(c_bool), value :: bra_neq_ket
        type(c_ptr), value :: slice
        integer(c_int) :: res
      end function
    end interface

    integer, dimension(2), optional, intent(in), target :: slice
    integer :: res
    logical,intent(in) :: bra_neq_ket
    type(c_ptr) :: slice_ptr

#ifdef _MOLCAS_MPP_
    if(dmrg_host_program_settings%myrank == 0)then
#endif
      if (present(slice)) then
          slice_ptr = c_loc(slice)
        else
          slice_ptr = c_null_ptr
      end if

      res = qcmaquis_interface_get_3rdm_elements_C(int(qcmaquis_param%L, c_int), logical(bra_neq_ket, c_bool), slice_ptr)
#ifdef _MOLCAS_MPP_
    endif
#endif
  end function qcmaquis_interface_get_3rdm_elements

  subroutine qcmaquis_interface_prepare_hirdm_template(filename, state, tpl, state_j)
    interface
      subroutine qcmaquis_interface_prepare_hirdm_template_C(filename, state, tpl, state_j) &
        bind(C,name='qcmaquis_interface_prepare_hirdm_template')
        import c_int, c_char
        character(len=1,kind=c_char), dimension(*), intent(in) :: filename
        integer(c_int), intent(in), value :: state
        integer(c_int), intent(in), value :: tpl ! enum
        integer(c_int), intent(in), value :: state_j
      end subroutine
    end interface

    character(len=*), intent(in) :: filename
    integer, intent(in) :: state ! bra state
    integer(c_int), intent(in) :: tpl ! enum that determines if we want 3-RDM or 4-RDM
    integer, intent(in), optional :: state_j ! ket state for 3-rdm
    integer(c_int) :: state_j_

#ifdef _MOLCAS_MPP_
        if(dmrg_host_program_settings%myrank == 0)then
#endif
      if (present(state_j)) then
        state_j_ = state_j
      else
        if (tpl .eq. TEMPLATE_TRANSITION_3RDM) then
          stop "Missing bra state missing for transition 3-RDM QCMaquis template"
        end if
        state_j_ = 0
      end if

      call qcmaquis_interface_prepare_hirdm_template_C(trim(filename)//c_null_char, &
                                                    int(state, c_int), tpl, state_j_)
#ifdef _MOLCAS_MPP_
    endif
#endif
  end subroutine qcmaquis_interface_prepare_hirdm_template

    !! do we really need this with modern Fortran?
    subroutine qcmaquis_interface_deinit

    if(allocated(dmrg_state%iroot))                   deallocate(dmrg_state%iroot)
    if(allocated(dmrg_state%weight))                  deallocate(dmrg_state%weight)
    if(allocated(dmrg_orbital_space%initial_occ))     deallocate(dmrg_orbital_space%initial_occ)
    if(allocated(dmrg_energy%dmrg_state_specific))    deallocate(dmrg_energy%dmrg_state_specific)
    if(allocated(dmrg_energy%num_sweeps))             deallocate(dmrg_energy%num_sweeps)
    if(allocated(dmrg_energy%max_truncW))             deallocate(dmrg_energy%max_truncW)
    if(allocated(dmrg_input%qcmaquis_input))          deallocate(dmrg_input%qcmaquis_input)
    if(allocated(dmrg_file%qcmaquis_checkpoint_file)) deallocate(dmrg_file%qcmaquis_checkpoint_file)

  end subroutine qcmaquis_interface_deinit


end module qcmaquis_interface

