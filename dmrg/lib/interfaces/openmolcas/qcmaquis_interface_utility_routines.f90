!!  dmrg-interface-utils: interface to the Maquis DMRG program for various
!!                        quantum-chemistry program packages.
!!  Copyright 2013-2018 Leon Freitag, Erik Hedegaard, Sebastian Keller,
!!                      Stefan Knecht, Yingjin Ma, Christopher Stein
!!                      and Markus Reiher
!!                      Laboratory for Physical Chemistry, ETH Zurich
!!
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

module qcmaquis_interface_utility_routines

  use qcmaquis_interface_cfg

contains

      SUBROUTINE pretty_print_util(AMATRX,ROWLOW,ROWHI,COLLOW,COLHI,ROWDIM,COLDIM,NCTL,LUPRI)
!.......................................................................
! Revised 15-Dec-1983 by Hans Jorgen Aa. Jensen.
!         16-Jun-1986 hjaaj ( removed Hollerith )
!
! OUTPUT PRINTS A REAL MATRIX IN FORMATTED FORM WITH NUMBERED ROWS
! AND COLUMNS.  THE INPUT IS AS FOLLOWS;
!
!        AMATRX(',').........MATRIX TO BE OUTPUT
!
!        ROWLOW..............ROW NUMBER AT WHICH OUTPUT IS TO BEGIN
!
!        ROWHI...............ROW NUMBER AT WHICH OUTPUT IS TO END
!
!        COLLOW..............COLUMN NUMBER AT WHICH OUTPUT IS TO BEGIN
!
!        COLHI...............COLUMN NUMBER AT WHICH OUTPUT IS TO END
!
!        ROWDIM..............ROW DIMENSION OF AMATRX(',')
!
!        COLDIM..............COLUMN DIMENSION OF AMATRX(',')
!
!        NCTL................CARRIAGE CONTROL FLAG; 1 FOR SINGLE SPACE
!                                                   2 FOR DOUBLE SPACE
!                                                   3 FOR TRIPLE SPACE
!                            hjaaj: negative for 132 col width
!
! THE PARAMETERS THAT FOLLOW MATRIX ARE ALL OF TYPE INTEGER.  THE
! PROGRAM IS SET Up TO HANDLE 5 COLUMNS/PAGE WITH A 1P,5D24.15 FORMAT
! FOR THE COLUMNS.  IF A DIFFERENT NUMBER OF COLUMNS IS REQUIRED,
! CHANGE FORMATS 1000 AND 2000, AND INITIALIZE KCOL WITH THE NEW NUMBER
! OF COLUMNS.
!
! AUTHOR;  NELSON H.F. BEEBE, QUANTUM THEORY PROJECT, UNIVERSITY OF
!          FLORIDA, GAINESVILLE
! REVISED; FEBRUARY 26, 1971
!
!.......................................................................
!
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER   ROWLOW,ROWHI,COLLOW,COLHI,ROWDIM,COLDIM,BEGIN,KCOL
      DIMENSION AMATRX(ROWDIM,COLDIM)
      CHARACTER*1 ASA(3), BLANK, CTL
      CHARACTER   PFMT*20, COLUMN*8
      ! LOGICAL, external :: IS_NAN
      PARAMETER (ZERO=0.D00, KCOLP=5, KCOLN=8)
      PARAMETER (FFMIN=1.D-3, FFMAX = 1.D3)
      DATA COLUMN/'Column  '/, BLANK/' '/, ASA/' ', '0', '-'/
!
      IF (ROWHI.LT.ROWLOW) GO TO 3
      IF (COLHI.LT.COLLOW) GO TO 3
!
      AMAX = ZERO
      N_NAN = 0
      DO J = COLLOW,COLHI
         DO I = ROWLOW,ROWHI
!           IF ( IS_NAN(AMATRX(I,J),AMATRX(I,J)) ) THEN
!              N_NAN = N_NAN + 1
!           ELSE
               AMAX = MAX( AMAX, ABS(AMATRX(I,J)) )
!           END IF
         END DO
     END DO
      IF (N_NAN .GT. 0) WRITE (LUPRI,'(/T6,A,I10,A)') 'WARNING: matrix contains',N_NAN,' NaN.'
      IF (AMAX <= 1.0d-20) THEN
         WRITE (LUPRI,'(/T6,A)') 'Zero matrix.'
         GO TO 3
      END IF
      IF (FFMIN .LE. AMAX .AND. AMAX .LE. FFMAX) THEN
!        use F output format
         PFMT = '(A1,I7,2X,8F18.11)'
         thrpri = 0.5D-10
      ELSE
!        use 1PE output format
         PFMT = '(A1,I7,2X,1P,8E15.6)'
         thrpri = 1.0D-10*AMAX
      END IF
!
      IF (NCTL .LT. 0) THEN
         KCOL = KCOLN
      ELSE
         KCOL = KCOLP
      END IF
      MCTL = ABS(NCTL)
      IF ((MCTL.LE.3).AND.(MCTL.GT.0)) THEN
         CTL = ASA(MCTL)
      ELSE
         CTL = BLANK
      END IF
!
      LAST = MIN(COLHI,COLLOW+KCOL-1)
      DO 2 BEGIN = COLLOW,COLHI,KCOL
         WRITE (LUPRI,1000) (COLUMN,I,I = BEGIN,LAST)
         DO 1 K = ROWLOW,ROWHI
            DO 4 I = BEGIN,LAST
               IF (abs(AMATRX(K,I)).gt.thrpri) GO TO 5
    4       CONTINUE
         GO TO 1
    5       WRITE (LUPRI,PFMT) CTL,K,(AMATRX(K,I), I = BEGIN,LAST)
    1    CONTINUE
    2    END DO
    LAST = MIN(LAST+KCOL,COLHI)
    3 WRITE(LUPRI,'(A)') '    ==== End of matrix output ===='
      RETURN
 1000 FORMAT (/10X,8(5X,A6,I4))

      end subroutine pretty_print_util

      LOGICAL FUNCTION IS_NAN(XA,XB)
!
!     May 2010, Hans Joergen Aa. Jensen
!     Purpose: IS_NAN(X,X) is true iff X is NAN
!
      REAL*8 XA, XB
      IS_NAN = XA .NE. XB
      END FUNCTION IS_NAN

      subroutine lower_to_upper(str)

        character(*), intent(inout) :: str

        integer                     :: i
        do i = 1, len(str)
          select case(str(i:i))
            case("a":"z")
              str(i:i) = achar(iachar(str(i:i))-32)
          end select
        end do
      end subroutine lower_to_upper


      subroutine find_qcmaquis_keyword(input_string,strdim,keyword,location)

        integer,            intent(in)                       :: strdim
        character(len=500), intent(in), dimension(*)         :: input_string
        character(len=500), intent(in)                       :: keyword
        integer,            intent(inout)                    :: location

        integer                                              :: i
        character(len=500)                                   :: string

        location = -1

        do i = 1, strdim,2
          string(1:500) = " "
          string        = trim(input_string(i))
          call lower_to_upper(string)
          if(trim(string) == trim(keyword))then
            location = i
            exit
          end if
        end do

      end subroutine find_qcmaquis_keyword

      ! Quick integer-to-string converter, not to mess with the conversion somewhere else.
      ! From https://stackoverflow.com/questions/1262695/convert-integers-to-strings-to-create-output-filenames-at-run-time
      character(len=30) function str(k)
          integer, intent(in) :: k
          write (str, *) k
          str = adjustl(str)
      end function str

      !! For Fiedler ordering: calculate the length of a string that would fit the Fiedler ordering
      !! i.e., for a given integer N, give the length of a string that fits numbers "1,2,...,N" with commas included
      integer*8 function fiedlerorder_length(L) result(res)
          implicit none
          integer, intent(in) :: L
          integer p,n,c ! temporary variables
          res = 0
          c = L
          do while (c.gt.0)
            ! Calculate the maximum power of 10 that does not exceed c
            p = int(log10(dble(c)))
            n = 10**p

            ! In total we will have (c-n+1) (p+1)-digit numbers, which will each,obviously,
            ! occupy p+1 digits plus a comma
            res = res + (c-n+1)*(p+2)

            ! subtract the numbers we already accounted for
            c = c - (l-n+1)
          end do

          res = res - 1 ! last number does not need a comma

      end function fiedlerorder_length

      ! ===========================================================================
      ! This subroutine generate the file name base on the prototype_name -Yingjin
      !  Input  : iroot
      !         : prototype_name
      !         : suffix
      !  Output : generated_name
      ! ===========================================================================

      subroutine file_name_generator(iroot, prototype_name, suffix, generated_name)
          implicit none
          integer*8, intent(in)                :: iroot
          character(len=*), intent(in)       :: prototype_name
          character(len=*), intent(in)       :: suffix
          character(len=2300), intent(inout) :: generated_name

          generated_name = ''

          if (len_trim(suffix).ne.0) then
            generated_name = trim(qcmaquis_param%currdir)//'/'//trim(qcmaquis_param%project_name)//'.'// &
                          trim(prototype_name)//trim(str(iroot))//trim(suffix)
          else ! should no longer be in use
            generated_name = trim(prototype_name)//trim(str(iroot))//"."//trim(str(iroot))
          endif

      end subroutine file_name_generator


! !**********************************************************************
! ! This subroutine is obsolete as FCIDUMP writing is achieved now directly in OpenMolcas
! ! and QCMaquis should no longer need FCIDUMP files anymore
! ! So it has to be removed in the future
! ! However, as of writing this, FCIDUMP files written in OpenMolcas have not been tested for
! ! compatibility with QCMaquis, so this subroutine remains here for compatibility
!
      subroutine qcmaquis_interface_fcidump(oneint,twoint,corenergy)
!     -----------------------------------------------------------------
!
!     purpose: write one- and two-electron integrals for a given active
!              space to disk in a formatted file.
!
!     filename: FCIDUMP
!     -----------------------------------------------------------------
!
      real*8 , intent(in), dimension(*) :: oneint
      real*8 , intent(in), dimension(*) :: twoint
      real*8 , intent(in)               :: corenergy
!     -----------------------------------------------------------------
      integer                           :: isym, ksym, lsym, jsym, ijsym, klsym
      integer                           :: i, j, k, l, ij, kl, ijkl
      integer                           :: ndummy,lenorbstring
      integer                           :: norbtot, noccend, noccendi, noccendj
      integer                           :: offset, offseti, offsetj, offsetk, offsetl
      character(len=30)                 :: form1
      character(len=30)                 :: form2
      character(len=30)                 :: form3
      character(len=5000)               :: orbstring
      integer, parameter                :: fcidump = 99
      real*8 , parameter                :: threshold = 1.0d-16
!     -----------------------------------------------------------------

!
!     define printing format | Using the same G21.12 as molpro
      form1="(G21.12, 4X, I6, I6, I6, I6)"
      form2="(A11,I3,A7,I2,A5,I2,A1)"
      form3="(G21.12, 4X, I6, I6, I6, I6)"

!     calculate total number of active orbitals
      norbtot = 0
      do ksym = 1, dmrg_symmetry%nirrep
        norbtot = norbtot + dmrg_orbital_space%nash(ksym)
      end do

!     Print header of FCIDUMP file using the MOLPRO format
      open(fcidump,file='FCIDUMP',status='replace',form='formatted',    &
           action='readwrite',position='rewind')

      write(fcidump,form2) ' &FCI NORB=', norbtot , ',NELEC=',          &
      dmrg_state%nactel, ',MS2=', qcmaquis_param%ms2, ','

      write(fcidump,"(A)",advance='no') '  ORBSYM='
      do isym = 1, dmrg_symmetry%nirrep
        if(dmrg_orbital_space%nash(isym) /= 0)then
          do i = 1, dmrg_orbital_space%nash(isym)
            write(fcidump,"(I1,A1)",advance='no') isym,','
          end do
        end if
      end do
      !> remove trailing ',' from FCIDUMP orbital symmetry string
      backspace(fcidump); read(fcidump,'(a)') orbstring; backspace(fcidump)
      lenorbstring = len_trim(orbstring); write(fcidump,'(a)',advance='no') orbstring(1:lenorbstring-1)

      write(fcidump,*)
      write(fcidump,"(A7,I1)") '  ISYM=',dmrg_state%irefsm
      write(fcidump,"(A5)") ' &END'

      if(dmrg_state%nactel > 1)then

        select case(dmrg_host_program_settings%dmrg_host_program)

        case ('dirac  ')

          if(dmrg_symmetry%nirrep > 1) stop 'symmetry handling not implemented for Dirac interface'

          print *, 'dump integrals in Dirac interface mode'
          !> next two-electron integrals
          !> integrals are sorted in (IJ|KL) order
          offset = 0

!         IJ KL
!         sum over all irreducible representations
          do isym = 1, dmrg_symmetry%nirrep
            if(dmrg_orbital_space%nash(isym) == 0) cycle
            do i = 1, dmrg_orbital_space%nash(isym)
              do j = 1, i
                do k = 1, i
                  do l = 1, k

                    !> check for redundant integrals
                    if(k == i .and. l > j) cycle
                    offset = (l-1)*(dmrg_orbital_space%nash(isym)**3)+(k-1)*(dmrg_orbital_space%nash(isym)**2)+&
                             (j-1)*dmrg_orbital_space%nash(isym)+i
                    !print *, 'offset for (ij|kl) ',i,j,k,l,' is ==> ',offset
                    if(dabs(twoint(offset)) < threshold) cycle
                    write(fcidump,form1) twoint(offset), i, j, k, l

                  end do ! do l
                end do ! do k
              end do ! do j
            end do ! do i
          end do ! do isym

        case default
          !> next two-electron integrals
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
                              if(dmrg_host_program_settings%dmrg_host_program == 'dalton ')then
                                if(k == i .and. l == i .and. j < i)then
                                   offset = offset + 1
                                   cycle
                                end if
                                if(k == i .and. j < l)then
                                   offset = offset + 1
                                   cycle
                                end if
                              else
                                if(k == i .and. l == i .and. j < i) cycle
                                if(k == i .and. j < l) cycle
                              end if
                            end if
                            if(dmrg_host_program_settings%dmrg_host_program == 'dalton ')then
                              ijkl   = offset
                              offset = offset + 1
                            else
                              ij   = max(i+offseti,j+offsetj)*(max(i+offseti,j+offsetj)-1)/2+min(i+offseti,j+offsetj)
                              kl   = max(k+offsetk,l+offsetl)*(max(k+offsetk,l+offsetl)-1)/2+min(k+offsetk,l+offsetl)
                              ijkl = max(ij,kl)*(max(ij,kl)-1)/2+min(ij,kl)
                            end if
                            if(dabs(twoint(ijkl)) < threshold)then
!                              write(fcidump,form1) 0.000000000000E-15, i+offseti, j+offsetj, k+offsetk, l+offsetl
                              cycle
                            else
                              write(fcidump,form1) twoint(ijkl), i+offseti, j+offsetj, k+offsetk, l+offsetl
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

        end select

      end if ! dmrg_state%nactel > 1

!     next step: one-electron integrals
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
              if (dmrg_host_program_settings%dmrg_host_program(1:7) == 'molcas ')then
                 if(dabs(oneint(offset)-(corenergy/dble(dmrg_state%nactel))) < threshold)then
                   cycle
                 else
                   ! subtract scaled inactive energy from diagonal elements
                   write(fcidump,form1) oneint(offset)-(corenergy/dble(dmrg_state%nactel)), &
                                                        i+ndummy, j+ndummy,0, 0
                 end if
              else
                 if(dabs(oneint(offset)) < threshold)then
                   cycle
                 else
                   write(fcidump,form1) oneint(offset), i+ndummy, j+ndummy,0, 0
                 end if
              end if
            else
              if(dabs(oneint(offset)) < threshold)then
                cycle
              else
                write(fcidump,form1) oneint(offset), i+ndummy, j+ndummy,0, 0
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
      write(fcidump,form3) corenergy , 0,0 ,0,0
!
      close(unit=fcidump,status='KEEP')

      end subroutine qcmaquis_interface_fcidump

      ! *********************************************************************
      subroutine print_dmrg_info(lupri,fmt2,switch,start_guess,nroots,thre)

        integer,            intent(in)    :: lupri
        integer,            intent(in)    :: switch
        integer,            intent(in)    :: nroots
        double precision,   intent(in)    :: thre
        character(len=8),   intent(in)    :: fmt2
        character(len=100), intent(inout) :: start_guess

        character(len=500)                :: mstates
        character(len=500)                :: sweeps
        character(len=500)                :: sweeps_tolerance
        character(len=500)                :: jcd_tolerance
        character(len=500)                :: svd_tolerance_initial
        character(len=500)                :: svd_tolerance_final
        character(len=500)                :: line
        integer                           :: i

        if(dmrg_host_program_settings%myrank == 0)then
          mstates               = '0'
          sweeps                = '0'
          svd_tolerance_initial = '1e-50'
          svd_tolerance_final   = ' '
          sweeps_tolerance      = ' '
          jcd_tolerance         = ' '

          write(      sweeps_tolerance,'(e10.3)') thre
          write(         jcd_tolerance,'(e10.3)') thre*0.001  ! same as molcas for Davidson
          write(   svd_tolerance_final,'(e10.3)') thre*0.001  !  in order to match Davidson

          do i = 1, size(dmrg_input%qcmaquis_input),2
            line(1:500) = dmrg_input%qcmaquis_input(i)(1:500)
            call lower_to_upper(line)
            if(trim(line) == 'TRUNCATION_INITIAL')then
              svd_tolerance_initial = trim(dmrg_input%qcmaquis_input(i+1))
            else if(trim(line) == 'TRUNCATION_FINAL')then
              svd_tolerance_final        = trim(dmrg_input%qcmaquis_input(i+1))
            else if(trim(line) == 'IETL_JCD_TOL')then
              jcd_tolerance         = trim(dmrg_input%qcmaquis_input(i+1))
            else if(trim(line) == 'CONV_THRESH')then
              sweeps_tolerance     = trim(dmrg_input%qcmaquis_input(i+1))
            else if(trim(line) == 'MAX_BOND_DIMENSION')then
              mstates              = trim(dmrg_input%qcmaquis_input(i+1))
            else if(trim(line) == 'NSWEEPS')then
              sweeps               = trim(dmrg_input%qcmaquis_input(i+1))
            end if
          end do

          if(trim(mstates) == '0') mstates = 'dynamically changing (according to sweep_bond_dimensions)'

          write(lupri,fmt2//'a,t45,5x,a)') 'Number of renormalized states           ', trim(mstates)
    !       if(dmrg_warmup%doCIDEAS)then
    !         write(lupri,fmt2//'a,t45,5x,a)') 'Start guess in warm-up sweep            ','CI-DEAS'
    !       else
    !         write(lupri,fmt2//'a,t45,5x,a)') 'Start guess in warm-up sweep            ', trim(start_guess)
    !       end if
          write(lupri,fmt2//'a,t45,5x,a)') '(Max) number of sweeps                  ', trim(sweeps)
          write(lupri,fmt2//'a,t45,5x,a)') 'Convergence threshold (sweep tolerance) ', trim(sweeps_tolerance)
          write(lupri,fmt2//'a,t45,5x,a)') 'Jacobi-Davidson threshold               ', trim(jcd_tolerance)
          write(lupri,fmt2//'a,t45,5x,a)') 'SVD truncation threshold (initial)      ', trim(svd_tolerance_initial)
          write(lupri,fmt2//'a,t45,5x,a)') 'SVD truncation threshold (final)        ', trim(svd_tolerance_final)

          !> output before optimization
    !       if (switch.eq.1) then
    !         if(dmrg_warmup%doFIEDLER) write(lupri,fmt2//'a,t45     )') 'Fiedler vector for orbital ordering     '
    !       endif
          ! switch==2 disabled for now as we don't write the ordering into files, and anyway the ordering can be only
          ! the same for all states
    !     else
    !       write(lupri,fmt2//'a,t45)') 'I am not master - no DMRG info print  '
        endif

    end subroutine print_dmrg_info


end module qcmaquis_interface_utility_routines

