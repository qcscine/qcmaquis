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

module qcmaquis_interface_caspt2
! Module for Fortran-C interoperability with the new (non-Python) QCMaquis interface

  use iso_c_binding
  use qcmaquis_interface_cfg
  implicit none
#ifdef _MOLCAS_MPP_
#include "mafdecls.fh"
#endif

  contains
  subroutine qcmaquis_interface_caspt2_init(epsa, nasht)
    interface
      subroutine qcmaquis_interface_caspt2_init_C(epsa, nasht) &
        bind(C,name='qcmaquis_interface_caspt2_init')
        import c_int, c_double
        real(c_double), dimension(*) :: epsa
        integer(c_int), intent(in), value :: nasht
      end subroutine
    end interface

    real(c_double), dimension(:) :: epsa
    integer, intent(in) :: nasht

    call qcmaquis_interface_caspt2_init_C(epsa, int(nasht, c_int))

  end subroutine qcmaquis_interface_caspt2_init
  end module qcmaquis_interface_caspt2
