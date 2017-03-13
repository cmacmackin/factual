!
!  chebyshev.f90
!  This file is part of FACTUAL.
!
!  Copyright 2017 Christopher MacMackin <cmacmackin@gmail.com>
!  
!  This program is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published
!  by the Free Software Foundation; either version 3 of the License,
!  or (at your option) any later version.
!  
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!  
!  You should have received a copy of the GNU General Public License
!  License along with this program; if not, write to the Free Software
!  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
!  MA 02110-1301, USA.
!  

module array_pointer_mod
  !* Author: Chris MacMackin
  !  Date: January 2017
  !  License: GPLv3
  !
  ! A module specifying transparent types which hold allocatable
  ! arrays. Field data types contain pointers to these. This allows
  ! the allocatable arrays to be moved and for the memory to be freed
  ! as soon as possible even in pure procedures.
  !
  use iso_fortran_env, only: r8 => real64
  implicit none
  private

  type, public :: array_1d
    !* Author: Chris MacMackin
    !  Date: January 2017
    !
    ! A transparent type containing a 1D allocatable array. Field data
    ! types can contain a pointer to this type. This allows the
    ! allocatable arrays to be moved and for the memory to be freed as
    ! soon as possible even in pure procedures.
    !
    real(r8), dimension(:), allocatable, public :: array
  end type array_1d

  type, public :: array_2d
    !* Author: Chris MacMackin
    !  Date: January 2017
    !
    ! A transparent type containing a 2D allocatable array. Field data
    ! types can contain a pointer to this type. This allows the
    ! allocatable arrays to be moved and for the memory to be freed as
    ! soon as possible even in pure procedures.
    !
    real(r8), dimension(:,:), allocatable, public :: array
  end type array_2d

end module array_pointer_mod
