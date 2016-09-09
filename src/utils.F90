!
!  utils.f90
!  This file is part of FACTUAL.
!
!  Copyright 2016 Christopher MacMackin <cmacmackin@gmail.com>
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

module utils_mod
  !* Author: Chris MacMackin
  !  Date: September 2016
  !  License: GPLv3
  !
  ! Provides useful helper functions for internal use with FACTUAL.
  !
  use iso_fortran_env, only: r8 => real64, stderr => error_unit
  use abstract_fields_mod, only: abstract_field
  implicit none
  private

  public :: is_nan, check_set_from_raw

contains

  elemental logical function is_nan(val)
    !* Author: Chris MacMackin
    !  Date: September 2016
    !
    ! Returns true if the argument is NaN, false otherwise. If using `gfortran`
    ! then the `isnan` extension is used. Otherwise, the value is tested for
    ! equality against itself.
    !
    real(r8), intent(in) :: val !! The value being checked to see if it is NaN
#ifdef __GFORTRAN__
    is_nan = isnan(val)
#else
    is_nan = val /= val
#endif
  end function is_nan

  subroutine check_set_from_raw(field,raw,provide_lower_bound,provide_upper_bound)
    !* Author: Chris MacMackin
    !  Date: September 2016
    !
    ! Checks to ensure that a raw array is the correct size to set the contents
    ! of a field object. If it is not then it will print an error message and
    ! stop the program.
    !
    class(abstract_field), intent(in) :: field
      !! The field being checked to see if it is compatible with the raw array.
    real(r8), dimension(:), intent(in) :: raw
      !! The raw data to be stored in this array.
    logical, dimension(:), optional, intent(in) :: provide_lower_bound
      !! Specifies whether raw data contains values at the lower
      !! boundary, for each dimension, with the index of the element
      !! corresponding to the dimension. Defaults to all `.true.`.
    logical, dimension(:), optional, intent(in) :: provide_upper_bound
      !! Specifies whether raw data contains values at the upper
      !! boundary, for each dimension, with the index of the element
      !! corresponding to the dimension. Defaults to all `.true.`.    
#ifdef DEBUG
    integer :: expected
    expected = field%raw_size(provide_lower_bound,provide_upper_bound)
    if (expected /= size(raw)) then
      write(stderr,*) 'array_vector_field: Error, setting from raw array of wrong size'
      write(stderr,"('    Needed: ',i0,', Actual: ',i0)") expected, size(raw)
#ifdef __GFORTRAN__
      call backtrace
#endif
      error stop
    end if
#endif
  end subroutine check_set_from_raw

end module utils_mod
