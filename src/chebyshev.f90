!
!  chebyshev.f90
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

module chebyshev_mod
  !* Author: Chris MacMackin
  !  Date: April 2016
  !  License: GPLv3
  !
  ! A module providing utility functions for Chebyshev pseudo-spectral
  ! differentiation.
  !
  use iso_fortran_env, only: r8 => real64
  use iso_c_binding
  use fftw3_mod
  use array_pointer_mod, only: array_1d
  implicit none
  private
  
  real(r8), parameter :: pi = 4.0_r8*atan(1.0_r8)
  
  type :: cached_points
    !* Author: Chris MacMackin
    !  Date: February 2017
    !
    ! A type to form a singly-linked list storing collocation points
    ! which have previously been calculated, as well as information on
    ! the parameters used in the calculation.
    !
    private
    type(array_1d), pointer      :: colloc_points => null()
    integer                      :: nodes
    real(r8)                     :: lower_bound, upper_bound
    type(cached_points), pointer :: next => null()
  end type cached_points

  type(cached_points), pointer :: cache_list => null()
    !! Contains pointers to sets of collocation points which have been
    !! created by the [[collocation_points]] function.

  public :: collocation_points, differentiate_1d

contains
  
  function collocation_points(nodes,lower_bound,upper_bound)
    !* Author: Chris MacMackin
    !  Date: April 2016
    !
    ! Returns a 1D array of Chebyshev collocation nodes. By default they
    ! will be on the domain [-1,1], but with optional arguments this
    ! can be changed.
    !
    integer, intent(in) :: nodes
      !! The number of collocation nodes to generate
    real(r8), optional, intent(in) :: lower_bound
      !! The position of the start of the domain. Default is -1.0.
    real(r8), optional, intent(in) :: upper_bound
      !! The position of the end of the domain. Default is 1.0.
    type(array_1d), pointer :: collocation_points
    integer :: i
    real(r8) :: upper, lower, factor
    type(cached_points), pointer :: list_location
    upper = 1.0_r8
    lower = -1.0_r8
    if (present(upper_bound)) upper = upper_bound
    if (present(lower_bound)) lower = lower_bound
    ! Check if collocation points have been calculated for these parameters before.
    list_location => cache_list
    do while (associated(list_location))
      !TODO: Convert this to relative difference check, rather than equality check.
      if ((list_location%nodes == nodes) .and. &
          (list_location%lower_bound == lower) .and. &
          (list_location%upper_bound == upper)) then
        collocation_points => list_location%colloc_points
        return
      else
        list_location => list_location%next
      end if
    end do
    ! If these collocation points have not already been calculated, do so now.
    allocate(collocation_points)
    allocate(collocation_points%array(nodes+1))
    factor = (upper - lower)/2.0_r8
    collocation_points%array = cos([(real(i,r8), i=0,nodes)] * (pi/real(nodes,r8)))
    collocation_points%array = factor*collocation_points%array
    collocation_points%array = (lower + 1.0_r8*factor) + collocation_points%array
    ! Cache this set of collocation points
    allocate(list_location)
    list_location%colloc_points => collocation_points
    list_location%nodes = nodes
    list_location%lower_bound = lower
    list_location%upper_bound = upper
    list_location%next => cache_list
    cache_list => list_location
  end function collocation_points

  subroutine differentiate_1d(field_data,xvals,order)
    !* Author: Chris MacMackin
    !  Date: April 2016
    !
    ! Computes the derivative (in place) of 1D data using a Chebyshev
    ! pseudo-spectral methods and a type-I discrete cosine transform,
    ! provided by the FFTW3 library. This is fastest for arrays of odd
    ! size.
    !
    ! @Note This routine works recursively for all derivatives with
    ! order greater than one. This is much simpler to implement than,
    ! but only about half as efficient as, directly computing any
    ! derivative. It also causes error to accumulate as higher
    ! derivatives are computed, which can not be avoided simply by using
    ! more nodes. This routine has been found to be somewhat unreliable
    ! for derivatives of order greater than three.
    !
    ! @Todo Add the ability to directly calculate second derivative.
    !
    real(c_double), dimension(:), intent(inout) :: field_data
      !! The data which is to be differentiated, along the direction of
      !! the array.
    real(c_double), dimension(size(field_data)), intent(in) :: xvals
      !! The location of each data point in `field_data`.
    integer, intent(in), optional :: order
      !! The order of the derivative to take. If absent, defaults to 1.
      !! Zero corresponds to no derivative and negative values are
      !! treated as zero.
    integer :: ord, order_, i
    integer(c_int) :: field_size
    real(c_double), dimension(:), pointer :: array1, array2
    real(r8) :: inverse_width, field_centre
    type(c_ptr) :: plan_dct, plan_dst, array_c1, array_c2
    field_size = size(field_data,kind=c_size_t)
    if (present(order)) then
      order_ = order
    else
      order_ = 1
    end if
    ord = order_
    inverse_width = 2.0_c_double/(xvals(1) - xvals(field_size))
    field_centre = xvals(1) - 1.0_r8/(inverse_width)
    if (ord > 0) then
      array_c1 = fftw_alloc_real(int(field_size,c_size_t))
      call c_f_pointer(array_c1, array1, [int(field_size,c_size_t)])
      array_c2 = fftw_alloc_real(int(field_size,c_size_t))
      call c_f_pointer(array_c2, array2, [int(field_size,c_size_t)])
    else
      return
    end if
    plan_dct = fftw_plan_r2r_1d(field_size, array1, array2, FFTW_REDFT00, FFTW_ESTIMATE)
    plan_dst = fftw_plan_r2r_1d(field_size-2, array1, array1, FFTW_RODFT00, FFTW_ESTIMATE)
    array1 = field_data
    do while(ord > 0)
      call fftw_execute_r2r(plan_dct,array1,array2)
      array2 = pi/real(field_size-1,c_double) * array2
      forall (i=1:field_size-1) array1(i) = -real(i,c_double)*array2(i+1)
      array1(field_size - 1) = array1(field_size - 1)*0.5_r8
      array1(field_size) = 0.0_c_double
      call fftw_execute_r2r(plan_dst,array1,array1)
      array1 = 0.5_c_double/pi * array1
      forall (i=2:field_size-1) array1(i) = -array1(i-1)/ &
          sqrt(1.0_c_double - ((xvals(i)-field_centre)*inverse_width)**2)
      array2 = 1.0_c_double/pi * array2
      forall (i=2:field_size) array2(i) = real(i-1,c_double)**2 * array2(i)
      array1(1) = sum(array2(2:field_size-1))
      array1(1) = array1(1) + 0.5_c_double*array2(field_size)
      forall (i=2:field_size) array2(i) = (-1._c_double)**(i) * array2(i)
      array1(field_size) = sum(array2(2:field_size-1))
      array1(field_size) = array1(field_size) + 0.5_c_double*array2(field_size)
      ord = ord - 1
    end do
    field_data = array1 * inverse_width**order_
    call fftw_free(array_c1)
    call fftw_free(array_c2)
    call fftw_destroy_plan(plan_dct)
    call fftw_destroy_plan(plan_dst)
  end subroutine differentiate_1d

end module chebyshev_mod
