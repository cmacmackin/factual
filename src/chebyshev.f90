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
  implicit none
  private
  
  real(r8), parameter :: pi = 4.0_r8*atan(1.0_r8)
  integer(c_int), parameter :: field_size_init = -565

  ! Cached values for reuse
  integer(c_int) :: field_size = field_size_init
  real(c_double), dimension(:), pointer :: array1, array2, array3
  type(c_ptr) :: plan_dct, plan_dst, array_c1, array_c2, array_c3

!logical, public :: verbose = .false.

  type :: cached_points
    !* Author: Chris MacMackin
    !  Date: February 2017
    !
    ! A type to form a singly-linked list storing collocation points
    ! which have previously been calculated, as well as information on
    ! the parameters used in the calculation.
    !
    private
    real(r8), dimension(:), pointer :: colloc_points => null()
    integer                         :: nodes
    real(r8)                        :: lower_bound, upper_bound
    type(cached_points), pointer    :: next => null()
  end type cached_points

  type(cached_points), pointer :: cache_list => null()
    !! Contains pointers to sets of collocation points which have been
    !! created by the [[collocation_points]] function.

  public :: collocation_points, differentiation_row, differentiate_1d, &
            integrate_1d

contains
  
  function collocation_points(nodes,lower_bound,upper_bound)
    !* Author: Chris MacMackin
    !  Date: April 2016
    !
    ! Returns a 1D array of Chebyshev collocation nodes. By default they
    ! will be on the domain [-1,1], but with optional arguments this
    ! can be changed.
    !
    integer, intent(in)             :: nodes
      !! The number of collocation nodes to generate
    real(r8), optional, intent(in)  :: lower_bound
      !! The position of the start of the domain. Default is -1.0.
    real(r8), optional, intent(in)  :: upper_bound
      !! The position of the end of the domain. Default is 1.0.
    real(r8), dimension(:), pointer :: collocation_points
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
    allocate(collocation_points(nodes+1))
    factor = (upper - lower)/2.0_r8
    collocation_points = cos([(real(i,r8), i=0,nodes)] * (pi/real(nodes,r8)))
    collocation_points = factor*collocation_points
    collocation_points = (lower + 1.0_r8*factor) + collocation_points
    ! Cache this set of collocation points
    allocate(list_location)
    list_location%colloc_points => collocation_points
    list_location%nodes = nodes
    list_location%lower_bound = lower
    list_location%upper_bound = upper
    list_location%next => cache_list
    cache_list => list_location
  end function collocation_points

  function differentiation_row(nodes,row,lower_bound,upper_bound) result(diff)
    integer, intent(in) :: nodes
      !! The number of collocation nodes
    integer, intent(in) :: row
      !! The point at which the derivative is being evaluated
    real(r8), optional, intent(in) :: lower_bound
      !! The position of the start of the domain. Default is -1.0.
    real(r8), optional, intent(in) :: upper_bound
      !! The position of the end of the domain. Default is 1.0.
    real(r8), dimension(nodes+1) :: diff
    integer  :: j
    real(r8) :: c_i, c_j
    real(r8), dimension(:), pointer :: collocs
    collocs => collocation_points(nodes,lower_bound,upper_bound)
    if (row == 1 .or. row == nodes + 1) then
      c_i = 2._r8
    else
      c_i = 1._r8
    end if
    do concurrent (j=1:nodes+1)
      if (j == row) then
        if (j == 1) then
          diff(j) = (2*nodes**2 + 1.0_r8)/6._r8
        else if (j == nodes + 1) then
          diff(j) = -(2*nodes**2 + 1.0_r8)/6._r8
        else
          diff(j) = -0.5_r8*collocs(j)/(1 - collocs(j)**2)
        end if
      else
        if (j == 1 .or. j == nodes + 1) then
          c_j = 2._r8
        else
          c_j = 1._r8
        end if
        diff(j) = c_i*(-1)**(row+j)/ &
                  (c_j*(collocs(row) - collocs(j)))
      end if
    end do
  end function differentiation_row

  subroutine setup_fftw3(field_size_new)
    !* Author: Chris MacMackin
    !  Date: September 2017
    !
    ! Allocates memory and computes plans for running DCT and DST with
    ! FFTW3. This is only recalculated if the size of arrays has
    ! changed, providing considerable savings.
    !
    integer(c_int), intent(in) :: field_size_new
      !! The size of the field to be integrated/differentiated.
    if (field_size_new /= field_size) then
      if (field_size /= field_size_init) then
        call fftw_free(array_c1)
        call fftw_free(array_c2)
        call fftw_free(array_c3)
        call fftw_destroy_plan(plan_dct)
        call fftw_destroy_plan(plan_dst)
      end if
      field_size = field_size_new
      array_c1 = fftw_alloc_real(int(field_size,c_size_t))
      call c_f_pointer(array_c1, array1, [int(field_size,c_size_t)])
      array_c2 = fftw_alloc_real(int(field_size,c_size_t))
      call c_f_pointer(array_c2, array2, [int(field_size,c_size_t)])
      array_c3 = fftw_alloc_real(int(field_size,c_size_t))
      call c_f_pointer(array_c3, array3, [int(field_size,c_size_t)])
      plan_dct = fftw_plan_r2r_1d(field_size, array1, array2, FFTW_REDFT00, &
                                  FFTW_PATIENT)
      plan_dst = fftw_plan_r2r_1d(field_size-2, array1, array3, FFTW_RODFT00, &
                                  FFTW_PATIENT)
    end if
  end subroutine setup_fftw3

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
    real(r8) :: inverse_width, field_centre, array1_lower
    if (present(order)) then
      order_ = order
    else
      order_ = 1
    end if
    ord = order_
    if (ord <= 0) return
    call setup_fftw3(size(field_data,kind=c_int))
    inverse_width = 2.0_c_double/(xvals(1) - xvals(field_size))
    field_centre = xvals(1) - 1.0_r8/(inverse_width)
    array1 = field_data
    do while(ord > 0)
      call fftw_execute_r2r(plan_dct,array1,array2)
      do concurrent (i=1:field_size-1)
        array2(i+1) = array2(i+1) / (field_size-1)
        array1(i) = -i*pi*array2(i+1)
        array2(i+1) = i**2 * array2(i+1)
      end do
      array1(field_size - 1) = array1(field_size - 1)*0.5_r8
      array1(field_size) = 0.0_c_double
      call fftw_execute_r2r(plan_dst,array1,array3)
      array1(1) = sum(array2(2:field_size-1)) + 0.5_c_double*array2(field_size)
 !     if (verbose) print*,sum(array2(2:field_size-1)),0.5_r8*array2(field_size),array1(1)
      do concurrent (i=2:field_size-1)
        array1(i) = -0.5_c_double/pi * array3(i-1)/ &
             sqrt(1.0_c_double - ((xvals(i)-field_centre)*inverse_width)**2)
        array2(i) = (-1)**(i) * array2(i)
      end do
      array2(field_size) = (-1)**(field_size) * array2(field_size)
      array1(field_size) = sum(array2(2:field_size-1)) + 0.5_c_double*array2(field_size)
!      if (verbose) print*,sum(array2(2:field_size-1)),0.5_r8*array2(field_size),array1(field_size)
      ord = ord - 1
    end do
    field_data = array1 * inverse_width**order_
  end subroutine differentiate_1d


  subroutine integrate_1d(field_data,xvals,boundary_point,boundary_val)
    !* Author: Chris MacMackin
    !  Date: September 2017
    !
    ! Computes the integral (in place) of 1D data using a Chebyshev
    ! pseudo-spectral methods and a type-I discrete cosine transform,
    ! provided by the FFTW3 library. This is fastest for arrays of odd
    ! size.
    !
    real(c_double), dimension(:), intent(inout) :: field_data
      !! The data which is to be integrated, along the direction of
      !! the array.
    real(c_double), dimension(size(field_data)), intent(in) :: xvals
      !! The location of each data point in `field_data`.
    integer, intent(in), optional :: boundary_point
      !! The array coordinate at which the value of the result is
      !! known. This is used in order to calculate the constant to be
      !! added to the integral. If not provided then the result will
      !! calculated such that the first Chebyshev mode (which is a
      !! constant) will have value 0. If the input is outside of the
      !! array bounds then it is ignored.
    real(c_double), intent(in), optional :: boundary_val
      !! The value which the boundary point should hold in the
      !! output. If `boundary_point` is not provided or outside of
      !! array bounds then this argument is ignored. If
      !! `boundary_point` is provided but this argument is note then
      !! it defaults to 0.
    integer :: i
    real(r8) :: width, field_centre, bval
    logical :: use_first_eq
    call setup_fftw3(size(field_data,kind=c_int))
    width = 0.5_c_double*(xvals(1) - xvals(field_size))
    field_centre = xvals(1) - width
    array1 = field_data
    do concurrent (i=2:field_size-1)
      array1(i-1) = -pi*width/real(field_size-1, c_double)*field_data(i)* &
           sqrt(1.0_c_double - ((xvals(i)-field_centre)/width)**2)
    end do
    call fftw_execute_r2r(plan_dst,array1,array3)
    do concurrent (i=2:field_size-1)
      array1(i) = -0.5_r8/(pi*(i-1)) * array3(i-1)
      array3(i-1) = (i-1)**2*array1(i)
    end do
    array1(1) = 0._r8
    ! For noisy input I find that the two equations I can use to find
    ! the highest order Chebyshev mode return different results. One
    ! of these comes from the equation for the derivative at the first
    ! collocation node and the other for the derivative at the
    ! last. This causes an errors when differentiating the result of
    ! this routine; the result of the differentiation at the end point
    ! whose equation is not used here is wrong. As a workaround for my
    ! purposes, I use whichever equation does not correspond to
    ! `boundary_point`, as in ISOFT the boundary point gets
    ! overwritten anyway.
    if (present(boundary_point)) then
      use_first_eq = boundary_point == field_size
    else
      use_first_eq = .true.
    end if
    if (use_first_eq) then
      array1(field_size) = 1._r8/(field_size-1)**2 * &
           (width*field_data(1) - 2*sum(array3(1:field_size-2)))
!    if (verbose) then
!      print*,field_data(1), field_data(field_size)
!      print*, width*(field_data(1)-field_data(field_size)), &
!              width*(field_data(1)+field_data(field_size)), &
!              4*sum(array3(1:field_size-2:2))
!    end if
    else
      array3(2:field_size:2) = -array3(2:field_size:2)
      array1(field_size) = real((-1)**field_size, r8)/(field_size-1)**2 * &
           (width*field_data(field_size) - 2_r8*sum(array3(1:field_size-2)))
    end if
    call fftw_execute_r2r(plan_dct,array1,array2)
    if (present(boundary_point)) then
      if (boundary_point < 1 .or. boundary_point > field_size) then
        field_data = array2
        return
      end if
      if (present(boundary_val)) then
        bval = boundary_val
      else
        bval = 0._c_double
      end if
      field_data = array2 + (bval - array2(boundary_point))
    else
      field_data = array2
    end if
  end subroutine integrate_1d

end module chebyshev_mod
