!
!  cheb2d_fields.f90
!  This file is part of FACTUAL.
!
!  Copyright 2016 Christopher MacMackin <cmacmackin@gmail.com>
!  
!  This program is free software; you can redistribute it and/or modify
!  it under the terms of the GNU Lesser General Public License as
!  published by the Free Software Foundation; either version 3 of the 
!  License, or (at your option) any later version.
!  
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU Lesser General Public License for more details.
!  
!  You should have received a copy of the GNU Lesser General Public
!  License along with this program; if not, write to the Free Software
!  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
!  MA 02110-1301, USA.
!  

module cheb1d_fields_mod
  !* Author: Chris MacMackin
  !  Date: March 2016
  !  License: LGPLv3
  !
  ! Provides concrete implementations of the [[scalar_field]] and
  ! [[vector_field]] types. This implementation is for 1D fields (with
  ! 3D vectors) only. It tracks the values of the fields at Chebyshev
  ! collocation nodes and uses a pseudospectral approach to 
  ! differentiate.
  !
  use iso_fortran_env, only: r8 => real64, stderr => error_unit
  use abstract_fields_mod
  use chebyshev_mod
  implicit none
  private

  public :: sin, cos, tan, asin, acos, atan, sinh, cosh, tanh, asinh, &
            acosh, atanh, log, log10, exp, abs, sqrt, minval, maxval

  interface check_compatible
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! Checks that the other field can be used in an operation with
    ! this one, stopping with an error message if not.
    !
    module procedure check_compatible_scalar
    module procedure check_compatible_vector
  end interface
  
  type, extends(scalar_field), public :: cheb1d_scalar_field
    !* Author: Chris MacMackin
    !  Date: March 2016
    !  License: LGPLv3
    !
    ! A concrete implementations of the [[scalar_field]] for 1D fields
    ! (with 3D vectors) only. It tracks the values of the fields at 
    ! Chebyshev collocation nodes and uses a pseudospectral approach to 
    ! differentiate. In addition to the type-bound procedures, all of
    ! the intrinsic mathematical procedures, with the exception of 
    ! Bessel functions and error functions.
    !
    private
    integer :: numpoints
      !! The number of datapoints used
    real(r8), dimension(1,2) :: extent = reshape([-1.0_r8, 1.0_r8],[1,2])
      !! The start and end values of the domain
    real(r8), dimension(:), allocatable :: colloc_points
      !! The location of the data-points
    real(r8), dimension(:), allocatable :: field_data
      !! The value of the scalar field at the data-points
  contains
    private
    procedure, public :: domain => cheb1d_scalar_domain
      !! Provides array with upper and lower limits of field's domain
    procedure, public :: dimensions => cheb1d_scalar_dimensions
      !! Returns number of dimensions of the field
    procedure, public :: raw_size => cheb1d_scalar_raw_size
      !! Provides the number of pieces of data needed to represent the
      !! field, i.e. the size of the array returned by `get_raw`.
    procedure, public :: raw => cheb1d_scalar_raw
      !! Returns array of data representing state of field. Can be
      !! useful for passing to nonlinear solvers.
    procedure, public :: resolution => cheb1d_scalar_resolution
      !! Returns array containing number of datapoints in each dimension.
    procedure, public :: set_from_raw => cheb1d_scalar_set_from_raw
      !! Assigns raw data, such as that produced by 
      !! [[cheb1d_scalar_field:raw]], to the field
    procedure :: field_multiply_field => cheb1d_scalar_sf_m_sf
      !! \({\rm field} \times {\rm field}\)
    procedure :: field_multiply_vecfield => cheb1d_scalar_sf_m_vf
      !! \({\rm field} \times {\rm \vec{field}}\)
    procedure, pass(rhs) :: real_multiply_field => cheb1d_scalar_r_m_sf
      !! \({\rm real}  \times {\rm field}\)
    procedure :: field_multiply_real => cheb1d_scalar_sf_m_r
      !! \({\rm field} \times {\rm real}\)
    procedure :: field_divide_field => cheb1d_scalar_sf_d_sf
      !! \(\frac{\rm field}{\rm field}\)
    procedure, pass(rhs) :: real_divide_field => cheb1d_scalar_r_d_sf
      !! \(\frac{\rm real}{\rm field}\)
    procedure :: field_divide_real => cheb1d_scalar_sf_d_r
      !! \(\frac{\rm field}{\rm real}\)
    procedure :: field_add_field => cheb1d_scalar_sf_a_sf
      !! \({\rm field} + {\rm field}\)
    procedure, pass(rhs) :: real_add_field => cheb1d_scalar_r_a_sf
      !! \({\rm real} + {\rm field}\)
    procedure :: field_add_real => cheb1d_scalar_sf_a_r
      !! \({\rm field} + {\rm real}\)
    procedure :: field_sub_field => cheb1d_scalar_sf_s_sf
      !! \({\rm field} - {\rm field}\)
    procedure, pass(rhs) :: real_sub_field => cheb1d_scalar_r_s_sf
      !! \({\rm real} - {\rm field}\)
    procedure :: field_sub_real => cheb1d_scalar_sf_s_r
      !! \({\rm field} - {\rm real}\)
    procedure :: field_pow_real => cheb1d_scalar_sf_p_r
      !! \({\rm field}^{\rm real}\)
    procedure :: field_pow_real4 => cheb1d_scalar_sf_p_r4
      !! \({\rm field}^{\rm real}\)
    procedure :: field_pow_int => cheb1d_scalar_sf_p_i
      !! \({\rm field}^{\rm int}\)
    procedure :: sin => cheb1d_scalar_sin
      !! \(\sin({\rm field})\)
    procedure :: cos => cheb1d_scalar_cos
      !! \(\cos({\rm field})\)
    procedure :: tan => cheb1d_scalar_tan
      !! \(\tan({\rm field})\)
    procedure :: asin => cheb1d_scalar_asin
      !! \(\sin^{-1}({\rm field})\)
    procedure :: acos => cheb1d_scalar_acos
      !! \(\cos^{-1}({\rm field})\)
    procedure :: atan => cheb1d_scalar_atan
      !! \(\tan^{-1}({\rm field})\)
    procedure :: sinh => cheb1d_scalar_sinh
      !! \(\sinh({\rm field})\)
    procedure :: cosh => cheb1d_scalar_cosh
      !! \(\cosh({\rm field})\)
    procedure :: tanh => cheb1d_scalar_tanh
      !! \(\tanh({\rm field})\)
    procedure :: asinh => cheb1d_scalar_asinh
      !! \(\sinh^{-1}({\rm field})\)
    procedure :: acosh => cheb1d_scalar_acosh
      !! \(\cosh^{-1}({\rm field})\)
    procedure :: atanh => cheb1d_scalar_atanh
      !! \(\tanh^{-1}({\rm field})\)
    procedure :: log => cheb1d_scalar_log
      !! \(\ln ({\rm field})\)
    procedure :: log10 => cheb1d_scalar_log10
      !! \(\log ({\rm field})\)
    procedure :: exp => cheb1d_scalar_exp
      !! \(\e^{\rm field})\)
    procedure :: abs => cheb1d_scalar_abs
      !! \(\|{\rm field}|\)
    procedure :: sqrt => cheb1d_scalar_sqrt
      !! \(\sqrt{\rm field}\)
    procedure :: minval => cheb1d_scalar_minval
      !! \(\min({\rm field})\)
    procedure :: maxval => cheb1d_scalar_maxval
      !! \(\max({\rm field})\)
    procedure, public :: d_dx => cheb1d_scalar_d_dx
      !! \(\frac{\partial^n}{\partial x_i^n}({\rm field})\)
    procedure :: gradient => cheb1d_scalar_gradient
      !! \(\nabla {\rm field}\)
    procedure :: laplacian => cheb1d_scalar_laplacian
      !! \(\nabla^2 {\rm field}\)
    procedure :: assign_field => cheb1d_scalar_assign
      !! \({\rm field} = {\rm field}\)
    procedure, public :: assign_meta_data => cheb1d_assign_meta_data
      !! Copies all data other than values stored in field from another
      !! field object to this one.
  end type
  
  interface cheb1d_scalar_field
    module procedure scalar_constructor
  end interface

contains

  function scalar_constructor(nodes,initializer,lower_bound, &
                              upper_bound) result(field)
    integer, intent(in) :: nodes
      !! The number of collocation nodes to use when modelling this
      !! field. This corresponds to resolution.
    procedure(scalar_init), optional :: initializer
      !! An elemental procedure taking which takes the position in the
      !! fields domain (an 8-byte real) as an argument and returns the
      !! fields value at that position. Default is for field to be zero
      !! everywhere.
    real(r8), intent(in), optional :: lower_bound
      !! The position of the lower bound of the field's domain. Default is -1.0.
    real(r8), intent(in), optional :: upper_bound
      !! The position of the upper bound of the field's domain. Default is 1.0.
    type(cheb1d_scalar_field) :: field
    integer :: i
    field%numpoints = nodes
    allocate(field%field_data(nodes+1))
    if (present(lower_bound)) field%extent(1,1) = lower_bound
    if (present(upper_bound)) field%extent(1,2) = upper_bound
    field%colloc_points = collocation_points(nodes,lower_bound,upper_bound)
    if (present(initializer)) then
      forall (i=1:nodes+1) field%field_data(i) = initializer(field%colloc_points(i))
    else
      field%field_data = 0.0_r8
    end if
  end function scalar_constructor

  pure function scalar_init(x) result(scalar)
    !* Author: Chris MacMackin
    !  Date: April 2016
    !
    ! Default function used to specify value held by a field at location
    ! `x`. It sets the field to be uniform zeros.
    !
    real(r8), intent(in) :: x 
      !! The position at which this function is evaluated
    real(r8) :: scalar !! The value of the field at this location
    scalar = 0.0_r8
  end function scalar_init

  pure function cheb1d_scalar_domain(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! Accessor for the upper and lower bounds of the domain.
    !
    class(cheb1d_scalar_field), intent(in) :: this
    real(r8), dimension(:,:), allocatable :: res
      !* A 2D array of shape \(n \times 2\), where n is the number of
      !  dimensions of the field. Each row contains the lower and the
      !  upper extent of the fields domain in the corresponding 
      !  dimension
    res = this%extent
  end function cheb1d_scalar_domain

  elemental function cheb1d_scalar_dimensions(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! Accessor for the number of dimensions.
    !
    class(cheb1d_scalar_field), intent(in) :: this
    integer :: res
    res = 1
  end function cheb1d_scalar_dimensions

  pure function cheb1d_scalar_raw_size(this,return_lower_bound, &
                                       return_upper_bound) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! Compute how many elements are in the raw representation of this
    ! field. This would be the number of data points, adjusted based on
    ! how boundary conditions are accounted for.
    !
    class(cheb1d_scalar_field), intent(in) :: this
    logical, dimension(:), optional, intent(in) :: return_lower_bound
      !! Specifies whether to return the values at the lower boundary
      !! for each dimension, with the index of the element
      !! corresponding to the dimension. Defaults to all `.true.`.
    logical, dimension(:), optional, intent(in) :: return_upper_bound
      !! Specifies whether to return the values at the upper boundary
      !! for each dimension, with the index of the element
      !! corresponding to the dimension. Defaults to all `.true.`.
    integer :: res
    res = this%numpoints + 1
    if (present(return_lower_bound)) then
      if (.not. return_lower_bound(1)) res = res - 1
    end if
    if (present(return_upper_bound)) then
      if (.not. return_upper_bound(1)) res = res - 1
    end if
  end function cheb1d_scalar_raw_size
  
  pure function cheb1d_scalar_raw(this,return_lower_bound, &
                                  return_upper_bound) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! Returns a representation of this field in the form of a 1D array
    ! of real numbers this allows, e.g. for manipulation by solver
    ! routines.
    !
    ! @BUG The returned value has length `this%raw_size()`, but
    ! a bug in gfortran 4.8 (fixed by version 5) caused the compiler
    ! to segfault if it was declared as such. As a workaround, it is
    ! allocatable isntead.
    !
    class(cheb1d_scalar_field), intent(in) :: this
    real(r8), dimension(:), allocatable :: res
      !! Array containing data needed to describe field
    logical, dimension(:), optional, intent(in) :: return_lower_bound
      !! Specifies whether to return the values at the lower boundary
      !! for each dimension, with the index of the element
      !! corresponding to the dimension. Defaults to all `.true.`.
    logical, dimension(:), optional, intent(in) :: return_upper_bound
      !! Specifies whether to return the values at the upper boundary
      !! for each dimension, with the index of the element
      !! corresponding to the dimension. Defaults to all `.true.`.
    integer :: lower, upper
    lower = 1
    upper = this%numpoints + 1
    ! Recall that Chebyshev points are in reverse order...
    if (present(return_lower_bound)) then
      if (.not. return_lower_bound(1)) upper = upper - 1
    end if
    if (present(return_upper_bound)) then
      if (.not. return_upper_bound(1)) lower = lower + 1
    end if
    if (allocated(this%field_data)) then
      res = this%field_data(lower:upper)
    else
      res = [0.0_r8]
    end if
  end function cheb1d_scalar_raw
  
  pure function cheb1d_scalar_resolution(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! Get the number of data points used to represent the field in
    ! its one dimension.
    !
    class(cheb1d_scalar_field), intent(in) :: this
    integer, dimension(:), allocatable :: res
      !! Array of length 1 specifying the number of data points.
    allocate(res(1))
    res(1) = this%numpoints
  end function cheb1d_scalar_resolution

  subroutine cheb1d_scalar_set_from_raw(this,raw,provide_lower_bound, &
                                        provide_upper_bound)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! Assigns raw data, such as that produced by 
    ! [[cheb1d_scalar_field:raw]], to the field. The routine will
    ! stop with an error if the array is the wrong length for a field
    ! of this resolution and with these boundary conditions. 
    !
    class(cheb1d_scalar_field), intent(inout) :: this
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
    integer :: lower, upper
#ifdef DEBUG
    integer :: expected
    expected = this%raw_size(provide_lower_bound,provide_upper_bound)
    if (expected /= size(raw)) then
      write(stderr,*) 'cheb1d_scalar_field: Error, setting from raw array of wrong size'
      write(stderr,"('    Needed: ',i,', Actual: ',i)") expected, size(raw)
#ifdef __GFORTRAN__
      call backtrace
#endif
      error stop
    end if
#endif
    lower = 1
    upper = this%numpoints + 1
    ! Recall that Chebyshev points are in reverse order...
    if (present(provide_lower_bound)) then
      if (.not. provide_lower_bound(1)) upper = upper - 1
    end if
    if (present(provide_upper_bound)) then
      if (.not. provide_upper_bound(1)) lower = lower + 1
    end if
    if (.not. allocated(this%field_data)) &
      allocate(this%field_data(this%numpoints + 1))
    this%field_data(lower:upper) = raw
  end subroutine cheb1d_scalar_set_from_raw
  
  function cheb1d_scalar_sf_m_sf(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm field} \times {\rm field}\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), intent(in) :: rhs
    class(scalar_field), allocatable :: res !! The restult of this operation
    type(cheb1d_scalar_field), allocatable :: local
#ifdef DEBUG
    call check_compatible(this, rhs)
#endif
    allocate(local)
    call local%assign_meta_data(this)
    select type(rhs)
    class is(cheb1d_scalar_field)
      allocate(local%field_data(size(this%field_data))) !Needed due to compiler bug
      local%field_data = this%field_data * rhs%field_data
    end select
    call move_alloc(local,res)
  end function cheb1d_scalar_sf_m_sf

  function cheb1d_scalar_sf_m_vf(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm field} \times {\rm \vec{field}}\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(vector_field), intent(in) :: rhs
    class(vector_field), allocatable :: res !! The result of this operation
  end function cheb1d_scalar_sf_m_vf

  function cheb1d_scalar_r_m_sf(lhs,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm real} \times {\rm field}\)
    !
    real(r8), intent(in) :: lhs
    class(cheb1d_scalar_field), intent(in) :: rhs
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: local
    allocate(local)
    allocate(local%field_data(size(rhs%field_data))) !Needed due to compiler bug
    local%field_data = lhs * rhs%field_data
    call move_alloc(local,res)
  end function cheb1d_scalar_r_m_sf

  function cheb1d_scalar_sf_m_r(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm field} \times {\rm real}\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    real(r8), intent(in) :: rhs
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: local
    allocate(local)
    allocate(local%field_data(size(this%field_data))) !Needed due to compiler bug
    local%field_data = this%field_data * rhs
    call move_alloc(local,res)
  end function cheb1d_scalar_sf_m_r
  
  function cheb1d_scalar_sf_d_sf(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm field} / {\rm field}\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), intent(in) :: rhs
    class(scalar_field), allocatable :: res !! The restult of this operation
    type(cheb1d_scalar_field), allocatable :: local
#ifdef DEBUG
    call check_compatible(this, rhs)
#endif
    allocate(local)
    call local%assign_meta_data(this)
    select type(rhs)
    class is(cheb1d_scalar_field)
      allocate(local%field_data(size(this%field_data))) !Needed due to compiler bug
      local%field_data = this%field_data / rhs%field_data
    end select
    call move_alloc(local,res)
  end function cheb1d_scalar_sf_d_sf

  function cheb1d_scalar_r_d_sf(lhs,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm real} / {\rm field}\)
    !
    real(r8), intent(in) :: lhs
    class(cheb1d_scalar_field), intent(in) :: rhs
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: local
    allocate(local)
    allocate(local%field_data(size(rhs%field_data))) !Needed due to compiler bug
    local%field_data = lhs / rhs%field_data
    call move_alloc(local,res)
  end function cheb1d_scalar_r_d_sf

  function cheb1d_scalar_sf_d_r(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm field} / {\rm real}\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    real(r8), intent(in) :: rhs
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: local
    allocate(local)
    allocate(local%field_data(size(this%field_data))) !Needed due to compiler bug
    local%field_data = this%field_data / rhs
    call move_alloc(local,res)
  end function cheb1d_scalar_sf_d_r
  
  function cheb1d_scalar_sf_s_sf(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm field} - {\rm field}\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), intent(in) :: rhs
    class(scalar_field), allocatable :: res !! The restult of this operation
    type(cheb1d_scalar_field), allocatable :: local
#ifdef DEBUG
    call check_compatible(this, rhs)
#endif
    allocate(local)
    call local%assign_meta_data(this)
    select type(rhs)
    class is(cheb1d_scalar_field)
      allocate(local%field_data(size(this%field_data))) !Needed due to compiler bug
      local%field_data = this%field_data - rhs%field_data
    end select
    call move_alloc(local,res)
  end function cheb1d_scalar_sf_s_sf

  function cheb1d_scalar_r_s_sf(lhs,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm real} - {\rm field}\)
    !
    real(r8), intent(in) :: lhs
    class(cheb1d_scalar_field), intent(in) :: rhs
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: local
    allocate(local)
    allocate(local%field_data(size(rhs%field_data))) !Needed due to compiler bug
    local%field_data = lhs - rhs%field_data
    call move_alloc(local,res)
  end function cheb1d_scalar_r_s_sf

  function cheb1d_scalar_sf_s_r(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm field} - {\rm real}\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    real(r8), intent(in) :: rhs
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: local
    allocate(local)
    allocate(local%field_data(size(this%field_data))) !Needed due to compiler bug
    local%field_data = this%field_data - rhs
    call move_alloc(local,res)
  end function cheb1d_scalar_sf_s_r
  
  function cheb1d_scalar_sf_a_sf(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm field} - {\rm field}\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), intent(in) :: rhs
    class(scalar_field), allocatable :: res !! The restult of this operation
    type(cheb1d_scalar_field), allocatable :: local
#ifdef DEBUG
    call check_compatible(this, rhs)
#endif
    allocate(local)
    call local%assign_meta_data(this)
    select type(rhs)
    class is(cheb1d_scalar_field)
      allocate(local%field_data(size(this%field_data))) !Needed due to compiler bug
      local%field_data = this%field_data + rhs%field_data
    end select
    call move_alloc(local,res)
  end function cheb1d_scalar_sf_a_sf

  function cheb1d_scalar_r_a_sf(lhs,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm real} - {\rm field}\)
    !
    real(r8), intent(in) :: lhs
    class(cheb1d_scalar_field), intent(in) :: rhs
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: local
    allocate(local)
    allocate(local%field_data(size(rhs%field_data))) !Needed due to compiler bug
    local%field_data = lhs + rhs%field_data
    call move_alloc(local,res)
  end function cheb1d_scalar_r_a_sf

  function cheb1d_scalar_sf_a_r(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm field} - {\rm real}\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    real(r8), intent(in) :: rhs
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: local
    allocate(local)
    allocate(local%field_data(size(this%field_data))) !Needed due to compiler bug
    local%field_data = this%field_data + rhs
    call move_alloc(local,res)
  end function cheb1d_scalar_sf_a_r

  function cheb1d_scalar_sf_p_r(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm field}^{\rm real}\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    real(r8), intent(in) :: rhs
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: local
    allocate(local)
    allocate(local%field_data(size(this%field_data))) !Needed due to compiler bug
    local%field_data = this%field_data ** rhs
    call move_alloc(local,res)
  end function cheb1d_scalar_sf_p_r

  function cheb1d_scalar_sf_p_r4(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm field}^{\rm real}\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    real, intent(in) :: rhs
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: local
    allocate(local)
    allocate(local%field_data(size(this%field_data))) !Needed due to compiler bug
    local%field_data = this%field_data * rhs
    call move_alloc(local,res)
  end function cheb1d_scalar_sf_p_r4

  function cheb1d_scalar_sf_p_i(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm field}^{\rm integer}\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    integer, intent(in) :: rhs
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: local
    allocate(local)
    allocate(local%field_data(size(this%field_data))) !Needed due to compiler bug
    local%field_data = this%field_data * rhs
    call move_alloc(local,res)
  end function cheb1d_scalar_sf_p_i

  pure function cheb1d_scalar_sin(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\sin({\rm field})\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: tmp
    allocate(tmp)
    if (allocated(this%field_data)) then
      allocate(tmp%field_data(size(this%field_data)))
      tmp%field_data = sin(this%field_data)
    end if
    call move_alloc(tmp, res)
  end function cheb1d_scalar_sin

  pure function cheb1d_scalar_cos(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\cos({\rm field})\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: tmp
    allocate(tmp)
    if (allocated(this%field_data)) then
      allocate(tmp%field_data(size(this%field_data)))
      tmp%field_data = cos(this%field_data)
    end if
    call move_alloc(tmp, res)
  end function cheb1d_scalar_cos

  pure function cheb1d_scalar_tan(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\tan({\rm field})\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: tmp
    allocate(tmp)
    if (allocated(this%field_data)) then
      allocate(tmp%field_data(size(this%field_data)))
      tmp%field_data = tan(this%field_data)
    end if
    call move_alloc(tmp, res)
  end function cheb1d_scalar_tan

  pure function cheb1d_scalar_asin(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\sin^{-1}({\rm field})\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: tmp
    allocate(tmp)
    if (allocated(this%field_data)) then
      allocate(tmp%field_data(size(this%field_data)))
      tmp%field_data = asin(this%field_data)
    end if
    call move_alloc(tmp, res)
  end function cheb1d_scalar_asin

  pure function cheb1d_scalar_acos(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\cos^{-1}({\rm field})\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: tmp
    allocate(tmp)
    if (allocated(this%field_data)) then
      allocate(tmp%field_data(size(this%field_data)))
      tmp%field_data = acos(this%field_data)
    end if
    call move_alloc(tmp, res)
  end function cheb1d_scalar_acos

  pure function cheb1d_scalar_atan(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\tan^{-1}({\rm field})\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: tmp
    allocate(tmp)
    if (allocated(this%field_data)) then
      allocate(tmp%field_data(size(this%field_data)))
      tmp%field_data = atan(this%field_data)
    end if
    call move_alloc(tmp, res)
  end function cheb1d_scalar_atan

  pure function cheb1d_scalar_sinh(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\sinh({\rm field})\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: tmp
    allocate(tmp)
    if (allocated(this%field_data)) then
      allocate(tmp%field_data(size(this%field_data)))
      tmp%field_data = sinh(this%field_data)
    end if
    call move_alloc(tmp, res)
  end function cheb1d_scalar_sinh

  pure function cheb1d_scalar_cosh(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\cosh({\rm field})\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: tmp
    allocate(tmp)
    if (allocated(this%field_data)) then
      allocate(tmp%field_data(size(this%field_data)))
      tmp%field_data = cosh(this%field_data)
    end if
    call move_alloc(tmp, res)
  end function cheb1d_scalar_cosh

  pure function cheb1d_scalar_tanh(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\tanh({\rm field})\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: tmp
    allocate(tmp)
    if (allocated(this%field_data)) then
      allocate(tmp%field_data(size(this%field_data)))
      tmp%field_data = tanh(this%field_data)
    end if
    call move_alloc(tmp, res)
  end function cheb1d_scalar_tanh

  pure function cheb1d_scalar_asinh(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\sinh^{-1}({\rm field})\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: tmp
    allocate(tmp)
    if (allocated(this%field_data)) then
      allocate(tmp%field_data(size(this%field_data)))
      tmp%field_data = asinh(this%field_data)
    end if
    call move_alloc(tmp, res)
  end function cheb1d_scalar_asinh

  pure function cheb1d_scalar_acosh(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\cosh^{-1}({\rm field})\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: tmp
    allocate(tmp)
    if (allocated(this%field_data)) then
      allocate(tmp%field_data(size(this%field_data)))
      tmp%field_data = acosh(this%field_data)
    end if
    call move_alloc(tmp, res)
  end function cheb1d_scalar_acosh

  pure function cheb1d_scalar_atanh(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\tanh^{-1}({\rm field})\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: tmp
    allocate(tmp)
    if (allocated(this%field_data)) then
      allocate(tmp%field_data(size(this%field_data)))
      tmp%field_data = atanh(this%field_data)
    end if
    call move_alloc(tmp, res)
  end function cheb1d_scalar_atanh

  pure function cheb1d_scalar_log(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\ln({\rm field})\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: tmp
    allocate(tmp)
    if (allocated(this%field_data)) then
      allocate(tmp%field_data(size(this%field_data)))
      tmp%field_data = log(this%field_data)
    end if
    call move_alloc(tmp, res)
  end function cheb1d_scalar_log

  pure function cheb1d_scalar_log10(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\log({\rm field})\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: tmp
    allocate(tmp)
    if (allocated(this%field_data)) then
      allocate(tmp%field_data(size(this%field_data)))
      tmp%field_data = log10(this%field_data)
    end if
    call move_alloc(tmp, res)
  end function cheb1d_scalar_log10

  pure function cheb1d_scalar_exp(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(e^{\rm field}\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: tmp
    allocate(tmp)
    if (allocated(this%field_data)) then
      allocate(tmp%field_data(size(this%field_data)))
      tmp%field_data = exp(this%field_data)
    end if
    call move_alloc(tmp, res)
  end function cheb1d_scalar_exp

  pure function cheb1d_scalar_abs(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(|{\rm field}|\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: tmp
    allocate(tmp)
    if (allocated(this%field_data)) then
      allocate(tmp%field_data(size(this%field_data)))
      tmp%field_data = abs(this%field_data)
    end if
    call move_alloc(tmp, res)
  end function cheb1d_scalar_abs

  pure function cheb1d_scalar_sqrt(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\sqrt{{\rm field}}\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: tmp
    allocate(tmp)
    if (allocated(this%field_data)) then
      allocate(tmp%field_data(size(this%field_data)))
      tmp%field_data = sqrt(this%field_data)
    end if
    call move_alloc(tmp, res)
  end function cheb1d_scalar_sqrt

  pure function cheb1d_scalar_minval(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\min({\rm field})\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    real(r8) :: res !! The result of this operation
    if (allocated(this%field_data)) then
      res = minval(this%field_data)
    else
      res = 0.0_r8
    end if
  end function cheb1d_scalar_minval

  pure function cheb1d_scalar_maxval(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\max({\rm field})\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    real(r8) :: res !! The result of this operation
    if (allocated(this%field_data)) then
      res = maxval(this%field_data)
    else
      res = 0.0_r8
    end if
  end function cheb1d_scalar_maxval
  
  pure function cheb1d_scalar_laplacian(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\nabla^2 {\rm field}\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: tmp
    allocate(tmp)
    if (allocated(this%field_data)) then
      allocate(tmp%field_data(size(this%field_data)))
    end if
    call move_alloc(tmp, res)
  end function cheb1d_scalar_laplacian
  
  pure function cheb1d_scalar_gradient(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\nabla{\rm field}\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(vector_field), allocatable :: res !! The result of this operation
  end function cheb1d_scalar_gradient
  
  elemental subroutine cheb1d_scalar_assign(this,rhs)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm field} = {\rm field}\)
    !
    class(cheb1d_scalar_field), intent(inout) :: this
    class(scalar_field), intent(in) :: rhs
    select type(rhs)
    class is(cheb1d_scalar_field)
      this%numpoints = rhs%numpoints
      this%extent = rhs%extent
      if (allocated(rhs%colloc_points)) then
        this%colloc_points = rhs%colloc_points
      else if (allocated(this%colloc_points)) then
        deallocate(this%colloc_points)
      end if
      if (allocated(rhs%field_data)) then
        this%field_data = rhs%field_data
      else if (allocated(this%field_data)) then
        deallocate(this%field_data)
      end if
    end select
  end subroutine cheb1d_scalar_assign

!~   subroutine cheb1d_scalar_set_boundaries(this,direction,lower, &
!~                                           free_bound,bound_val, &
!~                                           bound_deriv)
!~     !* Author: Chris MacMackin
!~     !  Date: March 2016
!~     !
!~     ! Sets boundary conditions and values. If a boundary is not
!~     ! explicitly set then it defaults to being free.
!~     !
!~ 
!~     class(cheb1d_scalar_field), intent(inout) :: this
!~       integer, intent(in) :: direction
!~         !! The number corresponding to the direction/dimension whose
!~         !! boundary condition is to be set
!~       logical, optional, intent(in) :: lower
!~         !! Sets lower boundary if true (default), upper if false
!~       logical, optional, intent(in) :: free_bound
!~         !! If true, makes this a free boundary. Any boundary values 
!~         !! passed will be ignored. Default is `.false.`.
!~       real(r8), optional, intent(in) :: bound_val
!~         !! Value of the field at the boundary. Default is 0.
!~       real(r8), optional, intent(in) :: bound_deriv
!~         !! Value of the first derivative of the field at the boundary.
!~         !! Default is 0.
!~   end subroutine cheb1d_scalar_set_boundaries

  function cheb1d_scalar_d_dx(this, dir, order) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\frac{\partial^{\rm order}}{\partial x_{\rm dir}^{\rm order}}{\rm field}\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    integer, intent(in) :: dir !! Direction in which to differentiate
    integer, optional, intent(in) :: order !! Order of the derivative, default = 1
    class(scalar_field), allocatable :: res
    type(cheb1d_scalar_field), allocatable :: tmp
    allocate(tmp)
    if (allocated(this%field_data)) then
      allocate(tmp%field_data(size(this%field_data)))
    end if
    call move_alloc(tmp, res)
  end function cheb1d_scalar_d_dx

  logical function cheb1d_scalar_same_bound(this,other) result(res)
    !* Author: Chris MacMackin
    !  Date: April 2016
    !
    ! Checks wither the boundary conditions for other field agree with
    ! those for this one, within tolerance.
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), intent(in) :: other
    res = .true.
  end function cheb1d_scalar_same_bound

  subroutine check_compatible_scalar(this,other)
    !* Author: Chris MacMackin
    !  Date: April 2016
    !
    ! Checks whether a field has the same type, boundaries, and
    ! resolution as this one, making it compatible for binary 
    ! operations. If incompatible, stops program with an error message.
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), intent(in) :: other
      !! The field being checked against this one
    character(:), allocatable :: err_message
    logical :: type_err, domain_err, res_err, has_message
    has_message = .false.
    select type(other)
    class is(cheb1d_scalar_field)
      type_err = .false.
      domain_err = any(this%extent /= other%extent)
      res_err = this%numpoints /= other%numpoints
    class default
      type_err = .true.
      domain_err = .false.
      res_err = .false.
    end select
    if (type_err) then
      err_message = 'incompatible types'
      has_message = .true.
    end if
    if (domain_err) then
      if (has_message) err_message = err_message//', '
      err_message = err_message//'different domains'
      has_message = .true.
    end if
    if (res_err) then
      if (has_message) err_message = err_message//', '
      err_message = 'different resolutions'
      has_message = .true.
    end if
    if (has_message) then
      write(stderr,*) 'cheb1d_scalar_field: Error, operation with incompatible fields'
      write(stderr,*) '    Following inconsistencies: '//err_message
#ifdef __GFORTRAN__
      call backtrace
#endif
      error stop
    end if
  end subroutine check_compatible_scalar

  subroutine check_compatible_vector(this,other)
    !* Author: Chris MacMackin
    !  Date: April 2016
    !
    ! Checks whether a vector field has the same boundaries, and
    ! resolution as this one and is implemented using the same numerics,
    ! making it compatible for binary operations. If incompatible, stops
    ! program with an error message.
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(vector_field), intent(in) :: other
      !! The field being checked against this one
    character(len=:), allocatable :: err_message
    logical :: type_err, domain_err, res_err, has_message
!~     has_message = .false.
!~     select type(other)
!~     case type is(cheb1d_vector_field)
!~       type_error = .true.
!~     case default
!~       type_error = .false.
!~     end select
!~     domain_err = (this%domain() \= other%domain())
!~     res_err = (this%resolution() \= other%resolution())
!~     if (type_err) then
!~       err_message = 'incompatible types'
!~       has_message = .true.
!~     end if
!~     if (domain_err) then
!~       if (has_message) err_message = err_message//', '
!~       err_message = err_message//'different domains'
!~       has_message = .true.
!~     end if
!~     if (res_err) then
!~       if (has_message) err_message = err_message//', '
!~       err_message = 'different resolutions'
!~       has_message = .true.
!~     end if
!~     if (has_message) then
!~       write(stderr,*) 'cheb1d_scalar_field: Error, operation with incompatible field'
!~       write(stderr,*) '    Following inconsistencies: '//err_message
!~ #ifdef __GFORTRAN__
!~       call backtrace
!~ #endif
!~       error stop
!~     end if
  end subroutine check_compatible_vector
  
  pure subroutine cheb1d_assign_meta_data(this, rhs)
    !* Author: Chris MacMackin
    !  Date: April 2016
    !
    ! Assigns all data in `rhs` to `this` except for the actual values
    ! of the field.
    !
    class(cheb1d_scalar_field), intent(inout) :: this
    class(abstract_field), intent(in) :: rhs
    select type(rhs)
    class is(cheb1d_scalar_field)
      this%numpoints = rhs%numpoints
      this%extent = rhs%extent
      if (allocated(rhs%colloc_points)) then
        if (.not. allocated(this%colloc_points)) &
          allocate(this%colloc_points(this%numpoints + 1))
        this%colloc_points = rhs%colloc_points
      else if (allocated(this%colloc_points)) then
        deallocate(this%colloc_points)
      end if
!~     class is(cheb1d_vector_field)
!~       this%numpoints = rhs%numpoints
!~       this%extent = rhs%extent
!~       this%colloc_points = rhs%colloc_points
    end select
  end subroutine cheb1d_assign_meta_data

end module cheb1d_fields_mod
