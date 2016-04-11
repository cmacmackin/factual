!
!  cheb2d_fields.f90
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

module cheb1d_fields_mod
  !* Author: Chris MacMackin
  !  Date: March 2016
  !  License: GPLv3
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

  type, extends(scalar_field), public :: cheb1d_scalar_field
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! A concrete implementations of the [[scalar_field]] for 1D fields.
    ! It tracks the values of the fields at Chebyshev collocation nodes
    ! and uses a pseudospectral approach to differentiate. In addition 
    ! to the type-bound procedures, all of the intrinsic mathematical
    ! procedures, with the exception of Bessel functions and error
    ! functions.
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
    procedure :: is_equal => cheb1d_scalar_is_equal
      !! Checks fields are equal within a tolerance
    procedure :: assign_field => cheb1d_scalar_assign
      !! \({\rm field} = {\rm field}\)
    procedure, public :: assign_meta_data => cheb1d_scalar_assign_meta_data
      !! Copies all data other than values stored in field from another
      !! field object to this one.
    procedure :: check_compatible => cheb1d_scalar_check_compatible
  end type
  
  interface cheb1d_scalar_field
    module procedure scalar_constructor
  end interface

  type, extends(vector_field), public :: cheb1d_vector_field
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! A concrete implementations of the [[vector_field]] for 1D fields.
    ! Although the field is 1D, it can handle vectors with arbitrary
    ! numbers of components. It tracks the values of the fields at 
    ! Chebyshev collocation nodes and uses a pseudospectral approach to 
    ! differentiate.
    !
    private
    integer :: numpoints
      !! The number of datapoints used
    real(r8), dimension(1,2) :: extent = reshape([-1.0_r8, 1.0_r8],[1,2])
      !! The start and end values of the domain
    real(r8), dimension(:), allocatable :: colloc_points
      !! The location of the data-points
    real(r8), dimension(:,:), allocatable :: field_data
      !! The value of the vector field at the data-points
    integer :: extra_dims = 0
      !! The number of vector components in addition to those parallel
      !! to the spatial directions represented in the grid.
  contains
    private
    procedure, public :: domain => cheb1d_vector_domain
      !! Provides array with upper and lower limits of field's domain
    procedure, public :: dimensions => cheb1d_vector_dimensions
      !! Returns number of dimensions of the field
    procedure, public :: raw_size => cheb1d_vector_raw_size
      !! Provides the number of pieces of data needed to represent the
      !! field, i.e. the size of the array returned by `get_raw`.
    procedure, public :: raw => cheb1d_vector_raw
      !! Returns array of data representing state of field. Can be
      !! useful for passing to nonlinear solvers.
    procedure, public :: resolution => cheb1d_vector_resolution
      !! Returns array containing number of datapoints in each dimension.
    procedure, public :: set_from_raw => cheb1d_vector_set_from_raw
      !! Assigns raw data, such as that produced by 
      !! [[cheb1d_vector_field:raw]], to the field
    procedure :: field_multiply_field => cheb1d_vector_vf_m_sf
      !! \({\rm field} \times {\rm field}\)
    procedure, pass(rhs) :: real_multiply_field => cheb1d_vector_r_m_vf
      !! \({\rm real}  \times {\rm field}\)
    procedure :: field_multiply_real => cheb1d_vector_vf_m_r
      !! \({\rm field} \times {\rm real}\)
    procedure :: field_divide_field => cheb1d_vector_vf_d_sf
      !! \(\frac{\rm field}{\rm field}\)
    procedure :: field_divide_real => cheb1d_vector_vf_d_r
      !! \(\frac{\rm field}{\rm real}\)
    procedure :: field_add_field => cheb1d_vector_vf_a_vf
      !! \({\rm field} + {\rm field}\)
    procedure, pass(rhs) :: real_add_field => cheb1d_vector_r_a_vf
      !! \({\rm real} + {\rm field}\)
    procedure :: field_add_real => cheb1d_vector_vf_a_r
      !! \({\rm field} + {\rm real}\)
    procedure :: field_sub_field => cheb1d_vector_vf_s_vf
      !! \({\rm field} - {\rm field}\)
    procedure, pass(rhs) :: real_sub_field => cheb1d_vector_r_s_vf
      !! \({\rm real} - {\rm field}\)
    procedure :: field_sub_real => cheb1d_vector_vf_s_r
      !! \({\rm field} - {\rm real}\)
    procedure, public :: norm => cheb1d_vector_norm
      !! \(\lVert {\rm \vec{field}} \rVert\)
    procedure, public :: component => cheb1d_vector_component
      !! Returns a scalar field containing the specified component of 
      !! the vector field
    procedure, public :: d_dx => cheb1d_vector_d_dx
      !! \(\frac{\partial^n}{\partial x_i^n}({\rm field})\)
    procedure :: divergence => cheb1d_vector_divergence
      !! \(\nabla\cdot {\rm field}\)
    procedure :: curl => cheb1d_vector_curl
      !! \(\nabla\times {\rm field}\)
    procedure :: laplacian => cheb1d_vector_laplacian
      !! \(\nabla^2 {\rm field}\)
    procedure :: dot_prod => cheb1d_vector_dot_prod
      !! \({\rm \vec{field}} \cdot {\rm \vec{field}}\)
    procedure :: cross_prod => cheb1d_vector_cross_prod
      !! \({\rm\vec{field}} \times {\rm\vec{field}}\)
    procedure :: assign_field => cheb1d_vector_assign
      !! \({\rm field} = {\rm field}\)
    procedure :: is_equal => cheb1d_vector_is_equal
      !! Checks fields are equal within a tolerance
    procedure, public :: assign_meta_data => cheb1d_vector_assign_meta_data
      !! Copies all data other than values stored in field from another
      !! field object to this one.
    procedure :: check_compatible => cheb1d_vector_check_compatible
  end type
  
  interface cheb1d_vector_field
    module procedure vector_constructor
  end interface
  
  abstract interface
    pure function scalar_init(x) result(scalar)
      !! Function used to specify value held by a scalar field at
      !! location `x`.
      import :: r8
      real(r8), intent(in) :: x 
        !! The position at which this function is evaluated
      real(r8) :: scalar
        !! The value of the field at this location
    end function scalar_init
  
    pure function vector_init(x) result(vector)
      !! Function used to specify value held by a vector field at
      !! location `x`.
      import :: r8
      real(r8), intent(in) :: x 
        !! The position at which this function is evaluated
      real(r8), dimension(:), allocatable :: vector
        !! The value of the field at this location
    end function vector_init
  end interface
  
contains

  !=====================================================================
  ! Scalar Field Methods
  !=====================================================================

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
      write(stderr,"('    Needed: ',i0,', Actual: ',i0)") expected, size(raw)
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
    call this%check_compatible(rhs)
#endif
    allocate(local)
    call local%assign_meta_data(this)
    select type(rhs)
    class is(cheb1d_scalar_field)
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
    type(cheb1d_vector_field), allocatable :: local
    integer :: i
#ifdef DEBUG
    call this%check_compatible(rhs)
#endif
    allocate(local)
    call local%assign_meta_data(rhs)
    select type(rhs)
    class is(cheb1d_vector_field)
      forall(i=1:this%numpoints+1) &
        local%field_data(i,:) = this%field_data(i) * rhs%field_data(i,:)
    end select
    call move_alloc(local,res)
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
    call local%assign_meta_data(rhs)
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
    call local%assign_meta_data(this)
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
    call this%check_compatible(rhs)
#endif
    allocate(local)
    call local%assign_meta_data(this)
    select type(rhs)
    class is(cheb1d_scalar_field)
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
    call local%assign_meta_data(rhs)
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
    call local%assign_meta_data(this)
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
    call this%check_compatible(rhs)
#endif
    allocate(local)
    call local%assign_meta_data(this)
    select type(rhs)
    class is(cheb1d_scalar_field)
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
    call local%assign_meta_data(rhs)
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
    call local%assign_meta_data(this)
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
    call this%check_compatible(rhs)
#endif
    allocate(local)
    call local%assign_meta_data(this)
    select type(rhs)
    class is(cheb1d_scalar_field)
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
    call local%assign_meta_data(rhs)
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
    call local%assign_meta_data(this)
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
    call local%assign_meta_data(this)
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
    call local%assign_meta_data(this)
    local%field_data = this%field_data ** rhs
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
    call local%assign_meta_data(this)
    local%field_data = this%field_data ** rhs
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
    type(cheb1d_scalar_field), allocatable :: local
    allocate(local)
    call local%assign_meta_data(this)
    if (allocated(this%field_data)) then
      local%field_data = sin(this%field_data)
    end if
    call move_alloc(local, res)
  end function cheb1d_scalar_sin

  pure function cheb1d_scalar_cos(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\cos({\rm field})\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: local
    allocate(local)
    call local%assign_meta_data(this)
    if (allocated(this%field_data)) then
      local%field_data = cos(this%field_data)
    end if
    call move_alloc(local, res)
  end function cheb1d_scalar_cos

  pure function cheb1d_scalar_tan(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\tan({\rm field})\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: local
    allocate(local)
    call local%assign_meta_data(this)
    if (allocated(this%field_data)) then
      local%field_data = tan(this%field_data)
    end if
    call move_alloc(local, res)
  end function cheb1d_scalar_tan

  pure function cheb1d_scalar_asin(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\sin^{-1}({\rm field})\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: local
    allocate(local)
    call local%assign_meta_data(this)
    if (allocated(this%field_data)) then
      local%field_data = asin(this%field_data)
    end if
    call move_alloc(local, res)
  end function cheb1d_scalar_asin

  pure function cheb1d_scalar_acos(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\cos^{-1}({\rm field})\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: local
    allocate(local)
    call local%assign_meta_data(this)
    if (allocated(this%field_data)) then
      local%field_data = acos(this%field_data)
    end if
    call move_alloc(local, res)
  end function cheb1d_scalar_acos

  pure function cheb1d_scalar_atan(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\tan^{-1}({\rm field})\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: local
    allocate(local)
    call local%assign_meta_data(this)
    if (allocated(this%field_data)) then
      local%field_data = atan(this%field_data)
    end if
    call move_alloc(local, res)
  end function cheb1d_scalar_atan

  pure function cheb1d_scalar_sinh(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\sinh({\rm field})\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: local
    allocate(local)
    call local%assign_meta_data(this)
    if (allocated(this%field_data)) then
      local%field_data = sinh(this%field_data)
    end if
    call move_alloc(local, res)
  end function cheb1d_scalar_sinh

  pure function cheb1d_scalar_cosh(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\cosh({\rm field})\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: local
    allocate(local)
    call local%assign_meta_data(this)
    if (allocated(this%field_data)) then
      local%field_data = cosh(this%field_data)
    end if
    call move_alloc(local, res)
  end function cheb1d_scalar_cosh

  pure function cheb1d_scalar_tanh(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\tanh({\rm field})\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: local
    allocate(local)
    call local%assign_meta_data(this)
    if (allocated(this%field_data)) then
      local%field_data = tanh(this%field_data)
    end if
    call move_alloc(local, res)
  end function cheb1d_scalar_tanh

  pure function cheb1d_scalar_asinh(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\sinh^{-1}({\rm field})\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: local
    allocate(local)
    call local%assign_meta_data(this)
    if (allocated(this%field_data)) then
      local%field_data = asinh(this%field_data)
    end if
    call move_alloc(local, res)
  end function cheb1d_scalar_asinh

  pure function cheb1d_scalar_acosh(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\cosh^{-1}({\rm field})\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: local
    allocate(local)
    call local%assign_meta_data(this)
    if (allocated(this%field_data)) then
      local%field_data = acosh(this%field_data)
    end if
    call move_alloc(local, res)
  end function cheb1d_scalar_acosh

  pure function cheb1d_scalar_atanh(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\tanh^{-1}({\rm field})\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: local
    allocate(local)
    call local%assign_meta_data(this)
    if (allocated(this%field_data)) then
      local%field_data = atanh(this%field_data)
    end if
    call move_alloc(local, res)
  end function cheb1d_scalar_atanh

  pure function cheb1d_scalar_log(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\ln({\rm field})\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: local
    allocate(local)
    call local%assign_meta_data(this)
    if (allocated(this%field_data)) then
      local%field_data = log(this%field_data)
    end if
    call move_alloc(local, res)
  end function cheb1d_scalar_log

  pure function cheb1d_scalar_log10(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\log({\rm field})\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: local
    allocate(local)
    call local%assign_meta_data(this)
    if (allocated(this%field_data)) then
      local%field_data = log10(this%field_data)
    end if
    call move_alloc(local, res)
  end function cheb1d_scalar_log10

  pure function cheb1d_scalar_exp(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(e^{\rm field}\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: local
    allocate(local)
    call local%assign_meta_data(this)
    if (allocated(this%field_data)) then
      local%field_data = exp(this%field_data)
    end if
    call move_alloc(local, res)
  end function cheb1d_scalar_exp

  pure function cheb1d_scalar_abs(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(|{\rm field}|\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: local
    allocate(local)
    call local%assign_meta_data(this)
    if (allocated(this%field_data)) then
      local%field_data = abs(this%field_data)
    end if
    call move_alloc(local, res)
  end function cheb1d_scalar_abs

  pure function cheb1d_scalar_sqrt(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\sqrt{{\rm field}}\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: local
    allocate(local)
    call local%assign_meta_data(this)
    if (allocated(this%field_data)) then
      local%field_data = sqrt(this%field_data)
    end if
    call move_alloc(local, res)
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
  
  function cheb1d_scalar_laplacian(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\nabla^2 {\rm field}\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: local
    allocate(local)
    call local%assign_meta_data(this)
    local = this%d_dx(1,2)
    call move_alloc(local, res)
  end function cheb1d_scalar_laplacian
  
  function cheb1d_scalar_gradient(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\nabla{\rm field}\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(vector_field), allocatable :: res !! The result of this operation
    type(cheb1d_vector_field), allocatable :: local
    allocate(local)
    call local%assign_meta_data(this)
    local%field_data(:,1) = this%field_data
    call differentiate_1d(local%field_data(:,1),local%colloc_points)
    call move_alloc(local,res)
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
      call this%assign_meta_data(rhs)
      if (allocated(rhs%field_data)) this%field_data = rhs%field_data
    end select
  end subroutine cheb1d_scalar_assign

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
    type(cheb1d_scalar_field), allocatable :: local
    allocate(local)
    call local%assign_meta_data(this)
    if (dir==1) then
      local%field_data = this%field_data
      call differentiate_1d(local%field_data,local%colloc_points,order)
    else
      local%field_data = 0.0
    end if    
    call move_alloc(local, res)
  end function cheb1d_scalar_d_dx

  logical function cheb1d_scalar_is_equal(this,rhs) result(iseq)
    !* Author: Chris MacMackin
    !  Date: April 2016
    !
    ! Evaluates whether two scalar fields are equal within a tolerance,
    ! specified by [[set_tol]].
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), intent(in) :: rhs
    integer :: i
    real(r8) :: normalization
    iseq = .true.
    select type(rhs)
    class is(cheb1d_scalar_field)
      if (this%numpoints /= rhs%numpoints) then
        iseq = .false.
        return
      end if
      do i=1,this%numpoints+1
        normalization = abs(this%field_data(i))
        if (normalization < tiny(normalization)) normalization = 1.0_r8
        iseq = iseq .and. ( ((this%field_data(i)-rhs%field_data(i))/normalization < &
                             get_tol()) .or. &
                            (isnan(this%field_data(i)).and.isnan(rhs%field_data(i))) )
        if (.not. iseq) return
      end do
    class default
      iseq = .false.
    end select
  end function cheb1d_scalar_is_equal

  subroutine cheb1d_scalar_check_compatible(this,other)
    !* Author: Chris MacMackin
    !  Date: April 2016
    !
    ! Checks whether a field has the same type, boundaries, and
    ! resolution as this one, making it compatible for binary 
    ! operations. If incompatible, stops program with an error message.
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(abstract_field), intent(in) :: other
      !! The field being checked against this one
    character(len=128) :: err_message
    logical :: type_err, domain_err, res_err, alloc_err, has_message
    has_message = .false.
    select type(other)
    class is(cheb1d_scalar_field)
      type_err = .false.
      domain_err = any(abs(this%extent - other%extent) > 1.e-15_r8)
      res_err = this%numpoints /= other%numpoints
      alloc_err = .not. (allocated(this%field_data).and.allocated(other%field_data))
    class is(cheb1d_vector_field)
      type_err = .false.
      domain_err = any(abs(this%extent - other%extent) > 1.e-15_r8)
      res_err = this%numpoints /= other%numpoints
      alloc_err = .not. (allocated(this%field_data).and.allocated(other%field_data))
    class default
      type_err = .true.
      domain_err = .false.
      res_err = .false.
      alloc_err = .false.
    end select
    err_message = ''
    if (type_err) then
      err_message = ' incompatible types'
      has_message = .true.
    end if
    if (domain_err) then
      if (has_message) err_message = trim(err_message)//','
      err_message = trim(err_message)//' different domains'
      has_message = .true.
    end if
    if (res_err) then
      if (has_message) err_message = trim(err_message)//','
      err_message = trim(err_message)//' different resolutions'
      has_message = .true.
    end if
    if (alloc_err) then
      if (has_message) err_message = trim(err_message)//','
      err_message = trim(err_message)//' uninitialized field'
      has_message = .true.
    end if
    if (has_message) then
      write(stderr,*) 'cheb1d_scalar_field: Error, operation with incompatible fields'
      write(stderr,*) '    Following inconsistencies:'//trim(err_message)
#ifdef __GFORTRAN__
      call backtrace
#endif
      error stop
    end if
  end subroutine cheb1d_scalar_check_compatible

  pure subroutine cheb1d_scalar_assign_meta_data(this, rhs, alloc)
    !* Author: Chris MacMackin
    !  Date: April 2016
    !
    ! Assigns all data in `rhs` to `this` except for the actual values
    ! of the field.
    !
    class(cheb1d_scalar_field), intent(inout) :: this
    class(abstract_field), intent(in) :: rhs
    logical, optional, intent(in) :: alloc
      !! If present and false, do not allocate the array of `this`.
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
      if (allocated(this%field_data)) deallocate(this%field_data)
      if (present(alloc)) then
        if (allocated(rhs%field_data) .and. alloc) &
          allocate(this%field_data(size(rhs%field_data))) !Needed due to compiler bug
      else
        allocate(this%field_data(size(rhs%field_data))) !Needed due to compiler bug
      end if
    class is(cheb1d_vector_field)
      this%numpoints = rhs%numpoints
      this%extent = rhs%extent
      if (allocated(rhs%colloc_points)) then
        if (.not. allocated(this%colloc_points)) &
          allocate(this%colloc_points(this%numpoints + 1))
        this%colloc_points = rhs%colloc_points
      else if (allocated(this%colloc_points)) then
        deallocate(this%colloc_points)
      end if
      if (allocated(this%field_data)) deallocate(this%field_data)
      if (present(alloc)) then
        if (allocated(rhs%field_data) .and. alloc) &
          allocate(this%field_data(size(rhs%field_data,1))) !Needed due to compiler bug
      else
        allocate(this%field_data(size(rhs%field_data,1))) !Needed due to compiler bug
      end if
    end select
  end subroutine cheb1d_scalar_assign_meta_data

  !=====================================================================
  ! Vector Field Methods
  !=====================================================================

  function vector_constructor(nodes,initializer,lower_bound, &
                              upper_bound, extra_dims) result(field)
    integer, intent(in) :: nodes
      !! The number of collocation nodes to use when modelling this
      !! field. This corresponds to resolution.
    procedure(vector_init), optional :: initializer
      !! An elemental procedure taking which takes the position in the
      !! fields domain (an 8-byte real) as an argument and returns the
      !! fields value at that position. Default is for field to be zero
      !! everywhere.
    real(r8), intent(in), optional :: lower_bound
      !! The position of the lower bound of the field's domain. Default is -1.0.
    real(r8), intent(in), optional :: upper_bound
      !! The position of the upper bound of the field's domain. Default is 1.0.
    integer, intent(in), optional :: extra_dims
      !! The number of vector components in addition to those parallel
      !! to the spatial directions represented in the grid. Ignored if 
      !! less than 0.
    type(cheb1d_vector_field) :: field
    real(r8), dimension(:), allocatable :: tmp
    integer :: i, j, dims, initializer_dims
    dims = 1
    if (present(extra_dims)) then
      if (extra_dims > 0) then
        dims = dims + extra_dims
        field%extra_dims = extra_dims
      end if
    end if
    field%numpoints = nodes
    if (present(lower_bound)) field%extent(1,1) = lower_bound
    if (present(upper_bound)) field%extent(1,2) = upper_bound
    field%colloc_points = collocation_points(nodes,lower_bound,upper_bound)
    allocate(field%field_data(nodes+1,dims))
    if (present(initializer)) then
      initializer_dims = size(initializer(field%extent(1,1)))
      if (initializer_dims < dims) then
        tmp = [(0.0, i=1,(1+dims-initializer_dims))]
        forall (i=1:nodes+1) &
          field%field_data(i,:) = [initializer(field%colloc_points(i)),&
                                   tmp]
      else if (initializer_dims > dims) then
        do i=1,nodes+1
          tmp = initializer(field%colloc_points(i))
          field%field_data(i,:) = tmp(1:dims)
        end do
      else
        forall (i=1:nodes+1) &
          field%field_data(i,:) = initializer(field%colloc_points(i))
      end if
    else
      field%field_data = 0.0_r8
    end if
  end function vector_constructor

  pure function cheb1d_vector_domain(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! Accessor for the upper and lower bounds of the domain.
    !
    class(cheb1d_vector_field), intent(in) :: this
    real(r8), dimension(:,:), allocatable :: res
      !* A 2D array of shape \(n \times 2\), where n is the number of
      !  dimensions of the field. Each row contains the lower and the
      !  upper extent of the fields domain in the corresponding 
      !  dimension
    res = this%extent
  end function cheb1d_vector_domain

  elemental function cheb1d_vector_dimensions(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! Accessor for the number of dimensions.
    !
    class(cheb1d_vector_field), intent(in) :: this
    integer :: res
    res = 1
  end function cheb1d_vector_dimensions

  pure function cheb1d_vector_raw_size(this,return_lower_bound, &
                                       return_upper_bound) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! Compute how many elements are in the raw representation of this
    ! field. This would be the number of data points, adjusted based on
    ! how boundary conditions are accounted for.
    !
    class(cheb1d_vector_field), intent(in) :: this
    logical, dimension(:), optional, intent(in) :: return_lower_bound
      !! Specifies whether to return the values at the lower boundary
      !! for each dimension, with the index of the element
      !! corresponding to the dimension. Defaults to all `.true.`.
    logical, dimension(:), optional, intent(in) :: return_upper_bound
      !! Specifies whether to return the values at the upper boundary
      !! for each dimension, with the index of the element
      !! corresponding to the dimension. Defaults to all `.true.`.
    integer :: res
    res = (this%numpoints + 1)
    if (present(return_lower_bound)) then
      if (.not. return_lower_bound(1)) res = res - 1
    end if
    if (present(return_upper_bound)) then
      if (.not. return_upper_bound(1)) res = res - 1
    end if
    res = res*(this%extra_dims + this%dimensions())
  end function cheb1d_vector_raw_size
  
  pure function cheb1d_vector_raw(this,return_lower_bound, &
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
    class(cheb1d_vector_field), intent(in) :: this
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
      res = pack(this%field_data(lower:upper,:),.true.)
    else
      res = [0.0_r8]
    end if
  end function cheb1d_vector_raw
  
  pure function cheb1d_vector_resolution(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! Get the number of data points used to represent the field in
    ! its one dimension.
    !
    class(cheb1d_vector_field), intent(in) :: this
    integer, dimension(:), allocatable :: res
      !! Array of length 1 specifying the number of data points.
    allocate(res(1))
    res(1) = this%numpoints
  end function cheb1d_vector_resolution

  subroutine cheb1d_vector_set_from_raw(this,raw,provide_lower_bound, &
                                        provide_upper_bound)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! Assigns raw data, such as that produced by 
    ! [[cheb1d_vector_field:raw]], to the field. The routine will
    ! stop with an error if the array is the wrong length for a field
    ! of this resolution and with these boundary conditions. 
    !
    class(cheb1d_vector_field), intent(inout) :: this
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
      write(stderr,*) 'cheb1d_vector_field: Error, setting from raw array of wrong size'
      write(stderr,"('    Needed: ',i0,', Actual: ',i0)") expected, size(raw)
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
      allocate(this%field_data(this%numpoints + 1,this%extra_dims + 1))
    this%field_data(lower:upper,:) = reshape(raw,[upper-lower+1,this%extra_dims + 1])
  end subroutine cheb1d_vector_set_from_raw
  
  function cheb1d_vector_vf_m_sf(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\vec{\rm field} \times {\rm field}\)
    !
    class(cheb1d_vector_field), intent(in) :: this
    class(scalar_field), intent(in) :: rhs
    class(vector_field), allocatable :: res !! The restult of this operation
    type(cheb1d_vector_field), allocatable :: local
    integer :: i
#ifdef DEBUG
    call this%check_compatible(rhs)
#endif
    allocate(local)
    call local%assign_meta_data(this)
    select type(rhs)
    class is(cheb1d_scalar_field)
      forall (i=1:this%extra_dims+1) &
        local%field_data(:,i) = this%field_data(:,i) * rhs%field_data
    end select
    call move_alloc(local,res)
  end function cheb1d_vector_vf_m_sf

  function cheb1d_vector_r_m_vf(lhs,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm real} \times \vec{\rm field}\)
    !
    real(r8), intent(in) :: lhs
    class(cheb1d_vector_field), intent(in) :: rhs
    class(vector_field), allocatable :: res !! The result of this operation
    type(cheb1d_vector_field), allocatable :: local
    allocate(local)
    call local%assign_meta_data(rhs)
    local%field_data = lhs * rhs%field_data
    call move_alloc(local,res)
  end function cheb1d_vector_r_m_vf

  function cheb1d_vector_vf_m_r(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\vec{\rm field} \times {\rm real}\)
    !
    class(cheb1d_vector_field), intent(in) :: this
    real(r8), intent(in) :: rhs
    class(vector_field), allocatable :: res !! The result of this operation
    type(cheb1d_vector_field), allocatable :: local
    allocate(local)
    call local%assign_meta_data(this)
    local%field_data = this%field_data * rhs
    call move_alloc(local,res)
  end function cheb1d_vector_vf_m_r
  
  function cheb1d_vector_vf_d_sf(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\vec{\rm field} / {\rm field}\)
    !
    class(cheb1d_vector_field), intent(in) :: this
    class(scalar_field), intent(in) :: rhs
    class(vector_field), allocatable :: res !! The restult of this operation
    type(cheb1d_vector_field), allocatable :: local
    integer :: i
#ifdef DEBUG
    call this%check_compatible(rhs)
#endif
    allocate(local)
    call local%assign_meta_data(this)
    select type(rhs)
    class is(cheb1d_scalar_field)
      forall (i=1:this%extra_dims+1) &
        local%field_data(:,i) = this%field_data(:,i) / rhs%field_data
    end select
    call move_alloc(local,res)
  end function cheb1d_vector_vf_d_sf

  function cheb1d_vector_vf_d_r(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\vec{\rm field} / {\rm real}\)
    !
    class(cheb1d_vector_field), intent(in) :: this
    real(r8), intent(in) :: rhs
    class(vector_field), allocatable :: res !! The result of this operation
    type(cheb1d_vector_field), allocatable :: local
    allocate(local)
    call local%assign_meta_data(this)
    local%field_data = this%field_data / rhs
    call move_alloc(local,res)
  end function cheb1d_vector_vf_d_r
  
  function cheb1d_vector_vf_s_vf(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\vec{\rm field} - \vec{\rm field}\)
    !
    class(cheb1d_vector_field), intent(in) :: this
    class(vector_field), intent(in) :: rhs
    class(vector_field), allocatable :: res !! The restult of this operation
    type(cheb1d_vector_field), allocatable :: local
    integer :: min_dims, max_dims
#ifdef DEBUG
    call this%check_compatible(rhs)
#endif
    allocate(local)
    select type(rhs)
    class is(cheb1d_vector_field)
      if (rhs%extra_dims > this%extra_dims) then
        call local%assign_meta_data(rhs)
      else
        call local%assign_meta_data(this)
      end if
      min_dims = min(this%extra_dims,rhs%extra_dims) + 1
      max_dims = max(this%extra_dims,rhs%extra_dims) + 1
      local%field_data(:,1:min_dims) = this%field_data(:,1:min_dims) &
                                     - rhs%field_data(:,1:min_dims)
      if (rhs%extra_dims > this%extra_dims) then
        local%field_data(:,min_dims+1:max_dims) = -rhs%field_data(:,min_dims+1:max_dims)
      else
        local%field_data(:,min_dims+1:max_dims) = this%field_data(:,min_dims+1:max_dims)
      end if
    end select
    call move_alloc(local,res)
  end function cheb1d_vector_vf_s_vf

  function cheb1d_vector_r_s_vf(lhs,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\vec{\rm real} - \vec{\rm field}\)
    !
    real(r8), dimension(:), intent(in) :: lhs
    class(cheb1d_vector_field), intent(in) :: rhs
    class(vector_field), allocatable :: res !! The result of this operation
    type(cheb1d_vector_field), allocatable :: local
    integer :: min_dims, max_dims, i
    allocate(local)
    call local%assign_meta_data(rhs,.false.)
    min_dims = min(rhs%extra_dims + 1, size(lhs))
    max_dims = max(rhs%extra_dims + 1, size(lhs))
    allocate(local%field_data(local%numpoints + 1, max_dims))
    forall(i=1:local%numpoints+1) &
      local%field_data(i,:min_dims) = lhs(:min_dims) - rhs%field_data(i,:min_dims)
    if (rhs%extra_dims+1 > size(lhs)) then
      local%field_data(:,min_dims+1:) = -rhs%field_data(:,min_dims+1:)
    else
      local%extra_dims = size(lhs) - 1
      forall(i=1:local%numpoints+1) &
        local%field_data(i,min_dims+1:) = lhs(min_dims+1:)
    end if
    call move_alloc(local,res)
  end function cheb1d_vector_r_s_vf

  function cheb1d_vector_vf_s_r(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\vec{\rm field} - \vec{\rm real}\)
    !
    class(cheb1d_vector_field), intent(in) :: this
    real(r8), dimension(:), intent(in) :: rhs
    class(vector_field), allocatable :: res !! The result of this operation
    type(cheb1d_vector_field), allocatable :: local
    integer :: min_dims, max_dims, i
    allocate(local)
    call local%assign_meta_data(this,.false.)
    min_dims = min(this%extra_dims + 1, size(rhs))
    max_dims = max(this%extra_dims + 1, size(rhs))
    allocate(local%field_data(local%numpoints + 1, max_dims))
    forall(i=1:local%numpoints+1) &
      local%field_data(i,:min_dims) = this%field_data(i,:min_dims) - rhs(:min_dims)
    if (this%extra_dims+1 > size(rhs)) then
      local%field_data(:,min_dims+1:) = this%field_data(:,min_dims+1:)
    else
      local%extra_dims = size(rhs) - 1
      forall(i=1:local%numpoints+1) &
        local%field_data(i,min_dims+1:) = -rhs(min_dims+1:)
    end if
    call move_alloc(local,res)
  end function cheb1d_vector_vf_s_r
  
  function cheb1d_vector_vf_a_vf(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\vec{\rm field} - \vec{\rm field}\)
    !
    class(cheb1d_vector_field), intent(in) :: this
    class(vector_field), intent(in) :: rhs
    class(vector_field), allocatable :: res !! The restult of this operation
    type(cheb1d_vector_field), allocatable :: local
    integer :: max_dims, min_dims
#ifdef DEBUG
    call this%check_compatible(rhs)
#endif
    allocate(local)
    select type(rhs)
    class is(cheb1d_vector_field)
      if (rhs%extra_dims > this%extra_dims) then
        call local%assign_meta_data(rhs)
      else
        call local%assign_meta_data(this)
      end if
      min_dims = min(this%extra_dims,rhs%extra_dims) + 1
      max_dims = max(this%extra_dims,rhs%extra_dims) + 1
      local%field_data(:,1:min_dims) = this%field_data(:,1:min_dims) &
                                     + rhs%field_data(:,1:min_dims)
      if (rhs%extra_dims > this%extra_dims) then
        local%field_data(:,min_dims+1:max_dims) = rhs%field_data(:,min_dims+1:max_dims)
      else
        local%field_data(:,min_dims+1:max_dims) = this%field_data(:,min_dims+1:max_dims)
      end if
    end select
    call move_alloc(local,res)
  end function cheb1d_vector_vf_a_vf

  function cheb1d_vector_r_a_vf(lhs,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\vec{\rm real} - \vec{\rm field}\)
    !
    real(r8), dimension(:), intent(in) :: lhs
    class(cheb1d_vector_field), intent(in) :: rhs
    class(vector_field), allocatable :: res !! The result of this operation
    type(cheb1d_vector_field), allocatable :: local
    integer :: min_dims, max_dims, i
    allocate(local)
    call local%assign_meta_data(rhs,.false.)
    min_dims = min(rhs%extra_dims + 1, size(lhs))
    max_dims = max(rhs%extra_dims + 1, size(lhs))
    allocate(local%field_data(local%numpoints + 1, max_dims))
    forall(i=1:local%numpoints+1) &
      local%field_data(i,:min_dims) = lhs(:min_dims) + rhs%field_data(i,:min_dims)
    if (rhs%extra_dims+1 > size(lhs)) then
      local%field_data(:,min_dims+1:) = rhs%field_data(:,min_dims+1:)
    else
      local%extra_dims = size(lhs) - 1
      forall(i=1:local%numpoints+1) &
        local%field_data(i,min_dims+1:) = lhs(min_dims+1:)
    end if
    call move_alloc(local,res)
  end function cheb1d_vector_r_a_vf

  function cheb1d_vector_vf_a_r(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\vec{\rm field} - \vec{\rm real}\)
    !
    class(cheb1d_vector_field), intent(in) :: this
    real(r8), dimension(:), intent(in) :: rhs
    class(vector_field), allocatable :: res !! The result of this operation
    type(cheb1d_vector_field), allocatable :: local
    integer :: min_dims, max_dims, i
    allocate(local)
    call local%assign_meta_data(this,.false.)
    min_dims = min(this%extra_dims + 1, size(rhs))
    max_dims = max(this%extra_dims + 1, size(rhs))
    allocate(local%field_data(local%numpoints + 1, max_dims))
    forall(i=1:local%numpoints+1) &
      local%field_data(i,:min_dims) = this%field_data(i,:min_dims) + rhs(:min_dims)
    if (this%extra_dims+1 > size(rhs)) then
      local%field_data(:,min_dims+1:) = this%field_data(:,min_dims+1:)
    else
      local%extra_dims = size(rhs) - 1
      forall(i=1:local%numpoints+1) &
        local%field_data(i,min_dims+1:) = rhs(min_dims+1:)
    end if
    call move_alloc(local,res)
  end function cheb1d_vector_vf_a_r

  elemental subroutine cheb1d_vector_assign(this,rhs)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\vec{\rm field} = \vec{\rm field}\)
    !
    class(cheb1d_vector_field), intent(inout) :: this
    class(vector_field), intent(in) :: rhs
    call this%assign_meta_data(rhs)
    select type(rhs)
    class is(cheb1d_vector_field)
      if (allocated(rhs%field_data)) this%field_data = rhs%field_data
    end select
  end subroutine cheb1d_vector_assign

  pure function cheb1d_vector_norm(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\lVert \vec{\rm field} \rVert\)
    !
    class(cheb1d_vector_field), intent(in) :: this
    class(scalar_field), allocatable :: res
    type(cheb1d_scalar_field), allocatable :: local
    allocate(local)
    call local%assign_meta_data(this)
    select type(this)
    class is(cheb1d_vector_field)
      local%field_data = sqrt(sum(this%field_data**2,2))
    end select
    call move_alloc(local,res)
  end function cheb1d_vector_norm

  function cheb1d_vector_component(this,comp) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! Scalar field containing specified component of the vector field
    !
    class(cheb1d_vector_field), intent(in) :: this
    integer, intent(in) :: comp
    class(scalar_field), allocatable :: res
    type(cheb1d_scalar_field), allocatable :: local
    allocate(local)
    call local%assign_meta_data(this)
    select type(this)
    class is(cheb1d_vector_field)
      if (comp <= this%extra_dims + 1) then
        local%field_data = this%field_data(:,comp)
      else
        local%field_data = 0.0_r8
      end if
    end select
    call move_alloc(local,res)
  end function cheb1d_vector_component

  function cheb1d_vector_d_dx(this, dir, order) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\frac{\partial^{\rm order}}{\partial x_{\rm dir}^{\rm order}}\vec{\rm field}\)
    !
    class(cheb1d_vector_field), intent(in) :: this
    integer, intent(in) :: dir !! Direction in which to differentiate
    integer, optional, intent(in) :: order !! Order of the derivative, default = 1
    class(vector_field), allocatable :: res
    type(cheb1d_vector_field), allocatable :: local
    integer :: i
    allocate(local)
    if (dir==1) then
      local = this
      do i=1,1+local%extra_dims
        call differentiate_1d(local%field_data(:,i),local%colloc_points,order)
      end do
    else
      call local%assign_meta_data(this,.false.)
      allocate(local%field_data(local%numpoints+1,1))
      local%field_data = 0.0_r8
      local%extra_dims = 0
    end if
    call move_alloc(local, res)
  end function cheb1d_vector_d_dx

  function cheb1d_vector_laplacian(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\nabla^2 \vec{\rm field}\)
    !
    class(cheb1d_vector_field), intent(in) :: this
    class(vector_field), allocatable :: res !! The result of this operation
    type(cheb1d_vector_field), allocatable :: local
    allocate(local)
    call local%assign_meta_data(this)
    call move_alloc(local, res)
  end function cheb1d_vector_laplacian
  
  function cheb1d_vector_divergence(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\nabla\cdot \vec{\rm field}\)
    !
    class(cheb1d_vector_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: local
    allocate(local)
    call local%assign_meta_data(this)
    local%field_data = this%field_data(:,1)
    call differentiate_1d(local%field_data,local%colloc_points)
    call move_alloc(local, res)
  end function cheb1d_vector_divergence
  
  function cheb1d_vector_curl(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\nabla\times \vec{\rm field}\)
    !
    class(cheb1d_vector_field), intent(in) :: this
    class(vector_field), allocatable :: res !! The result of this operation
    type(cheb1d_vector_field), allocatable :: local
    allocate(local)
    call local%assign_meta_data(this,.false.)
    local%extra_dims = 2
    allocate(local%field_data(local%numpoints+1,3))
    local%field_data(:,1) = 0.0_r8
    if (this%extra_dims >= 2) then
      local%field_data(:,2) = -this%field_data(:,3)
      call differentiate_1d(local%field_data(:,2),local%colloc_points)
    else
      local%field_data(:,2) = 0.0_r8
    end if
    if (this%extra_dims >= 1) then
      local%field_data(:,3) = this%field_data(:,2)
      call differentiate_1d(local%field_data(:,3),local%colloc_points)
    else
      local%field_data(:,3) = 0.0_r8
    end if
    call move_alloc(local, res)
  end function cheb1d_vector_curl

  function cheb1d_vector_cross_prod(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\vec{\rm field} \times \vec{\rm field}\)
    !
    ! The returned vector will always consist of three components,
    ! regardless of the number of components in the arguments.
    !
    class(cheb1d_vector_field), intent(in) :: this
    class(vector_field), intent(in) :: rhs
    class(vector_field), allocatable :: res !! The restult of this operation
    real(r8), dimension(3) :: vec1, vec2
    integer :: i, dims1, dims2
    type(cheb1d_vector_field), allocatable :: local
#ifdef DEBUG
    call this%check_compatible(rhs)
#endif
    allocate(local)
    call local%assign_meta_data(this,.false.)
    local%extra_dims = 2
    allocate(local%field_data(this%numpoints+1,3))
    vec1 = 0
    vec2 = 0
    select type(rhs)
    class is(cheb1d_vector_field)
      dims1 = min(3,this%dimensions()+this%extra_dims)
      dims2 = min(3,rhs%dimensions()+rhs%extra_dims)
      do i=1,this%numpoints+1
        vec1(:dims1) = this%field_data(i,:dims1)
        vec2(:dims2) = rhs%field_data(i,:dims2)
        local%field_data(i,:) = [vec1(2)*vec2(3) - vec2(2)*vec1(3), &
                                 vec1(3)*vec2(1) - vec2(3)*vec1(1), &
                                 vec1(1)*vec2(2) - vec2(1)*vec1(2)]
      end do
    end select
    call move_alloc(local,res)
  end function cheb1d_vector_cross_prod

  function cheb1d_vector_dot_prod(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\vec{\rm field} \cdot \vec{\rm field}\)
    !
    class(cheb1d_vector_field), intent(in) :: this
    class(vector_field), intent(in) :: rhs
    class(scalar_field), allocatable :: res !! The restult of this operation
    type(cheb1d_scalar_field), allocatable :: local
    integer :: min_dims
#ifdef DEBUG
    call this%check_compatible(rhs)
#endif
    allocate(local)
    call local%assign_meta_data(this)
    select type(rhs)
    class is(cheb1d_vector_field)
      min_dims = min(this%extra_dims,rhs%extra_dims) + 1
      local%field_data = &
        sum(this%field_data(:,1:min_dims)*rhs%field_data(:,1:min_dims),2)
    end select
    call move_alloc(local,res)
  end function cheb1d_vector_dot_prod

  logical function cheb1d_vector_is_equal(this,rhs) result(iseq)
    !* Author: Chris MacMackin
    !  Date: April 2016
    !
    ! Evaluates whether two vector fields are equal within a tolerance,
    ! specified by [[set_tol]].
    !
    class(cheb1d_vector_field), intent(in) :: this
    class(vector_field), intent(in) :: rhs
    integer :: i, dims
    real(r8) :: normalization
    iseq = .true.
    select type(rhs)
    class is(cheb1d_vector_field)
      if (this%numpoints/=rhs%numpoints) then
        iseq = .false.
        return
      end if
      if (this%extra_dims>rhs%extra_dims) then
        dims = rhs%extra_dims + 1
        iseq = all(this%field_data(:,this%dimensions()+dims:)==0.0_r8)
        if (.not. iseq) return
      else if (this%extra_dims<rhs%extra_dims) then
        dims = this%extra_dims + 1
        iseq = all(rhs%field_data(:,rhs%dimensions()+dims:)==0.0_r8)
        if (.not. iseq) return
      else
        dims = this%extra_dims + 1
      end if
      do i=1,this%numpoints+1
        normalization = norm2(this%field_data(i,:dims))
        if (normalization < tiny(normalization)) normalization = 1.0_r8
        iseq = iseq .and. ( (norm2(this%field_data(i,this%dimensions()+1:dims) &
             - rhs%field_data(i,rhs%dimensions()+1:dims))/normalization < get_tol()) &
             .or. (isnan(norm2(this%field_data(i,this%dimensions()+1:dims))) .and. &
                   isnan(norm2(rhs%field_data(i,rhs%dimensions()+1:dims)))) )
        if (.not. iseq) return
      end do
    class default
      iseq = .false.
    end select
  end function cheb1d_vector_is_equal
  
  subroutine cheb1d_vector_check_compatible(this,other)
    !* Author: Chris MacMackin
    !  Date: April 2016
    !
    ! Checks whether a field has the same type, boundaries, and
    ! resolution as this one, making it compatible for binary 
    ! operations. If incompatible, stops program with an error message.
    !
    class(cheb1d_vector_field), intent(in) :: this
    class(abstract_field), intent(in) :: other
      !! The field being checked against this one
    character(len=128) :: err_message
    logical :: type_err, domain_err, res_err, alloc_err, has_message
    has_message = .false.
    select type(other)
    class is(cheb1d_scalar_field)
      type_err = .false.
      domain_err = any(abs(this%extent - other%extent) > 1.e-15_r8)
      res_err = this%numpoints /= other%numpoints
      alloc_err = .not. (allocated(this%field_data).and.allocated(other%field_data))
    class is(cheb1d_vector_field)
      type_err = .false.
      domain_err = any(abs(this%extent - other%extent) > 1.e-15_r8)
      res_err = this%numpoints /= other%numpoints
      alloc_err = .not. (allocated(this%field_data).and.allocated(other%field_data))
    class default
      type_err = .true.
      domain_err = .false.
      res_err = .false.
      alloc_err = .false.
    end select
    err_message = ''
    if (type_err) then
      err_message = ' incompatible types'
      has_message = .true.
    end if
    if (domain_err) then
      if (has_message) err_message = trim(err_message)//','
      err_message = trim(err_message)//' different domains'
      has_message = .true.
    end if
    if (res_err) then
      if (has_message) err_message = trim(err_message)//','
      err_message = trim(err_message)//' different resolutions'
      has_message = .true.
    end if
    if (alloc_err) then
      if (has_message) err_message = trim(err_message)//','
      err_message = trim(err_message)//' uninitialized field'
      has_message = .true.
    end if
    if (has_message) then
      write(stderr,*) 'cheb1d_vector_field: Error, operation with incompatible fields'
      write(stderr,*) '    Following inconsistencies: '//err_message
#ifdef __GFORTRAN__
      call backtrace
#endif
      error stop
    end if
  end subroutine cheb1d_vector_check_compatible

  pure subroutine cheb1d_vector_assign_meta_data(this, rhs, alloc)
    !* Author: Chris MacMackin
    !  Date: April 2016
    !
    ! Assigns all data in `rhs` to `this` except for the actual values
    ! of the field.
    !
    class(cheb1d_vector_field), intent(inout) :: this
    class(abstract_field), intent(in) :: rhs
    logical, optional, intent(in) :: alloc
      !! If present and false, do not allocate the array of `this`.
   select type(rhs)
    class is(cheb1d_scalar_field)
      this%numpoints = rhs%numpoints
      this%extent = rhs%extent
      this%extra_dims = 0
      if (allocated(rhs%colloc_points)) then
        if (.not. allocated(this%colloc_points)) &
          allocate(this%colloc_points(this%numpoints + 1))
        this%colloc_points = rhs%colloc_points
      else if (allocated(this%colloc_points)) then
        deallocate(this%colloc_points)
      end if
      if (allocated(this%field_data)) deallocate(this%field_data)
      if (present(alloc)) then
      !Needed due to compiler bug:
        if (allocated(rhs%field_data) .and. alloc) &
          allocate(this%field_data(size(rhs%field_data),1))
      else
        allocate(this%field_data(size(rhs%field_data),1))
      end if
    class is(cheb1d_vector_field)
      this%numpoints = rhs%numpoints
      this%extent = rhs%extent
      this%extra_dims = rhs%extra_dims
      if (allocated(rhs%colloc_points)) then
        if (.not. allocated(this%colloc_points)) &
          allocate(this%colloc_points(this%numpoints + 1))
        this%colloc_points = rhs%colloc_points
      else if (allocated(this%colloc_points)) then
        deallocate(this%colloc_points)
      end if
      if (allocated(this%field_data)) deallocate(this%field_data)
      if (present(alloc)) then
        !Needed due to compiler bug:
        if (allocated(rhs%field_data) .and. alloc) &
          allocate(this%field_data(size(rhs%field_data,1),size(rhs%field_data,2)))
      else
        allocate(this%field_data(size(rhs%field_data,1),size(rhs%field_data,2)))
      end if
    end select
  end subroutine cheb1d_vector_assign_meta_data

end module cheb1d_fields_mod
