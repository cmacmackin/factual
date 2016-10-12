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
  use array_fields_mod, only: array_scalar_field, array_vector_field
  use chebyshev_mod
  implicit none
  private

  public :: sin, cos, tan, asin, acos, atan, sinh, cosh, tanh, asinh, &
            acosh, atanh, log, log10, exp, abs, sqrt, minval, maxval

  type, extends(array_scalar_field), public :: cheb1d_scalar_field
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
    real(r8), dimension(:,:), allocatable :: extent
      !! The start and end values of the domain
    real(r8), dimension(:), allocatable :: colloc_points
      !! The location of the data-points
  contains
    private
    procedure, public :: dimensions => cheb1d_scalar_dimensions
      !! Returns the number of dimensions of the field
    procedure, public :: domain => cheb1d_scalar_domain
      !! Provides array with upper and lower limits of field's domain
    procedure, public :: resolution => cheb1d_scalar_resolution
      !! Returns array containing number of datapoints in each
      !! dimension.
    procedure :: assign_subtype_meta_data => cheb1d_scalar_assign_meta
      !! Copies all data unique to this field type from another field
      !! object to this one.
    procedure :: check_compatible => cheb1d_scalar_check_compatible
      !! Tests whether two fields are suitable for binary operations
      !! together
    procedure :: array_dx => cheb1d_scalar_array_dx
      !! Takes the derivative of the field using a 1-D array of data
      !! passed as an argument.
    procedure, public :: allocate_scalar_field => cheb1d_scalar_allocate_scalar
      !! Allocates a scalar field with concrete type [[cheb1d_scalar_field]]
    procedure, public :: allocate_vector_field => cheb1d_scalar_allocate_vector
      !! Allocates a vector field with concrete type [[cheb1d_vector_field]]
  end type cheb1d_scalar_field
  
  interface cheb1d_scalar_field
    module procedure scalar_constructor
  end interface

  type, extends(array_vector_field), public :: cheb1d_vector_field
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
    real(r8), dimension(:,:), allocatable :: extent
      !! The start and end values of the domain
    real(r8), dimension(:), allocatable :: colloc_points
      !! The location of the data-points
  contains
    private
    procedure, public :: dimensions => cheb1d_vector_dimensions
      !! Returns number of dimensions of the field
    procedure, public :: domain => cheb1d_vector_domain
      !! Provides array with upper and lower limits of field's domain
    procedure, public :: resolution => cheb1d_vector_resolution
      !! Returns array containing number of datapoints in each dimension.
    procedure :: assign_subtype_meta_data => cheb1d_vector_assign_meta
      !! Copies all data other than values stored in field from another
      !! field object to this one.
    procedure :: check_compatible => cheb1d_vector_check_compatible
      !! Tests whether two fields are suitable for binary operations
      !! together
    procedure :: array_component_dx => cheb1d_vector_component_array_dx
      !! Takes the derivative of particular vector component of the
      !! field, using a 1-D array of data passed to it.
    procedure :: array_dx => cheb1d_vector_array_dx
      !! Takes the derivative of the field using a 1-D array of data
      !! passed as an argument.
    procedure, public :: allocate_scalar_field => cheb1d_vector_allocate_scalar
      !! Allocates a scalar field with concrete type [[cheb1d_scalar_field]]
    procedure, public :: allocate_vector_field => cheb1d_vector_allocate_vector
      !! Allocates a vector field with concrete type [[cheb1d_vector_field]]
  end type cheb1d_vector_field
  
  interface cheb1d_vector_field
    module procedure vector_constructor
  end interface
  
  abstract interface
    pure function scalar_init(x) result(scalar)
      !! Function used to specify value held by a scalar field at
      !! location `x`.
      import :: r8
      real(r8), dimension(:), intent(in) :: x 
        !! The position at which this function is evaluated
      real(r8) :: scalar
        !! The value of the field at this location
    end function scalar_init
  
    pure function vector_init(x) result(vector)
      !! Function used to specify value held by a vector field at
      !! location `x`.
      import :: r8
      real(r8), dimension(:), intent(in) :: x 
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
      forall (i=1:nodes+1) field%field_data(i) = initializer([field%colloc_points(i)])
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
    res(1) = this%elements()
  end function cheb1d_scalar_resolution

  function cheb1d_scalar_array_dx(this, data_array, dir, order) result(res)
    !* Author: Chris MacMackin
    !  Date: October 2016
    !
    ! \(\frac{\partial^{\rm order}}{\partial x_{\rm dir}^{\rm order}}{\rm field}\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    integer, intent(in) :: dir
      !! Direction in which to differentiate
    real(r8), dimension(:), intent(in)     :: data_array
      !! An array holding the datapoints for this field, identical
      !! in layout to that stored the field itself.
    integer, optional, intent(in) :: order
      !! Order of the derivative, default = 1
    real(r8), dimension(:), allocatable   :: res
      !! The spatial derivative of order `order` taken in direction `dir`
    if (dir==1) then
      res = data_array
      call differentiate_1d(res,this%colloc_points,order)
    else
      res = 0.0
    end if
  end function cheb1d_scalar_d_dx

  pure subroutine cheb1d_scalar_check_compatible(this,other)
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
      res_err = this%elements() /= other%elements()
      alloc_err = .not. (allocated(this%field_data).and.allocated(other%field_data))
    class is(cheb1d_vector_field)
      type_err = .false.
      domain_err = any(abs(this%extent - other%extent) > 1.e-15_r8)
      res_err = this%elements() /= other%elements()
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
      error stop('cheb1d_scalar_field: Error, operation with incompatible '//&
                 'fields')
    end if
  end subroutine cheb1d_scalar_check_compatible

  pure subroutine cheb1d_scalar_assign_meta(this, rhs, alloc)
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
      this%extent = rhs%extent
      if (allocated(rhs%colloc_points)) then
        if (.not. allocated(this%colloc_points)) &
          allocate(this%colloc_points(this%elements() + 1))
        this%colloc_points = rhs%colloc_points
      else if (allocated(this%colloc_points)) then
        deallocate(this%colloc_points)
      end if
    class is(cheb1d_vector_field)
      this%extent = rhs%extent
      if (allocated(rhs%colloc_points)) then
        if (.not. allocated(this%colloc_points)) &
          allocate(this%colloc_points(this%elements() + 1))
        this%colloc_points = rhs%colloc_points
      else if (allocated(this%colloc_points)) then
        deallocate(this%colloc_points)
      end if
    end select
  end subroutine cheb1d_scalar_assign_meta

  !=====================================================================
  ! Vector Field Methods
  !=====================================================================

  function vector_constructor(nodes,initializer,lower_bound, &
                              upper_bound, extra_dims) result(field)
    integer, intent(in) :: nodes
      !! The number of collocation nodes to use when modelling this
      !! field. This corresponds to resolution.
    procedure(vector_init), optional :: initializer
      !! An elemental procedure which takes the position in the
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
    integer :: i, dims, initializer_dims
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
      initializer_dims = size(initializer([field%extent(1,1)]))
      if (initializer_dims < dims) then
        tmp = [(0.0, i=1,(1+dims-initializer_dims))]
        forall (i=1:nodes+1) &
          field%field_data(i,:) = [initializer([field%colloc_points(i)]),&
                                   tmp]
      else if (initializer_dims > dims) then
        do i=1,nodes+1
          tmp = initializer([field%colloc_points(i)])
          field%field_data(i,:) = tmp(1:dims)
        end do
      else
        forall (i=1:nodes+1) &
          field%field_data(i,:) = initializer([field%colloc_points(i)])
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
    res(1) = this%elements()
  end function cheb1d_vector_resolution
  
  function cheb1d_vector_array_dx(this, data_array, dir, order) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\frac{\partial^{\rm order}}{\partial x_{\rm dir}^{\rm order}}\vec{\rm field}\)
    !
    class(cheb1d_vector_field), intent(in) :: this
      real(r8), dimension(:,:), intent(in)  :: data_array
        !! An array holding the datapoints for the vectors in this
        !! field, with identical in layout to the storage in the field
        !! itself. Each row represents a different spatial location,
        !! while each column represents a different component of the
        !! vector.
    integer, intent(in) :: dir
      !! Direction in which to differentiate
    integer, optional, intent(in) :: order
      !! Order of the derivative, default = 1
    real(r8), dimension(:,:), allocatable :: res
      !! The spatial derivative of order `order` taken in direction `dir`
    integer :: i
    if (dir==1) then
      res = data_array
      do i=1,local%vector_dimensions()
        call differentiate_1d(res(:,i),this%colloc_points,order)
      end do
    else
      res = 0.0_r8
    end if
  end function cheb1d_vector_array_dx

  function cheb1d_vector_component_array_dx(this, data_array, dir, order) &
                                                              result(res)
    !* Author: Chris MacMackin
    !  Date: October 2016
    !
    ! \(\frac{\partial^{\rm order}}{\partial x_{\rm dir}^{\rm order}}{\rm field}\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    integer, intent(in) :: dir
      !! Direction in which to differentiate
    real(r8), dimension(:), intent(in)     :: data_array
      !! An array holding the datapoints for this field, identical
      !! in layout to that stored the field itself.
    integer, optional, intent(in) :: order
      !! Order of the derivative, default = 1
    real(r8), dimension(:), allocatable   :: res
      !! The spatial derivative of order `order` taken in direction `dir`
    if (dir==1) then
      res = data_array
      call differentiate_1d(res,this%colloc_points,order)
    else
      res = 0.0
    end if
  end function cheb1d_vector_component_array_dx

  pure subroutine cheb1d_vector_check_compatible(this,other)
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
      res_err = this%elements() /= other%elements()
      alloc_err = .not. (allocated(this%field_data).and.allocated(other%field_data))
    class is(cheb1d_vector_field)
      type_err = .false.
      domain_err = any(abs(this%extent - other%extent) > 1.e-15_r8)
      res_err = this%elements() /= other%elements()
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
      error stop('cheb1d_vector_field: Error, operation with incompatible '//&
                 'fields')
    end if
  end subroutine cheb1d_vector_check_compatible

  pure subroutine cheb1d_vector_assign_meta(this, rhs, alloc)
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
      this%extent = rhs%extent
      if (allocated(rhs%colloc_points)) then
        if (.not. allocated(this%colloc_points)) &
          allocate(this%colloc_points(this%elements()))
        this%colloc_points = rhs%colloc_points
      else if (allocated(this%colloc_points)) then
        deallocate(this%colloc_points)
      end if
    class is(cheb1d_vector_field)
      this%extent = rhs%extent
      if (allocated(rhs%colloc_points)) then
        if (.not. allocated(this%colloc_points)) &
          allocate(this%colloc_points(this%elements()))
        this%colloc_points = rhs%colloc_points
      else if (allocated(this%colloc_points)) then
        deallocate(this%colloc_points)
      end if
    end select
  end subroutine cheb1d_vector_assign_meta

end module cheb1d_fields_mod
