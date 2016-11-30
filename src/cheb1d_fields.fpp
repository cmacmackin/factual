!
!  cheb1d_fields.f90
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

! Make procedures non-pure for debugging, so that messages can be
! printed to the screen.
#:if defined('DEBUG')
#define pure 
#define elemental 
#:endif

#:include 'fypp_utils.fpp'

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
  use iso_fortran_env, only: r8 => real64
  use abstract_fields_mod
  use uniform_fields_mod, only: uniform_scalar_field, &
                                uniform_vector_field
  use array_fields_mod, only: array_scalar_field, array_vector_field, &
                              scalar_init, vector_init
  use chebyshev_mod
  use h5lt, only: hid_t, size_t, h5ltmake_dataset_double_f, &
                  h5ltset_attribute_string_f, h5ltset_attribute_int_f
  implicit none
  private

  character(len=19), parameter, public :: hdf_scalar_name = 'cheb1d_scalar_field'
  character(len=19), parameter, public :: hdf_vector_name = 'cheb1d_vector_field'

$:public_unary()
  public :: minval
  public :: maxval
  
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
    real(r8), dimension(2)              :: extent
      !! The start and end values of the domain
    real(r8), dimension(:), allocatable :: colloc_points
      !! The location of the data-points
    logical                             :: differentiable = .true.
      !! Indicates whether this is a complete Chebyshev field which
      !! can be differentiated using the Psuedo-spectral method.
  contains
    private
    procedure, public :: dimensions => cheb1d_scalar_dimensions
      !! Returns the number of dimensions of the field
    procedure, public :: domain => cheb1d_scalar_domain
      !! Provides array with upper and lower limits of field's domain
    procedure, public :: resolution => cheb1d_scalar_resolution
    !! Returns array containing number of datapoints in each
      !! dimension.
    procedure :: raw_slices => cheb1d_scalar_raw_slices
      !! Returns an array of integers used for getting the correct
      !! data for a given raw representation of the field.
    procedure, public :: id_to_position => cheb1d_scalar_id_to_position
      !! Given the ID number of a location in the field, returns the
      !! coordinates of that position
    procedure :: assign_subtype_meta_data => cheb1d_scalar_assign_meta
      !! Copies all data unique to this field type from another field
      !! object to this one.
    procedure :: check_subtype_compatible => cheb1d_scalar_check_compatible
      !! Tests whether two fields are suitable for binary operations
      !! together
    procedure :: array_dx => cheb1d_scalar_array_dx
      !! Takes the derivative of the field using a 1-D array of data
      !! passed as an argument.
    procedure, public :: allocate_scalar_field => cheb1d_scalar_allocate_scalar
      !! Allocates a scalar field with concrete type [[cheb1d_scalar_field]]
    procedure, public :: allocate_vector_field => cheb1d_scalar_allocate_vector
      !! Allocates a vector field with concrete type [[cheb1d_vector_field]]
    procedure, public :: subtype_boundary => cheb1d_scalar_bound
      !! Performs whatever operations are needed on the field to get
      !! a boundary field, including returning the slices needed to
      !! extract the appropriate data.
    procedure, public :: write_hdf => cheb1d_scalar_write_hdf
      !! Write field data to a new dataset in an HDF5 file.
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
    real(r8), dimension(2)              :: extent
      !! The start and end values of the domain
    real(r8), dimension(:), allocatable :: colloc_points
      !! The location of the data-points
    logical                             :: differentiable = .true.
      !! Indicates whether this is a complete Chebyshev field which
      !! can be differentiated using the Psuedo-spectral method.
  contains
    private
    procedure, public :: dimensions => cheb1d_vector_dimensions
      !! Returns number of dimensions of the field
    procedure, public :: domain => cheb1d_vector_domain
      !! Provides array with upper and lower limits of field's domain
    procedure, public :: resolution => cheb1d_vector_resolution
      !! Returns array containing number of datapoints in each dimension.
    procedure :: raw_slices => cheb1d_vector_raw_slices
      !! Returns an array of integers used for getting the correct
      !! data for a given raw representation of the field.
    procedure, public :: id_to_position => cheb1d_vector_id_to_position
      !! Given the ID number of a location in the field, returns the
      !! coordinates of that position
    procedure :: assign_subtype_meta_data => cheb1d_vector_assign_meta
      !! Copies all data other than values stored in field from another
      !! field object to this one.
    procedure :: check_subtype_compatible => cheb1d_vector_check_compatible
      !! Tests whether two fields are suitable for binary operations
      !! together
    procedure :: array_dx => cheb1d_vector_array_dx
      !! Takes the derivative of particular vector component of the
      !! field, using a 1-D array of data passed to it.
    procedure, public :: allocate_scalar_field => cheb1d_vector_allocate_scalar
      !! Allocates a scalar field with concrete type [[cheb1d_scalar_field]]
    procedure, public :: allocate_vector_field => cheb1d_vector_allocate_vector
      !! Allocates a vector field with concrete type [[cheb1d_vector_field]]
    procedure, public :: subtype_boundary => cheb1d_vector_bound
      !! Performs whatever operations are needed on the field to get
      !! a boundary field, including returning the slices needed to
      !! extract the appropriate data.
    procedure, public :: write_hdf => cheb1d_vector_write_hdf
      !! Write field data to a new dataset in an HDF5 file.
  end type cheb1d_vector_field
  
  interface cheb1d_vector_field
    module procedure vector_constructor
  end interface cheb1d_vector_field
 
contains

  !=====================================================================
  ! Scalar Field Methods
  !=====================================================================

  function scalar_constructor(numpoints,initializer,lower_bound, &
                              upper_bound) result(field)
    integer, intent(in) :: numpoints
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
    if (present(lower_bound)) then
      field%extent(1) = lower_bound
    else
      field%extent(1) = -1.0_r8
    end if
    if (present(upper_bound)) then
      field%extent(2) = upper_bound
    else
      field%extent(2) = 1.0_r8
    end if
    field%colloc_points = collocation_points(numpoints-1,lower_bound,upper_bound)
    field = array_scalar_field(field,numpoints,initializer)
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
    allocate(res(1,2))
    res(1,:) = this%extent
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

  pure function cheb1d_scalar_raw_slices(this,exclude_lower_bound,exclude_upper_bound) &
                                                                        result(slices)
    class(cheb1d_scalar_field), intent(in)      :: this
    integer, dimension(:), optional, intent(in) :: exclude_lower_bound
      !! Specifies how many layers of data points should be excluded
      !! from the result at the lower boundary for each
      !! dimension. The number in element `n` of the array indicates
      !! how many layers of cells at the lower boundary normal to
      !! dimension `n` will be ignored. Defaults to 0 for all.
    integer, dimension(:), optional, intent(in) :: exclude_upper_bound
      !! Specifies how many layers of data points should be excluded
      !! from the result at the upper boundary for each
      !! dimension. The number in element `n` of the array indicates
      !! how many layers of cells at the upper boundary normal to
      !! dimension `n` will be ignored. Defaults to 0 for all.
    integer, dimension(:,:), allocatable        :: slices
      !! An array containing array slice data which can be used to
      !! construct the raw representation of a field, with the given
      !! boundary conditions. The form of the array is
      !! ```
      !! slices(1,i) = start_index
      !! slices(2,i) = end_index
      !! slices(3,i) = stride
      !! ```
    allocate(slices(3,1))
    ! Remember that Chebyshev collocation nodes end up in reverse order
    if (present(exclude_upper_bound)) then
      slices(1,1) = 1 + exclude_upper_bound(1)
    else
      slices(1,1) = 1
    end if
    if (present(exclude_lower_bound)) then
      slices(2,1) = this%elements() - exclude_lower_bound(1)
    else
      slices(2,1) = this%elements()
    end if
    slices(3,1) = 1
  end function cheb1d_scalar_raw_slices

  pure function cheb1d_scalar_id_to_position(this, id) result(pos)
    !* Author: Chris MacMackin
    !  Date: October 2016
    !
    ! From an ID number (representing the index of the element in this
    ! field's 1-D storage array), the coordinates of a location in the
    ! field with that ID are returned.
    !
    class(cheb1d_scalar_field), intent(in) :: this
    integer, intent(in)                    :: id
      !! The ID number for some location in the field
    real(r8), dimension(:), allocatable    :: pos
      !! The coordinates for this location in the field
    allocate(pos(1))
    if (id <= this%elements()) then
      pos(1) = this%colloc_points(id)
    else
      error stop ('cheb1d_scalar_field: Invalid ID number provided.')
    end if
  end function cheb1d_scalar_id_to_position

  function cheb1d_scalar_array_dx(this, data_array, dir, order) result(res)
    !* Author: Chris MacMackin
    !  Date: October 2016
    !
    ! \(\frac{\partial^{\rm order}}{\partial x_{\rm dir}^{\rm order}}{\rm field}\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    real(r8), dimension(:), intent(in)     :: data_array
      !! An array holding the datapoints for this field, identical
      !! in layout to that stored the field itself.
    integer, intent(in) :: dir
      !! Direction in which to differentiate
    integer, optional, intent(in) :: order
      !! Order of the derivative, default = 1
    real(r8), dimension(:), allocatable   :: res
      !! The spatial derivative of order `order` taken in direction `dir`
    integer :: i
    if (dir==1) then
#:if defined('DEBUG')
      if (.not. this%differentiable) &
        error stop ('Trying to take derivative of undifferentiable instance of '// &
                   '`cheb1d_scalar_field`.')
#:endif
      res = data_array
      call differentiate_1d(res,this%colloc_points,order)
    else
      allocate(res(size(data_array)))
      res = 0.0_r8
    end if
  end function cheb1d_scalar_array_dx

  pure subroutine cheb1d_scalar_check_compatible(this,other)
    !* Author: Chris MacMackin
    !  Date: April 2016
    !
    ! Checks whether a field has the same type, and domain as this
    ! one, making it compatible for binary operations. If
    ! incompatible, stops program with an error message.
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(abstract_field), intent(in) :: other
      !! The field being checked against this one
    character(len=70), parameter :: err_message = 'cheb1d_scalar_field: '//&
         'Error, operation with incompatible fields due to '
    select type(other)
    class is(cheb1d_scalar_field)
      if (any(abs(this%extent - other%extent) > 1.e-15_r8)) &
           error stop (err_message//'different domains.')
    class is(cheb1d_vector_field)
      if (any(abs(this%extent - other%extent) > 1.e-15_r8)) &
           error stop (err_message//'different domains.')
    class is(uniform_scalar_field)
      continue
    class is(uniform_vector_field)
      continue
    class default
      error stop (err_message//'incompatible types.')
    end select
  end subroutine cheb1d_scalar_check_compatible

  pure subroutine cheb1d_scalar_assign_meta(this, rhs)
    !* Author: Chris MacMackin
    !  Date: April 2016
    !
    ! Assigns all data in `rhs` to `this` except for the actual values
    ! of the field.
    !
    class(cheb1d_scalar_field), intent(inout) :: this
    class(abstract_field), intent(in) :: rhs
      !! If present and false, do not allocate the array of `this`.
    select type(rhs)
    class is(cheb1d_scalar_field)
      this%extent = rhs%extent
      if (allocated(rhs%colloc_points)) then
      !  if (.not. allocated(this%colloc_points)) &
       !   allocate(this%colloc_points(this%elements() + 1))
        this%colloc_points = rhs%colloc_points
      else if (allocated(this%colloc_points)) then
        deallocate(this%colloc_points)
      end if
    class is(cheb1d_vector_field)
      this%extent = rhs%extent
      if (allocated(rhs%colloc_points)) then
        !if (.not. allocated(this%colloc_points)) &
         ! allocate(this%colloc_points(this%elements() + 1))
        this%colloc_points = rhs%colloc_points
      else if (allocated(this%colloc_points)) then
        deallocate(this%colloc_points)
      end if
    end select
  end subroutine cheb1d_scalar_assign_meta

  pure subroutine cheb1d_scalar_allocate_scalar(this, new_field)
    !* Author: Chris MacMackin
    !  Date: October 2016
    !
    ! Allocates an abstract [[scalar_field(type)]] to have a type
    ! compatible for operations on this field and to be the type
    ! returned by this field's methods which produce scalars.
    !
    class(cheb1d_scalar_field), intent(in)          :: this
    class(scalar_field), allocatable, intent(inout) :: new_field
      !! A field which, upon return, is allocated to be of the same
      !! concrete type as scalar fields produced by `this`.
    if (allocated(new_field)) deallocate(new_field)
    allocate(cheb1d_scalar_field :: new_field)
  end subroutine cheb1d_scalar_allocate_scalar

  pure subroutine cheb1d_scalar_allocate_vector(this, new_field)
    !* Author: Chris MacMackin
    !  Date: October 2016
    !
    ! Allocates an abstract [[vector_field(type)]] to have a type
    ! compatible for operations on this field and to be the type
    ! returned by this field's methods which produce vectors.
    !
    class(cheb1d_scalar_field), intent(in)          :: this
    class(vector_field), allocatable, intent(inout) :: new_field
      !! A field which, upon return, is allocated to be of the same
      !! concrete type as vector fields produced by `this`.
    if (allocated(new_field)) deallocate(new_field)
    allocate(cheb1d_vector_field :: new_field)
  end subroutine cheb1d_scalar_allocate_vector

  pure subroutine cheb1d_scalar_bound(this,src,boundary,depth,slices)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Alters the meta-data of the field to be appropriate for the
    ! specified boundary. Also computes the slices needed to access
    ! the field contents for the specified boundary.
    !
    class(cheb1d_scalar_field), intent(out)           :: this
    class(array_scalar_field), intent(in)             :: src
      !! The field for which the boundary data is to be provided.
    integer, intent(in)                               :: boundary
      !! Specifies which boundary is to be returned. The boundary
      !! will be the one normal to dimension of number
      !! `abs(boundary)`. If the argument is negative, then the
      !! lower boundary is returned. If positive, then the upper
      !! boundary is returned.
    integer, intent(in)                               :: depth
      !! The number of layers of data-points to return at the
      !! specified boundary.
    integer, dimension(:,:), allocatable, intent(out) :: slices
      !! An array which, on return, contains array slice data which
      !! can be used to construct extract the data for the given field
      !! boundary. The form of the array is
      !! ```
      !! slices(1,i) = start_index
      !! slices(2,i) = end_index
      !! slices(3,i) = stride
      !! ```
    integer :: length
    select type(src)
    class is(cheb1d_scalar_field)
      length = size(src%colloc_points)
      select case(boundary)
      case(-1)
         
        this%colloc_points = src%colloc_points(length+1-depth:)
        this%extent = [src%colloc_points(length+1-depth), &
                      src%colloc_points(length)]
        this%differentiable = .false.
        allocate(slices(3,1))
        slices(:,1) = [length+1-depth,length,1]
      case(1)
        this%colloc_points = src%colloc_points(1:depth)
        this%extent = [src%colloc_points(depth), src%colloc_points(1)]
        this%differentiable = .false.
        allocate(slices(3,1))
        slices(:,1) = [1,depth,1]
      case default
        this%colloc_points = src%colloc_points
        this%extent = src%extent
        allocate(slices(3,1))
        slices(:,1) = [1,length,1]
      end select
    class default
      error stop ('Trying to compute boundary metadata for type other '// &
                 'than `cheb1d_scalar_field`.')
    end select
  end subroutine cheb1d_scalar_bound

  subroutine cheb1d_scalar_field_write_hdf(this, hdf_id, dataset_name, &
                                           error)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Writes the contents of the field to a dataset in an HDF
    ! file. The dataset will have an attribute specifying the name of
    ! the field type and an attribute with a zero value, indicating
    ! that this is not a vector field.
    !
    class(cheb1d_scalar_field), intent(in) :: this
    integer(hid_t), intent(in)             :: hdf_id
      !! The identifier for the HDF file/group in which the field
      !! data is to be written.
    character(len=*), intent(in)           :: dataset_name
      !! The name of the dataset to be created in the HDF file
      !! containing this field's data.
    integer, intent(out)                   :: error
      !! An error code which, upon succesful completion of the
      !! routine, is 0. Otherwise, contains the error code returned
      !! by the HDF library.
    error = 0
    call this%write_hdf_array(hdf_id, dataset_name, [this%numpoints], error)
    if (error /= 0) return
    call h5ltset_attribute_string_f(hdf_id, dataset_name, hdf_field_type_attr, &
                                    hdf_scalar_type, error)
    if (error /= 0) return
    call h5ltset_attribute_int_f(hdf_id, dataset_name, hdf_vector_attr, [0], &
                                 1_size_t, error)
    if (error /= 0) return
    call h5ltset_attribute_double_f(hdf_id, dataset_name, hdf_grid_attr//'1', &
                                    int(this%colloc_points,size_t), error)
  end subroutine cheb1d_scalar_field_write_hdf


  !=====================================================================
  ! Vector Field Methods
  !=====================================================================


  function vector_constructor(numpoints,initializer,lower_bound, &
                              upper_bound, extra_dims) result(field)
    integer, intent(in) :: numpoints
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
    integer :: dims
    dims = 1
    if (present(extra_dims)) then
      if (extra_dims > 0) then
        dims = dims + extra_dims
      end if
    end if
    if (present(lower_bound)) then
       field%extent(1) = lower_bound
    else
       field%extent(1) = -1.0_r8
    end if
    if (present(upper_bound)) then
       field%extent(2) = upper_bound
    else
       field%extent(2) = 1.0_r8
    end if
    field%colloc_points = collocation_points(numpoints-1,lower_bound,upper_bound)
    field = array_vector_field(field,numpoints,dims,initializer)
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
    allocate(res(1,2))
    res(1,:) = this%extent
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

  pure function cheb1d_vector_raw_slices(this,exclude_lower_bound,exclude_upper_bound) &
                                                                        result(slices)
    class(cheb1d_vector_field), intent(in)      :: this
    integer, dimension(:), optional, intent(in) :: exclude_lower_bound
      !! Specifies how many layers of data points should be excluded
      !! from the result at the lower boundary for each
      !! dimension. The number in element `n` of the array indicates
      !! how many layers of cells at the lower boundary normal to
      !! dimension `n` will be ignored. Defaults to 0 for all.
    integer, dimension(:), optional, intent(in) :: exclude_upper_bound
      !! Specifies how many layers of data points should be excluded
      !! from the result at the upper boundary for each
      !! dimension. The number in element `n` of the array indicates
      !! how many layers of cells at the upper boundary normal to
      !! dimension `n` will be ignored. Defaults to 0 for all.
    integer, dimension(:,:), allocatable        :: slices
      !! An array containing array slice data which can be used to
      !! construct the raw representation of a field, with the given
      !! boundary conditions. The form of the array is
      !! ```
      !! slices(1,i) = start_index
      !! slices(2,i) = end_index
      !! slices(3,i) = stride
      !! ```
    allocate(slices(3,1))
    ! Remember that Chebyshev collocation nodes end up in reverse order
    if (present(exclude_upper_bound)) then
      slices(1,1) = 1 + exclude_upper_bound(1)
    else
      slices(1,1) = 1
    end if
    if (present(exclude_lower_bound)) then
      slices(2,1) = this%elements() - exclude_lower_bound(1)
    else
      slices(2,1) = this%elements()
    end if
    slices(3,1) = 1
  end function cheb1d_vector_raw_slices
  
  pure function cheb1d_vector_id_to_position(this, id) result(pos)
    !* Author: Chris MacMackin
    !  Date: October 2016
    !
    ! From an ID number (representing the index of the element in this
    ! field's 1-D storage array), the coordinates of a location in the
    ! field with that ID are returned.
    !
    class(cheb1d_vector_field), intent(in) :: this
    integer, intent(in)                    :: id
    !! The ID number for some location in the field
    real(r8), dimension(:), allocatable    :: pos
    !! The coordinates for this location in the field
    allocate(pos(1))
    if (id <= this%elements()) then
      pos(1) = this%colloc_points(id)
    else
      error stop ('cheb1d_vector_field: Invalid ID number provided.')
    end if
  end function cheb1d_vector_id_to_position

  function cheb1d_vector_array_dx(this, data_array, dir, order) &
                                                              result(res)
    !* Author: Chris MacMackin
    !  Date: October 2016
    !
    ! \(\frac{\partial^{\rm order}}{\partial x_{\rm dir}^{\rm order}}{\rm field}\)
    !
    class(cheb1d_vector_field), intent(in) :: this
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
#:if defined('DEBUG')
      if (.not. this%differentiable) &
        error stop ('Trying to take derivative of undifferentiable instance of '// &
                   '`cheb1d_vector_field`.')
#:endif
      res = data_array
      call differentiate_1d(res,this%colloc_points,order)
    else
      allocate(res, mold=data_array)
      res = 0.0
    end if
  end function cheb1d_vector_array_dx

  pure subroutine cheb1d_vector_check_compatible(this,other)
    !* Author: Chris MacMackin
    !  Date: April 2016
    !
    ! Checks whether a field has the same type, and domain as this
    ! one, making it compatible for binary operations. If
    ! incompatible, stops program with an error message.
    !
    class(cheb1d_vector_field), intent(in) :: this
    class(abstract_field), intent(in) :: other
      !! The field being checked against this one
    character(len=70), parameter :: err_message = 'cheb1d_vector_field: '//&
         'Error, operation with incompatible fields due to '
    select type(other)
    class is(cheb1d_scalar_field)
      if (any(abs(this%extent - other%extent) > 1.e-15_r8)) &
           error stop (err_message//'different domains.')
    class is(cheb1d_vector_field)
      if (any(abs(this%extent - other%extent) > 1.e-15_r8)) &
           error stop (err_message//'different domains.')
    class is(uniform_scalar_field)
      continue
    class is(uniform_vector_field)
      continue
    class default
      error stop (err_message//'incompatible types.')
    end select
  end subroutine cheb1d_vector_check_compatible

  pure subroutine cheb1d_vector_assign_meta(this, rhs)
    !* Author: Chris MacMackin
    !  Date: April 2016
    !
    ! Assigns all data in `rhs` to `this` except for the actual values
    ! of the field.
    !
    class(cheb1d_vector_field), intent(inout) :: this
    class(abstract_field), intent(in) :: rhs
      !! If present and false, do not allocate the array of `this`.
    select type(rhs)
    class is(cheb1d_scalar_field)
      this%extent = rhs%extent
      if (allocated(rhs%colloc_points)) then
        !if (.not. allocated(this%colloc_points)) &
        !  allocate(this%colloc_points(this%elements()))
        this%colloc_points = rhs%colloc_points
      else if (allocated(this%colloc_points)) then
        deallocate(this%colloc_points)
      end if
    class is(cheb1d_vector_field)
      this%extent = rhs%extent
      if (allocated(rhs%colloc_points)) then
        !if (.not. allocated(this%colloc_points)) &
        !  allocate(this%colloc_points(this%elements()))
        this%colloc_points = rhs%colloc_points
      else if (allocated(this%colloc_points)) then
        deallocate(this%colloc_points)
      end if
    end select
  end subroutine cheb1d_vector_assign_meta

  pure subroutine cheb1d_vector_allocate_scalar(this, new_field)
    !* Author: Chris MacMackin
    !  Date: October 2016
    !
    ! Allocates an abstract [[scalar_field(type)]] to have a type
    ! compatible for operations on this field and to be the type
    ! returned by this field's methods which produce scalars.
    !
    class(cheb1d_vector_field), intent(in)          :: this
    class(scalar_field), allocatable, intent(inout) :: new_field
      !! A field which, upon return, is allocated to be of the same
      !! concrete type as scalar fields produced by `this`.
    if (allocated(new_field)) deallocate(new_field)
    allocate(cheb1d_scalar_field :: new_field)
  end subroutine cheb1d_vector_allocate_scalar

  pure subroutine cheb1d_vector_allocate_vector(this, new_field)
    !* Author: Chris MacMackin
    !  Date: October 2016
    !
    ! Allocates an abstract [[vector_field(type)]] to have a type
    ! compatible for operations on this field and to be the type
    ! returned by this field's methods which produce vectors.
    !
    class(cheb1d_vector_field), intent(in)          :: this
    class(vector_field), allocatable, intent(inout) :: new_field
      !! A field which, upon return, is allocated to be of the same
      !! concrete type as vector fields produced by `this`.
    if (allocated(new_field)) deallocate(new_field)
    allocate(cheb1d_vector_field :: new_field)
  end subroutine cheb1d_vector_allocate_vector

  pure subroutine cheb1d_vector_bound(this,src,boundary,depth,slices)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Alters the meta-data of the field to be appropriate for the
    ! specified boundary. Also computes the slices needed to access
    ! the field contents for the specified boundary.
    !
    class(cheb1d_vector_field), intent(out)           :: this
    class(array_vector_field), intent(in)             :: src
      !! The field for which the boundary data is to be provided.
    integer, intent(in)                               :: boundary
      !! Specifies which boundary is to be returned. The boundary
      !! will be the one normal to dimension of number
      !! `abs(boundary)`. If the argument is negative, then the
      !! lower boundary is returned. If positive, then the upper
      !! boundary is returned.
    integer, intent(in)                               :: depth
      !! The number of layers of data-points to return at the
      !! specified boundary.
    integer, dimension(:,:), allocatable, intent(out) :: slices
      !! An array which, on return, contains array slice data which
      !! can be used to construct extract the data for the given field
      !! boundary. The form of the array is
      !! ```
      !! slices(1,i) = start_index
      !! slices(2,i) = end_index
      !! slices(3,i) = stride
      !! ```
    integer :: length
    select type(src)
    class is(cheb1d_vector_field)
      length = size(src%colloc_points)
      select case(boundary)
      case(-1)
        this%colloc_points = src%colloc_points(length+1-depth:)
        this%extent = [src%colloc_points(length+1-depth), &
                      src%colloc_points(length)]
        this%differentiable = .false.
        allocate(slices(3,1))
        slices(:,1) = [length+1-depth,length,1]
      case(1)
        this%colloc_points = src%colloc_points(1:depth)
        this%extent = [src%colloc_points(depth), src%colloc_points(1)]
        this%differentiable = .false.
        allocate(slices(3,1))
        slices(:,1) = [1,depth,1]
      case default
        this%colloc_points = src%colloc_points
        this%extent = src%extent
        allocate(slices(3,1))
        slices(:,1) = [1,length,1]
      end select
    class default
      error stop ('Trying to compute boundary metadata for type other '// &
                 'than `cheb1d_vector_field`.')
    end select
  end subroutine cheb1d_vector_bound
  
  subroutine cheb1d_vector_field_write_hdf(this, hdf_id, dataset_name, &
                                           error)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Writes the contents of the field to a dataset in an HDF
    ! file. The dataset will have an attribute specifying the name of
    ! the field type and an attribute with a non-zero value,
    ! indicating that this is a vector field.
    !
    class(cheb1d_vector_field), intent(in) :: this
    integer(hid_t), intent(in)             :: hdf_id
      !! The identifier for the HDF file/group in which the field
      !! data is to be written.
    character(len=*), intent(in)           :: dataset_name
      !! The name of the dataset to be created in the HDF file
      !! containing this field's data.
    integer, intent(out)                   :: error
      !! An error code which, upon succesful completion of the
      !! routine, is 0. Otherwise, contains the error code returned
      !! by the HDF library.
    error = 0
    call this%write_hdf_array(hdf_id, dataset_name, [this%numpoints], error)
    if (error /= 0) return
    call h5ltset_attribute_string_f(hdf_id, dataset_name, hdf_field_type_attr, &
                                    hdf_vector_type, error)
    if (error /= 0) return
    call h5ltset_attribute_int_f(hdf_id, dataset_name, hdf_vector_attr, [1], &
                                 1_size_t, error)
    if (error /= 0) return
    call h5ltset_attribute_double_f(hdf_id, dataset_name, hdf_grid_attr//'1', &
                                    int(this%colloc_points,size_t), error)
  end subroutine cheb1d_vector_field_write_hdf

end module cheb1d_fields_mod
