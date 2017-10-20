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
  use iso_fortran_env, only: r8 => real64, stderr => error_unit
  use utils_mod, only: grid_to_spacing, open_or_create_hdf_group, &
                       create_hdf_dset, linear_interp
  use abstract_fields_mod
  use scalar_pool_mod, only: scalar_pool
  use vector_pool_mod, only: vector_pool
  use uniform_fields_mod, only: uniform_scalar_field, &
                                uniform_vector_field
  use array_fields_mod, only: array_scalar_field, array_vector_field, &
                              scalar_init, vector_init
  use chebyshev_mod
  use h5lt, only: hid_t, size_t, hsize_t, h5ltmake_dataset_double_f,     &
                  h5ltset_attribute_string_f, h5ltset_attribute_int_f,   &
                  h5ltread_dataset_double_f, h5ltget_attribute_string_f, &
                  h5ltget_dataset_info_f
  use hdf5, only: h5gclose_f
  implicit none
  private

  character(len=19), parameter, public :: hdf_scalar_name = 'cheb1d_scalar_field'
  character(len=19), parameter, public :: hdf_vector_name = 'cheb1d_vector_field'
  character(len=34), parameter, public :: hdf_grid_template = &
                                          '("cheb1d-",i0.5,"-",f0.7,"-",f0.7)'

$:public_unary()
  public :: minval
  public :: maxval

  logical :: initialised = .false.
  
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
    real(r8), dimension(2)          :: extent
      !! The start and end values of the domain
    real(r8), dimension(:), pointer :: colloc_points => null()
      !! The location of the data-points
    logical                         :: differentiable = .true.
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
    procedure, public :: set_boundary => cheb1d_scalar_set_bound
      !! Sets the specified boundary to hold the same values as the
      !! passed field argument.
    procedure, public :: read_hdf => cheb1d_scalar_read_hdf
      !! Read field data from a dataset in an HDF5 file.
    procedure, public :: write_hdf => cheb1d_scalar_write_hdf
      !! Write field data to a new dataset in an HDF5 file.
    procedure, public :: grid_spacing => cheb1d_scalar_grid_spacing
      !! Returns a vector field where the values are the size of cells
      !! in the grid of the field at each point, in each
      !! direction. This can be useful, e.g., when calculating time
      !! steps for integration.
    procedure :: force_finalise => cheb1d_scalar_force_finalise
      !! Frees the data array for this field, in order to reduce the
      !! volume of any memory leaks.
    final :: cheb1d_scalar_finalise
      !! Deallocates the contents of the field. The bulk of the
      !! field's memory should be deallocated with the `clean_temp`
      !! method, but that is unable to deallocate the pointer to the
      !! data array--only the array itself. In compilers which support
      !! finalisation, this method eliminates the small memory leak
      !! from the pointer.
    procedure :: array_interpolate => cheb1d_scalar_interp
      !! Interpolates a value in the field, using a 1-D array of data
      !! passed to it.
  end type cheb1d_scalar_field
  
  interface cheb1d_scalar_field
    module procedure scalar_constructor
  end interface

  type(scalar_pool) :: scalars

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
    real(r8), dimension(2)          :: extent
      !! The start and end values of the domain
    real(r8), dimension(:), pointer :: colloc_points => null()
      !! The location of the data-points
    logical                         :: differentiable = .true.
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
    procedure, public :: set_boundary => cheb1d_vector_set_bound
      !! Sets the specified boundary to hold the same values as the
      !! passed field argument.
    procedure, public :: read_hdf => cheb1d_vector_read_hdf
      !! Read field data from a dataset in an HDF5 file.
    procedure, public :: write_hdf => cheb1d_vector_write_hdf
      !! Write field data to a new dataset in an HDF5 file.
    procedure, public :: grid_spacing => cheb1d_vector_grid_spacing
      !! Returns a vector field where the values are the size of cells
      !! in the grid of the field at each point, in each
      !! direction. This can be useful, e.g., when calculating time
      !! steps for integration.
    procedure :: force_finalise => cheb1d_vector_force_finalise
      !! Frees the data array for this field, in order to reduce the
      !! volume of any memory leaks.
    final :: cheb1d_vector_finalise
      !! Deallocates the contents of the field. The bulk of the
      !! field's memory should be deallocated with the `clean_temp`
      !! method, but that is unable to deallocate the pointer to the
      !! data array--only the array itself. In compilers which support
      !! finalisation, this method eliminates the small memory leak
      !! from the pointer..
    procedure :: array_interpolate => cheb1d_vector_interp
      !! Interpolates a value in the field, using a 2-D array of data
      !! passed to it.
  end type cheb1d_vector_field
  
  interface cheb1d_vector_field
    module procedure vector_constructor
  end interface cheb1d_vector_field

  type(vector_pool) :: vectors
  
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
      !! An impure elemental procedure taking which takes the position in the
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
    field%colloc_points => collocation_points(numpoints-1,lower_bound,upper_bound)
    field = array_scalar_field(field,numpoints,initializer)
  end function scalar_constructor

  function cheb1d_scalar_domain(this) result(res)
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
    call this%guard_temp()
    allocate(res(1,2))
    res(1,:) = this%extent
    call this%clean_temp()
  end function cheb1d_scalar_domain

  impure elemental function cheb1d_scalar_dimensions(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! Accessor for the number of dimensions.
    !
    class(cheb1d_scalar_field), intent(in) :: this
    integer :: res
    call this%guard_temp()
    res = 1
    call this%clean_temp()
  end function cheb1d_scalar_dimensions

  function cheb1d_scalar_resolution(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! Get the number of data points used to represent the field in
    ! its one dimension.
    !
    class(cheb1d_scalar_field), intent(in) :: this
    integer, dimension(:), allocatable :: res
      !! Array of length 1 specifying the number of data points.
    call this%guard_temp()
    allocate(res(1))
    res(1) = this%elements()
    call this%clean_temp()
  end function cheb1d_scalar_resolution

  function cheb1d_scalar_raw_slices(this,exclude_lower_bound,exclude_upper_bound) &
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
    call this%guard_temp()
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
    call this%clean_temp()
  end function cheb1d_scalar_raw_slices

  function cheb1d_scalar_id_to_position(this, id) result(pos)
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
    call this%guard_temp()
    allocate(pos(1))
    if (id <= this%elements()) then
      pos(1) = this%colloc_points(id)
    else
      error stop ('cheb1d_scalar_field: Invalid ID number provided.')
    end if
    call this%clean_temp()
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
    call this%guard_temp()
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
    call this%clean_temp()
  end function cheb1d_scalar_array_dx

  subroutine cheb1d_scalar_check_compatible(this,other)
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
    call this%guard_temp(); call other%guard_temp()
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
    call this%clean_temp(); call other%clean_temp()
  end subroutine cheb1d_scalar_check_compatible

  subroutine cheb1d_scalar_assign_meta(this, rhs)
    !* Author: Chris MacMackin
    !  Date: April 2016
    !
    ! Assigns all data in `rhs` to `this` except for the actual values
    ! of the field.
    !
    class(cheb1d_scalar_field), intent(inout) :: this
    class(abstract_field), intent(in) :: rhs
    call rhs%guard_temp()
    if (associated(this%colloc_points)) then
      if (this%differentiable) then
        nullify(this%colloc_points)
      else
        deallocate(this%colloc_points)
      end if
    end if
    select type(rhs)
    class is(cheb1d_scalar_field)
      this%extent = rhs%extent
      if (associated(rhs%colloc_points)) then
        if (rhs%differentiable) then
          this%colloc_points => rhs%colloc_points
        else
          allocate(this%colloc_points(size(rhs%colloc_points)))
          this%colloc_points = rhs%colloc_points
        end if
      end if
      this%differentiable = rhs%differentiable
    class is(cheb1d_vector_field)
      this%extent = rhs%extent
      if (associated(rhs%colloc_points)) then
        if (rhs%differentiable) then
          this%colloc_points => rhs%colloc_points
        else
          allocate(this%colloc_points(size(rhs%colloc_points)))
          this%colloc_points = rhs%colloc_points
        end if
      end if
      this%differentiable = rhs%differentiable
    end select
    call rhs%clean_temp()
  end subroutine cheb1d_scalar_assign_meta

  subroutine cheb1d_scalar_allocate_scalar(this, new_field)
    !* Author: Chris MacMackin
    !  Date: October 2016
    !
    ! Allocates an abstract [[scalar_field(type)]] to have a type
    ! compatible for operations on this field and to be the type
    ! returned by this field's methods which produce scalars.
    !
    class(cheb1d_scalar_field), intent(in)      :: this
    class(scalar_field), pointer, intent(inout) :: new_field
      !! A field which, upon return, is allocated to be of the same
      !! concrete type as scalar fields produced by `this`.
    call this%guard_temp()
    if (.not. initialised) then
      scalars = scalar_pool(object_pool_size, this)
      block
        type(cheb1d_vector_field) :: vf
        vectors = vector_pool(object_pool_size, vf)
      end block
      initialised = .true.
    end if
    new_field => scalars%acquire()
    select type(new_field)
    class is(cheb1d_scalar_field)
      if (associated(new_field%colloc_points)) then
        if (new_field%differentiable) then
          nullify(new_field%colloc_points)
        else
          deallocate(new_field%colloc_points)
        end if
      end if
    end select
    call this%clean_temp()
  end subroutine cheb1d_scalar_allocate_scalar

  subroutine cheb1d_scalar_allocate_vector(this, new_field)
    !* Author: Chris MacMackin
    !  Date: October 2016
    !
    ! Allocates an abstract [[vector_field(type)]] to have a type
    ! compatible for operations on this field and to be the type
    ! returned by this field's methods which produce vectors.
    !
    class(cheb1d_scalar_field), intent(in)      :: this
    class(vector_field), pointer, intent(inout) :: new_field
      !! A field which, upon return, is allocated to be of the same
      !! concrete type as scalar fields produced by `this`.
    call this%guard_temp()
    if (.not. initialised) then
      scalars = scalar_pool(object_pool_size, this)
      block
        type(cheb1d_vector_field) :: vf
        vectors = vector_pool(object_pool_size, vf)
      end block
      initialised = .true.
    end if
    new_field => vectors%acquire()
    select type(new_field)
    class is(cheb1d_vector_field)
      if (associated(new_field%colloc_points)) then
        if (new_field%differentiable) then
          nullify(new_field%colloc_points)
        else
          deallocate(new_field%colloc_points)
        end if
      end if
    end select
    call this%clean_temp()
  end subroutine cheb1d_scalar_allocate_vector

  subroutine cheb1d_scalar_bound(this,src,boundary,depth,slices)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Alters the meta-data of the field to be appropriate for the
    ! specified boundary. Also computes the slices needed to access
    ! the field contents for the specified boundary.
    !
    class(cheb1d_scalar_field), intent(inout)         :: this
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
    call src%guard_temp()
    select type(src)
    class is(cheb1d_scalar_field)
      length = size(src%colloc_points)
      select case(boundary)
      case(-1)
        allocate(this%colloc_points(depth))
        this%colloc_points = src%colloc_points(length+1-depth:)
        this%extent = [src%colloc_points(length+1-depth), &
                      src%colloc_points(length)]
        this%differentiable = .false.
        allocate(slices(3,1))
        slices(:,1) = [length+1-depth,length,1]
      case(1)
        allocate(this%colloc_points(depth))
        this%colloc_points = src%colloc_points(1:depth)
        this%extent = [src%colloc_points(depth), src%colloc_points(1)]
        this%differentiable = .false.
        allocate(slices(3,1))
        slices(:,1) = [1,depth,1]
      case default
        if (src%differentiable) then
          this%colloc_points => src%colloc_points
        else
          if (associated(this%colloc_points)) then
            if (this%differentiable) then
              nullify(this%colloc_points)
            else
              deallocate(this%colloc_points)
            end if
          end if
          this%colloc_points = src%colloc_points
        end if
        this%extent = src%extent
        allocate(slices(3,1))
        slices(:,1) = [1,length,1]
      end select
    class default
      error stop ('Trying to compute boundary metadata for type other '// &
                 'than `cheb1d_scalar_field`.')
    end select
    call src%clean_temp()
  end subroutine cheb1d_scalar_bound

  subroutine cheb1d_scalar_set_bound(this,boundary,depth,boundary_field)
    !* Author: Chris MacMackin
    !  Date: March 2017
    !
    ! Sets the value of this field at the specified boundary to be the
    ! same as those in the passed field.
    !
    class(cheb1d_scalar_field), intent(inout) :: this
    integer, intent(in) :: boundary
      !! Specifies which boundary is to be set. The boundary
      !! will be the one normal to dimension of number
      !! `abs(boundary)`. If the argument is negative, then the
      !! lower boundary is set. If positive, then the upper
      !! boundary is set. If 0, then the whole field is
      !! set.
    integer, intent(in) :: depth
      !! The number of layers of data-points to set at the
      !! specified boundary.
    class(scalar_field), intent(in) :: boundary_field
      !! A field, of the same type as `this` and with the same
      !! resolution, number of dimensions etc., but containing only
      !! the points within the specified number of layers of cells
      !! adjecent to the specified boundary. Alternatively, it may
      !! be a [[uniform_scalar_field]].
    integer :: length, i
    real(r8) :: val
#:if defined('DEBUG')
    integer, dimension(:), allocatable :: res
#:endif
    call this%guard_temp(); call boundary_field%guard_temp()
    length = this%elements()
    select type(boundary_field)
    class is(cheb1d_scalar_field)
#:if defined('DEBUG')
      res = boundary_field%resolution()
      if (abs(boundary) >= lbound(res,1) .and. abs(boundary) <= ubound(res,1)) then
        if (res(abs(boundary)) /= depth) then
          error stop ('Boundary array not deep enough.')
        end if
      end if
#:endif
      select case(boundary)
      case(-1)
#:if defined('DEBUG')
        if (abs(this%colloc_points(length+1-depth) - boundary_field%extent(1)) &
            > 1e-12+r8) error stop ('Domain mismatch.')
        if (abs(this%colloc_points(length) - boundary_field%extent(2)) &
            > 1e-12+r8) error stop ('Domain mismatch.')
#:endif
        do i = 1, depth
          call this%set_element(length-depth+i, boundary_field%get_element(i))
        end do
      case(1)
#:if defined('DEBUG')
        if (abs(this%colloc_points(1) - boundary_field%extent(1)) &
            > 1e-12+r8) error stop ('Domain mismatch.')
        if (abs(this%colloc_points(depth) - boundary_field%extent(2)) &
            > 1e-12+r8) error stop ('Domain mismatch.')
#:endif
        do i = 1, depth
          call this%set_element(i, boundary_field%get_element(i))
        end do
      case default
#:if defined('DEBUG')
        if (abs(this%colloc_points(1) - boundary_field%extent(1)) &
            > 1e-12+r8) error stop ('Domain mismatch.')
        if (abs(this%colloc_points(length) - boundary_field%extent(2)) &
            > 1e-12+r8) error stop ('Domain mismatch.')
        if (length /= boundary_field%elements()) error stop ('Resolution mismatch.')
#:endif        
        this = boundary_field
      end select
    class is(uniform_scalar_field)
      val = boundary_field%get_value()
      select case(boundary)
      case(-1)
        do i = 1, depth
          call this%set_element(length-depth+i, val)
        end do
      case(1)
        do i = 1, depth
          call this%set_element(i, val)
        end do
      case default
        this = boundary_field
      end select      
    class default
      error stop ('Trying to set boundary from type other than '// &
                  '`cheb1d_scalar_field` or `uniform_scalar_field`.')
    end select
    call this%clean_temp(); call boundary_field%clean_temp()
  end subroutine cheb1d_scalar_set_bound

  subroutine cheb1d_scalar_read_hdf(this, hdf_id, dataset_name, error)
    !* Author: Chris MacMackin
    !  Date: April 2017
    !
    ! Reads the contents of the field from a dataset in an HDF
    ! file. The dataset should have an attribute specifying the name of
    ! the field type.
    !
    ! @Note It is assumed a differentiable field is being read (i.e.,
    ! one which isn't just a region near the boundary).
    !
    class(cheb1d_scalar_field), intent(inout) :: this
    integer(hid_t), intent(in)                :: hdf_id
      !! The identifier for the HDF file/group from which the field
      !! data is to be read.
    character(len=*), intent(in)              :: dataset_name
      !! The name of the dataset in the HDF file containing this
      !! field's data.
    integer, intent(out)                      :: error
      !! An error code which, upon succesful completion of the
      !! routine, is 0. Otherwise, contains the error code returned
      !! by the HDF library.
    integer :: type_class
    integer(size_t) :: type_size
    integer(hsize_t), dimension(1) :: dims
    character(len=50) :: string
    real(r8), dimension(:), allocatable :: grid
    call this%guard_temp()
    error = 0
    call h5ltget_dataset_info_f(hdf_id, dataset_name, dims, type_class, &
                                type_size, error)
    if (error < 0) then
      write(stderr,*) 'Error occurred when reading HDF dataset "'// &
                       dataset_name//'".'
      call this%clean_temp()
      return
    end if
    call h5ltget_attribute_string_f(hdf_id, dataset_name, hdf_field_type_attr, &
                                    string, error)
    if (error < 0 .or. trim(string) /= hdf_scalar_name) then
      write(stderr,*) 'HDF dataset "'//dataset_name//'" not '// &
                      'produced by cheb1d_scalar_field type.'
      error = -1
      call this%clean_temp()
      return
    end if
    call this%read_hdf_array(hdf_id, dataset_name, dims, error)
    call h5ltget_attribute_string_f(hdf_id, dataset_name, hdf_grid_attr//'1', &
                                    string, error)
    allocate(grid(this%elements()))
    call h5ltread_dataset_double_f(hdf_id, string, grid, dims, error)
    this%extent(1) = grid(this%elements())
    this%extent(2) = grid(1)
    deallocate(grid)
    this%colloc_points => collocation_points(this%elements()-1, this%extent(1), &
                                             this%extent(2))
    this%differentiable = .true.
    call this%clean_temp()
  end subroutine cheb1d_scalar_read_hdf

  subroutine cheb1d_scalar_write_hdf(this, hdf_id, dataset_name, error)
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
    integer(hid_t) :: group_id
    character(len=50) :: grid_name
    call this%guard_temp()
    error = 0
    call this%write_hdf_array(hdf_id, dataset_name, [int(this%elements(), &
                              hsize_t)], error)
    call h5ltset_attribute_string_f(hdf_id, dataset_name, hdf_field_type_attr, &
                                    hdf_scalar_name, error)
    call h5ltset_attribute_int_f(hdf_id, dataset_name, hdf_vector_attr, [0], &
                                 1_size_t, error)
    ! Create reference to grid
    write(grid_name,hdf_grid_template) size(this%colloc_points), &
                                       this%extent(1), this%extent(2)
    call open_or_create_hdf_group(hdf_id, hdf_grid_group, group_id, error)
    call create_hdf_dset(group_id, trim(grid_name), 1,   &
                         [int(this%elements(),hsize_t)], &
                         this%colloc_points, error)
    call h5ltset_attribute_string_f(hdf_id, dataset_name, hdf_grid_attr//'1', &
                                    hdf_grid_group//'/'//trim(grid_name), error)
    call h5gclose_f(group_id, error)
    call this%clean_temp()
  end subroutine cheb1d_scalar_write_hdf

  function cheb1d_scalar_grid_spacing(this) result(grid)
    !* Author: Chris MacMackin
    !  Date: December 2016
    !
    ! Returns a field containing the grid spacing at each point.
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(vector_field), pointer           :: grid
    class(vector_field), pointer           :: tmp
      !! A field where the values indicate the grid spacing that
      !! point. Each vector dimension representes the spacing of the
      !! grid in that direction.
    call this%guard_temp()
    call this%allocate_vector_field(grid)
    call grid%assign_meta_data(this)
    select type(grid)
    class is(array_vector_field)
      tmp => array_vector_field(grid, grid_to_spacing(this%colloc_points))
    end select
    grid => tmp
    call this%clean_temp()
  end function cheb1d_scalar_grid_spacing

  subroutine cheb1d_scalar_force_finalise(this)
    !* Author: Chris MacMackin
    !  Date: January 2017
    !
    ! Deallocates the field data for this object.
    !
    class(cheb1d_scalar_field), intent(in) :: this
    if (this%get_pool_id() /= non_pool_id) call scalars%release(this%get_pool_id())
  end subroutine cheb1d_scalar_force_finalise

  elemental subroutine cheb1d_scalar_finalise(this)
    !* Author: Chris MacMackin
    !  Date: January 2017
    !
    ! Deallocates the pointer to the field data for this object.
    ! 
    ! @Warning This only works for scalar field objects. For some
    ! reason gfortran segfualts when an elemental or higher rank
    ! finalisation routine is used.
    !
    type(cheb1d_scalar_field), intent(inout) :: this
    if (associated(this%colloc_points)) then
      call this%finalise()
      if (this%differentiable) then
        nullify(this%colloc_points)
      else
        deallocate(this%colloc_points)
      end if
    end if
  end subroutine cheb1d_scalar_finalise

  function cheb1d_scalar_interp(this, data_array, location) result(val)
    !* Author: Chris MacMackin
    !  Date: July 2017
    !
    ! Interpolates the value of the field at the specified location.
    !
    class(cheb1d_scalar_field), intent(in) :: this
    real(r8), dimension(:), intent(in)     :: data_array
      !! An array holding the datapoints for this field, identical
      !! in layout to that stored the field itself.
    real(r8), dimension(:), intent(in)     :: location
      !! The location at which to calculate the interpolated value.
    real(r8)                               :: val
    call this%guard_temp()
    val = linear_interp(location(1), data_array, this%colloc_points)
    call this%clean_temp()
  end function cheb1d_scalar_interp


  !=====================================================================
  ! Vector Field Method
  !=====================================================================


  function vector_constructor(numpoints,initializer,lower_bound, &
                              upper_bound, extra_dims) result(field)
    integer, intent(in) :: numpoints
      !! The number of collocation nodes to use when modelling this
      !! field. This corresponds to resolution.
    procedure(vector_init), optional :: initializer
      !! An impure elemental procedure which takes the position in the
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
    field%colloc_points => collocation_points(numpoints-1,lower_bound,upper_bound)
    field = array_vector_field(field,numpoints,dims,initializer)
  end function vector_constructor

  function cheb1d_vector_domain(this) result(res)
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
    call this%guard_temp()
    allocate(res(1,2))
    res(1,:) = this%extent
    call this%clean_temp()
  end function cheb1d_vector_domain

  impure elemental function cheb1d_vector_dimensions(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! Accessor for the number of dimensions.
    !
    class(cheb1d_vector_field), intent(in) :: this
    integer :: res
    call this%guard_temp()
    res = 1
    call this%clean_temp()
  end function cheb1d_vector_dimensions

  function cheb1d_vector_resolution(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! Get the number of data points used to represent the field in
    ! its one dimension.
    !
    class(cheb1d_vector_field), intent(in) :: this
    integer, dimension(:), allocatable :: res
      !! Array of length 1 specifying the number of data points.
    call this%guard_temp()
    allocate(res(1))
    res(1) = this%elements()
    call this%clean_temp()
  end function cheb1d_vector_resolution

  function cheb1d_vector_raw_slices(this,exclude_lower_bound,exclude_upper_bound) &
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
    call this%guard_temp()
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
    call this%clean_temp()
  end function cheb1d_vector_raw_slices
  
  function cheb1d_vector_id_to_position(this, id) result(pos)
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
    call this%guard_temp()
    allocate(pos(1))
    if (id <= this%elements()) then
      pos(1) = this%colloc_points(id)
    else
      error stop ('cheb1d_vector_field: Invalid ID number provided.')
    end if
    call this%clean_temp()
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
    call this%guard_temp()
    if (dir==1) then
#:if defined('DEBUG')
      if (.not. this%differentiable) &
        error stop ('Trying to take derivative of undifferentiable instance of '// &
                   '`cheb1d_vector_field`.')
#:endif
      res = data_array
      call differentiate_1d(res,this%colloc_points,order)
    else
      allocate(res(size(data_array)), mold=data_array)
      res = 0.0
    end if
    call this%clean_temp()
  end function cheb1d_vector_array_dx

  subroutine cheb1d_vector_check_compatible(this,other)
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
    call this%guard_temp(); call other%guard_temp()
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
    call this%clean_temp(); call other%clean_temp()
  end subroutine cheb1d_vector_check_compatible

  subroutine cheb1d_vector_assign_meta(this, rhs)
    !* Author: Chris MacMackin
    !  Date: April 2016
    !
    ! Assigns all data in `rhs` to `this` except for the actual values
    ! of the field.
    !
    class(cheb1d_vector_field), intent(inout) :: this
    class(abstract_field), intent(in) :: rhs
      !! If present and false, do not allocate the array of `this`.
    call rhs%guard_temp()
    if (associated(this%colloc_points)) then
      if (this%differentiable) then
        nullify(this%colloc_points)
      else
        deallocate(this%colloc_points)
      end if
    end if
    select type(rhs)
    class is(cheb1d_scalar_field)
      this%extent = rhs%extent
      if (associated(rhs%colloc_points)) then
        if (rhs%differentiable) then
          this%colloc_points => rhs%colloc_points
        else
          allocate(this%colloc_points(size(rhs%colloc_points)))
          this%colloc_points = rhs%colloc_points
        end if
      end if
      this%differentiable = rhs%differentiable
    class is(cheb1d_vector_field)
      this%extent = rhs%extent
      if (associated(rhs%colloc_points)) then
        if (rhs%differentiable) then
          this%colloc_points => rhs%colloc_points
        else
          allocate(this%colloc_points(size(rhs%colloc_points)))
          this%colloc_points = rhs%colloc_points
        end if
      end if
      this%differentiable = rhs%differentiable
    end select
    call rhs%clean_temp()
  end subroutine cheb1d_vector_assign_meta

  subroutine cheb1d_vector_allocate_scalar(this, new_field)
    !* Author: Chris MacMackin
    !  Date: October 2016
    !
    ! Allocates an abstract [[scalar_field(type)]] to have a type
    ! compatible for operations on this field and to be the type
    ! returned by this field's methods which produce scalars.
    !
    class(cheb1d_vector_field), intent(in)      :: this
    class(scalar_field), pointer, intent(inout) :: new_field
      !! A field which, upon return, is allocated to be of the same
      !! concrete type as scalar fields produced by `this`.
    call this%guard_temp()
    if (.not. initialised) then
      vectors = vector_pool(object_pool_size, this)
      block
        type(cheb1d_scalar_field) :: sf
        scalars = scalar_pool(object_pool_size, sf)
      end block
      initialised = .true.
    end if
    new_field => scalars%acquire()
    select type(new_field)
    class is(cheb1d_scalar_field)
      if (associated(new_field%colloc_points)) then
        if (new_field%differentiable) then
          nullify(new_field%colloc_points)
        else
          deallocate(new_field%colloc_points)
        end if
      end if
    end select
    call this%clean_temp()
  end subroutine cheb1d_vector_allocate_scalar

  subroutine cheb1d_vector_allocate_vector(this, new_field)
    !* Author: Chris MacMackin
    !  Date: October 2016
    !
    ! Allocates an abstract [[vector_field(type)]] to have a type
    ! compatible for operations on this field and to be the type
    ! returned by this field's methods which produce vectors.
    !
    class(cheb1d_vector_field), intent(in)      :: this
    class(vector_field), pointer, intent(inout) :: new_field
      !! A field which, upon return, is allocated to be of the same
      !! concrete type as vector fields produced by `this`.
    call this%guard_temp()
    if (.not. initialised) then
      vectors = vector_pool(object_pool_size, this)
      block
        type(cheb1d_scalar_field) :: sf
        scalars = scalar_pool(object_pool_size, sf)
      end block
      initialised = .true.
    end if
    new_field => vectors%acquire()
    select type(new_field)
    class is(cheb1d_vector_field)
      if (associated(new_field%colloc_points)) then
        if (new_field%differentiable) then
          nullify(new_field%colloc_points)
        else
          deallocate(new_field%colloc_points)
        end if
      end if
    end select
    call this%clean_temp()
  end subroutine cheb1d_vector_allocate_vector

  subroutine cheb1d_vector_bound(this,src,boundary,depth,slices)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Alters the meta-data of the field to be appropriate for the
    ! specified boundary. Also computes the slices needed to access
    ! the field contents for the specified boundary.
    !
    class(cheb1d_vector_field), intent(inout)         :: this
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
    call src%guard_temp()
    select type(src)
    class is(cheb1d_vector_field)
      length = size(src%colloc_points)
      select case(boundary)
      case(-1)
        allocate(this%colloc_points(depth))
        this%colloc_points = src%colloc_points(length+1-depth:)
        this%extent = [src%colloc_points(length+1-depth), &
                      src%colloc_points(length)]
        this%differentiable = .false.
        allocate(slices(3,1))
        slices(:,1) = [length+1-depth,length,1]
      case(1)
        allocate(this%colloc_points(depth))
        this%colloc_points = src%colloc_points(1:depth)
        this%extent = [src%colloc_points(depth), &
                      src%colloc_points(1)]
        this%differentiable = .false.
        allocate(slices(3,1))
        slices(:,1) = [1,depth,1]
      case default
        if (this%differentiable) then
          this%colloc_points => src%colloc_points
        else
          if (associated(this%colloc_points)) then
            if (this%differentiable) then
              nullify(this%colloc_points)
            else
              deallocate(this%colloc_points)
            end if
          end if
          this%colloc_points = src%colloc_points
        end if
        this%extent = src%extent
        allocate(slices(3,1))
        slices(:,1) = [1,length,1]
      end select
    class default
      error stop ('Trying to compute boundary metadata for type other '// &
                 'than `cheb1d_vector_field`.')
    end select
    call src%clean_temp()
  end subroutine cheb1d_vector_bound

  subroutine cheb1d_vector_set_bound(this,boundary,depth,boundary_field)
    !* Author: Chris MacMackin
    !  Date: March 2017
    !
    ! Sets the value of this field at the specified boundary to be the
    ! same as those in the passed field.
    !
    class(cheb1d_vector_field), intent(inout) :: this
    integer, intent(in) :: boundary
      !! Specifies which boundary is to be set. The boundary
      !! will be the one normal to dimension of number
      !! `abs(boundary)`. If the argument is negative, then the
      !! lower boundary is set. If positive, then the upper
      !! boundary is set. If 0, then the whole field is
      !! set.
    integer, intent(in) :: depth
      !! The number of layers of data-points to set at the
      !! specified boundary.
    class(vector_field), intent(in) :: boundary_field
      !! A field, of the same type as `this` and with the same
      !! resolution, number of dimensions etc., but containing only
      !! the points within the specified number of layers of cells
      !! adjecent to the specified boundary. Alternatively, it may
      !! be a [[uniform_vector_field]].
    integer :: length, i, n
    real(r8), dimension(:), allocatable :: val
#:if defined('DEBUG')
    integer, dimension(:), allocatable :: res
#:endif
    call this%guard_temp(); call boundary_field%guard_temp()
    length = this%elements()
    n = this%vector_dimensions()
    select type(boundary_field)
    class is(cheb1d_vector_field)
#:if defined('DEBUG')
      res = boundary_field%resolution()
      if (abs(boundary) >= lbound(res,1) .and. abs(boundary) <= ubound(res,1)) then
        if (res(abs(boundary)) /= depth) then
          error stop ('Boundary array not deep enough.')
        end if
      end if
#:endif
      select case(boundary)
      case(-1)
#:if defined('DEBUG')
        if (abs(this%colloc_points(length+1-depth) - boundary_field%extent(1)) &
            > 1e-12+r8) error stop ('Domain mismatch.')
        if (abs(this%colloc_points(length) - boundary_field%extent(2)) &
            > 1e-12+r8) error stop ('Domain mismatch.')
#:endif
        do i = 1, depth
          val = boundary_field%get_element(i)
          call this%set_element(length-depth+i, val(1:n))
        end do
      case(1)
#:if defined('DEBUG')
        if (abs(this%colloc_points(1) - boundary_field%extent(1)) &
            > 1e-12+r8) error stop ('Domain mismatch.')
        if (abs(this%colloc_points(depth) - boundary_field%extent(2)) &
            > 1e-12+r8) error stop ('Domain mismatch.')
#:endif
        do i = 1, depth
          val = boundary_field%get_element(i)
          call this%set_element(i, val(1:n))
        end do
      case default
#:if defined('DEBUG')
        if (abs(this%colloc_points(1) - boundary_field%extent(1)) &
            > 1e-12+r8) error stop ('Domain mismatch.')
        if (abs(this%colloc_points(length) - boundary_field%extent(2)) &
            > 1e-12+r8) error stop ('Domain mismatch.')
        if (length /= boundary_field%elements()) error stop ('Resolution mismatch.')
#:endif
        this = boundary_field
      end select
    class is(uniform_vector_field)
      val = boundary_field%get_value()
      select case(boundary)
      case(-1)
        do i = 1, depth
          call this%set_element(length-depth+i, val(1:n))
        end do
      case(1)
        do i = 1, depth
          call this%set_element(i, val(1:n))
        end do
      case default
        this = boundary_field
      end select      
    class default
      error stop ('Trying to set boundary from type other than '// &
                  '`cheb1d_vector_field` or `uniform_vector_field`.')
    end select
    call this%clean_temp(); call boundary_field%clean_temp()
  end subroutine cheb1d_vector_set_bound

  subroutine cheb1d_vector_read_hdf(this, hdf_id, dataset_name, error)
    !* Author: Chris MacMackin
    !  Date: April 2017
    !
    ! Reads the contents of the field from a dataset in an HDF
    ! file. The dataset should have an attribute specifying the name of
    ! the field type.
    !
    ! @Note It is assumed a differentiable field is being read (i.e.,
    ! one which isn't just a region near the boundary).
    !
    class(cheb1d_vector_field), intent(inout) :: this
    integer(hid_t), intent(in)                :: hdf_id
      !! The identifier for the HDF file/group from which the field
      !! data is to be read.
    character(len=*), intent(in)              :: dataset_name
      !! The name of the dataset in the HDF file containing this
      !! field's data.
    integer, intent(out)                      :: error
      !! An error code which, upon succesful completion of the
      !! routine, is 0. Otherwise, contains the error code returned
      !! by the HDF library.
    integer :: type_class
    integer(size_t) :: type_size
    integer(hsize_t), dimension(2) :: dims
    character(len=50) :: string
    real(r8), dimension(:), allocatable :: grid
    call this%guard_temp()
    error = 0
    call h5ltget_dataset_info_f(hdf_id, dataset_name, dims, type_class, &
                                type_size, error)
    if (error < 0) then
      write(stderr,*) 'Error occurred when reading HDF dataset "'// &
                       dataset_name//'".'
      call this%clean_temp()
      return
    end if
    call h5ltget_attribute_string_f(hdf_id, dataset_name, hdf_field_type_attr, &
                                    string, error)
    if (error < 0 .or. trim(string) /= hdf_vector_name) then
      write(stderr,*) 'HDF dataset "'//dataset_name//'" not '// &
                      'produced by cheb1d_vector_field type.'
      error = -1
      call this%clean_temp()
      return
    end if
    call this%read_hdf_array(hdf_id, dataset_name, dims, error)
    call h5ltget_attribute_string_f(hdf_id, dataset_name, hdf_grid_attr//'1', &
                                    string, error)
    allocate(grid(this%elements()))
    call h5ltread_dataset_double_f(hdf_id, string, grid, dims, error)
    this%extent(1) = grid(this%elements())
    this%extent(2) = grid(1)
    deallocate(grid)
    this%colloc_points => collocation_points(this%elements()-1, this%extent(1), &
                                             this%extent(2))
    this%differentiable = .true.
    call this%clean_temp()
  end subroutine cheb1d_vector_read_hdf

  subroutine cheb1d_vector_write_hdf(this, hdf_id, dataset_name, error)
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
    integer(hid_t) :: group_id
    character(len=50) :: grid_name
    call this%guard_temp()
    error = 0
    call this%write_hdf_array(hdf_id, dataset_name, [int(this%elements(), &
                              hsize_t)], error)
    if (error /= 0) return
    call h5ltset_attribute_string_f(hdf_id, dataset_name, hdf_field_type_attr, &
                                    hdf_vector_name, error)
    if (error /= 0) return
    call h5ltset_attribute_int_f(hdf_id, dataset_name, hdf_vector_attr, [1], &
                                 1_size_t, error)
    if (error /= 0) return
    ! Create reference to grid
    write(grid_name,hdf_grid_template) size(this%colloc_points), &
                                       this%extent(1), this%extent(2)
    call open_or_create_hdf_group(hdf_id, hdf_grid_group, group_id, error)
    call create_hdf_dset(group_id, trim(grid_name), 1,   &
                         [int(this%elements(),hsize_t)], &
                         this%colloc_points, error)
    call h5ltset_attribute_string_f(hdf_id, dataset_name, hdf_grid_attr//'1', &
                                    hdf_grid_group//'/'//trim(grid_name), error)
    call h5gclose_f(group_id, error)
    call this%clean_temp()
  end subroutine cheb1d_vector_write_hdf

  function cheb1d_vector_grid_spacing(this) result(grid)
    !* Author: Chris MacMackin
    !  Date: December 2016
    !
    ! Returns a field containing the grid spacing at each point.
    !
    class(cheb1d_vector_field), intent(in) :: this
    class(vector_field), pointer           :: grid
      !! A field where the values indicate the grid spacing that
      !! point. Each vector dimension representes the spacing of the
      !! grid in that direction.
    class(vector_field), pointer :: tmp
    call this%guard_temp()
    call this%allocate_vector_field(grid)
    call grid%assign_meta_data(this)
    select type(grid)
    class is(array_vector_field)
      tmp => array_vector_field(grid, grid_to_spacing(this%colloc_points))
    end select
    grid => tmp
    call this%clean_temp()
  end function cheb1d_vector_grid_spacing

  subroutine cheb1d_vector_force_finalise(this)
    !* Author: Chris MacMackin
    !  Date: January 2017
    !
    ! Deallocates the field data for this object.
    !
    class(cheb1d_vector_field), intent(in) :: this
    if (this%get_pool_id() /= non_pool_id) call vectors%release(this%get_pool_id())
  end subroutine cheb1d_vector_force_finalise

  elemental subroutine cheb1d_vector_finalise(this)
    !* Author: Chris MacMackin
    !  Date: January 2017
    !
    ! Deallocates the field data for this object.
    !
    ! @Warning This only works for scalar field objects. For some
    ! reason gfortran segfualts when an elemental or higher rank
    ! finalisation routine is used.
    !
    type(cheb1d_vector_field), intent(inout) :: this
    if (associated(this%colloc_points)) then
      call this%finalise()
      if (this%differentiable) then
        nullify(this%colloc_points)
      else
        deallocate(this%colloc_points)
      end if
    end if
  end subroutine cheb1d_vector_finalise

  function cheb1d_vector_interp(this, data_array, location) result(val)
    !* Author: Chris MacMackin
    !  Date: July 2017
    !
    ! Interpolates the value of the field at the specified location.
    !
    class(cheb1d_vector_field), intent(in) :: this
    real(r8), dimension(:,:), intent(in)   :: data_array
      !! An array holding the datapoints for this field, identical
      !! in layout to that stored the field itself.
    real(r8), dimension(:), intent(in)     :: location
      !! The location at which to calculate the interpolated value.
    real(r8), dimension(:), allocatable    :: val
    integer :: i
    call this%guard_temp()
    do i = 1, size(data_array, 2)
      val(i) = linear_interp(location(1), data_array(:,i), this%colloc_points)
    end do
    call this%clean_temp()
  end function cheb1d_vector_interp

end module cheb1d_fields_mod
