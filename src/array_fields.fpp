!
!  array_fields.f90
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

module array_fields_mod
  !* Author: Chris MacMackin
  !  Date: September 2016
  !  License: GPLv3
  !
  ! Provides abstract extensions of the fileds in [[abstract_fields_mod]]
  ! which store field data in an array. These fields implement all of the
  ! operations except for the differentiation. Currently only fields using
  ! Cartesian coordinates are supported.
  !
  use iso_fortran_env, only: r8 => real64, stderr => error_unit
  use abstract_fields_mod, only: abstract_field, scalar_field, vector_field, &
                                 get_tol
  use utils_mod, only: is_nan, check_set_from_raw, elements_in_slice
  implicit none
  private

  type, extends(scalar_field), abstract, public :: array_scalar_field
    !* Author: Chris MacMackin
    !  Date: September 2016
    !
    ! An abstract data type representing a mathematical field of
    ! scalar values, with the field contents stored in a 1-D
    ! array. All of the operators are implemented except for
    ! differentiation. Currently only fields using Cartesian
    ! coordinates are supported by the vector and calculus routines.
    !
    ! Note that, when performing an operation on two fields, an error
    ! will occur if one of the fields is of the incorrect type or if
    ! the sizes of the fields do not match.
    !
    ! @TODO Add option to interpolate values of one of the fields when
    ! using a binary operator, if they are not otherwise compatible.
    ! Similarly, have lower-dimensional fields act as though they are 
    ! uniform on all higher dimensions when being operated on with a 
    ! higher dimensional field.
    !
    private
    integer :: numpoints
      !! The number of datapoints used
    real(r8), dimension(:), allocatable :: field_data
      !! The value of the scalar field at the data-points. Each
      !! element represents the value at a different location.
  contains
    private
    procedure, non_overridable, public :: elements => array_scalar_elements
      !! Specifies the number of individual data points present in this field. 
    procedure, public :: raw_size => array_scalar_raw_size
      !! Provides the number of pieces of data needed to represent the
      !! field, i.e. the size of the array returned by `raw`.
    procedure, public :: raw => array_scalar_raw
      !! Returns array of data representing state of field. Can be
      !! useful for passing to nonlinear solvers
    procedure, public :: set_from_raw => array_scalar_set_from_raw
      !! Assigns raw data, such as that produced by 
      !! [[array_scalar_field:raw]], to the field
    procedure(sf_raw_slices), deferred :: raw_slices
      !! Returns an array of integers used for getting the correct
      !! data for a given raw representation of the field.
    procedure :: field_multiply_field => array_scalar_sf_m_sf
      !! \({\rm field} \times {\rm field}\)
    procedure :: field_multiply_vecfield => array_scalar_sf_m_vf
      !! \({\rm field} \times {\rm \vec{field}}\)
    procedure, pass(rhs) :: real_multiply_field => array_scalar_r_m_sf
      !! \({\rm real}  \times {\rm field}\)
    procedure :: field_multiply_real => array_scalar_sf_m_r
      !! \({\rm field} \times {\rm real}\)
    procedure :: field_divide_field => array_scalar_sf_d_sf
      !! \(\frac{\rm field}{\rm field}\)
    procedure, pass(rhs) :: real_divide_field => array_scalar_r_d_sf
      !! \(\frac{\rm real}{\rm field}\)
    procedure :: field_divide_real => array_scalar_sf_d_r
      !! \(\frac{\rm field}{\rm real}\)
    procedure :: field_add_field => array_scalar_sf_a_sf
      !! \({\rm field} + {\rm field}\)
    procedure, pass(rhs) :: real_add_field => array_scalar_r_a_sf
      !! \({\rm real} + {\rm field}\)
    procedure :: field_add_real => array_scalar_sf_a_r
      !! \({\rm field} + {\rm real}\)
    procedure :: field_sub_field => array_scalar_sf_s_sf
      !! \({\rm field} - {\rm field}\)
    procedure, pass(rhs) :: real_sub_field => array_scalar_r_s_sf
      !! \({\rm real} - {\rm field}\)
    procedure :: field_sub_real => array_scalar_sf_s_r
      !! \({\rm field} - {\rm real}\)
    procedure :: field_pow_real => array_scalar_sf_p_r
      !! \({\rm field}^{\rm real}\)
    procedure :: field_pow_real4 => array_scalar_sf_p_r4
      !! \({\rm field}^{\rm real}\)
    procedure :: field_pow_int => array_scalar_sf_p_i
      !! \({\rm field}^{\rm int}\)
#:for FUNC, TEX in UNARY_FUNCTIONS
    $:unary_binding(FUNC, TEX, 'array_scalar')
#:endfor
    procedure :: minval => array_scalar_minval
      !! \(\min({\rm field})\)
    procedure :: maxval => array_scalar_maxval
      !! \(\max({\rm field})\)
    procedure, public :: d_dx => array_scalar_d_dx
      !! \(\frac{\partial^n}{\partial x_i^n}({\rm field})\)
    procedure :: gradient => array_scalar_gradient
      !! \(\nabla {\rm field}\)
    procedure :: laplacian => array_scalar_laplacian
      !! \(\nabla^2 {\rm field}\)
    procedure :: is_equal => array_scalar_is_equal
      !! Checks fields are equal within a tolerance
    procedure :: assign_field => array_scalar_assign
      !! \({\rm field} = {\rm field}\)
    procedure, public :: assign_meta_data => array_scalar_assign_meta_data
      !! Copies all data other than values stored in field from another
      !! field object to this one.
    procedure(sf_meta), deferred :: assign_subtype_meta_data
      !! Copies all data stored in a subtype of [[array_scalar_field(type)]]
      !! from another field object to this one.
    procedure, non_overridable :: check_compatible => array_scalar_compatible
      !! Tests whether two fields are suitable for binary operations together
    procedure(sf_compatible), deferred :: check_subtype_compatible
      !! Tests whether two fields are suitable for binary operations
      !! together, checking that any properties of subtypes of
      !! [[array_scalar_field(type)]] are compatible.
    procedure(sf_scalar_dx), deferred :: array_dx
      !! Takes the derivative of the scalar field using a 1-D array of
      !! data passed to it.
  end type array_scalar_field

  interface array_scalar_field
    module procedure array_scalar_constructor
  end interface array_scalar_field

  abstract interface
    pure function sf_raw_slices(this,return_lower_bound,return_upper_bound) &
                                                             result(slices)
      import array_scalar_field
      class(array_scalar_field), intent(in)       :: this
      logical, dimension(:), optional, intent(in) :: return_lower_bound
        !! Specifies whether to return the values at the lower boundary
        !! for each dimension, with the index of the element
        !! corresponding to the dimension. Defaults to all `.true.`.
      logical, dimension(:), optional, intent(in) :: return_upper_bound
        !! Specifies whether to return the values at the upper boundary
        !! for each dimension, with the index of the element
        !! corresponding to the dimension. Defaults to all `.true.`.
      integer, dimension(:,:), allocatable        :: slices
        !! An array containing array slice data which can be used to
        !! construct the raw representation of a field, with the given
        !! boundary conditions. The form of the array is
        !! ```
        !! slices(1,i) = start_index
        !! slices(2,i) = end_index
        !! slices(3,i) = stride
        !! ```
    end function sf_raw_slices
     
    pure subroutine sf_compatible(this,other)
      !* Author: Chris MacMackin
      !  Date: April 2016
      !
      ! Checks whether a field has the same type, boundaries, and
      ! resolution as this one, making it compatible for binary 
      ! operations. If incompatible, stops program with an error message.
      !
      import :: array_scalar_field
      import :: abstract_field
      class(array_scalar_field), intent(in) :: this
      class(abstract_field), intent(in) :: other
        !! The field being checked against this one
    end subroutine sf_compatible

    pure subroutine sf_meta(this, rhs)
      import :: array_scalar_field
      import :: abstract_field
      class(array_scalar_field), intent(inout) :: this
      class(abstract_field), intent(in) :: rhs
        !! The field whose metadata (data which are not field vales)
        !! is to be copied
    end subroutine sf_meta

    function sf_scalar_dx(this, data_array, dir, order) result(res)
      import :: array_scalar_field
      import :: r8
      class(array_scalar_field), intent(in) :: this
      real(r8), dimension(:), intent(in)    :: data_array
        !! An array holding the datapoints for this field, identical
        !! in layout to that stored the field itself.
      integer, intent(in)                   :: dir
        !! Direction in which to differentiate
      integer, intent(in), optional         :: order
        !! Order of the derivative, default = 1
      real(r8), dimension(:), allocatable   :: res
        !! The spatial derivative of order `order` taken in direction `dir`
    end function sf_scalar_dx

    pure function scalar_init(x) result(scalar)
      !! Function used to specify value held by a scalar field at
      !! location `x`.
      import :: r8
      real(r8), dimension(:), intent(in) :: x 
      !! The position at which this function is evaluated
      real(r8) :: scalar
      !! The value of the field at this location
    end function scalar_init
  end interface


  type, extends(vector_field), abstract, public :: array_vector_field
    !* Author: Chris MacMackin
    !  Date: September 2016
    !
    ! An abstract data type representing a mathematical field of vector
    ! values, with the field contents stored in a 1-D array. All of the
    ! non-calculus operators are implemented.
    !
    ! Note that, when performing an operation on two fields, an error
    ! will occur if one of the fields is of the incorrect type or if
    ! the sizes of the fields do not match.
    !
    ! @TODO Add option to interpolate values of one of the fields when
    ! using a binary operator, if they are not otherwise compatible.
    ! Similarly, have lower-dimensional fields act as though they are 
    ! uniform on all higher dimensions when being operated on with a 
    ! higher dimensional field.
    !
    private
    integer :: numpoints
      !! The number of datapoints used
    real(r8), dimension(:,:), allocatable :: field_data
      !! The value of the vector field at the data-points. Each row
      !! represents a different spatial location, while each column
      !! represents a different component of the vector.
    integer :: vector_dims = 0
      !! The number of vector components
  contains
    private
    procedure, non_overridable, public :: &
         vector_dimensions => array_vector_vector_dimensions
      !! Returns dimension of the vectors in the field
    procedure, public :: elements => array_vector_elements
      !! Specifies the number of individual data points present in this field. 
    procedure, public :: raw_size => array_vector_raw_size
      !! Provides the number of pieces of data needed to represent the
      !! field, i.e. the size of the array returned by `raw`.
    procedure, public :: raw => array_vector_raw
      !! Returns array of data representing state of field. Can be
      !! useful for passing to nonlinear solvers.
    procedure, public :: set_from_raw => array_vector_set_from_raw
      !! Assigns raw data, such as that produced by 
      !! [[array_vector_field:raw]], to the field
    procedure(vf_raw_slices), deferred :: raw_slices
      !! Returns an array of integers used for getting the correct
      !! data for a given raw representation of the field.
    procedure :: field_multiply_field => array_vector_vf_m_sf
      !! \({\rm field} \times {\rm field}\)
    procedure, pass(rhs) :: real_multiply_field => array_vector_r_m_vf
      !! \({\rm real}  \times {\rm field}\)
    procedure :: field_multiply_real => array_vector_vf_m_r
      !! \({\rm field} \times {\rm real}\)
    procedure :: field_divide_field => array_vector_vf_d_sf
      !! \(\frac{\rm field}{\rm field}\)
    procedure :: field_divide_real => array_vector_vf_d_r
      !! \(\frac{\rm field}{\rm real}\)
    procedure :: field_add_field => array_vector_vf_a_vf
      !! \({\rm field} + {\rm field}\)
    procedure, pass(rhs) :: real_add_field => array_vector_r_a_vf
      !! \({\rm real} + {\rm field}\)
    procedure :: field_add_real => array_vector_vf_a_r
      !! \({\rm field} + {\rm real}\)
    procedure :: field_sub_field => array_vector_vf_s_vf
      !! \({\rm field} - {\rm field}\)
    procedure, pass(rhs) :: real_sub_field => array_vector_r_s_vf
      !! \({\rm real} - {\rm field}\)
    procedure :: field_sub_real => array_vector_vf_s_r
      !! \({\rm field} - {\rm real}\)
    procedure, public :: norm => array_vector_norm
      !! \(\lVert {\rm \vec{field}} \rVert\)
    procedure, public :: component => array_vector_component
      !! Returns a scalar field containing the specified component of 
      !! the vector field
    procedure, public :: d_dx => array_vector_d_dx
      !! \(\frac{\partial^n}{\partial x_i^n}({\rm \vec{field}})\)
    procedure, public :: component_d_dx => array_vector_component_d_dx
      !! \(\frac{\partial^n}{\partial x_i^n}({\rm field_j})\)
    procedure :: divergence => array_vector_divergence
      !! \(\nabla\cdot {\rm field}\)
    procedure :: curl => array_vector_curl
      !! \(\nabla\times {\rm field}\)
    procedure :: laplacian => array_vector_laplacian
      !! \(\nabla^2 {\rm field}\)
    procedure :: dot_prod => array_vector_dot_prod
      !! \({\rm \vec{field}} \cdot {\rm \vec{field}}\)
    procedure :: cross_prod => array_vector_cross_prod
      !! \({\rm\vec{field}} \times {\rm\vec{field}}\)
    procedure :: assign_field => array_vector_assign
      !! \({\rm field} = {\rm field}\)
    procedure :: is_equal => array_vector_is_equal
      !! Checks fields are equal within a tolerance
    procedure, public :: assign_meta_data => array_vector_assign_meta_data
      !! Copies all data other than values stored in field from another
      !! field object to this one.
    procedure(vf_meta), deferred :: assign_subtype_meta_data
      !! Copies all data stored in a subtype of [[array_vector_field(type)]]
      !! from another field object to this one.
    procedure, non_overridable :: check_compatible => array_vector_compatible
      !! Tests whether two fields are suitable for binary operations together
    procedure(vf_compatible), deferred :: check_subtype_compatible
      !! Tests whether two fields are suitable for binary operations
      !! together, checking that any properties of subtypes of
      !! [[array_vector_field(type)]] are compatible.
    procedure(vf_scalar_dx), deferred :: array_dx
      !! Takes the derivative of particular vector component of the
      !! field, using a 1-D array of data passed to it.
  end type array_vector_field

  interface array_vector_field
    module procedure array_vector_constructor
  end interface array_vector_field

  abstract interface
    pure function vf_raw_slices(this,return_lower_bound,return_upper_bound) &
                                                             result(slices)
      import array_vector_field
      class(array_vector_field), intent(in)       :: this
      logical, dimension(:), optional, intent(in) :: return_lower_bound
        !! Specifies whether to return the values at the lower boundary
        !! for each dimension, with the index of the element
        !! corresponding to the dimension. Defaults to all `.true.`.
      logical, dimension(:), optional, intent(in) :: return_upper_bound
        !! Specifies whether to return the values at the upper boundary
        !! for each dimension, with the index of the element
        !! corresponding to the dimension. Defaults to all `.true.`.
      integer, dimension(:,:), allocatable        :: slices
        !! An array containing array slice data which can be used to
        !! construct the raw representation of a field, with the given
        !! boundary conditions. The form of the array is
        !! ```
        !! slices(1,i) = start_index
        !! slices(2,i) = end_index
        !! slices(3,i) = stride
        !! ```
    end function vf_raw_slices

    pure subroutine vf_compatible(this,other)
      !* Author: Chris MacMackin
      !  Date: September 2016
      !
      ! Checks whether a field has the same type, boundaries, and
      ! resolution as this one, making it compatible for binary 
      ! operations. If incompatible, stops program with an error message.
      !
      import :: array_vector_field
      import :: abstract_field
      class(array_vector_field), intent(in) :: this
      class(abstract_field), intent(in) :: other
        !! The field being checked against this one
    end subroutine vf_compatible  

    pure subroutine vf_meta(this, rhs)
      import :: array_vector_field
      import :: abstract_field
      class(array_vector_field), intent(inout) :: this
      class(abstract_field), intent(in) :: rhs
        !! The field whose metadata (data which are not field vales)
        !! is to be copied
    end subroutine vf_meta

    function vf_scalar_dx(this, data_array, dir, order) result(res)
      import :: array_vector_field
      import :: r8
      class(array_vector_field), intent(in) :: this
      real(r8), dimension(:), intent(in)    :: data_array
        !! An array holding the datapoints for a component of the
        !! vectors in this field, with identical in layout to the
        !! storage in the field itself.
      integer, intent(in)                   :: dir
        !! Direction in which to differentiate
      integer, intent(in), optional         :: order
        !! Order of the derivative, default = 1
      real(r8), dimension(:), allocatable   :: res
        !! The spatial derivative of order `order` taken in direction `dir`
    end function vf_scalar_dx

    function vf_vector_dx(this, data_array, dir, order) result(res)
      import :: array_vector_field
      import :: r8
      class(array_vector_field), intent(in) :: this
      real(r8), dimension(:,:), intent(in)  :: data_array
        !! An array holding the datapoints for the vectors in this
        !! field, with identical in layout to the storage in the field
        !! itself. Each row represents a different spatial location,
        !! while each column represents a different component of the
        !! vector.
      integer, intent(in)                   :: dir
        !! Direction in which to differentiate
      integer, intent(in), optional         :: order
        !! Order of the derivative, default = 1
      real(r8), dimension(:,:), allocatable :: res
        !! The spatial derivative of order `order` taken in direction `dir`
    end function vf_vector_dx
    
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

  public :: scalar_init, vector_init

contains


  !=====================================================================
  ! Scalar Field Methods
  !=====================================================================

  function array_scalar_constructor(template,numpoints, &
       initializer) result(this)
    !* Author: Chris MacMackin
    !  Date: October 2016
    !
    ! Creates a new field with the same concrete type as the template
    ! argument. The array of values will be allocated an initiated.
    !
    class(array_scalar_field), intent(in)  :: template
      !! A scalar field object which will act as a mold for the concrete
      !! type of the returned type, also copying over any metadata.
    integer, intent(in)                    :: numpoints
      !! The number of data points needed in the array when modelling this
      !! field.
    procedure(scalar_init), optional       :: initializer
      !! An elemental procedure taking which takes the position in the
      !! fields domain (an 8-byte real) as an argument and returns the
      !! fields value at that position. Default is for field to be zero
      !! everywhere.
    class(array_scalar_field), allocatable :: this
      !! A scalar field initiated based on the arguments to this function.
    integer :: i
    allocate(this, source=template)
    if (allocated(this%field_data)) deallocate(this%field_data)
    allocate(this%field_data(numpoints))
    this%numpoints = numpoints
    if (present(initializer)) then
      do i = 1, numpoints
        this%field_data(i) = initializer(this%id_to_position(i))
      end do
    end if
  end function array_scalar_constructor

  pure function array_scalar_elements(this) result(elements)
    !* Author: Chris MacMackin
    !  Date: October 2016
    !
    ! Gives the number of individual data points present in the field.
    !
    class(array_scalar_field), intent(in) :: this
    integer :: elements
    elements = this%numpoints
  end function array_scalar_elements

  pure function array_scalar_raw_size(this,return_lower_bound, &
                                       return_upper_bound) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! Compute how many elements are in the raw representation of this
    ! field. This would be the number of data points, adjusted based on
    ! how boundary conditions are accounted for.
    !
    class(array_scalar_field), intent(in)       :: this
    logical, dimension(:), optional, intent(in) :: return_lower_bound
      !! Specifies whether to return the values at the lower boundary
      !! for each dimension, with the index of the element
      !! corresponding to the dimension. Defaults to all `.true.`.
    logical, dimension(:), optional, intent(in) :: return_upper_bound
      !! Specifies whether to return the values at the upper boundary
      !! for each dimension, with the index of the element
      !! corresponding to the dimension. Defaults to all `.true.`.
    integer :: res
    integer, dimension(:,:), allocatable :: slices
    integer :: i
    if (allocated(this%field_data)) then
      slices = this%raw_slices(return_lower_bound,return_upper_bound)
      res = sum(elements_in_slice(slices(1,:),slices(2,:),slices(3,:)))
    else
      res = 0
    end if
  end function array_scalar_raw_size
  
  pure function array_scalar_raw(this,return_lower_bound, &
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
    class(array_scalar_field), intent(in) :: this
    logical, dimension(:), optional, intent(in) :: return_lower_bound
      !! Specifies whether to return the values at the lower boundary
      !! for each dimension, with the index of the element
      !! corresponding to the dimension. Defaults to all `.true.`.
    logical, dimension(:), optional, intent(in) :: return_upper_bound
      !! Specifies whether to return the values at the upper boundary
      !! for each dimension, with the index of the element
      !! corresponding to the dimension. Defaults to all `.true.`.
    real(r8), dimension(:), allocatable :: res
      !! Array containing data needed to describe field
    integer, dimension(:,:), allocatable :: slices
    integer :: i
    slices = this%raw_slices(return_lower_bound,return_upper_bound)
    if (allocated(this%field_data)) then
      res = [(this%field_data(slices(1,i):slices(2,i):slices(3,i)), i=1, &
             size(slices,2))]
    else
      allocate(res(0))
    end if
  end function array_scalar_raw
  
  pure subroutine array_scalar_set_from_raw(this,raw,provide_lower_bound, &
                                            provide_upper_bound)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! Assigns raw data, such as that produced by 
    ! [[array_scalar_field:raw]], to the field. The routine will
    ! stop with an error if the array is the wrong length for a field
    ! of this resolution and with these boundary conditions. 
    !
    class(array_scalar_field), intent(inout) :: this
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
    integer, dimension(:,:), allocatable :: slices
    integer, dimension(:), allocatable :: counts
    integer :: i, start, finish
    call check_set_from_raw(this,raw,provide_lower_bound,provide_upper_bound)
    slices = this%raw_slices(provide_lower_bound,provide_upper_bound)
    counts = elements_in_slice(slices(1,:),slices(2,:),slices(3,:))
    if (.not. allocated(this%field_data)) then
       allocate(this%field_data(maxval(slices(2,:))))
    end if
    do concurrent (i = 1:size(counts))
      start = sum(counts(1:i-1)) + 1
      finish = start + counts(i) - 1 
      this%field_data(slices(1,i):slices(2,i):slices(3,i)) = raw(start:finish)
    end do
  end subroutine array_scalar_set_from_raw
  
  pure function array_scalar_sf_m_sf(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm field} \times {\rm field}\)
    !
    class(array_scalar_field), intent(in) :: this
    class(scalar_field), intent(in) :: rhs
    class(scalar_field), allocatable :: res !! The restult of this operation
    class(array_scalar_field), allocatable :: local
#:if defined('DEBUG')
    call this%check_compatible(rhs)
#:endif
    allocate(local, mold=this)
    call local%assign_meta_data(this)
    select type(rhs)
    class is(array_scalar_field)
      local%field_data = this%field_data * rhs%field_data
    end select
    call move_alloc(local,res)
  end function array_scalar_sf_m_sf

  pure function array_scalar_sf_m_vf(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm field} \times {\rm \vec{field}}\)
    !
    class(array_scalar_field), intent(in) :: this
    class(vector_field), intent(in) :: rhs
    class(vector_field), allocatable :: res !! The result of this operation
    class(array_vector_field), allocatable :: local
    integer :: i
#:if defined('DEBUG')
    call this%check_compatible(rhs)
#:endif
    select type(rhs)
    class is(array_vector_field)
      allocate(local, mold=rhs)
      call local%assign_meta_data(rhs)
      do concurrent (i=1:this%numpoints)
        local%field_data(i,:) = this%field_data(i) * rhs%field_data(i,:)
      end do
    end select
    call move_alloc(local,res)
  end function array_scalar_sf_m_vf

  pure function array_scalar_r_m_sf(lhs,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm real} \times {\rm field}\)
    !
    real(r8), intent(in) :: lhs
    class(array_scalar_field), intent(in) :: rhs
    class(scalar_field), allocatable :: res !! The result of this operation
    class(array_scalar_field), allocatable :: local
    allocate(local, mold=rhs)
    call local%assign_meta_data(rhs)
    local%field_data = lhs * rhs%field_data
    call move_alloc(local,res)
  end function array_scalar_r_m_sf

  pure function array_scalar_sf_m_r(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm field} \times {\rm real}\)
    !
    class(array_scalar_field), intent(in) :: this
    real(r8), intent(in) :: rhs
    class(scalar_field), allocatable :: res !! The result of this operation
    class(array_scalar_field), allocatable :: local
    allocate(local, mold=this)
    call local%assign_meta_data(this)
    local%field_data = this%field_data * rhs
    call move_alloc(local,res)
  end function array_scalar_sf_m_r
  
  pure function array_scalar_sf_d_sf(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm field} / {\rm field}\)
    !
    class(array_scalar_field), intent(in) :: this
    class(scalar_field), intent(in) :: rhs
    class(scalar_field), allocatable :: res !! The restult of this operation
    class(array_scalar_field), allocatable :: local
#:if defined('DEBUG')
    call this%check_compatible(rhs)
#:endif
    allocate(local, mold=this)
    call local%assign_meta_data(this)
    select type(rhs)
    class is(array_scalar_field)
      local%field_data = this%field_data / rhs%field_data
    end select
    call move_alloc(local,res)
  end function array_scalar_sf_d_sf

  pure function array_scalar_r_d_sf(lhs,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm real} / {\rm field}\)
    !
    real(r8), intent(in) :: lhs
    class(array_scalar_field), intent(in) :: rhs
    class(scalar_field), allocatable :: res !! The result of this operation
    class(array_scalar_field), allocatable :: local
    allocate(local, mold=rhs)
    call local%assign_meta_data(rhs)
    local%field_data = lhs / rhs%field_data
    call move_alloc(local,res)
  end function array_scalar_r_d_sf

  pure function array_scalar_sf_d_r(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm field} / {\rm real}\)
    !
    class(array_scalar_field), intent(in) :: this
    real(r8), intent(in) :: rhs
    class(scalar_field), allocatable :: res !! The result of this operation
    class(array_scalar_field), allocatable :: local
    allocate(local, mold=this)
    call local%assign_meta_data(this)
    local%field_data = this%field_data / rhs
    call move_alloc(local,res)
  end function array_scalar_sf_d_r
  
  pure function array_scalar_sf_s_sf(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm field} - {\rm field}\)
    !
    class(array_scalar_field), intent(in) :: this
    class(scalar_field), intent(in) :: rhs
    class(scalar_field), allocatable :: res !! The restult of this operation
    class(array_scalar_field), allocatable :: local
#:if defined('DEBUG')
    call this%check_compatible(rhs)
#:endif
    allocate(local, mold=this)
    call local%assign_meta_data(this)
    select type(rhs)
    class is(array_scalar_field)
      local%field_data = this%field_data - rhs%field_data
    end select
    call move_alloc(local,res)
  end function array_scalar_sf_s_sf

  pure function array_scalar_r_s_sf(lhs,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm real} - {\rm field}\)
    !
    real(r8), intent(in) :: lhs
    class(array_scalar_field), intent(in) :: rhs
    class(scalar_field), allocatable :: res !! The result of this operation
    class(array_scalar_field), allocatable :: local
    allocate(local, mold=rhs)
    call local%assign_meta_data(rhs)
    local%field_data = lhs - rhs%field_data
    call move_alloc(local,res)
  end function array_scalar_r_s_sf

  pure function array_scalar_sf_s_r(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm field} - {\rm real}\)
    !
    class(array_scalar_field), intent(in) :: this
    real(r8), intent(in) :: rhs
    class(scalar_field), allocatable :: res !! The result of this operation
    class(array_scalar_field), allocatable :: local
    allocate(local, mold=this)
    call local%assign_meta_data(this)
    local%field_data = this%field_data - rhs
    call move_alloc(local,res)
  end function array_scalar_sf_s_r
  
  pure function array_scalar_sf_a_sf(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm field} - {\rm field}\)
    !
    class(array_scalar_field), intent(in) :: this
    class(scalar_field), intent(in) :: rhs
    class(scalar_field), allocatable :: res !! The restult of this operation
    class(array_scalar_field), allocatable :: local
#:if defined('DEBUG')
    call this%check_compatible(rhs)
#:endif
    allocate(local, mold=this)
    call local%assign_meta_data(this)
    select type(rhs)
    class is(array_scalar_field)
      local%field_data = this%field_data + rhs%field_data
    end select
    call move_alloc(local,res)
  end function array_scalar_sf_a_sf

  pure function array_scalar_r_a_sf(lhs,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm real} - {\rm field}\)
    !
    real(r8), intent(in) :: lhs
    class(array_scalar_field), intent(in) :: rhs
    class(scalar_field), allocatable :: res !! The result of this operation
    class(array_scalar_field), allocatable :: local
    allocate(local, mold=rhs)
    call local%assign_meta_data(rhs)
    local%field_data = lhs + rhs%field_data
    call move_alloc(local,res)
  end function array_scalar_r_a_sf

  pure function array_scalar_sf_a_r(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm field} - {\rm real}\)
    !
    class(array_scalar_field), intent(in) :: this
    real(r8), intent(in) :: rhs
    class(scalar_field), allocatable :: res !! The result of this operation
    class(array_scalar_field), allocatable :: local
    allocate(local, mold=this)
    call local%assign_meta_data(this)
    local%field_data = this%field_data + rhs
    call move_alloc(local,res)
  end function array_scalar_sf_a_r

  pure function array_scalar_sf_p_r(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm field}^{\rm real}\)
    !
    class(array_scalar_field), intent(in) :: this
    real(r8), intent(in) :: rhs
    class(scalar_field), allocatable :: res !! The result of this operation
    class(array_scalar_field), allocatable :: local
    allocate(local, mold=this)
    call local%assign_meta_data(this)
    local%field_data = this%field_data ** rhs
    call move_alloc(local,res)
  end function array_scalar_sf_p_r

  pure function array_scalar_sf_p_r4(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm field}^{\rm real}\)
    !
    class(array_scalar_field), intent(in) :: this
    real, intent(in) :: rhs
    class(scalar_field), allocatable :: res !! The result of this operation
    class(array_scalar_field), allocatable :: local
    allocate(local, mold=this)
    call local%assign_meta_data(this)
    local%field_data = this%field_data ** rhs
    call move_alloc(local,res)
  end function array_scalar_sf_p_r4

  pure function array_scalar_sf_p_i(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm field}^{\rm integer}\)
    !
    class(array_scalar_field), intent(in) :: this
    integer, intent(in) :: rhs
    class(scalar_field), allocatable :: res !! The result of this operation
    class(array_scalar_field), allocatable :: local
    allocate(local, mold=this)
    call local%assign_meta_data(this)
    local%field_data = this%field_data ** rhs
    call move_alloc(local,res)
  end function array_scalar_sf_p_i

#:for FUNC, TEX in UNARY_FUNCTIONS
  $:unary_func(FUNC, TEX, 'array_scalar')

#:endfor

  pure function array_scalar_minval(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\min({\rm field})\)
    !
    class(array_scalar_field), intent(in) :: this
    real(r8) :: res !! The result of this operation
    if (allocated(this%field_data)) then
      res = minval(this%field_data)
    else
      res = 0.0_r8
    end if
  end function array_scalar_minval

  pure function array_scalar_maxval(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\max({\rm field})\)
    !
    class(array_scalar_field), intent(in) :: this
    real(r8) :: res !! The result of this operation
    if (allocated(this%field_data)) then
      res = maxval(this%field_data)
    else
      res = 0.0_r8
    end if
  end function array_scalar_maxval

  function array_scalar_d_dx(this, dir, order) result(res)
    !* Author: Chris MacMackin
    !  Date: October 2016
    !
    ! \(\frac{\partial^{\rm order}}{\partial x_{\rm dir}^{\rm order}}{\rm field}\)
    !
    class(array_scalar_field), intent(in) :: this
    integer, intent(in) :: dir !! Direction in which to differentiate
    integer, optional, intent(in) :: order !! Order of the derivative, default = 1
    class(scalar_field), allocatable :: res
    class(array_scalar_field), allocatable :: local
    allocate(local, mold=this)
    call local%assign_meta_data(this)
    local%field_data = this%array_dx(this%field_data, dir, order)
    call move_alloc(local, res)
  end function array_scalar_d_dx

  function array_scalar_laplacian(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\nabla^2 {\rm field}\)
    !
    class(array_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    integer :: i
    call this%allocate_scalar_field(res)
    call res%assign_meta_data(this)
    select type(res)
    class is(array_scalar_field)
      res%field_data = this%array_dx(this%field_data,1,2)
      do i = 2, this%dimensions()
        res%field_data = res%field_data + this%array_dx(this%field_data,i,2)
      end do
    class default
      error stop('Non-array_scalar_field type allocated by '//&
                 '`allocate_scalar_field` routine.')
    end select
  end function array_scalar_laplacian
  
  function array_scalar_gradient(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\nabla{\rm field}\)
    !
    class(array_scalar_field), intent(in) :: this
    class(vector_field), allocatable :: res !! The result of this operation
    class(array_vector_field), allocatable :: local
    integer :: i
    call this%allocate_vector_field(res)
    call res%assign_meta_data(this, .false.)
    select type(res)
    class is(array_vector_field)
      allocate(res%field_data(size(this%field_data),this%dimensions()))
      do i = 1, this%dimensions()
        res%field_data(:,i) = this%array_dx(this%field_data,i,1)
      end do
    class default
      error stop('Non-array_vector_field type allocated by '//&
                 '`allocate_vector_field` routine.')
    end select
  end function array_scalar_gradient
  
  elemental subroutine array_scalar_assign(this,rhs)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm field} = {\rm field}\)
    !
    class(array_scalar_field), intent(inout) :: this
    class(scalar_field), intent(in) :: rhs
    select type(rhs)
    class is(array_scalar_field)
      call this%assign_meta_data(rhs)
      if (allocated(rhs%field_data)) this%field_data = rhs%field_data
    end select
  end subroutine array_scalar_assign


  pure logical function array_scalar_is_equal(this,rhs) result(iseq)
    !* Author: Chris MacMackin
    !  Date: April 2016
    !
    ! Evaluates whether two scalar fields are equal within a tolerance,
    ! specified by [[set_tol]].
    !
    class(array_scalar_field), intent(in) :: this
    class(scalar_field), intent(in) :: rhs
    integer :: i
    real(r8) :: normalization
    iseq = .true.
    select type(rhs)
    class is(array_scalar_field)
      if (this%numpoints /= rhs%numpoints) then
        iseq = .false.
        return
      end if
      do i=1,this%numpoints
        normalization = abs(this%field_data(i))
        if (normalization < tiny(normalization)) normalization = 1.0_r8
        iseq = iseq .and.( ((this%field_data(i)-rhs%field_data(i))/normalization < &
                             get_tol()) .or. (is_nan(this%field_data(i)).and. &
                                              is_nan(rhs%field_data(i))) )
        if (.not. iseq) return
      end do
    class default
      iseq = .false.
    end select
  end function array_scalar_is_equal

  pure subroutine array_scalar_assign_meta_data(this, rhs, alloc)
    !* Author: Chris MacMackin
    !  Date: April 2016
    !
    ! Assigns all data in `rhs` to `this` except for the actual values
    ! of the field.
    !
    class(array_scalar_field), intent(inout) :: this
    class(abstract_field), intent(in) :: rhs
    logical, optional, intent(in) :: alloc
      !! If present and false, do not allocate the array of `this`.
    call this%assign_subtype_meta_data(rhs)
    select type(rhs)
    class is(array_scalar_field)
      this%numpoints = rhs%numpoints
      if (allocated(this%field_data)) deallocate(this%field_data)
      if (present(alloc)) then
        if (allocated(rhs%field_data) .and. alloc) &
          allocate(this%field_data(size(rhs%field_data))) !Needed due to compiler bug
      else
        allocate(this%field_data(size(rhs%field_data))) !Needed due to compiler bug
      end if
    class is(array_vector_field)
      this%numpoints = rhs%numpoints
      if (allocated(this%field_data)) deallocate(this%field_data)
      if (present(alloc)) then
        if (allocated(rhs%field_data) .and. alloc) &
          allocate(this%field_data(size(rhs%field_data,1))) !Needed due to compiler bug
      else
        allocate(this%field_data(size(rhs%field_data,1))) !Needed due to compiler bug
      end if
    end select
  end subroutine array_scalar_assign_meta_data

  pure subroutine array_scalar_compatible(this,other)
    !* Author: Chris MacMackin
    !  Date: October 2016
    !
    ! Checks whether a field has the same type, boundaries, and
    ! resolution as this one, making it compatible for binary 
    ! operations. If incompatible, stops program with an error message.
    !
    class(array_scalar_field), intent(in) :: this
    class(abstract_field), intent(in) :: other
      !! The field being checked against this one
    character(len=69), parameter :: err_message = 'array_scalar_field: '//&
         'Error, operation with incompatible fields due to'//new_line('a')
    call this%check_subtype_compatible(other)
    select type(other)
    class is(array_scalar_field)
      if (this%numpoints /= other%numpoints) &
           error stop(err_message//'    different resolutions.')
      if (.not.(allocated(this%field_data).and.allocated(other%field_data))) &
           error stop(err_message//'    uninitialised fields.')
    class is(array_vector_field)
      if (this%numpoints /= other%numpoints) &
           error stop(err_message//'    different resolutions.')
      if (.not.(allocated(this%field_data).and.allocated(other%field_data))) &
           error stop(err_message//'    uninitialised fields.')
    class default
      error stop(err_message//'    incompatible types.')
    end select
  end subroutine array_scalar_compatible


  !=====================================================================
  ! Vector Field Methods
  !=====================================================================

  function array_vector_constructor(template,numpoints,vector_dims, &
                                    initializer) result(this)
    !* Author: Chris MacMackin
    !  Date: October 2016
    !
    ! Creates a new field with the same concrete type as the template
    ! argument. The array of values will be allocated an initiated.
    !
    class(array_vector_field), intent(in)  :: template
    !! A scalar field object which will act as a mold for the concrete
    !! type of the returned type.
    integer, intent(in)                    :: numpoints
    !! The number of data points needed in the array when modelling this
    !! field.
    integer, intent(in)                    :: vector_dims
    !! The number of components of vectors in this field
    procedure(vector_init), optional       :: initializer
    !! An elemental procedure taking which takes the position in the
    !! fields domain (an 8-byte real) as an argument and returns the
    !! fields value at that position. Default is for field to be zero
    !! everywhere.
    class(array_vector_field), allocatable :: this
    !! A scalar field initiated based on the arguments to this function.
    integer :: i
    allocate(this, source=template)
    if (allocated(this%field_data)) deallocate(this%field_data)
    allocate(this%field_data(numpoints,vector_dims))
    this%numpoints = numpoints
    this%vector_dims = vector_dims
    if (present(initializer)) then
      do i = 1, numpoints
        this%field_data(i,:) = initializer(this%id_to_position(i))
      end do
    end if
  end function array_vector_constructor

  elemental function array_vector_vector_dimensions(this) result(dims)
    !* Author: Chris MacMackin
    !  Date: October 2016
    !
    ! Returns the number of dimensions/components are in the vectors
    ! of this field.
    !
    class(array_vector_field), intent(in) :: this
    integer :: dims !! Number of vector components
    dims = this%vector_dims
  end function array_vector_vector_dimensions

  pure function array_vector_elements(this) result(elements)
    !* Author: Chris MacMackin
    !  Date: October 2016
    !
    ! Gives the number of individual data points present in the field.
    !
    class(array_vector_field), intent(in) :: this
    integer :: elements
    elements = this%numpoints
  end function array_vector_elements

  pure function array_vector_raw_size(this,return_lower_bound, &
                                       return_upper_bound) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! Compute how many elements are in the raw representation of this
    ! field. This would be the number of data points, adjusted based on
    ! how boundary conditions are accounted for.
    !
    class(array_vector_field), intent(in)       :: this
    logical, dimension(:), optional, intent(in) :: return_lower_bound
      !! Specifies whether to return the values at the lower boundary
      !! for each dimension, with the index of the element
      !! corresponding to the dimension. Defaults to all `.true.`.
    logical, dimension(:), optional, intent(in) :: return_upper_bound
      !! Specifies whether to return the values at the upper boundary
      !! for each dimension, with the index of the element
      !! corresponding to the dimension. Defaults to all `.true.`.
    integer :: res
    integer, dimension(:,:), allocatable :: slices
    integer :: i
    if (allocated(this%field_data)) then
      slices = this%raw_slices(return_lower_bound,return_upper_bound)
      res = sum(elements_in_slice(slices(1,:),slices(2,:),slices(3,:))) * &
            this%vector_dims
    else
      res = 0
    end if
  end function array_vector_raw_size
  
  pure function array_vector_raw(this,return_lower_bound, &
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
    class(array_vector_field), intent(in) :: this
    logical, dimension(:), optional, intent(in) :: return_lower_bound
      !! Specifies whether to return the values at the lower boundary
      !! for each dimension, with the index of the element
      !! corresponding to the dimension. Defaults to all `.true.`.
    logical, dimension(:), optional, intent(in) :: return_upper_bound
      !! Specifies whether to return the values at the upper boundary
      !! for each dimension, with the index of the element
      !! corresponding to the dimension. Defaults to all `.true.`.
    real(r8), dimension(:), allocatable :: res
      !! Array containing data needed to describe field
    integer, dimension(:,:), allocatable :: slices
    integer :: i
    slices = this%raw_slices(return_lower_bound,return_upper_bound)
    if (allocated(this%field_data)) then
      res = [(this%field_data(slices(1,i):slices(2,i):slices(3,i),:), i=1, &
             size(slices,2))]
    else
      allocate(res(0))
    end if
  end function array_vector_raw
  
  pure subroutine array_vector_set_from_raw(this,raw,provide_lower_bound, &
                                            provide_upper_bound)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! Assigns raw data, such as that produced by 
    ! [[array_scalar_field:raw]], to the field. The routine will
    ! stop with an error if the array is the wrong length for a field
    ! of this resolution and with these boundary conditions. 
    !
    class(array_vector_field), intent(inout) :: this
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
    integer, dimension(:,:), allocatable :: slices
    integer, dimension(:), allocatable :: counts
    integer :: i, start, finish
    call check_set_from_raw(this,raw,provide_lower_bound,provide_upper_bound)
    slices = this%raw_slices(provide_lower_bound,provide_upper_bound)
    counts = elements_in_slice(slices(1,:),slices(2,:),slices(3,:)) * this%vector_dims
    if (.not. allocated(this%field_data)) then
       allocate(this%field_data(maxval(slices(2,:)),this%vector_dims))
    end if
    do concurrent (i = 1:size(counts))
      start = sum(counts(1:i-1)) + 1
      finish = start + counts(i) - 1
      this%field_data(slices(1,i):slices(2,i):slices(3,i),:) = reshape(raw(start:finish), &
           [counts(i)/this%vector_dims,this%vector_dims])
    end do
  end subroutine array_vector_set_from_raw

  pure function array_vector_vf_m_sf(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\vec{\rm field} \times {\rm field}\)
    !
    class(array_vector_field), intent(in) :: this
    class(scalar_field), intent(in) :: rhs
    class(vector_field), allocatable :: res !! The restult of this operation
    class(array_vector_field), allocatable :: local
    integer :: i
#:if defined('DEBUG')
    call this%check_compatible(rhs)
#:endif
    allocate(local, mold=this)
    call local%assign_meta_data(this)
    select type(rhs)
    class is(array_scalar_field)
      do concurrent (i=1:this%vector_dims)
        local%field_data(:,i) = this%field_data(:,i) * rhs%field_data
      end do
    end select
    call move_alloc(local,res)
  end function array_vector_vf_m_sf

  pure function array_vector_r_m_vf(lhs,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm real} \times \vec{\rm field}\)
    !
    real(r8), intent(in) :: lhs
    class(array_vector_field), intent(in) :: rhs
    class(vector_field), allocatable :: res !! The result of this operation
    class(array_vector_field), allocatable :: local
    allocate(local, mold=rhs)
    call local%assign_meta_data(rhs)
    local%field_data = lhs * rhs%field_data
    call move_alloc(local,res)
  end function array_vector_r_m_vf

  pure function array_vector_vf_m_r(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\vec{\rm field} \times {\rm real}\)
    !
    class(array_vector_field), intent(in) :: this
    real(r8), intent(in) :: rhs
    class(vector_field), allocatable :: res !! The result of this operation
    class(array_vector_field), allocatable :: local
    allocate(local, mold=this)
    call local%assign_meta_data(this)
    local%field_data = this%field_data * rhs
    call move_alloc(local,res)
  end function array_vector_vf_m_r
  
  pure function array_vector_vf_d_sf(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\vec{\rm field} / {\rm field}\)
    !
    class(array_vector_field), intent(in) :: this
    class(scalar_field), intent(in) :: rhs
    class(vector_field), allocatable :: res !! The restult of this operation
    class(array_vector_field), allocatable :: local
    integer :: i
#:if defined('DEBUG')
    call this%check_compatible(rhs)
#:endif
    allocate(local, mold=this)
    call local%assign_meta_data(this)
    select type(rhs)
    class is(array_scalar_field)
      do concurrent (i=1:this%vector_dims)
        local%field_data(:,i) = this%field_data(:,i) / rhs%field_data
      end do
    end select
    call move_alloc(local,res)
  end function array_vector_vf_d_sf

  pure function array_vector_vf_d_r(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\vec{\rm field} / {\rm real}\)
    !
    class(array_vector_field), intent(in) :: this
    real(r8), intent(in) :: rhs
    class(vector_field), allocatable :: res !! The result of this operation
    class(array_vector_field), allocatable :: local
    allocate(local, mold=this)
    call local%assign_meta_data(this)
    local%field_data = this%field_data / rhs
    call move_alloc(local,res)
  end function array_vector_vf_d_r
  
  pure function array_vector_vf_s_vf(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\vec{\rm field} - \vec{\rm field}\)
    !
    class(array_vector_field), intent(in) :: this
    class(vector_field), intent(in) :: rhs
    class(vector_field), allocatable :: res !! The restult of this operation
    class(array_vector_field), allocatable :: local
    integer :: min_dims, max_dims
#:if defined('DEBUG')
    call this%check_compatible(rhs)
#:endif
    allocate(local, mold=this)
    select type(rhs)
    class is(array_vector_field)
      if (rhs%vector_dims > this%vector_dims) then
        call local%assign_meta_data(rhs)
      else
        call local%assign_meta_data(this)
      end if
      min_dims = min(this%vector_dims,rhs%vector_dims)
      max_dims = max(this%vector_dims,rhs%vector_dims)
      local%field_data(:,1:min_dims) = this%field_data(:,1:min_dims) &
                                     - rhs%field_data(:,1:min_dims)
      if (rhs%vector_dims > this%vector_dims) then
        local%field_data(:,min_dims+1:max_dims) = -rhs%field_data(:,min_dims+1:max_dims)
      else
        local%field_data(:,min_dims+1:max_dims) = this%field_data(:,min_dims+1:max_dims)
      end if
    end select
    call move_alloc(local,res)
  end function array_vector_vf_s_vf

  pure function array_vector_r_s_vf(lhs,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\vec{\rm real} - \vec{\rm field}\)
    !
    real(r8), dimension(:), intent(in) :: lhs
    class(array_vector_field), intent(in) :: rhs
    class(vector_field), allocatable :: res !! The result of this operation
    class(array_vector_field), allocatable :: local
    integer :: min_dims, max_dims, i
    allocate(local, mold=rhs)
    call local%assign_meta_data(rhs,.false.)
    min_dims = min(rhs%vector_dims, size(lhs))
    max_dims = max(rhs%vector_dims, size(lhs))
    allocate(local%field_data(local%numpoints, max_dims))
    do concurrent (i=1:local%numpoints)
      local%field_data(i,:min_dims) = lhs(:min_dims) - rhs%field_data(i,:min_dims)
    end do
    if (rhs%vector_dims > size(lhs)) then
      local%field_data(:,min_dims+1:) = -rhs%field_data(:,min_dims+1:)
    else
      local%vector_dims = size(lhs)
      do concurrent (i=1:local%numpoints)
        local%field_data(i,min_dims+1:) = lhs(min_dims+1:)
      end do
    end if
    call move_alloc(local,res)
  end function array_vector_r_s_vf

  pure function array_vector_vf_s_r(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\vec{\rm field} - \vec{\rm real}\)
    !
    class(array_vector_field), intent(in) :: this
    real(r8), dimension(:), intent(in) :: rhs
    class(vector_field), allocatable :: res !! The result of this operation
    class(array_vector_field), allocatable :: local
    integer :: min_dims, max_dims, i
    allocate(local, mold=this)
    call local%assign_meta_data(this,.false.)
    min_dims = min(this%vector_dims, size(rhs))
    max_dims = max(this%vector_dims, size(rhs))
    allocate(local%field_data(local%numpoints, max_dims))
    do concurrent (i=1:local%numpoints)
      local%field_data(i,:min_dims) = this%field_data(i,:min_dims) - rhs(:min_dims)
    end do
    if (this%vector_dims > size(rhs)) then
      local%field_data(:,min_dims+1:) = this%field_data(:,min_dims+1:)
    else
      local%vector_dims = size(rhs)
      do concurrent (i=1:local%numpoints)
        local%field_data(i,min_dims+1:) = -rhs(min_dims+1:)
      end do
    end if
    call move_alloc(local,res)
  end function array_vector_vf_s_r
  
  pure function array_vector_vf_a_vf(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\vec{\rm field} - \vec{\rm field}\)
    !
    class(array_vector_field), intent(in) :: this
    class(vector_field), intent(in) :: rhs
    class(vector_field), allocatable :: res !! The restult of this operation
    class(array_vector_field), allocatable :: local
    integer :: max_dims, min_dims
#:if defined('DEBUG')
    call this%check_compatible(rhs)
#:endif
    allocate(local, mold=this)
    select type(rhs)
    class is(array_vector_field)
      if (rhs%vector_dims > this%vector_dims) then
        call local%assign_meta_data(rhs)
      else
        call local%assign_meta_data(this)
      end if
      min_dims = min(this%vector_dims,rhs%vector_dims)
      max_dims = max(this%vector_dims,rhs%vector_dims)
      local%field_data(:,1:min_dims) = this%field_data(:,1:min_dims) &
                                     + rhs%field_data(:,1:min_dims)
      if (rhs%vector_dims > this%vector_dims) then
        local%field_data(:,min_dims+1:max_dims) = rhs%field_data(:,min_dims+1:max_dims)
      else
        local%field_data(:,min_dims+1:max_dims) = this%field_data(:,min_dims+1:max_dims)
      end if
    end select
    call move_alloc(local,res)
  end function array_vector_vf_a_vf

  pure function array_vector_r_a_vf(lhs,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\vec{\rm real} - \vec{\rm field}\)
    !
    real(r8), dimension(:), intent(in) :: lhs
    class(array_vector_field), intent(in) :: rhs
    class(vector_field), allocatable :: res !! The result of this operation
    class(array_vector_field), allocatable :: local
    integer :: min_dims, max_dims, i
    allocate(local, mold=rhs)
    call local%assign_meta_data(rhs,.false.)
    min_dims = min(rhs%vector_dims, size(lhs))
    max_dims = max(rhs%vector_dims, size(lhs))
    allocate(local%field_data(local%numpoints, max_dims))
    do concurrent (i=1:local%numpoints)
      local%field_data(i,:min_dims) = lhs(:min_dims) + rhs%field_data(i,:min_dims)
    end do
    if (rhs%vector_dims > size(lhs)) then
      local%field_data(:,min_dims+1:) = rhs%field_data(:,min_dims+1:)
    else
      local%vector_dims = size(lhs)
      do concurrent (i=1:local%numpoints)
        local%field_data(i,min_dims+1:) = lhs(min_dims+1:)
      end do
    end if
    call move_alloc(local,res)
  end function array_vector_r_a_vf

  pure function array_vector_vf_a_r(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\vec{\rm field} - \vec{\rm real}\)
    !
    class(array_vector_field), intent(in) :: this
    real(r8), dimension(:), intent(in) :: rhs
    class(vector_field), allocatable :: res !! The result of this operation
    class(array_vector_field), allocatable :: local
    integer :: min_dims, max_dims, i
    allocate(local, mold=this)
    call local%assign_meta_data(this,.false.)
    min_dims = min(this%vector_dims, size(rhs))
    max_dims = max(this%vector_dims, size(rhs))
    allocate(local%field_data(local%numpoints, max_dims))
    do concurrent (i=1:local%numpoints)
      local%field_data(i,:min_dims) = this%field_data(i,:min_dims) + rhs(:min_dims)
    end do
    if (this%vector_dims > size(rhs)) then
      local%field_data(:,min_dims+1:) = this%field_data(:,min_dims+1:)
    else
      local%vector_dims = size(rhs)
      do concurrent (i=1:local%numpoints) 
        local%field_data(i,min_dims+1:) = rhs(min_dims+1:)
      end do
    end if
    call move_alloc(local,res)
  end function array_vector_vf_a_r

  elemental subroutine array_vector_assign(this,rhs)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\vec{\rm field} = \vec{\rm field}\)
    !
    class(array_vector_field), intent(inout) :: this
    class(vector_field), intent(in) :: rhs
    call this%assign_meta_data(rhs)
    select type(rhs)
    class is(array_vector_field)
      if (allocated(rhs%field_data)) this%field_data = rhs%field_data
    end select
  end subroutine array_vector_assign

  pure function array_vector_norm(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\lVert \vec{\rm field} \rVert\)
    !
    class(array_vector_field), intent(in) :: this
    class(scalar_field), allocatable :: res
    class(array_scalar_field), allocatable :: local
    call this%allocate_scalar_field(res)
    call res%assign_meta_data(this)
    select type(res)
    class is(array_scalar_field)
      res%field_data = sqrt(sum(this%field_data**2,2))
    class default
      error stop('Non-array_scalar_field type allocated by '//&
                 '`allocate_scalar_field` routine.')
    end select
  end function array_vector_norm

  pure function array_vector_component(this,comp) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! Returns a calar field containing specified component of the vector field
    !
    class(array_vector_field), intent(in) :: this
    integer, intent(in) :: comp
    class(scalar_field), allocatable :: res
    call this%allocate_scalar_field(res)
    call res%assign_meta_data(this)
    select type(res)
    class is(array_scalar_field)
      if (comp <= this%vector_dims) then
        res%field_data = this%field_data(:,comp)
      else
        res%field_data = 0.0_r8
      end if
    class default
      error stop('Non-array_scalar_field type allocated by '//&
                 '`allocate_scalar_field` routine.')
    end select
  end function array_vector_component

  function array_vector_d_dx(this, dir, order) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\frac{\partial^{\rm order}}{\partial x_{\rm dir}^{\rm order}}\vec{\rm field}\)
    !
    class(array_vector_field), intent(in) :: this
    integer, intent(in) :: dir !! Direction in which to differentiate
    integer, optional, intent(in) :: order !! Order of the derivative, default = 1
    class(vector_field), allocatable :: res
    class(array_vector_field), allocatable :: local
    integer :: i
    allocate(local, mold=this)
    allocate(local%field_data, mold=this%field_data)
    call local%assign_meta_data(this)
    do i = 1, this%vector_dims
      local%field_data(:,i) = this%array_dx(this%field_data(:,i), &
                                                     dir, order)
    end do
    call move_alloc(local, res)
  end function array_vector_d_dx

  function array_vector_component_d_dx(this, dir, component, order) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\frac{\partial^{\rm order}}{\partial x_{\rm dir}^{\rm order}}{\rm field_{component}}\)
    !
    class(array_vector_field), intent(in) :: this
    integer, intent(in) :: dir !! Direction in which to differentiate
    integer, intent(in) :: component !! Which component of the vector is being differentiated
    integer, optional, intent(in) :: order !! Order of the derivative, default = 1
    class(scalar_field), allocatable :: res
    call this%allocate_scalar_field(res)
    call res%assign_meta_data(this)
    select type(res)
    class is(array_scalar_field)
      res%field_data(:) = this%array_dx(this%field_data(:,component), &
                                                 dir, order)
    class default
      error stop('Non-array_scalar_field type allocated by '//&
                 '`allocate_scalar_field` routine.')
    end select
  end function array_vector_component_d_dx

  function array_vector_laplacian(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\nabla^2 \vec{\rm field}\)
    !
    class(array_vector_field), intent(in) :: this
    class(vector_field), allocatable :: res !! The result of this operation
    integer :: i, j
    call this%allocate_vector_field(res)
    call res%assign_meta_data(this, .false.)
    select type(res)
    class is(array_vector_field)
      allocate(res%field_data(size(this%field_data,1),this%vector_dims))
      do i = 1, this%vector_dims
        res%field_data(:,i) = this%array_dx(this%field_data(:,i),1,2)
        do j = 2, this%dimensions()
          res%field_data(:,i) = res%field_data(:,i) + &
                              this%array_dx(this%field_data(:,i),j,2)
        end do
      end do
    class default
      error stop('Non-array_vector_field type allocated by '//&
                 '`allocate_vector_field` routine.')
    end select
    call this%allocate_vector_field(res)
    call res%assign_meta_data(this)
  end function array_vector_laplacian
  
  function array_vector_divergence(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\nabla\cdot \vec{\rm field}\)
    !
    class(array_vector_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    integer :: i
    call this%allocate_scalar_field(res)
    call res%assign_meta_data(this)
    select type(res)
    class is(array_scalar_field)
      res%field_data = this%array_dx(this%field_data(:,1),1,1)
      do i = 2, this%dimensions()
        res%field_data = res%field_data + &
                         this%array_dx(this%field_data(:,i),i,1)
      end do
    class default
      error stop('Non-array_scalar_field type allocated by '//&
                 '`allocate_scalar_field` routine.')
    end select
  end function array_vector_divergence
  
  function array_vector_curl(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\nabla\times \vec{\rm field}\)
    !
    class(array_vector_field), intent(in) :: this
    class(vector_field), allocatable :: res !! The result of this operation
    logical, dimension(3) :: been_set
    integer :: i
    been_set = .false.
    call this%allocate_vector_field(res)
    call res%assign_meta_data(this,.false.)
    select type(res)
    class is(array_vector_field)
      res%vector_dims = 3
      allocate(res%field_data(res%numpoints,3))
      if (this%dimensions() >= 3) then
        res%field_data(:,2) = this%array_dx(this%field_data(:,1),3)
        been_set(2) = .true.
      end if
      if (this%dimensions() >= 2) then
        res%field_data(:,3) = -this%array_dx(this%field_data(:,1),2)
        been_set(3) = .true.
      end if
      if (this%vector_dims >= 2) then
        if (this%dimensions() >= 3) then
          res%field_data(:,1) = -this%array_dx(this%field_data(:,2),3)
          been_set(1) = .true.
        end if
        if (been_set(3)) then
          res%field_data(:,3) = res%field_data(:,3) &
                             + this%array_dx(this%field_data(:,2),1)
        else
          res%field_data(:,3) = this%array_dx(this%field_data(:,2),1)
          been_set(3) = .true.
        end if
      end if
      if (this%vector_dims >= 3) then
        if (this%dimensions() >= 2) then
          if (been_set(1)) then
            res%field_data(:,1) = res%field_data(:,1) &
                               + this%array_dx(this%field_data(:,3),2)
          else
            res%field_data(:,1) = this%array_dx(this%field_data(:,3),2)
            been_set(1) = .true.
          end if
        end if
        if (been_set(2)) then
          res%field_data(:,2) = res%field_data(:,2) &
                             - this%array_dx(this%field_data(:,3),1)
        else
          res%field_data(:,2) = -this%array_dx(this%field_data(:,3),1)
          been_set(2) = .true.
        end if
      end if
      do i = 1, 3
        if (.not. been_set(i)) res%field_data(:,i) = 0.0_r8
      end do
    class default
      error stop('Non-array_scalar_field type allocated by '//&
                 'allocate_scalar_field routine.')
    end select
  end function array_vector_curl

  pure function array_vector_cross_prod(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\vec{\rm field} \times \vec{\rm field}\)
    !
    ! The returned vector will always consist of three components,
    ! regardless of the number of components in the arguments.
    !
    class(array_vector_field), intent(in) :: this
    class(vector_field), intent(in) :: rhs
    class(vector_field), allocatable :: res !! The restult of this operation
    real(r8), dimension(3) :: vec1, vec2
    integer :: i, dims1, dims2
    class(array_vector_field), allocatable :: local
#:if defined('DEBUG')
    call this%check_compatible(rhs)
#:endif
    allocate(local, mold=this)
    call local%assign_meta_data(this,.false.)
    local%vector_dims = 3
    allocate(local%field_data(this%numpoints,3))
    vec1 = 0
    vec2 = 0
    select type(rhs)
    class is(array_vector_field)
      dims1 = min(3,this%vector_dims)
      dims2 = min(3,rhs%vector_dims)
      do i=1,this%numpoints
        vec1(:dims1) = this%field_data(i,:dims1)
        vec2(:dims2) = rhs%field_data(i,:dims2)
        local%field_data(i,:) = [vec1(2)*vec2(3) - vec2(2)*vec1(3), &
                                 vec1(3)*vec2(1) - vec2(3)*vec1(1), &
                                 vec1(1)*vec2(2) - vec2(1)*vec1(2)]
      end do
    end select
    call move_alloc(local,res)
  end function array_vector_cross_prod

  pure function array_vector_dot_prod(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\vec{\rm field} \cdot \vec{\rm field}\)
    !
    class(array_vector_field), intent(in) :: this
    class(vector_field), intent(in) :: rhs
    class(scalar_field), allocatable :: res !! The restult of this operation
    integer :: min_dims
#:if defined('DEBUG')
    call this%check_compatible(rhs)
#:endif
    call this%allocate_scalar_field(res)
    call res%assign_meta_data(this)
    select type(res)
    class is(array_scalar_field)
      select type(rhs)
      class is(array_vector_field)
        min_dims = min(this%vector_dims,rhs%vector_dims)
        res%field_data = &
          sum(this%field_data(:,1:min_dims)*rhs%field_data(:,1:min_dims),2)
        min_dims = min(this%vector_dims,rhs%vector_dims)
        res%field_data = &
          sum(this%field_data(:,1:min_dims)*rhs%field_data(:,1:min_dims),2)
      end select
    class default
      error stop('Non-array_scalar_field type allocated by '//&
                 'allocate_scalar_field routine.')
    end select
  end function array_vector_dot_prod

  pure logical function array_vector_is_equal(this,rhs) result(iseq)
    !* Author: Chris MacMackin
    !  Date: April 2016
    !
    ! Evaluates whether two vector fields are equal within a tolerance,
    ! specified by [[set_tol]].
    !
    class(array_vector_field), intent(in) :: this
    class(vector_field), intent(in) :: rhs
    integer :: i, dims
    real(r8) :: normalization
    iseq = .true.
    select type(rhs)
    class is(array_vector_field)
      if (this%numpoints/=rhs%numpoints) then
        iseq = .false.
        return
      end if
      if (this%vector_dims > rhs%vector_dims) then
        dims = rhs%vector_dims
        iseq = all(this%field_data(:,dims:)==0.0_r8)
        if (.not. iseq) return
      else if (this%vector_dims < rhs%vector_dims) then
        dims = this%vector_dims
        iseq = all(rhs%field_data(:,dims:)==0.0_r8)
        if (.not. iseq) return
      else
        dims = this%vector_dims
      end if
      do i=1,this%numpoints
        normalization = norm2(this%field_data(i,:dims))
        if (normalization < tiny(normalization)) normalization = 1.0_r8
        iseq = iseq .and. ( (norm2(this%field_data(i,this%dimensions()+1:dims) &
             - rhs%field_data(i,rhs%dimensions()+1:dims))/normalization < get_tol()) &
             .or. (is_nan(norm2(this%field_data(i,this%dimensions()+1:dims))).and. &
                   is_nan(norm2(rhs%field_data(i,rhs%dimensions()+1:dims)))) )
        if (.not. iseq) return
      end do
    class default
      iseq = .false.
    end select
  end function array_vector_is_equal

  pure subroutine array_vector_assign_meta_data(this, rhs, alloc)
    !* Author: Chris MacMackin
    !  Date: April 2016
    !
    ! Assigns all data in `rhs` to `this` except for the actual values
    ! of the field.
    !
    class(array_vector_field), intent(inout) :: this
    class(abstract_field), intent(in) :: rhs
    logical, optional, intent(in) :: alloc
      !! If present and false, do not allocate the array of `this`.
    call this%assign_subtype_meta_data(rhs)
    select type(rhs)
    class is(array_scalar_field)
      this%numpoints = rhs%numpoints
      this%vector_dims = this%dimensions()
      if (allocated(this%field_data)) deallocate(this%field_data)
      if (present(alloc)) then
      !Needed due to compiler bug:
        if (allocated(rhs%field_data) .and. alloc) &
          allocate(this%field_data(size(rhs%field_data),1))
      else
        allocate(this%field_data(size(rhs%field_data),1))
      end if
    class is(array_vector_field)
      this%numpoints = rhs%numpoints
      this%vector_dims = rhs%vector_dims
      if (allocated(this%field_data)) deallocate(this%field_data)
      if (present(alloc)) then
        !Needed due to compiler bug:
        if (allocated(rhs%field_data) .and. alloc) &
          allocate(this%field_data(size(rhs%field_data,1),size(rhs%field_data,2)))
      else
        allocate(this%field_data(size(rhs%field_data,1),size(rhs%field_data,2)))
      end if
    end select
  end subroutine array_vector_assign_meta_data

  pure subroutine array_vector_compatible(this,other)
    !* Author: Chris MacMackin
    !  Date: October 2016
    !
    ! Checks whether a field has the same type, boundaries, and
    ! resolution as this one, making it compatible for binary 
    ! operations. If incompatible, stops program with an error message.
    !
    class(array_vector_field), intent(in) :: this
    class(abstract_field), intent(in) :: other
      !! The field being checked against this one
    character(len=69), parameter :: err_message = 'array_vector_field: '//&
         'Error, operation with incompatible fields due to'//new_line('a')
    call this%check_subtype_compatible(other)
    select type(other)
    class is(array_scalar_field)
      if (this%numpoints /= other%numpoints) &
           error stop(err_message//'    different resolutions.')
      if (.not.(allocated(this%field_data).and.allocated(other%field_data))) &
           error stop(err_message//'    uninitialised fields.')
    class is(array_vector_field)
      if (this%numpoints /= other%numpoints) &
           error stop(err_message//'    different resolutions.')
      if (.not.(allocated(this%field_data).and.allocated(other%field_data))) &
           error stop(err_message//'    uninitialised fields.')
    class default
      error stop(err_message//'    incompatible types.')
    end select
  end subroutine array_vector_compatible

end module array_fields_mod