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
  use uniform_fields_mod, only: uniform_scalar_field, uniform_vector_field
  use utils_mod, only: is_nan, check_set_from_raw, elements_in_slice
  use h5lt, only: hid_t, hsize_t, h5ltmake_dataset_double_f, &
                  h5ltread_dataset_double_f
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
    integer                             :: numpoints
      !! The number of datapoints used
    real(r8), dimension(:), allocatable :: field_data
      !! The value of the scalar field at the data-points. Each
      !! element represents the value at a different location.
    logical                             :: has_deriv = .false.
      !! Whether this field has a derivative being used for automatic
      !! differentiation.
    real(r8), dimension(:), allocatable :: deriv_data
      !! The value of the differential of the scalar field at the
      !! data-points. Each element represents the value at a different
      !! location.
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
    procedure, pass(rhs) :: real_array_multiply_field => array_scalar_vr_m_sf
      !! \(\vec{\rm real} \times {\rm field}\)
    procedure :: field_multiply_real => array_scalar_sf_m_r
      !! \({\rm field} \times {\rm real}\)
    procedure :: field_multiply_real_array => array_scalar_sf_m_vr
      !! \({\rm field} \times \vec{\rm real}\)
    procedure :: field_divide_field => array_scalar_sf_d_sf
      !! \(\frac{\rm field}{\rm field}\)
    procedure, pass(rhs) :: real_divide_field => array_scalar_r_d_sf
      !! \(\frac{\rm real}{\rm field}\)
    procedure, pass(rhs) :: real_array_divide_field => array_scalar_vr_d_sf
      !! \(\frac{\vec{\rm real}}{\rm field}\)
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
    procedure, non_overridable, public :: read_hdf_array => array_scalar_read_hdf
      !! Reads the array of field data from an HDF file
    procedure, non_overridable, public :: write_hdf_array => array_scalar_write_hdf
      !! Writes the array of field data to an HDF file
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
    procedure(sf_array_interp), deferred :: array_interpolate
      !! Interpolates a value in the field, using a 1-D array of data
      !! passed to it.
    procedure, public :: get_element => array_scalar_get_element
      !! Returns one of the constituent values of the field, i.e. the 
      !! field's value at a particular location.
    procedure, public :: set_element => array_scalar_set_element
      !! Sets one of the constituent values of the field, i.e. the 
      !! field's value at a particular location.
    procedure, public :: get_boundary => array_scalar_get_bound
      !! Returns a field of the same type, containing only the
      !! specified ammount of data at the specified boundary.
    procedure, public :: set_derivative => array_scalar_set_deriv
      !! Sets a derivative value for this field, which gets propagated
      !! through operations using automatic differentiation (tangent
      !! mode).
    procedure, public :: get_derivative => array_scalar_get_deriv
      !! Provides the derivative value for this field. This was either
      !! set by the user or calculated using automatic
      !! differentiation. If no derivative information is available
      !! for this object then a field of zeros is returned.
    procedure, non_overridable :: allocate_deriv => array_scalar_alloc_deriv
      !! Decides whether to allocate the array containing the
      !! derivative.
    procedure(sf_bound), deferred :: subtype_boundary
      !! Performs whatever operations are needed on the subtype to get
      !! a boundary field, including returning the slices needed to
      !! extract the appropriate data.
    procedure, non_overridable :: is_allocated => array_scalar_is_allocated
      !! Indicates whether both the array pointer _and_ the array it
      !! contains are allocated.
    procedure, public :: finalise => array_scalar_finalise
      !! Deallocates the contents of the field. The bulk of the
      !! field's memory should be deallocated with the `clean_temp`
      !! method, but that is unable to deallocate the pointer to the
      !! data array--only the array itself. In compilers which support
      !! finalisation, this method eliminates the small memory leak
      !! from the pointer.
    procedure, public :: interpolate => array_scalar_interp
      !! Interpolates the value of the field at the specified
      !! location.
  end type array_scalar_field

  interface array_scalar_field
    module procedure array_scalar_constructor1
    module procedure array_scalar_constructor2
  end interface array_scalar_field

  abstract interface
    function sf_raw_slices(this,exclude_lower_bound,exclude_upper_bound) &
                                                             result(slices)
      import :: array_scalar_field
      class(array_scalar_field), intent(in)       :: this
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
    end function sf_raw_slices
     
    subroutine sf_compatible(this,other)
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

    subroutine sf_bound(this,src,boundary,depth,slices)
      import :: array_scalar_field
      class(array_scalar_field), intent(inout)          :: this
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
        !! An array containing array slice data which can be used to
        !! construct extract the data for the given field
        !! boundary. The form of the array is
        !! ```
        !! slices(1,i) = start_index
        !! slices(2,i) = end_index
        !! slices(3,i) = stride
        !! ```
    end subroutine sf_bound

    subroutine sf_meta(this, rhs)
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

    function sf_array_interp(this, data_array, location) result(val)
      import :: array_scalar_field
      import :: r8
      class(array_scalar_field), intent(in) :: this
      real(r8), dimension(:), intent(in)    :: data_array
        !! An array holding the datapoints for this field, identical
        !! in layout to that stored the field itself.
      real(r8), dimension(:), intent(in)    :: location
        !! The location at which to calculate the interpolated value.
      real(r8)                              :: val
        !! The spatial derivative of order `order` taken in direction `dir`
    end function sf_array_interp

    function scalar_init(x) result(scalar)
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
    integer                               :: numpoints
      !! The number of datapoints used
    real(r8), dimension(:,:), allocatable :: field_data
      !! The value of the vector field at the data-points. Each row
      !! represents a different spatial location, while each column
      !! represents a different component of the vector.
    integer                               :: vector_dims = 0
      !! The number of vector components
    logical                               :: has_deriv = .false.
      !! Whether this field has a derivative being used for automatic
      !! differentiation.
    real(r8), dimension(:,:), allocatable :: deriv_data
      !! The value of the differential of the scalar field at the
      !! data-points. Each row represents a different spatial
      !! location, while each column represents a different component
      !! of the vector.
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
    procedure :: field_dot_field => array_vector_vf_dot_vf
      !! \({\rm \vec{field}} \cdot {\rm \vec{field}}\)
    procedure :: field_dot_real => array_vector_vf_dot_vr
      !! \({\rm \vec{field}} \cdot {\rm \vec{real}}\)
    procedure, pass(rhs) :: real_dot_field => array_vector_vr_dot_vf
      !! \({\rm \vec{real}} \cdot {\rm \vec{field}}\)
    procedure :: field_cross_field => array_vector_vf_cross_vf
      !! \({\rm\vec{field}} \times {\rm\vec{field}}\)
    procedure :: field_cross_real => array_vector_vf_cross_vr
      !! \({\rm\vec{field}} \times {\rm\vec{real}}\)
    procedure, pass(rhs) :: real_cross_field => array_vector_vr_cross_vf
      !! \({\rm\vec{real}} \times {\rm\vec{field}}\)
    procedure :: assign_field => array_vector_assign
      !! \({\rm field} = {\rm field}\)
    procedure :: assign_scalar_fields => array_vector_assign_scalar
      !! \({\rm \vec{field}} = [{\rm field1, field2, \ldots}]\)
    procedure :: is_equal => array_vector_is_equal
      !! Checks fields are equal within a tolerance
    procedure, public :: assign_meta_data => array_vector_assign_meta_data
      !! Copies all data other than values stored in field from another
      !! field object to this one.
    procedure, non_overridable, public :: read_hdf_array => array_vector_read_hdf
      !! Reads the array of field data from an HDF file
    procedure, non_overridable, public :: write_hdf_array => array_vector_write_hdf
      !! Writes the array of field data to an HDF file
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
    procedure(vf_array_interp), deferred :: array_interpolate
      !! Interpolates a value in the field, using a 2-D array of data
      !! passed to it.
    procedure :: get_element_vector => array_vector_get_element_vec
      !! Returns ones of the constituent vectors of the field, i.e. the 
      !! field's value at a particular location.
    procedure :: get_element_component => array_vector_get_element_comp
      !! Returns one of the components of a constituent vector of the 
      !! field, i.e. the component of the field's value at a particular 
      !! location.
    procedure :: set_element_vector => array_vector_set_element_vec
      !! Sets ones of the constituent vectors of the field, i.e. the 
      !! field's value at a particular location.
    procedure :: set_element_component => array_vector_set_element_comp
      !! Sets one of the components of a constituent vector of the 
      !! field, i.e. the component of the field's value at a particular 
      !! location.
    procedure, public :: get_boundary => array_vector_get_bound
      !! Returns a field of the same type, containing only the
      !! specified ammount of data at the specified boundary.
    procedure, public :: set_derivative => array_vector_set_deriv
      !! Sets a derivative value for this field, which gets propagated
      !! through operations using automatic differentiation (tangent
      !! mode).
    procedure, public :: get_derivative => array_vector_get_deriv
      !! Provides the derivative value for this field. This was either
      !! set by the user or calculated using automatic
      !! differentiation. If no derivative information is available
      !! for this object then a field of zeros is returned.
    procedure, non_overridable :: allocate_deriv => array_vector_alloc_deriv
      !! Decides whether to allocate the array containing the
      !! derivative.
    procedure(vf_bound), deferred :: subtype_boundary
      !! Performs whatever operations are needed by the subtype to get
      !! a boundary field, including returning the slices needed to
      !! extract the appropriate data.
    procedure, non_overridable :: is_allocated => array_vector_is_allocated
      !! Indicates whether both the array pointer _and_ the array it
      !! contains are allocated.
    procedure, public :: finalise => array_vector_finalise
      !! Deallocates the contents of the field. The bulk of the
      !! field's memory should be deallocated with the `clean_temp`
      !! method, but that is unable to deallocate the pointer to the
      !! data array--only the array itself. In compilers which support
      !! finalisation, this method eliminates the small memory leak
      !! from the pointer.
    procedure, public :: interpolate => array_vector_interp
      !! Interpolates the value of the field at the specified
      !! location.
  end type array_vector_field

  interface array_vector_field
    module procedure array_vector_constructor1
    module procedure array_vector_constructor2
 end interface array_vector_field

  abstract interface
    function vf_raw_slices(this,exclude_lower_bound,exclude_upper_bound) &
                                                             result(slices)
      import array_vector_field
      class(array_vector_field), intent(in)       :: this
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
    end function vf_raw_slices

    subroutine vf_bound(this,src,boundary,depth,slices)
      import :: array_vector_field
      class(array_vector_field), intent(inout)          :: this
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
        !! An array containing array slice data which can be used to
        !! construct extract the data for the given field
        !! boundary. The form of the array is
        !! ```
        !! slices(1,i) = start_index
        !! slices(2,i) = end_index
        !! slices(3,i) = stride
        !! ```
    end subroutine vf_bound

    subroutine vf_compatible(this,other)
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

    subroutine vf_meta(this, rhs)
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
    
    function vf_array_interp(this, data_array, location) result(val)
      import :: array_vector_field
      import :: r8
      class(array_vector_field), intent(in) :: this
      real(r8), dimension(:,:), intent(in)  :: data_array
        !! An array holding the datapoints for this field, identical
        !! in layout to that stored the field itself.
      real(r8), dimension(:), intent(in)    :: location
        !! The location at which to calculate the interpolated value.
      real(r8), dimension(:), allocatable   :: val
        !! The spatial derivative of order `order` taken in direction `dir`
    end function vf_array_interp

    function vector_init(x) result(vector)
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

  interface allocate_like
    !! Helper routines to allocate one array to have a 
    !! similar shape to another.
    module procedure ::  allocate_like_1d_1d
    module procedure ::  allocate_like_1d_2d
    module procedure ::  allocate_like_2d_1d
    module procedure ::  allocate_like_2d_2d
    module procedure ::  allocate_like_2d_2d_1d
  end interface

contains

  subroutine allocate_like_1d_1d(array, mold, alloc)
    !* Author: Chris MacMackin
    !  Date: March 2018
    !
    ! Allocates an array to have the same shape as an existing one.
    !
    real(r8), allocatable, dimension(:), intent(inout) :: array
      !! The array to be allocated
    real(r8), allocatable, dimension(:), intent(in) :: mold
      !! The array whose shape is used to allocate the first argument
    logical, optional, intent(in) :: alloc
      !! If present and false, do not allocate the array
    logical :: al
    integer :: n
    if (present(alloc)) then
      al = alloc
    else
      al = .true.
    end if
    if (al .and. allocated(mold)) then
      n = size(mold)
      if (allocated(array)) then
        if (size(array) /= n) then
          deallocate(array)
          allocate(array(n))
        end if
      else
        allocate(array(n))
      end if
    else if (allocated(array)) then
       deallocate(array)
    end if
  end subroutine allocate_like_1d_1d

  subroutine allocate_like_1d_2d(array, mold, alloc)
    !* Author: Chris MacMackin
    !  Date: March 2018
    !
    ! Allocates an array to have the same shape as an existing one.
    !
    real(r8), allocatable, dimension(:), intent(inout) :: array
      !! The array to be allocated
    real(r8), allocatable, dimension(:,:), intent(in) :: mold
      !! The array whose shape is used to allocate the first argument
    logical, optional, intent(in) :: alloc
      !! If present and false, do not allocate the array
    logical :: al
    integer :: n
    if (present(alloc)) then
      al = alloc
    else
      al = .true.
    end if
    if (al .and. allocated(mold)) then
      n = size(mold, 1)
      if (allocated(array)) then
        if (size(array) /= n) then
          deallocate(array)
          allocate(array(n))
        end if
      else
        allocate(array(n))
      end if
    else if (allocated(array)) then
       deallocate(array)
    end if
  end subroutine allocate_like_1d_2d

  subroutine allocate_like_2d_1d(array, mold, n2, alloc)
    !* Author: Chris MacMackin
    !  Date: March 2018
    !
    ! Allocates an array to have the same shape as an existing one.
    !
    real(r8), allocatable, dimension(:,:), intent(inout) :: array
      !! The array to be allocated
    real(r8), allocatable, dimension(:), intent(in) :: mold
      !! The array whose shape is used to allocate the first argument
    integer, optional, intent(in) :: n2
      !! Desired size of array in the second dimension. Deafualt is 1.
    logical, optional, intent(in) :: alloc
      !! If present and false, do not allocate the array
    logical :: al
    integer :: n0, n1
    if (present(alloc)) then
      al = alloc
    else
      al = .true.
    end if
    if (al .and. allocated(mold)) then
      n0 = size(mold, 1)
      if (present(n2)) then
        n1 = n2
      else
        n1 = 1
      end if
      if (allocated(array)) then
        if (size(array, 1) /= n0 .or. size(array, 2) /= n1) then
          deallocate(array)
          allocate(array(n0, n1))
        end if
      else
        allocate(array(n0, n1))
      end if
    else if (allocated(array)) then
       deallocate(array)
    end if
  end subroutine allocate_like_2d_1d

  subroutine allocate_like_2d_2d(array, mold, mold2, alloc)
    !* Author: Chris MacMackin
    !  Date: March 2018
    !
    ! Allocates an array to have the same shape as an existing one.
    !
    real(r8), allocatable, dimension(:,:), intent(inout) :: array
      !! The array to be allocated
    real(r8), allocatable, dimension(:,:), intent(in) :: mold
      !! The array whose shape is used to allocate the first argument
    real(r8), optional, allocatable, dimension(:,:), intent(in) :: mold2
      !! A second array which, if present, will be compared with that
      !! of `mold` to decide the size to allocate the second
      !! dimenions. This will be the larger of the second dimensions
      !! in `mold1` and `mold2`
    logical, optional, intent(in) :: alloc
      !! If present and false, do not allocate the array
    logical :: al
    integer :: n0, n1
    if (present(alloc)) then
      al = alloc
    else
      al = .true.
    end if
    if (al .and. allocated(mold)) then
      n0 = size(mold, 1)
      n1 = size(mold, 2)
      if (present(mold2)) n1 = max(n1, size(mold2, 2))
      if (allocated(array)) then
        if (size(array, 1) /= n0 .or. size(array, 2) /= n1) then
          deallocate(array)
          allocate(array(n0, n1))
        end if
      else
        allocate(array(n0, n1))
      end if
    else if (allocated(array)) then
       deallocate(array)
    end if
  end subroutine allocate_like_2d_2d

  subroutine allocate_like_2d_2d_1d(array, mold, mold2, alloc)
    !* Author: Chris MacMackin
    !  Date: March 2018
    !
    ! Allocates an array to have the same shape as an existing one.
    !
    real(r8), allocatable, dimension(:,:), intent(inout) :: array
      !! The array to be allocated
    real(r8), allocatable, dimension(:,:), intent(in) :: mold
      !! The array whose shape is used to allocate the first argument
    real(r8), dimension(:), intent(in) :: mold2
      !! A second array which will be compared with that
      !! of `mold` to decide the size to allocate the second
      !! dimenions. This will be the larger of the second dimension
      !! of `mold1` and the first of `mold2`
    logical, optional, intent(in) :: alloc
      !! If present and false, do not allocate the array
    logical :: al
    integer :: n0, n1
    if (present(alloc)) then
      al = alloc
    else
      al = .true.
    end if
    if (al .and. allocated(mold)) then
      n0 = size(mold, 1)
      n1 = max(size(mold, 2), size(mold2))
      if (allocated(array)) then
        if (size(array, 1) /= n0 .or. size(array, 2) /= n1) then
          deallocate(array)
          allocate(array(n0, n1))
        end if
      else
        allocate(array(n0, n1))
      end if
    else if (allocated(array)) then
       deallocate(array)
    end if
  end subroutine allocate_like_2d_2d_1d

  !=====================================================================
  ! Scalar Field Methods
  !=====================================================================

  function array_scalar_constructor1(template,numpoints, &
       initializer) result(this)
    !* Author: Chris MacMackin
    !  Date: October 2016
    !
    ! Creates a new field with the same concrete type as the template
    ! argument. The array of values will be allocated and initiated.
    !
    class(array_scalar_field), intent(in)  :: template
      !! A scalar field object which will act as a mold for the concrete
      !! type of the returned type, also copying over any metadata.
    integer, intent(in)                    :: numpoints
      !! The number of data points needed in the array when modelling this
      !! field.
    procedure(scalar_init), optional       :: initializer
      !! An impure elemental procedure taking which takes the position in the
      !! fields domain (an 8-byte real) as an argument and returns the
      !! fields value at that position. Default is for field to be zero
      !! everywhere.
    class(scalar_field), pointer           :: this
      !! A scalar field initiated based on the arguments to this function.
    integer :: i
    call template%guard_temp()
    call template%allocate_scalar_field(this)
    call this%unset_temp()
    select type(this)
    class is(array_scalar_field)
      this%has_deriv = .false.
      call this%assign_subtype_meta_data(template)
      if (allocated(this%field_data)) then
        if (size(this%field_data) /= numpoints) then
          deallocate(this%field_data)
          allocate(this%field_data(numpoints))
        end if
      else
        allocate(this%field_data(numpoints))
      end if
      this%numpoints = numpoints
      if (present(initializer)) then
        do i = 1, numpoints
          this%field_data(i) = initializer(this%id_to_position(i))
        end do
      else
        this%field_data = 0
      end if
    class default
      error stop ('Non-array_scalar_field type allocated by '//&
                  '`allocate_scalar_field` routine.')
    end select
    call template%clean_temp()
    call this%set_temp()
  end function array_scalar_constructor1

  function array_scalar_constructor2(template, array) result(this)
    !* Author: Chris MacMackin
    !  Date: December 2016
    !
    ! Creates a new field with the same concrete type as the template
    ! argument. The field values will be copied from the provided array.
    !
    class(array_scalar_field), intent(in)  :: template
      !! A scalar field object which will act as a mold for the concrete
      !! type of the returned type, also copying over any metadata.
    real(r8), dimension(:), intent(in)     :: array
      !! An array containing the values which this field will be
      !! initialised with.
    class(scalar_field), pointer           :: this
      !! A scalar field initiated based on the arguments to this function.
    call template%guard_temp()
    call template%allocate_scalar_field(this)
    select type(this)
    class is(array_scalar_field)
      this%has_deriv = .false.
      call this%assign_subtype_meta_data(template)
      if (.not. allocated(this%field_data)) then
        allocate(this%field_data(size(array)))
      else
        if (size(this%field_data) /= size(array)) then
          deallocate(this%field_data)
          allocate(this%field_data(size(array)))
        end if
      end if
      this%numpoints = size(array)
      this%field_data = array
    class default
      error stop ('Non-array_scalar_field type allocated by '//&
                  '`allocate_scalar_field` routine.')
    end select
    call template%clean_temp()
  end function array_scalar_constructor2

  function array_scalar_elements(this) result(elements)
    !* Author: Chris MacMackin
    !  Date: October 2016
    !
    ! Gives the number of individual data points present in the field.
    !
    class(array_scalar_field), intent(in) :: this
    integer :: elements
    call this%guard_temp()
    elements = this%numpoints
    call this%clean_temp()
  end function array_scalar_elements

  function array_scalar_raw_size(this,exclude_lower_bound, &
                                      exclude_upper_bound) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! Compute how many elements are in the raw representation of this
    ! field. This would be the number of data points, adjusted based on
    ! how boundary conditions are accounted for.
    !
    class(array_scalar_field), intent(in)       :: this
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
    integer :: res
    integer, dimension(:,:), allocatable :: slices
    integer :: i
    call this%guard_temp()
    if (this%is_allocated()) then
      slices = this%raw_slices(exclude_lower_bound,exclude_upper_bound)
      res = sum(elements_in_slice(slices(1,:),slices(2,:),slices(3,:)))
    else
      res = 0
    end if
    call this%clean_temp()
  end function array_scalar_raw_size
  
  function array_scalar_raw(this,exclude_lower_bound, &
                                 exclude_upper_bound) result(res)
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
    real(r8), dimension(:), allocatable :: res
      !! Array containing data needed to describe field
    integer, dimension(:,:), allocatable :: slices
    integer :: i
    call this%guard_temp()
    slices = this%raw_slices(exclude_lower_bound,exclude_upper_bound)
    if (this%is_allocated()) then
      res = [(this%field_data(slices(1,i):slices(2,i):slices(3,i)), i=1, &
             size(slices,2))]
    else
      allocate(res(0))
    end if
    call this%clean_temp()
  end function array_scalar_raw
  
  subroutine array_scalar_set_from_raw(this,raw,provide_lower_bound, &
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
    integer, dimension(:), optional, intent(in) :: provide_lower_bound
      !! Specifies how many layers of data points are excluded
      !! from the raw data at the lower boundary for each
      !! dimension. The number in element `n` of the array indicates
      !! how many layers of cells at the lower boundary normal to
      !! dimension `n` are missed. Defaults to 0 for all.
    integer, dimension(:), optional, intent(in) :: provide_upper_bound
      !! Specifies how many layers of data points are excluded
      !! from the raw data at the upper boundary for each
      !! dimension. The number in element `n` of the array indicates
      !! how many layers of cells at the upper boundary normal to
      !! dimension `n` are missed. Defaults to 0 for all.
    integer, dimension(:,:), allocatable :: slices
    integer, dimension(:), allocatable :: counts
    integer :: i, start, finish
    call this%guard_temp()
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
    call this%clean_temp()
  end subroutine array_scalar_set_from_raw
  
  function array_scalar_sf_m_sf(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm field} \times {\rm field}\)
    !
    class(array_scalar_field), intent(in) :: this
    class(scalar_field), intent(in) :: rhs
    class(scalar_field), pointer :: res !! The restult of this operation
    call this%guard_temp(); call rhs%guard_temp()
#:if defined('DEBUG')
    call this%check_compatible(rhs)
#:endif
    call this%allocate_scalar_field(res)
    call res%assign_meta_data(this)
    select type (res)
    class is(array_scalar_field)
      call res%allocate_deriv(this, rhs)
      select type(rhs)
      class is(array_scalar_field)
        res%field_data = this%field_data * rhs%field_data
        if (this%has_deriv .and. rhs%has_deriv) then
          res%deriv_data = this%field_data*rhs%deriv_data + this%deriv_data*rhs%field_data
        else if (this%has_deriv) then
          res%deriv_data = this%deriv_data*rhs%field_data
        else if (rhs%has_deriv) then
          res%deriv_data = this%field_data*rhs%deriv_data
        end if
      class is(uniform_scalar_field)
        res%field_data = this%field_data * rhs%get_value()
        if (this%has_deriv) then
          res%deriv_data = rhs%get_value()*this%deriv_data
        end if
      end select
    class default
      error stop ('Non-array_scalar_field type allocated by '//&
                  '`allocate_scalar_field` routine.')
    end select
    call this%clean_temp(); call rhs%clean_temp()
  end function array_scalar_sf_m_sf

  function array_scalar_sf_m_vf(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm field} \times {\rm \vec{field}}\)
    !
    class(array_scalar_field), intent(in) :: this
    class(vector_field), intent(in) :: rhs
    class(vector_field), pointer :: res !! The result of this operation
    integer :: i
    call this%guard_temp(); call rhs%guard_temp()
#:if defined('DEBUG')
    call this%check_compatible(rhs)
#:endif
    call this%allocate_vector_field(res)
    select type(res)
    class is(array_vector_field)
      select type(rhs)
      class is(array_vector_field)
        call res%assign_meta_data(rhs)
        call res%allocate_deriv(this, rhs)
        if (this%has_deriv .and. rhs%has_deriv) then
          do concurrent (i=1:this%numpoints)
            res%field_data(i,:) = this%field_data(i) * rhs%field_data(i,:)
            res%deriv_data(i,:) = this%field_data(i) * rhs%deriv_data(i,:) + &
                                  this%deriv_data(i) * rhs%field_data(i,:)
          end do
        else if (this%has_deriv) then
          do concurrent (i=1:this%numpoints)
            res%field_data(i,:) = this%field_data(i) * rhs%field_data(i,:)
            res%deriv_data(i,:) = this%deriv_data(i) * rhs%field_data(i,:)
         end do
        else if (rhs%has_deriv) then
          do concurrent (i=1:this%numpoints)
            res%field_data(i,:) = this%field_data(i) * rhs%field_data(i,:)
            res%deriv_data(i,:) = this%field_data(i) * rhs%deriv_data(i,:)
          end do
        else
          do concurrent (i=1:this%numpoints)
            res%field_data(i,:) = this%field_data(i) * rhs%field_data(i,:)
          end do
        end if
      class is(uniform_vector_field)
        call res%assign_meta_data(rhs, .false.)
        allocate(res%field_data(this%numpoints,rhs%vector_dimensions()))
        call res%allocate_deriv(this, rhs)
        res%vector_dims = size(rhs%get_value())
        call res%allocate_deriv(this, rhs)
        if (this%has_deriv) then
          do i = 1, this%numpoints
            res%field_data(i,:) = this%field_data(i) * rhs%get_value()
            res%deriv_data(i,:) = this%deriv_data(i) * rhs%get_value()
          end do
        else
          do i = 1, this%numpoints
            res%field_data(i,:) = this%field_data(i) * rhs%get_value()
          end do
        end if
      end select
    class default
      error stop ('Non-array_vector_field type allocated by '//&
                  '`allocate_vector_field` routine.')
    end select
    call this%clean_temp(); call rhs%clean_temp()
  end function array_scalar_sf_m_vf

  function array_scalar_r_m_sf(lhs,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm real} \times {\rm field}\)
    !
    real(r8), intent(in) :: lhs
    class(array_scalar_field), intent(in) :: rhs
    class(scalar_field), pointer :: res !! The result of this operation
    call rhs%guard_temp()
    call rhs%allocate_scalar_field(res)
    select type(res)
    class is(array_scalar_field)
      call res%assign_meta_data(rhs)
      res%field_data = lhs * rhs%field_data
      if (res%has_deriv) res%deriv_data = lhs*rhs%deriv_data
    class default
      error stop ('Non-array_scalar_field type allocated by '//&
                  '`allocate_scalar_field` routine.')
    end select
    call rhs%clean_temp()
  end function array_scalar_r_m_sf

  function array_scalar_vr_m_sf(lhs,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \(\vec{\rm real} \times {\rm field}\)
    !
    real(r8), dimension(:), intent(in) :: lhs
    class(array_scalar_field), intent(in) :: rhs
    class(vector_field), pointer :: res !! The result of this operation
    integer :: i
    call rhs%guard_temp()
    call rhs%allocate_vector_field(res)
    select type(res)
    class is(array_vector_field)
      call res%assign_meta_data(rhs, .false.)
      call allocate_like(res%field_data, rhs%field_data, size(lhs))
      call res%allocate_deriv(rhs)
      res%vector_dims = size(lhs)
      if (res%has_deriv) then
        do concurrent (i=1:rhs%numpoints)
          res%field_data(i,:) = rhs%field_data(i) * lhs
          res%deriv_data(i,:) = rhs%deriv_data(i) * lhs
        end do
      else
        do concurrent (i=1:rhs%numpoints)
          res%field_data(i,:) = rhs%field_data(i) * lhs
        end do
      end if
    class default
      error stop ('Non-array_vector_field type allocated by '//&
                  '`allocate_vector_field` routine.')
    end select
    call rhs%clean_temp()
  end function array_scalar_vr_m_sf

  function array_scalar_sf_m_r(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm field} \times {\rm real}\)
    !
    class(array_scalar_field), intent(in) :: this
    real(r8), intent(in) :: rhs
    class(scalar_field), pointer :: res !! The result of this operation
    call this%guard_temp()
    call this%allocate_scalar_field(res)
    select type(res)
    class is(array_scalar_field)
      call res%assign_meta_data(this)
      call res%allocate_deriv(this)
      res%field_data = this%field_data * rhs
      if (res%has_deriv) res%deriv_data = this%deriv_data * rhs
    class default
      error stop ('Non-array_scalar_field type allocated by '//&
                  '`allocate_scalar_field` routine.')
    end select
    call this%clean_temp()
  end function array_scalar_sf_m_r

  function array_scalar_sf_m_vr(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \({\rm field} \times \vec{\rm real}\)
    !
    class(array_scalar_field), intent(in) :: this
    real(r8), dimension(:), intent(in) :: rhs
    class(vector_field), pointer :: res !! The result of this operation
    integer :: i
    call this%guard_temp()
    call this%allocate_vector_field(res)
    select type(res)
    class is(array_vector_field)
      call res%assign_meta_data(this, .false.)
      allocate(res%field_data(this%numpoints,size(rhs)))
      call res%allocate_deriv(this)
      res%vector_dims = size(rhs)
      if (res%has_deriv) then
        do concurrent (i=1:this%numpoints)
          res%field_data(i,:) = this%field_data(i) * rhs
          res%deriv_data(i,:) = this%deriv_data(i) * rhs
        end do
      else
        do concurrent (i=1:this%numpoints)
          res%field_data(i,:) = this%field_data(i) * rhs
        end do
      end if
    class default
      error stop ('Non-array_vector_field type allocated by '//&
                  '`allocate_vector_field` routine.')
    end select
    call this%clean_temp()
  end function array_scalar_sf_m_vr
  
  function array_scalar_sf_d_sf(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm field} / {\rm field}\)
    !
    class(array_scalar_field), intent(in) :: this
    class(scalar_field), intent(in) :: rhs
    class(scalar_field), pointer :: res !! The restult of this operation
    call this%guard_temp(); call rhs%guard_temp()
#:if defined('DEBUG')
    call this%check_compatible(rhs)
#:endif
    call this%allocate_scalar_field(res)
    call res%assign_meta_data(this)
    select type(res)
    class is(array_scalar_field)
      call res%allocate_deriv(this, rhs)
      select type(rhs)
      class is(array_scalar_field)
        res%field_data = this%field_data / rhs%field_data
        if (this%has_deriv .and. rhs%has_deriv) then
          res%field_data = (this%deriv_data*rhs%field_data - &
                            this%field_data*rhs%deriv_data)/rhs%field_data**2
        else if (this%has_deriv) then
          res%field_data = this%deriv_data/rhs%field_data
        else if (rhs%has_deriv) then
          res%field_data = -this%field_data*rhs%deriv_data/rhs%field_data**2
        end if
      class is(uniform_scalar_field)
        res%field_data = this%field_data / rhs%get_value()
        if (res%has_deriv) then
          res%deriv_data = this%deriv_data / rhs%get_value()
        end if
      end select
    class default
      error stop ('Non-array_scalar_field type allocated by '//&
                  '`allocate_scalar_field` routine.')
    end select
    call this%clean_temp(); call rhs%clean_temp()
  end function array_scalar_sf_d_sf

  function array_scalar_r_d_sf(lhs,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm real} / {\rm field}\)
    !
    real(r8), intent(in) :: lhs
    class(array_scalar_field), intent(in) :: rhs
    class(scalar_field), pointer :: res !! The result of this operation
    call rhs%guard_temp()
    call rhs%allocate_scalar_field(res)
    select type(res)
    class is(array_scalar_field)
      call res%assign_meta_data(rhs)
      call res%allocate_deriv(rhs)
      res%field_data = lhs / rhs%field_data
      if (res%has_deriv) then
        res%field_data = -lhs*rhs%deriv_data/rhs%field_data**2
      end if
    class default
      error stop ('Non-array_scalar_field type allocated by '//&
                  '`allocate_scalar_field` routine.')
    end select      
    call rhs%clean_temp()
  end function array_scalar_r_d_sf

  function array_scalar_vr_d_sf(lhs,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \(\vec{\rm real} / {\rm field}\)
    !
    real(r8), dimension(:), intent(in) :: lhs
    class(array_scalar_field), intent(in) :: rhs
    class(vector_field), pointer :: res !! The result of this operation
    integer :: i
    call rhs%guard_temp()
    call rhs%allocate_vector_field(res)
    select type(res)
    class is(array_vector_field)
      call res%assign_meta_data(rhs, .false.)
      res%vector_dims = size(lhs)
      call allocate_like(res%field_data, rhs%field_data, res%vector_dims)
      call res%allocate_deriv(rhs)
      if (res%has_deriv) then
        do concurrent (i=1:rhs%numpoints)
          res%field_data(i,:) = lhs / rhs%field_data(i)
          res%field_data(i,:) = -lhs*rhs%deriv_data(i)/rhs%field_data(i)**2
        end do
      else
        do concurrent (i=1:rhs%numpoints)
          res%field_data(i,:) = lhs / rhs%field_data(i)
        end do
      end if
    class default
      error stop ('Non-array_scalar_field type allocated by '//&
                  '`allocate_scalar_field` routine.')
    end select
    call rhs%clean_temp()
  end function array_scalar_vr_d_sf

  function array_scalar_sf_d_r(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm field} / {\rm real}\)
    !
    class(array_scalar_field), intent(in) :: this
    real(r8), intent(in) :: rhs
    class(scalar_field), pointer :: res !! The result of this operation
    call this%guard_temp()
    call this%allocate_scalar_field(res)
    call res%assign_meta_data(this)
    select type(res)
    class is(array_scalar_field)
      call res%allocate_deriv(this)
      res%field_data = this%field_data / rhs
      if (res%has_deriv) then
        res%deriv_data = this%deriv_data / rhs
      end if
    class default
      error stop ('Non-array_scalar_field type allocated by '//&
                  '`allocate_scalar_field` routine.')
    end select
    call this%clean_temp()
  end function array_scalar_sf_d_r
  
  function array_scalar_sf_s_sf(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm field} - {\rm field}\)
    !
    class(array_scalar_field), intent(in) :: this
    class(scalar_field), intent(in) :: rhs
    class(scalar_field), pointer :: res !! The restult of this operation
    call this%guard_temp(); call rhs%guard_temp()
#:if defined('DEBUG')
    call this%check_compatible(rhs)
#:endif
    call this%allocate_scalar_field(res)
    call res%assign_meta_data(this)
    select type(res)
    class is(array_scalar_field)
      call res%allocate_deriv(this, rhs)
      select type(rhs)
      class is(array_scalar_field)
        res%field_data = this%field_data - rhs%field_data
        if (this%has_deriv .and. rhs%has_deriv) then
          res%deriv_data = this%deriv_data - rhs%deriv_data
        else if (this%has_deriv) then
          res%deriv_data = this%deriv_data
        else if (rhs%has_deriv) then
          res%deriv_data = -rhs%deriv_data
        end if
      class is(uniform_scalar_field)
        res%field_data = this%field_data - rhs%get_value()
        if (res%has_deriv) then
          res%deriv_data = this%deriv_data
        end if
      end select
    class default
      error stop ('Non-array_scalar_field type allocated by '//&
                  '`allocate_scalar_field` routine.')
    end select
    call this%clean_temp(); call rhs%clean_temp()
  end function array_scalar_sf_s_sf

  function array_scalar_r_s_sf(lhs,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm real} - {\rm field}\)
    !
    real(r8), intent(in) :: lhs
    class(array_scalar_field), intent(in) :: rhs
    class(scalar_field), pointer :: res !! The result of this operation
    call rhs%guard_temp()
    call rhs%allocate_scalar_field(res)
    call res%assign_meta_data(rhs)
    select type(res)
    class is(array_scalar_field)
      call res%allocate_deriv(rhs)
      res%field_data = lhs - rhs%field_data
      if (res%has_deriv) then
        res%deriv_data = -rhs%deriv_data
      end if
    class default
      error stop ('Non-array_scalar_field type allocated by '//&
                  '`allocate_scalar_field` routine.')
    end select
    call rhs%clean_temp()
  end function array_scalar_r_s_sf

  function array_scalar_sf_s_r(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm field} - {\rm real}\)
    !
    class(array_scalar_field), intent(in) :: this
    real(r8), intent(in) :: rhs
    class(scalar_field), pointer :: res !! The result of this operation
    call this%guard_temp()
    call this%allocate_scalar_field(res)
    call res%assign_meta_data(this)
    select type(res)
    class is(array_scalar_field)
      call res%allocate_deriv(this)
      res%field_data = this%field_data - rhs
      if (res%has_deriv) then
        res%deriv_data = this%deriv_data
      end if
    class default
      error stop ('Non-array_scalar_field type allocated by '//&
                  '`allocate_scalar_field` routine.')
    end select    
    call this%clean_temp()
  end function array_scalar_sf_s_r
  
  function array_scalar_sf_a_sf(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm field} + {\rm field}\)
    !
    class(array_scalar_field), intent(in) :: this
    class(scalar_field), intent(in) :: rhs
    class(scalar_field), pointer :: res !! The restult of this operation
    call this%guard_temp(); call rhs%guard_temp()
#:if defined('DEBUG')
    call this%check_compatible(rhs)
#:endif
    call this%allocate_scalar_field(res)
    call res%assign_meta_data(this)
    select type(res)
    class is(array_scalar_field)
      call res%allocate_deriv(this, rhs)
      select type(rhs)
      class is(array_scalar_field)
        res%field_data = this%field_data + rhs%field_data
        if (this%has_deriv .and. rhs%has_deriv) then
          res%deriv_data = this%deriv_data + rhs%deriv_data
        else if (this%has_deriv) then
          res%deriv_data = this%deriv_data
        else if (rhs%has_deriv) then
          res%deriv_data = rhs%deriv_data
        end if
      class is(uniform_scalar_field)
        res%field_data = this%field_data + rhs%get_value()
        if (res%has_deriv) then
          res%deriv_data = this%deriv_data
        end if
      end select
    class default
      error stop ('Non-array_scalar_field type allocated by '//&
                  '`allocate_scalar_field` routine.')
    end select    
    call this%clean_temp(); call rhs%clean_temp()
  end function array_scalar_sf_a_sf

  function array_scalar_r_a_sf(lhs,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm real} + {\rm field}\)
    !
    real(r8), intent(in) :: lhs
    class(array_scalar_field), intent(in) :: rhs
    class(scalar_field), pointer :: res !! The result of this operation
    call rhs%guard_temp()
    call rhs%allocate_scalar_field(res)
    call res%assign_meta_data(rhs)
    select type(res)
    class is(array_scalar_field)
      call res%allocate_deriv(rhs)
      res%field_data = lhs + rhs%field_data
      if (res%has_deriv) then
        res%deriv_data = rhs%deriv_data
      end if
    class default
      error stop ('Non-array_scalar_field type allocated by '//&
                  '`allocate_scalar_field` routine.')
    end select    
    call rhs%clean_temp()
  end function array_scalar_r_a_sf

  function array_scalar_sf_a_r(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm field} + {\rm real}\)
    !
    class(array_scalar_field), intent(in) :: this
    real(r8), intent(in) :: rhs
    class(scalar_field), pointer :: res !! The result of this operation
    call this%guard_temp()
    call this%allocate_scalar_field(res)
    call res%assign_meta_data(this)
    select type(res)
    class is(array_scalar_field)
      call res%allocate_deriv(this)
      res%field_data = this%field_data + rhs
      if (res%has_deriv) then
        res%deriv_data = this%deriv_data
      end if
    class default
      error stop ('Non-array_scalar_field type allocated by '//&
                  '`allocate_scalar_field` routine.')
    end select    
    call this%clean_temp()
  end function array_scalar_sf_a_r

  function array_scalar_sf_p_r(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm field}^{\rm real}\)
    !
    class(array_scalar_field), intent(in) :: this
    real(r8), intent(in) :: rhs
    class(scalar_field), pointer :: res !! The result of this operation
    call this%guard_temp()
    call this%allocate_scalar_field(res)
    call res%assign_meta_data(this)
    select type(res)
    class is(array_scalar_field)
      call res%allocate_deriv(this)
      res%field_data = this%field_data ** rhs
      if (res%has_deriv) then
        res%deriv_data = rhs*this%field_data ** (rhs-1._r8) * this%deriv_data
      end if
    class default
      error stop ('Non-array_scalar_field type allocated by '//&
                  '`allocate_scalar_field` routine.')
    end select    
    call this%clean_temp()
  end function array_scalar_sf_p_r

  function array_scalar_sf_p_r4(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm field}^{\rm real}\)
    !
    class(array_scalar_field), intent(in) :: this
    real, intent(in) :: rhs
    class(scalar_field), pointer :: res !! The result of this operation
    call this%guard_temp()
    call this%allocate_scalar_field(res)
    call res%assign_meta_data(this)
    select type(res)
    class is(array_scalar_field)
      call res%allocate_deriv(this)
      res%field_data = this%field_data ** rhs
      if (res%has_deriv) then
        res%deriv_data = rhs*this%field_data ** (rhs-1.) * this%deriv_data
      end if
    class default
      error stop ('Non-array_scalar_field type allocated by '//&
                  '`allocate_scalar_field` routine.')
    end select    
    call this%clean_temp()
  end function array_scalar_sf_p_r4

  function array_scalar_sf_p_i(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm field}^{\rm integer}\)
    !
    class(array_scalar_field), intent(in) :: this
    integer, intent(in) :: rhs
    class(scalar_field), pointer :: res !! The result of this operation
    call this%guard_temp()
    call this%allocate_scalar_field(res)
    call res%assign_meta_data(this)
    select type(res)
    class is(array_scalar_field)
      call res%allocate_deriv(this)
      res%field_data = this%field_data ** rhs
      if (res%has_deriv) then
        res%deriv_data = rhs*this%field_data ** (rhs-1) * this%deriv_data
      end if
    class default
      error stop ('Non-array_scalar_field type allocated by '//&
                  '`allocate_scalar_field` routine.')
    end select    
    call this%clean_temp()
  end function array_scalar_sf_p_i


  function array_scalar_sin(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\sin({\rm field})\)
    !
    class(array_scalar_field), intent(in) :: this
    class(scalar_field), pointer :: res !! The result of this operation
    call this%guard_temp()
    call this%allocate_scalar_field(res)
    call res%assign_meta_data(this)
    select type(res)
    class is(array_scalar_field)
      call res%allocate_deriv(this)
      if (allocated(this%field_data)) then
        res%field_data = sin(this%field_data)
      end if
      if (res%has_deriv) then
        res%deriv_data = this%deriv_data*cos(this%field_data)
      end if
    class default
      error stop ('Non-array_scalar_field type allocated by '//&
                  '`allocate_scalar_field` routine.')
    end select    
    call this%clean_temp()
  end function array_scalar_sin

  function array_scalar_cos(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\cos({\rm field})\)
    !
    class(array_scalar_field), intent(in) :: this
    class(scalar_field), pointer :: res !! The result of this operation
    call this%guard_temp()
    call this%allocate_scalar_field(res)
    call res%assign_meta_data(this)
    select type(res)
    class is(array_scalar_field)
      call res%allocate_deriv(this)
      if (allocated(this%field_data)) then
        res%field_data = cos(this%field_data)
      end if
      if (res%has_deriv) then
        res%deriv_data = -this%deriv_data*sin(this%field_data)
      end if
    class default
      error stop ('Non-array_scalar_field type allocated by '//&
                  '`allocate_scalar_field` routine.')
    end select    
    call this%clean_temp()
  end function array_scalar_cos

  function array_scalar_tan(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\tan({\rm field})\)
    !
    class(array_scalar_field), intent(in) :: this
    class(scalar_field), pointer :: res !! The result of this operation
    call this%guard_temp()
    call this%allocate_scalar_field(res)
    call res%assign_meta_data(this)
    select type(res)
    class is(array_scalar_field)
      call res%allocate_deriv(this)
      if (allocated(this%field_data)) then
        res%field_data = tan(this%field_data)
      end if
      if (res%has_deriv) then
        res%deriv_data = this%deriv_data/cos(this%field_data)**2
      end if
    class default
      error stop ('Non-array_scalar_field type allocated by '//&
                  '`allocate_scalar_field` routine.')
    end select    
    call this%clean_temp()
  end function array_scalar_tan

  function array_scalar_asin(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\sin^{-1}({\rm field})\)
    !
    class(array_scalar_field), intent(in) :: this
    class(scalar_field), pointer :: res !! The result of this operation
    call this%guard_temp()
    call this%allocate_scalar_field(res)
    call res%assign_meta_data(this)
    select type(res)
    class is(array_scalar_field)
      call res%allocate_deriv(this)
      if (allocated(this%field_data)) then
        res%field_data = asin(this%field_data)
      end if
      if (res%has_deriv) then
        res%deriv_data = this%deriv_data/sqrt(1._r8 - this%field_data**2)
      end if
    class default
      error stop ('Non-array_scalar_field type allocated by '//&
                  '`allocate_scalar_field` routine.')
    end select    
    call this%clean_temp()
  end function array_scalar_asin

  function array_scalar_acos(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\cos^{-1}({\rm field})\)
    !
    class(array_scalar_field), intent(in) :: this
    class(scalar_field), pointer :: res !! The result of this operation
    call this%guard_temp()
    call this%allocate_scalar_field(res)
    call res%assign_meta_data(this)
    select type(res)
    class is(array_scalar_field)
      call res%allocate_deriv(this)
      if (allocated(this%field_data)) then
        res%field_data = acos(this%field_data)
      end if
      if (res%has_deriv) then
        res%deriv_data = -this%deriv_data/sqrt(1._r8 - this%field_data**2)
      end if
    class default
      error stop ('Non-array_scalar_field type allocated by '//&
                  '`allocate_scalar_field` routine.')
    end select    
    call this%clean_temp()
  end function array_scalar_acos

  function array_scalar_atan(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\tan^{-1}({\rm field})\)
    !
    class(array_scalar_field), intent(in) :: this
    class(scalar_field), pointer :: res !! The result of this operation
    call this%guard_temp()
    call this%allocate_scalar_field(res)
    call res%assign_meta_data(this)
    select type(res)
    class is(array_scalar_field)
      call res%allocate_deriv(this)
      if (allocated(this%field_data)) then
        res%field_data = atan(this%field_data)
      end if
      if (res%has_deriv) then
        res%deriv_data = this%deriv_data/(1._r8 + this%field_data**2)
      end if
    class default
      error stop ('Non-array_scalar_field type allocated by '//&
                  '`allocate_scalar_field` routine.')
    end select    
    call this%clean_temp()
  end function array_scalar_atan

  function array_scalar_sinh(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\sinh({\rm field})\)
    !
    class(array_scalar_field), intent(in) :: this
    class(scalar_field), pointer :: res !! The result of this operation
    call this%guard_temp()
    call this%allocate_scalar_field(res)
    call res%assign_meta_data(this)
    select type(res)
    class is(array_scalar_field)
      call res%allocate_deriv(this)
      if (allocated(this%field_data)) then
        res%field_data = sinh(this%field_data)
      end if
      if (res%has_deriv) then
        res%deriv_data = this%deriv_data*cosh(this%field_data)
      end if
    class default
      error stop ('Non-array_scalar_field type allocated by '//&
                  '`allocate_scalar_field` routine.')
    end select    
    call this%clean_temp()
  end function array_scalar_sinh

  function array_scalar_cosh(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\cosh({\rm field})\)
    !
    class(array_scalar_field), intent(in) :: this
    class(scalar_field), pointer :: res !! The result of this operation
    call this%guard_temp()
    call this%allocate_scalar_field(res)
    call res%assign_meta_data(this)
    select type(res)
    class is(array_scalar_field)
      call res%allocate_deriv(this)
      if (allocated(this%field_data)) then
        res%field_data = cosh(this%field_data)
      end if
      if (res%has_deriv) then
        res%deriv_data = this%deriv_data*sinh(this%field_data)
      end if
    class default
      error stop ('Non-array_scalar_field type allocated by '//&
                  '`allocate_scalar_field` routine.')
    end select    
    call this%clean_temp()
  end function array_scalar_cosh

  function array_scalar_tanh(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\tanh({\rm field})\)
    !
    class(array_scalar_field), intent(in) :: this
    class(scalar_field), pointer :: res !! The result of this operation
    call this%guard_temp()
    call this%allocate_scalar_field(res)
    call res%assign_meta_data(this)
    select type(res)
    class is(array_scalar_field)
      call res%allocate_deriv(this)
      if (allocated(this%field_data)) then
        res%field_data = tanh(this%field_data)
      end if
      if (res%has_deriv) then
        res%deriv_data = this%deriv_data/cosh(this%field_data)**2
      end if
    class default
      error stop ('Non-array_scalar_field type allocated by '//&
                  '`allocate_scalar_field` routine.')
    end select    
    call this%clean_temp()
  end function array_scalar_tanh

  function array_scalar_asinh(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\sinh^{-1}({\rm field})\)
    !
    class(array_scalar_field), intent(in) :: this
    class(scalar_field), pointer :: res !! The result of this operation
    call this%guard_temp()
    call this%allocate_scalar_field(res)
    call res%assign_meta_data(this)
    select type(res)
    class is(array_scalar_field)
      call res%allocate_deriv(this)
      if (allocated(this%field_data)) then
        res%field_data = asinh(this%field_data)
      end if
      if (res%has_deriv) then
        res%deriv_data = this%deriv_data/sqrt(this%deriv_data**2 + 1._r8)
      end if
    class default
      error stop ('Non-array_scalar_field type allocated by '//&
                  '`allocate_scalar_field` routine.')
    end select    
    call this%clean_temp()
  end function array_scalar_asinh

  function array_scalar_acosh(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\cosh^{-1}({\rm field})\)
    !
    class(array_scalar_field), intent(in) :: this
    class(scalar_field), pointer :: res !! The result of this operation
    call this%guard_temp()
    call this%allocate_scalar_field(res)
    call res%assign_meta_data(this)
    select type(res)
    class is(array_scalar_field)
      call res%allocate_deriv(this)
      if (allocated(this%field_data)) then
        res%field_data = acosh(this%field_data)
      end if
      if (res%has_deriv) then
        res%deriv_data = this%deriv_data/sqrt(this%deriv_data**2 - 1._r8)
      end if
    class default
      error stop ('Non-array_scalar_field type allocated by '//&
                  '`allocate_scalar_field` routine.')
    end select    
    call this%clean_temp()
  end function array_scalar_acosh

  function array_scalar_atanh(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\tanh^{-1}({\rm field})\)
    !
    class(array_scalar_field), intent(in) :: this
    class(scalar_field), pointer :: res !! The result of this operation
    call this%guard_temp()
    call this%allocate_scalar_field(res)
    call res%assign_meta_data(this)
    select type(res)
    class is(array_scalar_field)
      call res%allocate_deriv(this)
      if (allocated(this%field_data)) then
        res%field_data = atanh(this%field_data)
      end if
      if (res%has_deriv) then
        res%deriv_data = this%deriv_data/(1-this%field_data**2)
      end if
    class default
      error stop ('Non-array_scalar_field type allocated by '//&
                  '`allocate_scalar_field` routine.')
    end select    
    call this%clean_temp()
  end function array_scalar_atanh

  function array_scalar_log(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\ln({\rm field})\)
    !
    class(array_scalar_field), intent(in) :: this
    class(scalar_field), pointer :: res !! The result of this operation
    call this%guard_temp()
    call this%allocate_scalar_field(res)
    call res%assign_meta_data(this)
    select type(res)
    class is(array_scalar_field)
      call res%allocate_deriv(this)
      if (allocated(this%field_data)) then
        res%field_data = log(this%field_data)
      end if
      if (res%has_deriv) then
        res%deriv_data = this%deriv_data/this%field_data
      end if
    class default
      error stop ('Non-array_scalar_field type allocated by '//&
                  '`allocate_scalar_field` routine.')
    end select    
    call this%clean_temp()
  end function array_scalar_log

  function array_scalar_log10(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\log({\rm field})\)
    !
    class(array_scalar_field), intent(in) :: this
    class(scalar_field), pointer :: res !! The result of this operation
    call this%guard_temp()
    call this%allocate_scalar_field(res)
    call res%assign_meta_data(this)
    select type(res)
    class is(array_scalar_field)
      call res%allocate_deriv(this)
      if (allocated(this%field_data)) then
        res%field_data = log10(this%field_data)
      end if
      if (res%has_deriv) then
        res%deriv_data = this%deriv_data/(this%field_data*log(10._r8))
      end if
    class default
      error stop ('Non-array_scalar_field type allocated by '//&
                  '`allocate_scalar_field` routine.')
    end select    
    call this%clean_temp()
  end function array_scalar_log10

  function array_scalar_exp(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(e^({\rm field})\)
    !
    class(array_scalar_field), intent(in) :: this
    class(scalar_field), pointer :: res !! The result of this operation
    call this%guard_temp()
    call this%allocate_scalar_field(res)
    call res%assign_meta_data(this)
    select type(res)
    class is(array_scalar_field)
      call res%allocate_deriv(this)
      if (allocated(this%field_data)) then
        res%field_data = exp(this%field_data)
      end if
      if (res%has_deriv) then
        res%deriv_data = this%deriv_data*res%field_data
      end if
    class default
      error stop ('Non-array_scalar_field type allocated by '//&
                  '`allocate_scalar_field` routine.')
    end select    
    call this%clean_temp()
  end function array_scalar_exp

  function array_scalar_abs(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\abs({\rm field})\)
    !
    class(array_scalar_field), intent(in) :: this
    class(scalar_field), pointer :: res !! The result of this operation
    call this%guard_temp()
    call this%allocate_scalar_field(res)
    call res%assign_meta_data(this)
    select type(res)
    class is(array_scalar_field)
      call res%allocate_deriv(this)
      if (allocated(this%field_data)) then
        res%field_data = abs(this%field_data)
      end if
      if (res%has_deriv) then
        res%deriv_data = this%deriv_data*sign(1._r8,this%field_data)
      end if
    class default
      error stop ('Non-array_scalar_field type allocated by '//&
                  '`allocate_scalar_field` routine.')
    end select    
    call this%clean_temp()
  end function array_scalar_abs

  function array_scalar_sqrt(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\sqrt({\rm field})\)
    !
    class(array_scalar_field), intent(in) :: this
    class(scalar_field), pointer :: res !! The result of this operation
    call this%guard_temp()
    call this%allocate_scalar_field(res)
    call res%assign_meta_data(this)
    select type(res)
    class is(array_scalar_field)
      call res%allocate_deriv(this)
      if (allocated(this%field_data)) then
        res%field_data = sqrt(this%field_data)
      end if
      if (res%has_deriv) then
        res%deriv_data = this%deriv_data*0.5_r8/res%field_data
      end if
    class default
      error stop ('Non-array_scalar_field type allocated by '//&
                  '`allocate_scalar_field` routine.')
    end select    
    call this%clean_temp()
  end function array_scalar_sqrt

  function array_scalar_minval(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\min({\rm field})\)
    !
    class(array_scalar_field), intent(in) :: this
    real(r8) :: res !! The result of this operation
    call this%guard_temp()
    if (allocated(this%field_data)) then
      res = minval(this%field_data)
    else
      res = 0.0_r8
    end if
    call this%clean_temp()
  end function array_scalar_minval

  function array_scalar_maxval(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\max({\rm field})\)
    !
    class(array_scalar_field), intent(in) :: this
    real(r8) :: res !! The result of this operation
    call this%guard_temp()
    if (allocated(this%field_data)) then
      res = maxval(this%field_data)
    else
      res = 0.0_r8
    end if
    call this%clean_temp()
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
    class(scalar_field), pointer :: res
    call this%guard_temp()
    call this%allocate_scalar_field(res)
    call res%assign_meta_data(this)
    select type(res)
    class is(array_scalar_field)
      call res%allocate_deriv(this)
      res%field_data = this%array_dx(this%field_data, dir, order)
      if (res%has_deriv) then
        res%deriv_data = this%array_dx(this%deriv_data, dir, order)
      end if
    class default
      error stop ('Non-array_scalar_field type allocated by '//&
                  '`allocate_scalar_field` routine.')
    end select    
    call this%clean_temp()
  end function array_scalar_d_dx

  function array_scalar_laplacian(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\nabla^2 {\rm field}\)
    !
    class(array_scalar_field), intent(in) :: this
    class(scalar_field), pointer :: res !! The result of this operation
    integer :: i
    call this%guard_temp()
    call this%allocate_scalar_field(res)
    select type(res)
    class is(array_scalar_field)
      call res%assign_meta_data(this)
      call res%allocate_deriv(this)
      res%field_data = this%array_dx(this%field_data,1,2)
      if (res%has_deriv) then
        res%deriv_data = this%array_dx(this%deriv_data,1,2)
      end if
      do i = 2, this%dimensions()
        res%field_data = res%field_data + this%array_dx(this%field_data,i,2)
        if (res%has_deriv) then
          res%deriv_data = res%deriv_data + this%array_dx(this%deriv_data,i,2)
        end if
      end do
    class default
      error stop ('Non-array_scalar_field type allocated by '//&
                 '`allocate_scalar_field` routine.')
    end select
    call this%clean_temp()
  end function array_scalar_laplacian
  
  function array_scalar_gradient(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\nabla{\rm field}\)
    !
    class(array_scalar_field), intent(in) :: this
    class(vector_field), pointer :: res !! The result of this operation
    integer :: i
    call this%guard_temp()
    call this%allocate_vector_field(res)
    select type(res)
    class is(array_vector_field)
      call res%assign_meta_data(this, .false.)
      call allocate_like(res%field_data, this%field_data, this%dimensions())
      call res%allocate_deriv(this)
      do i = 1, this%dimensions()
      call res%allocate_deriv(this)
        res%field_data(:,i) = this%array_dx(this%field_data,i,1)
        if (res%has_deriv) then
          res%deriv_data(:,i) = this%array_dx(this%deriv_data,i,1)
        end if
      end do
    class default
      error stop ('Non-array_vector_field type allocated by '//&
                 '`allocate_vector_field` routine.')
    end select
    call this%clean_temp()
  end function array_scalar_gradient
  
  impure elemental subroutine array_scalar_assign(this,rhs)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm field} = {\rm field}\)
    !
    class(array_scalar_field), intent(inout) :: this
    class(scalar_field), intent(in) :: rhs
    call rhs%guard_temp()
    select type(rhs)
    class is(array_scalar_field)
      call this%assign_meta_data(rhs, .not. rhs%memory_reusable())
      if (allocated(rhs%field_data)) then
        if (rhs%memory_reusable()) then
          if (allocated(this%field_data)) deallocate(this%field_data)
          call move_alloc(rhs%field_data, this%field_data)
          if (rhs%has_deriv) then
            if (allocated(this%deriv_data)) deallocate(this%deriv_data)
            call move_alloc(rhs%deriv_data, this%deriv_data)
          end if
        else
          this%field_data = rhs%field_data
          if (this%has_deriv) then 
            call this%allocate_deriv(rhs)
            this%deriv_data = rhs%deriv_data
          end if
        end if
      end if
    class is(uniform_scalar_field)
      if (allocated(this%field_data)) then
        this%field_data = rhs%get_value()
      else
        error stop ('Trying to assign `uniform_scalar_field` to '// &
                   '`array_scalar_field` with unallocated contents.')
      end if
    end select
    call rhs%clean_temp()
  end subroutine array_scalar_assign

  logical function array_scalar_is_equal(this,rhs) result(iseq)
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
    call this%guard_temp(); call rhs%guard_temp()
    iseq = .true.
    select type(rhs)
    class is(array_scalar_field)
      if (this%numpoints /= rhs%numpoints) then
        iseq = .false.
        call this%clean_temp(); call rhs%clean_temp()
        return
      end if
      do i=1,this%numpoints
        normalization = abs(this%field_data(i))
        if (normalization < get_tol()) normalization = 1.0_r8
        iseq = iseq .and.( ((this%field_data(i)-rhs%field_data(i)) &
                            /normalization < get_tol()) .or. &
                           (is_nan(this%field_data(i)).and. &
                            is_nan(rhs%field_data(i))) )
        if (.not. iseq) then
          call this%clean_temp(); call rhs%clean_temp()
          return
        end if
      end do
    class is(uniform_scalar_field)
      do i=1,this%numpoints
        normalization = abs(this%field_data(i))
        if (normalization < get_tol()) normalization = 1.0_r8
        iseq = iseq .and.( ((this%field_data(i)-rhs%get_value())/normalization < &
                             get_tol()) .or. (is_nan(this%field_data(i)).and. &
                                              is_nan(rhs%get_value())) )
        if (.not. iseq) then
          call this%clean_temp(); call rhs%clean_temp()
          return
        end if
      end do
    class default
      iseq = .false.
    end select
    call this%clean_temp(); call rhs%clean_temp()
  end function array_scalar_is_equal

  subroutine array_scalar_assign_meta_data(this, rhs, alloc)
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
    call rhs%guard_temp()
    call this%assign_subtype_meta_data(rhs)
    select type(rhs)
    class is(array_scalar_field)
      this%numpoints = rhs%numpoints
      call allocate_like(this%field_data, rhs%field_data, alloc)
      this%has_deriv = rhs%has_deriv
    class is(array_vector_field)
      this%numpoints = rhs%numpoints
      call allocate_like(this%field_data, rhs%field_data, alloc)
      this%has_deriv = rhs%has_deriv
    class default
      this%has_deriv = .false.
    end select
    call rhs%clean_temp()
  end subroutine array_scalar_assign_meta_data

  subroutine array_scalar_read_hdf(this, hdf_id, dataset_name, dims, &
                                   error)
    !* Author: Chris MacMackin
    !  Date: April 2017
    !
    ! Reads the contents of the field from a dataset in an HDF file.
    !
    class(array_scalar_field), intent(inout)   :: this
    integer(hid_t), intent(in)                 :: hdf_id
      !! The identifier for the HDF file/group from which the field
      !! data is to be read
    character(len=*), intent(in)               :: dataset_name
      !! The name of the dataset to be read from the HDF file
    integer(hsize_t), dimension(:), intent(in) :: dims
      !! The dimensions of the data set to be read in
    integer, intent(out)                       :: error
      !! An error code which, upon succesful completion of the
      !! routine, is 0. Otherwise, contains the error code returned
      !! by the HDF library.
    call this%guard_temp()
    error = 0
    this%numpoints = int(product(dims))
    if (.not. allocated(this%field_data)) then
      allocate(this%field_data(this%numpoints))
    else if (size(this%field_data) /= this%numpoints) then
      deallocate(this%field_data)
      allocate(this%field_data(this%numpoints))
    end if
    call h5ltread_dataset_double_f(hdf_id, dataset_name,        &
                                   this%field_data, dims, &
                                   error)
    call this%clean_temp()
  end subroutine array_scalar_read_hdf

  subroutine array_scalar_write_hdf(this, hdf_id, dataset_name, dims, &
                                    error)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Writes the contents of the field to a dataset in an HDF file.
    !
    class(array_scalar_field), intent(in)      :: this
    integer(hid_t), intent(in)                 :: hdf_id
      !! The identifier for the HDF file/group in which the field
      !! data is to be written.
    character(len=*), intent(in)               :: dataset_name
      !! The name of the dataset to be created in the HDF file
      !! containing this field's data.
    integer(hsize_t), dimension(:), intent(in) :: dims
      !! The number of elements in each dimension of the fields
    integer, intent(out)                       :: error
      !! An error code which, upon succesful completion of the
      !! routine, is 0. Otherwise, contains the error code returned
      !! by the HDF library.
    call this%guard_temp()
    error = 0
#:if defined('DEBUG')
    if (product(dims) /= this%numpoints) &
      error stop ('Dimensions in requested output array do not match '// &
                  'number of data-points in the field.')
#:endif
    call h5ltmake_dataset_double_f(hdf_id, dataset_name, size(dims), &
                                   dims, this%field_data, error)
    call this%clean_temp()
  end subroutine array_scalar_write_hdf

  subroutine array_scalar_compatible(this,other)
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
    call this%guard_temp(); call other%guard_temp()
    call this%check_subtype_compatible(other)
    select type(other)
    class is(array_scalar_field)
      if (this%numpoints /= other%numpoints) &
           error stop (err_message//'    different resolutions.')
      if (.not.(allocated(this%field_data).and.allocated(other%field_data))) &
           error stop (err_message//'    uninitialised fields.')
    class is(array_vector_field)
      if (this%numpoints /= other%numpoints) &
           error stop (err_message//'    different resolutions.')
      if (.not.(allocated(this%field_data).and.allocated(other%field_data))) &
           error stop (err_message//'    uninitialised fields.')
    class is(uniform_scalar_field)
      if (.not.(allocated(this%field_data))) &
           error stop (err_message//'    uninitialised fields.')
    class is(uniform_vector_field)
      if (.not.(allocated(this%field_data))) &
           error stop (err_message//'    uninitialised fields.')
    class default
      error stop (err_message//'    incompatible types.')
    end select
    call this%clean_temp(); call other%clean_temp()
  end subroutine array_scalar_compatible

  function array_scalar_get_element(this,element) result(val)
    !* Author: Chris MacMackin
    !  Date: January 2017
    !
    ! Returns an element of the field corresponding to the provided ID 
    ! number.
    !
    class(array_scalar_field), intent(in) :: this
    integer, intent(in) :: element
      !! The ID number of the field element to be returned
    real(r8) :: val
      !! The value of the field corresponding to the specified ID
    call this%guard_temp()
    val = this%field_data(element)
    call this%clean_temp()
  end function array_scalar_get_element

  subroutine array_scalar_set_element(this,element,val)
    !* Author: Chris MacMackin
    !  Date: January 2017
    !
    ! Sets the element of the field corresponding to the given ID to
    ! the given value.
    !
    class(array_scalar_field), intent(inout) :: this
    integer, intent(in) :: element
      !! The ID number of the field element to be set
    real(r8), intent(in) :: val
      !! The new value the field element is to be set to
    call this%guard_temp()
    this%field_data(element) = val
    call this%clean_temp()
  end subroutine array_scalar_set_element

  function array_scalar_get_bound(this,boundary,depth) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Returns a field containing the specified boundary
    ! information.
    !
    class(array_scalar_field), intent(in) :: this
    integer, intent(in) :: boundary
      !! Specifies which boundary is to be returned. The boundary will
      !! be the one normal to dimension of number `abs(boundary)`. If
      !! the argument is negative, then the lower boundary is
      !! returned. If positive, then the upper boundary is
      !! returned. If 0 then the whole field is returned.
    integer, intent(in) :: depth
      !! The number of layers of data-points to return at the
      !! specified boundary.
    class(scalar_field), pointer :: res
      !! A field, of the same type as `this` and with the same
      !! resolution, number of dimensions etc., but containing only
      !! the points within the specified number of layers of cells
      !! adjecent to the specified boundary.
    integer, dimension(:,:), allocatable :: slices
    integer :: i
    call this%guard_temp()
    if (boundary == 0) then
      allocate(res, mold=this)
      res = this
        call this%clean_temp()
      return
    end if
    call this%allocate_scalar_field(res)
    select type(res)
    class is(array_scalar_field)
      call res%subtype_boundary(this,boundary,depth,slices)
      res%numpoints = sum(elements_in_slice(slices(1,:), slices(2,:), slices(3,:)))
      if (.not. allocated(res%field_data)) then
        allocate(res%field_data(res%numpoints))
      else if (size(res%field_data) /= res%numpoints) then
        deallocate(res%field_data)
        allocate(res%field_data(res%numpoints))
      end if
      res%field_data = &
           [(this%field_data(slices(1,i):slices(2,i):slices(3,i)), &
            i=1,size(slices,2))]
      res%has_deriv = this%has_deriv
      if (res%has_deriv) then
        if (.not. allocated(res%deriv_data)) then
          allocate(res%deriv_data(res%numpoints))
        else if (size(res%deriv_data) /= res%numpoints) then
          deallocate(res%deriv_data)
          allocate(res%deriv_data(res%numpoints))
        end if
        res%deriv_data = &
             [(this%deriv_data(slices(1,i):slices(2,i):slices(3,i)), &
              i=1,size(slices,2))]
      end if
    class default
      error stop ('Non-array_scalar_field type allocated by '//&
                 '`allocate_scalar_field` routine.')

    end select
    call this%clean_temp()
  end function array_scalar_get_bound

  subroutine array_scalar_set_deriv(this, deriv)
    !* Author: Chris MacMackin
    !  Date: March 2018
    !
    ! Sets a differential/derivative value for the field, which will
    ! be propagated through operations using automatic differentiation
    ! (tangent mode).
    !
    class(array_scalar_field), intent(inout) :: this
    class(scalar_field), intent(in)          :: deriv
      !! The derivative (i.e. differential) value of this field
    call this%guard_temp(); call deriv%guard_temp()
#:if defined('DEBUG')
    call this%check_compatible(deriv)
#:endif
    select type(deriv)
    class is(array_scalar_field)
      call allocate_like(this%deriv_data, deriv%field_data)
      this%deriv_data = deriv%field_data
    class is(uniform_scalar_field)
      if (allocated(this%field_data)) then
        call allocate_like(this%deriv_data, this%field_data)
        this%deriv_data = deriv%get_value()
      end if
    class default
      error stop ('Non-array_scalar_field and non-uniform_scalar_field '//&
                  'type passed as derivative value.')
    end select
    this%has_deriv = .true.
    call this%clean_temp(); call deriv%clean_temp()
  end subroutine array_scalar_set_deriv

  function array_scalar_get_deriv(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2018
    !
    ! Gets the differential/derivative value for this field. If no
    ! such information is available, returns a field with zero values.
    !
    class(array_scalar_field), intent(in) :: this
    class(scalar_field), pointer    :: res
    call this%guard_temp()
    call this%allocate_scalar_field(res)
    call res%assign_meta_data(this)
    select type (res)
    class is(array_scalar_field)
      if (allocated(this%deriv_data)) then
        res%field_data = this%deriv_data
      else
        res%field_data = 0._r8 * this%field_data
      end if   
      res%has_deriv = .false.
    class default
      error stop ('Non-array_scalar_field type allocated by '//&
                  '`allocate_scalar_field` routine.')
    end select
    call this%clean_temp()
  end function array_scalar_get_deriv
 
  subroutine array_scalar_alloc_deriv(this, field1, field2)
    !* Author: Chris MacMackin
    !  Date: March 2018
    !
    ! Helper function which decides whether a derivative is going to
    ! be calculated and allocates the array if necessary.
    !
    class(array_scalar_field), intent(inout)    :: this
    class(abstract_field), intent(in)           :: field1
      !! The first field which is checked to see if it has a
      !! derivative.
    class(abstract_field), optional, intent(in) :: field2
      !! The optional second field which is checked to see if it has a
      !! derivative.
    logical :: alloc
    select type(field1)
    class is(array_scalar_field)
      alloc = field1%has_deriv
      if (present(field2) .and. .not. alloc) then
        select type(field2)
        class is(array_scalar_field)
          alloc = field2%has_deriv
        class is(array_vector_field)
          alloc = field2%has_deriv
        end select
      end if
    class is(array_vector_field)
      alloc = field1%has_deriv
      if (present(field2) .and. .not. alloc) then
        select type(field2)
        class is(array_scalar_field)
          alloc = field2%has_deriv
        class is(array_vector_field)
          alloc = field2%has_deriv
        end select
      end if
    class default
      alloc = .false.
    end select
    if (alloc) then
      call allocate_like(this%deriv_data, this%field_data)
      this%has_deriv = allocated(this%deriv_data)
    else
      this%has_deriv = .false.
    end if
  end subroutine array_scalar_alloc_deriv

  elemental subroutine array_scalar_finalise(this)
    !* Author: Chris MacMackin
    !  Date: January 2017
    !
    ! Deallocates the field pointer for this object and removes the
    ! 'temporary' setting.
    !
    class(array_scalar_field), intent(inout) :: this
    if (allocated(this%field_data)) deallocate(this%field_data)
    call this%unset_temp()
  end subroutine array_scalar_finalise

  function array_scalar_is_allocated(this) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Returns a logical indicating whether both the array pointer
    ! _and_ the array which it holds are allocated.
    !
    class(array_scalar_field), intent(in) :: this
    logical :: res
    call this%guard_temp()
    res = allocated(this%field_data)
    call this%clean_temp()
  end function array_scalar_is_allocated

  function array_scalar_interp(this, location) result(val)
    !* Author: Chris MacMackin
    !  Date: July 2017
    !
    ! Interpolates the value of the field at the specified location.
    !
    class(array_scalar_field), intent(in) :: this
    real(r8), dimension(:), intent(in)    :: location
      !! The location at which to calculate the interpolated value.
    real(r8)                              :: val
    call this%guard_temp()
    val = this%array_interpolate(this%field_data, location)
    call this%clean_temp()
  end function array_scalar_interp


  !=====================================================================
  ! Vector Field Methods
  !=====================================================================

  function array_vector_constructor1(template,numpoints,vector_dims, &
                                    initializer) result(this)
    !* Author: Chris MacMackin
    !  Date: October 2016
    !
    ! Creates a new field with the same concrete type as the template
    ! argument. The array of values will be allocated and initiated.
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
      !! An impure elemental procedure taking which takes the position in the
      !! fields domain (an 8-byte real) as an argument and returns the
      !! fields value at that position. Default is for field to be zero
      !! everywhere.
    class(vector_field), pointer           :: this
      !! A scalar field initiated based on the arguments to this function.
    integer :: i
    call template%guard_temp()
    call template%allocate_vector_field(this)
    call this%unset_temp()
    select type(this)
    class is(array_vector_field)
      this%has_deriv = .false.
      call this%assign_subtype_meta_data(template)
      if (allocated(this%field_data)) deallocate(this%field_data)
      allocate(this%field_data(numpoints,vector_dims))
      this%numpoints = numpoints
      this%vector_dims = vector_dims
      if (present(initializer)) then
        do i = 1, numpoints
          this%field_data(i,:) = initializer(this%id_to_position(i))
        end do
      else
        this%field_data = 0.0_r8
      end if
    class default
      error stop ('Non-array_vector_field type allocated by '//&
                  '`allocate_vector_field` routine.')
    end select
    call template%clean_temp()
    call this%set_temp()
  end function array_vector_constructor1

  function array_vector_constructor2(template,array) result(this)
    !* Author: Chris MacMackin
    !  Date: December 2016
    !
    ! Creates a new field with the same concrete type as the template
    ! argument. The array of values will be allocated and initiated.
    !
    class(array_vector_field), intent(in)  :: template
      !! A scalar field object which will act as a mold for the concrete
      !! type of the returned type and which provides the meta-data.
    real(r8), dimension(:,:), intent(in)   :: array
      !! An array containing the values which this field will be
      !! initialised with.
    class(vector_field), pointer           :: this
      !! A scalar field initiated based on the arguments to this function.
    integer :: i
    call template%guard_temp()
    call template%allocate_vector_field(this)
    select type(this)
    class is(array_vector_field)
      this%has_deriv = .false.
      call this%assign_subtype_meta_data(template)
      if (.not. allocated(this%field_data)) then
        allocate(this%field_data(size(array,1),size(array,2)))
      else
        if (any(shape(this%field_data) /= shape(array))) then
          deallocate(this%field_data)
          allocate(this%field_data(size(array,1),size(array,2)))
        end if
      end if
      this%numpoints = size(array,1)
      this%vector_dims = size(array,2)
      this%field_data = array
    class default
      error stop ('Non-array_vector_field type allocated by '//&
                  '`allocate_vector_field` routine.')
    end select
    call template%clean_temp()
  end function array_vector_constructor2

  impure elemental function array_vector_vector_dimensions(this) result(dims)
    !* Author: Chris MacMackin
    !  Date: October 2016
    !
    ! Returns the number of dimensions/components are in the vectors
    ! of this field.
    !
    class(array_vector_field), intent(in) :: this
    integer :: dims !! Number of vector components
    call this%guard_temp()
    dims = this%vector_dims
    call this%clean_temp()
  end function array_vector_vector_dimensions

  function array_vector_elements(this) result(elements)
    !* Author: Chris MacMackin
    !  Date: October 2016
    !
    ! Gives the number of individual data points present in the field.
    !
    class(array_vector_field), intent(in) :: this
    integer :: elements
    call this%guard_temp()
    elements = this%numpoints
    call this%clean_temp()
  end function array_vector_elements

  function array_vector_raw_size(this,exclude_lower_bound, &
                                       exclude_upper_bound) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! Compute how many elements are in the raw representation of this
    ! field. This would be the number of data points, adjusted based on
    ! how boundary conditions are accounted for.
    !
    class(array_vector_field), intent(in)       :: this
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
    integer :: res
    integer, dimension(:,:), allocatable :: slices
    integer :: i
    call this%guard_temp()
    if (this%is_allocated()) then
      slices = this%raw_slices(exclude_lower_bound,exclude_upper_bound)
      res = sum(elements_in_slice(slices(1,:),slices(2,:),slices(3,:))) * &
            this%vector_dims
    else
      res = 0
    end if
    call this%clean_temp()
  end function array_vector_raw_size
  
  function array_vector_raw(this,exclude_lower_bound, &
                                  exclude_upper_bound) result(res)
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
    real(r8), dimension(:), allocatable :: res
      !! Array containing data needed to describe field
    integer, dimension(:,:), allocatable :: slices
    integer :: i
    call this%guard_temp()
    slices = this%raw_slices(exclude_lower_bound,exclude_upper_bound)
    if (this%is_allocated()) then
      res = [(this%field_data(slices(1,i):slices(2,i):slices(3,i),:), i=1, &
             size(slices,2))]
    else
      allocate(res(0))
    end if
    call this%clean_temp()
  end function array_vector_raw
  
  subroutine array_vector_set_from_raw(this,raw,provide_lower_bound, &
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
    integer, dimension(:), optional, intent(in) :: provide_lower_bound
      !! Specifies how many layers of data points are excluded
      !! from the raw data at the lower boundary for each
      !! dimension. The number in element `n` of the array indicates
      !! how many layers of cells at the lower boundary normal to
      !! dimension `n` are missed. Defaults to 0 for all.
    integer, dimension(:), optional, intent(in) :: provide_upper_bound
      !! Specifies how many layers of data points are excluded
      !! from the raw data at the upper boundary for each
      !! dimension. The number in element `n` of the array indicates
      !! how many layers of cells at the upper boundary normal to
      !! dimension `n` are missed. Defaults to 0 for all.
    integer, dimension(:,:), allocatable :: slices
    integer, dimension(:), allocatable :: counts
    integer :: i, start, finish
    call this%guard_temp()
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
    call this%clean_temp()
  end subroutine array_vector_set_from_raw

  function array_vector_vf_m_sf(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\vec{\rm field} \times {\rm field}\)
    !
    class(array_vector_field), intent(in) :: this
    class(scalar_field), intent(in) :: rhs
    class(vector_field), pointer :: res !! The restult of this operation
    integer :: i
    call this%guard_temp(); call rhs%guard_temp()
#:if defined('DEBUG')
    call this%check_compatible(rhs)
#:endif
    call this%allocate_vector_field(res)
    call res%assign_meta_data(this)
    select type(res)
    class is(array_vector_field)
      call res%allocate_deriv(this, rhs)
      select type(rhs)
      class is(array_scalar_field)
        if (this%has_deriv .and. rhs%has_deriv) then
          do concurrent(i=1:this%vector_dims)
            res%field_data(:,i) = this%field_data(:,i) * rhs%field_data
            res%deriv_data(:,i) = this%field_data(:,i) * rhs%deriv_data + &
                                  this%deriv_data(:,i) * rhs%field_data
          end do
        else if (this%has_deriv) then
          do concurrent(i=1:this%vector_dims)
            res%field_data(:,i) = this%field_data(:,i) * rhs%field_data
            res%deriv_data(:,i) = this%deriv_data(:,i) * rhs%field_data
          end do
        else if (rhs%has_deriv) then
          do concurrent(i=1:this%vector_dims)
            res%field_data(:,i) = this%field_data(:,i) * rhs%field_data
            res%deriv_data(:,i) = this%field_data(:,i) * rhs%deriv_data
          end do
        else
          do concurrent(i=1:this%vector_dims)
            res%field_data(:,i) = this%field_data(:,i) * rhs%field_data
          end do
        end if
      class is(uniform_scalar_field)
        call res%allocate_deriv(this)
        res%field_data = this%field_data * rhs%get_value()
        if (this%has_deriv) then
          res%deriv_data = this%deriv_data * rhs%get_value()
        end if
      end select
    class default
      error stop ('Non-array_vector_field type allocated by '//&
                  '`allocate_vector_field` routine.')
    end select
    call this%clean_temp(); call rhs%clean_temp()
  end function array_vector_vf_m_sf

  function array_vector_r_m_vf(lhs,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm real} \times \vec{\rm field}\)
    !
    real(r8), intent(in) :: lhs
    class(array_vector_field), intent(in) :: rhs
    class(vector_field), pointer :: res !! The result of this operation
    call rhs%guard_temp()
    call rhs%allocate_vector_field(res)
    call res%assign_meta_data(rhs)
    select type(res)
    class is(array_vector_field)
      call res%allocate_deriv(rhs)
      res%field_data = lhs * rhs%field_data
      if (res%has_deriv) res%deriv_data = lhs*rhs%deriv_data 
    class default
      error stop ('Non-array_vector_field type allocated by '//&
                  '`allocate_vector_field` routine.')
    end select
    call rhs%clean_temp()
  end function array_vector_r_m_vf

  function array_vector_vf_m_r(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\vec{\rm field} \times {\rm real}\)
    !
    class(array_vector_field), intent(in) :: this
    real(r8), intent(in) :: rhs
    class(vector_field), pointer :: res !! The result of this operation
    call this%guard_temp()
    call this%allocate_vector_field(res)
    call res%assign_meta_data(this)
    select type(res)
    class is(array_vector_field)
      call res%allocate_deriv(this)
      res%field_data = this%field_data * rhs
      if (res%has_deriv) res%deriv_data = rhs*this%deriv_data 
    class default
      error stop ('Non-array_vector_field type allocated by '//&
                  '`allocate_vector_field` routine.')
    end select
    call this%clean_temp()
  end function array_vector_vf_m_r
  
  function array_vector_vf_d_sf(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\vec{\rm field} / {\rm field}\)
    !
    class(array_vector_field), intent(in) :: this
    class(scalar_field), intent(in) :: rhs
    class(vector_field), pointer :: res !! The restult of this operation
    integer :: i
    call this%guard_temp(); call rhs%guard_temp()
#:if defined('DEBUG')
    call this%check_compatible(rhs)
#:endif
    call this%allocate_vector_field(res)
    call res%assign_meta_data(this)
    select type(res)
    class is(array_vector_field)
      select type(rhs)
      class is(array_scalar_field)
        call res%allocate_deriv(this, rhs)
        if (this%has_deriv .and. rhs%has_deriv) then
          do concurrent(i=1:this%vector_dims)
            res%field_data(:,i) = this%field_data(:,i) / rhs%field_data
            res%deriv_data(:,i) = (this%deriv_data(:,i)*rhs%field_data - &
                                   this%field_data(:,i)*rhs%deriv_data)/ &
                                  rhs%field_data**2
          end do
        else if (this%has_deriv) then
          do concurrent(i=1:this%vector_dims)
            res%field_data(:,i) = this%field_data(:,i) / rhs%field_data
            res%field_data(:,i) = this%deriv_data(:,i)/rhs%field_data
          end do
        else if (rhs%has_deriv) then
          do concurrent(i=1:this%vector_dims)
            res%field_data(:,i) = this%field_data(:,i) / rhs%field_data
            res%field_data(:,i) = -this%field_data(:,i)*rhs%deriv_data/rhs%field_data**2
          end do
        else
          do concurrent(i=1:this%vector_dims)
            res%field_data(:,i) = this%field_data(:,i) / rhs%field_data
          end do
        end if
      class is(uniform_scalar_field)
        call res%allocate_deriv(this)
        res%field_data = this%field_data / rhs%get_value()
        if (res%has_deriv) res%deriv_data = this%deriv_data / rhs%get_value()
      end select
    class default
      error stop ('Non-array_vector_field type allocated by '//&
                  '`allocate_vector_field` routine.')
    end select
    call this%clean_temp(); call rhs%clean_temp()
  end function array_vector_vf_d_sf

  function array_vector_vf_d_r(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\vec{\rm field} / {\rm real}\)
    !
    class(array_vector_field), intent(in) :: this
    real(r8), intent(in) :: rhs
    class(vector_field), pointer :: res !! The result of this operation
    call this%guard_temp()
    call this%allocate_vector_field(res)
    call res%assign_meta_data(this)
    select type(res)
    class is(array_vector_field)
      call res%allocate_deriv(this)
      res%field_data = this%field_data / rhs
      if (res%has_deriv) then
        res%deriv_data = this%deriv_data / rhs
      end if
    class default
      error stop ('Non-array_vector_field type allocated by '//&
                  '`allocate_vector_field` routine.')
    end select
    call this%clean_temp()
  end function array_vector_vf_d_r
  
  function array_vector_vf_s_vf(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\vec{\rm field} - \vec{\rm field}\)
    !
    class(array_vector_field), intent(in) :: this
    class(vector_field), intent(in) :: rhs
    class(vector_field), pointer :: res !! The restult of this operation
    integer :: i, min_dims, max_dims
    real(r8), dimension(:), allocatable :: vec
    call this%guard_temp(); call rhs%guard_temp()
#:if defined('DEBUG')
    call this%check_compatible(rhs)
#:endif
    call this%allocate_vector_field(res)
    select type(res)
    class is(array_vector_field)
      select type(rhs)
      class is(array_vector_field)
        min_dims = min(this%vector_dims,rhs%vector_dimensions())
        max_dims = max(this%vector_dims,rhs%vector_dimensions())
        call res%assign_meta_data(this, .false.)
        call allocate_like(res%field_data, this%field_data, rhs%field_data)
        call res%allocate_deriv(this, rhs)
        res%vector_dims = max_dims
        res%field_data(:,1:min_dims) = this%field_data(:,1:min_dims) &
                                       - rhs%field_data(:,1:min_dims)
        if (rhs%vector_dims > this%vector_dims) then
          res%field_data(:,min_dims+1:max_dims) = -rhs%field_data(:,min_dims+1:max_dims)
        else
          res%field_data(:,min_dims+1:max_dims) = this%field_data(:,min_dims+1:max_dims)
        end if
        if (this%has_deriv .and. rhs%has_deriv) then
          res%deriv_data(:,1:min_dims) = this%deriv_data(:,1:min_dims) &
                                         - rhs%deriv_data(:,1:min_dims)
          if (rhs%vector_dims > this%vector_dims) then
            res%deriv_data(:,min_dims+1:max_dims) = -rhs%deriv_data(:,min_dims+1:max_dims)
          else
            res%deriv_data(:,min_dims+1:max_dims) = this%deriv_data(:,min_dims+1:max_dims)
          end if
        else if (this%has_deriv) then
          res%deriv_data(:,1:min_dims) = this%deriv_data(:,1:min_dims)
          if (rhs%vector_dims > this%vector_dims) then
            res%deriv_data(:,min_dims+1:max_dims) = 0._r8
          else
            res%deriv_data(:,min_dims+1:max_dims) = this%deriv_data(:,min_dims+1:max_dims)
          end if
        else if (rhs%has_deriv) then
          res%deriv_data(:,1:min_dims) = -rhs%deriv_data(:,1:min_dims)
          if (rhs%vector_dims > this%vector_dims) then
            res%deriv_data(:,min_dims+1:max_dims) = -rhs%deriv_data(:,min_dims+1:max_dims)
          else
            res%deriv_data(:,min_dims+1:max_dims) = 0._r8
          end if
        end if
      class is(uniform_vector_field)
        min_dims = min(this%vector_dims,rhs%vector_dimensions())
        max_dims = max(this%vector_dims,rhs%vector_dimensions())
        vec = rhs%get_value()
        call res%assign_meta_data(this, .false.)
        call allocate_like(res%field_data, this%field_data, vec)
        call res%allocate_deriv(this)
        res%vector_dims = max_dims
        if (res%has_deriv) then
          do concurrent(i=1:min_dims)
            res%field_data(:,i) = this%field_data(:,i) - vec(i)
            res%deriv_data(:,i) = this%deriv_data(:,i) 
          end do
        else
          do concurrent(i=1:min_dims)
            res%field_data(:,i) = this%field_data(:,i) - vec(i)
          end do
        end if
        if (rhs%vector_dimensions() > this%vector_dims) then
          if (res%has_deriv) then
            do concurrent(i=min_dims+1:max_dims)
              res%field_data(:,i) = -vec(i)
              res%deriv_data(:,i) = 0._r8
            end do
          else
            do concurrent(i=min_dims+1:max_dims)
              res%field_data(:,i) = -vec(i)
            end do
          end if
        else
          if (res%has_deriv) then
            do concurrent(i=min_dims+1:max_dims)
              res%field_data(:,i) = this%field_data(:,i)
              res%deriv_data(:,i) = this%deriv_data(:,i) 
            end do
          else
            do concurrent(i=min_dims+1:max_dims)
              res%field_data(:,i) = this%field_data(:,i)
            end do
          end if
        end if
      end select
    class default
      error stop ('Non-array_vector_field type allocated by '//&
                  '`allocate_vector_field` routine.')
    end select
    call this%clean_temp(); call rhs%clean_temp()
  end function array_vector_vf_s_vf

  function array_vector_r_s_vf(lhs,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\vec{\rm real} - \vec{\rm field}\)
    !
    real(r8), dimension(:), intent(in) :: lhs
    class(array_vector_field), intent(in) :: rhs
    class(vector_field), pointer :: res !! The result of this operation
    integer :: min_dims, max_dims, i
    call rhs%guard_temp()
    call rhs%allocate_vector_field(res)
    call res%assign_meta_data(rhs,.false.)
    select type(res)
    class is(array_vector_field)
      call res%allocate_deriv(rhs)
      min_dims = min(rhs%vector_dims, size(lhs))
      max_dims = max(rhs%vector_dims, size(lhs))
      call allocate_like(res%field_data, rhs%field_data, lhs)
      if (res%has_deriv) then
        do concurrent (i=1:res%numpoints)
          res%field_data(i,:min_dims) = lhs(:min_dims) - rhs%field_data(i,:min_dims)
          res%deriv_data(i,:min_dims) = -rhs%deriv_data(i,:min_dims) 
        end do
      else
        do concurrent (i=1:res%numpoints)
          res%field_data(i,:min_dims) = lhs(:min_dims) - rhs%field_data(i,:min_dims)
        end do
      end if
      if (rhs%vector_dims > size(lhs)) then
        res%field_data(:,min_dims+1:) = -rhs%field_data(:,min_dims+1:)
        if (res%has_deriv) res%deriv_data(:,min_dims+1:) = -rhs%deriv_data(:,min_dims+1:)
      else
        res%vector_dims = size(lhs)
        if (res%has_deriv) then
          do concurrent (i=1:res%numpoints)
            res%field_data(i,min_dims+1:) = lhs(min_dims+1:)
            res%deriv_data(i,min_dims+1:) = 0._r8 
          end do
        else
          do concurrent (i=1:res%numpoints)
            res%field_data(i,min_dims+1:) = lhs(min_dims+1:)
          end do
        end if
      end if
    class default
      error stop ('Non-array_vector_field type allocated by '//&
                  '`allocate_vector_field` routine.')
    end select
    call rhs%clean_temp()
  end function array_vector_r_s_vf

  function array_vector_vf_s_r(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\vec{\rm field} - \vec{\rm real}\)
    !
    class(array_vector_field), intent(in) :: this
    real(r8), dimension(:), intent(in) :: rhs
    class(vector_field), pointer :: res !! The result of this operation
    integer :: min_dims, max_dims, i
    call this%guard_temp()
    call this%allocate_vector_field(res)
    call res%assign_meta_data(this,.false.)
    min_dims = min(this%vector_dims, size(rhs))
    max_dims = max(this%vector_dims, size(rhs))
    select type(res)
    class is(array_vector_field)
      call allocate_like(res%field_data, this%field_data, rhs)
      call res%allocate_deriv(this)
      if (res%has_deriv) then
        do concurrent (i=1:res%numpoints)
          res%field_data(i,:min_dims) = this%field_data(i,:min_dims) - rhs(:min_dims)
          res%deriv_data(i,:min_dims) = this%deriv_data(i,:min_dims)
        end do
      else
        do concurrent (i=1:res%numpoints)
          res%field_data(i,:min_dims) = this%field_data(i,:min_dims) - rhs(:min_dims)
        end do
      end if
      if (this%vector_dims > size(rhs)) then
        res%field_data(:,min_dims+1:) = this%field_data(:,min_dims+1:)
        if (res%has_deriv) res%deriv_data(:,min_dims+1:) = this%deriv_data(:,min_dims+1:)
      else
        res%vector_dims = size(rhs)
        if (res%has_deriv) then
          do concurrent (i=1:res%numpoints)
            res%field_data(i,min_dims+1:) = -rhs(min_dims+1:)
            res%deriv_data(i,min_dims+1:) = 0._r8
          end do
        else
          do concurrent (i=1:res%numpoints)
            res%field_data(i,min_dims+1:) = -rhs(min_dims+1:)
          end do
        end if
      end if
    class default
      error stop ('Non-array_vector_field type allocated by '//&
                  '`allocate_vector_field` routine.')
    end select
    call this%clean_temp()
  end function array_vector_vf_s_r
  
  function array_vector_vf_a_vf(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\vec{\rm field} + \vec{\rm field}\)
    !
    class(array_vector_field), intent(in) :: this
    class(vector_field), intent(in) :: rhs
    class(vector_field), pointer :: res !! The restult of this operation
    integer :: i, max_dims, min_dims
    real(r8), dimension(:), allocatable :: vec
    call this%guard_temp(); call rhs%guard_temp()
#:if defined('DEBUG')
    call this%check_compatible(rhs)
#:endif
    call this%allocate_vector_field(res)
    select type(res)
    class is(array_vector_field)
      select type(rhs)
      class is(array_vector_field)
        min_dims = min(this%vector_dims,rhs%vector_dimensions())
        max_dims = max(this%vector_dims,rhs%vector_dimensions())
        call res%assign_meta_data(rhs, .false.)
        call allocate_like(res%field_data, this%field_data, rhs%field_data)
        call res%allocate_deriv(this, rhs)
        res%vector_dims = max_dims
        res%field_data(:,1:min_dims) = this%field_data(:,1:min_dims) &
                                       + rhs%field_data(:,1:min_dims)
        if (rhs%vector_dims > this%vector_dims) then
          res%field_data(:,min_dims+1:max_dims) = rhs%field_data(:,min_dims+1:max_dims)
        else
          res%field_data(:,min_dims+1:max_dims) = this%field_data(:,min_dims+1:max_dims)
        end if
        if (this%has_deriv .and. rhs%has_deriv) then
          res%deriv_data(:,1:min_dims) = this%deriv_data(:,1:min_dims) &
                                         + rhs%deriv_data(:,1:min_dims)
          if (rhs%vector_dims > this%vector_dims) then
            res%deriv_data(:,min_dims+1:max_dims) = rhs%deriv_data(:,min_dims+1:max_dims)
          else
            res%deriv_data(:,min_dims+1:max_dims) = this%deriv_data(:,min_dims+1:max_dims)
          end if
        else if (this%has_deriv) then
          res%deriv_data(:,1:min_dims) = this%deriv_data(:,1:min_dims)
          if (rhs%vector_dims > this%vector_dims) then
            res%deriv_data(:,min_dims+1:max_dims) = 0._r8
          else
            res%deriv_data(:,min_dims+1:max_dims) = this%deriv_data(:,min_dims+1:max_dims)
          end if
        else if (rhs%has_deriv) then
          res%deriv_data(:,1:min_dims) = rhs%deriv_data(:,1:min_dims)
          if (rhs%vector_dims > this%vector_dims) then
            res%deriv_data(:,min_dims+1:max_dims) = rhs%deriv_data(:,min_dims+1:max_dims)
          else
            res%deriv_data(:,min_dims+1:max_dims) = 0._r8
          end if
        end if
      class is(uniform_vector_field)
        min_dims = min(this%vector_dims,rhs%vector_dimensions())
        max_dims = max(this%vector_dims,rhs%vector_dimensions())
        vec = rhs%get_value()
        call res%assign_meta_data(rhs, .false.)
        call allocate_like(res%field_data, this%field_data, vec)
        call res%allocate_deriv(this)
        res%vector_dims = max_dims
        if (res%has_deriv) then
          do concurrent(i=1:min_dims)
            res%field_data(:,i) = this%field_data(:,i) + vec(i)
            res%field_data(:,i) = this%deriv_data(:,i)
          end do
        else
          do concurrent(i=1:min_dims)
            res%field_data(:,i) = this%field_data(:,i) + vec(i)
          end do
        end if
        if (rhs%vector_dimensions() > this%vector_dims) then
          if (res%has_deriv) then
            do concurrent(i=min_dims+1:max_dims)
              res%field_data(:,i) = vec(i)
              res%field_data(:,i) = 0._r8
            end do
          else
            do concurrent(i=min_dims+1:max_dims)
              res%field_data(:,i) = vec(i)
            end do
          end if
        else
          if (res%has_deriv) then
            do concurrent(i=min_dims+1:max_dims)
              res%field_data(:,i) = this%field_data(:,i)
              res%deriv_data(:,i) = this%deriv_data(:,i)
            end do
          else
            do concurrent(i=min_dims+1:max_dims)
              res%field_data(:,i) = this%field_data(:,i)
            end do
          end if
        end if
      end select
    class default
      error stop ('Non-array_vector_field type allocated by '//&
                  '`allocate_vector_field` routine.')
    end select
    call this%clean_temp(); call rhs%clean_temp()
  end function array_vector_vf_a_vf

  function array_vector_r_a_vf(lhs,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\vec{\rm real} + \vec{\rm field}\)
    !
    real(r8), dimension(:), intent(in) :: lhs
    class(array_vector_field), intent(in) :: rhs
    class(vector_field), pointer :: res !! The result of this operation
    integer :: min_dims, max_dims, i
    call rhs%guard_temp()
    call rhs%allocate_vector_field(res)
    call res%assign_meta_data(rhs,.false.)
    min_dims = min(rhs%vector_dims, size(lhs))
    max_dims = max(rhs%vector_dims, size(lhs))
    select type(res)
    class is(array_vector_field)
      call allocate_like(res%field_data, rhs%field_data, lhs)
      call res%allocate_deriv(rhs)
      if (res%has_deriv) then
        do concurrent (i=1:res%numpoints)
          res%field_data(i,:min_dims) = lhs(:min_dims) + rhs%field_data(i,:min_dims)
          res%deriv_data(i,:min_dims) = rhs%field_data(i,:min_dims)
        end do
      else
        do concurrent (i=1:res%numpoints)
          res%field_data(i,:min_dims) = lhs(:min_dims) + rhs%field_data(i,:min_dims)
        end do
      end if
      if (rhs%vector_dims > size(lhs)) then
        res%field_data(:,min_dims+1:) = rhs%field_data(:,min_dims+1:)
        if (res%has_deriv) res%deriv_data(:,min_dims+1:) = rhs%deriv_data(:,min_dims+1:)
      else
        res%vector_dims = size(lhs)
        if (res%has_deriv) then
          do concurrent (i=1:res%numpoints)
            res%field_data(i,min_dims+1:) = lhs(min_dims+1:)
            res%deriv_data(i,min_dims+1:) = 0._r8
          end do
        else
          do concurrent (i=1:res%numpoints)
            res%field_data(i,min_dims+1:) = lhs(min_dims+1:)
          end do
        end if
      end if
    class default
      error stop ('Non-array_vector_field type allocated by '//&
                  '`allocate_vector_field` routine.')
    end select
    call rhs%clean_temp()
  end function array_vector_r_a_vf

  function array_vector_vf_a_r(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\vec{\rm field} + \vec{\rm real}\)
    !
    class(array_vector_field), intent(in) :: this
    real(r8), dimension(:), intent(in) :: rhs
    class(vector_field), pointer :: res !! The result of this operation
    integer :: min_dims, max_dims, i
    call this%guard_temp()
    call this%allocate_vector_field(res)
    call res%assign_meta_data(this,.false.)
    min_dims = min(this%vector_dims, size(rhs))
    max_dims = max(this%vector_dims, size(rhs))
    select type(res)
    class is(array_vector_field)
      call allocate_like(res%field_data, this%field_data, rhs)
      call res%allocate_deriv(this)
      if (res%has_deriv) then
        do concurrent (i=1:res%numpoints)
          res%field_data(i,:min_dims) = this%field_data(i,:min_dims) + rhs(:min_dims)
          res%deriv_data(i,:min_dims) = this%deriv_data(i,:min_dims)
        end do
      else
        do concurrent (i=1:res%numpoints)
          res%field_data(i,:min_dims) = this%field_data(i,:min_dims) + rhs(:min_dims)
        end do
      end if
      if (this%vector_dims > size(rhs)) then
        res%field_data(:,min_dims+1:) = this%field_data(:,min_dims+1:)
        if (res%has_deriv) res%deriv_data(:,min_dims+1:) = this%deriv_data(:,min_dims+1:)
      else
        res%vector_dims = size(rhs)
        if (res%has_deriv) then
          do concurrent (i=1:res%numpoints) 
            res%field_data(i,min_dims+1:) = rhs(min_dims+1:)
            res%deriv_data(i,min_dims+1:) = 0._r8
          end do
        else
          do concurrent (i=1:res%numpoints) 
            res%field_data(i,min_dims+1:) = rhs(min_dims+1:)
          end do
        end if
      end if
    class default
      error stop ('Non-array_vector_field type allocated by '//&
                  '`allocate_vector_field` routine.')
    end select
    call this%clean_temp()
  end function array_vector_vf_a_r

  impure elemental subroutine array_vector_assign(this,rhs)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\vec{\rm field} = \vec{\rm field}\)
    !
    class(array_vector_field), intent(inout) :: this
    class(vector_field), intent(in) :: rhs
    integer :: i
    real(r8), dimension(:), allocatable :: tmp
    call rhs%guard_temp()
    call this%assign_meta_data(rhs, .not. this%memory_reusable())
    select type(rhs)
    class is(array_vector_field)
      if (allocated(rhs%field_data)) then
        if (rhs%memory_reusable()) then
          if (allocated(this%field_data)) deallocate(this%field_data)
          call move_alloc(rhs%field_data, this%field_data)
          if (rhs%has_deriv) then
            if (allocated(this%deriv_data)) deallocate(this%deriv_data)
            call move_alloc(rhs%deriv_data, this%deriv_data)
          end if
        else
          this%field_data = rhs%field_data
          if (this%has_deriv) then 
            call this%allocate_deriv(rhs)
            this%deriv_data = rhs%deriv_data
          end if
        end if
      end if
    class is(uniform_vector_field)
      if (allocated(this%field_data)) then
        tmp = rhs%get_value()
        do concurrent (i=1:this%vector_dims)
          this%field_data(:,i) = tmp(i)
        end do
      else
        error stop ('Trying to assign `uniform_vector_field` to '// &
                   '`array_vector_field` with unallocated contents.')
      end if
    class default
      error stop ('Incompatible field assignment.')
    end select
    call rhs%clean_temp()
  end subroutine array_vector_assign

  subroutine array_vector_assign_scalar(this,rhs)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \(\vec{\rm field} = [{\rm field1, field2, \ldots}]\)
    !
    class(array_vector_field), intent(inout)      :: this
    class(scalar_field), dimension(:), intent(in) :: rhs
    integer :: i
    call rhs%guard_temp()
    call this%assign_meta_data(rhs(1), .false.)
    select type(rhs)
    class is(array_scalar_field)
      this%has_deriv = any(rhs%has_deriv)
      this%vector_dims = size(rhs)
      this%numpoints = rhs(1)%numpoints
      if (allocated(this%field_data)) then
        if (size(this%field_data,1) /= this%numpoints .or. &
            size(this%field_data,2) /= this%vector_dims) then
          deallocate(this%field_data)
        end if
      end if
      if (.not. allocated(this%field_data)) then
        allocate(this%field_data(this%numpoints,this%vector_dims))
      end if
      if (this%has_deriv) call allocate_like(this%deriv_data, this%field_data)
      do concurrent (i=1:this%vector_dims)
        this%field_data(:,i) = rhs(i)%field_data
        if (rhs(i)%has_deriv) then
          this%deriv_data(:,i) = rhs(i)%deriv_data
        else if (this%has_deriv) then
          this%deriv_data(:,i) = 0._r8
        end if
      end do
    class is(uniform_scalar_field)
      if (allocated(this%field_data)) then
        do i = 1, this%vector_dims
          this%field_data(:,i) = rhs(i)%get_value()
        end do
      else
        error stop ('Trying to assign `uniform_scalar_field` to '// &
                   '`array_vector_field` with unallocated contents.')
      end if
    end select
	call rhs%clean_temp()
  end subroutine array_vector_assign_scalar

  function array_vector_norm(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\lVert \vec{\rm field} \rVert\)
    !
    class(array_vector_field), intent(in) :: this
    class(scalar_field), pointer :: res
    call this%guard_temp()
    call this%allocate_scalar_field(res)
    select type(res)
    class is(array_scalar_field)
      call res%assign_meta_data(this)
      call res%allocate_deriv(this)
      res%field_data = sqrt(sum(this%field_data**2,2))
      if (res%has_deriv) then
        res%deriv_data =sum(this%field_data*this%deriv_data,2)/res%field_data
      end if
    class default
      error stop ('Non-array_scalar_field type allocated by '//&
                 '`allocate_scalar_field` routine.')
    end select
    call this%clean_temp()
  end function array_vector_norm

  function array_vector_component(this,comp) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! Returns a calar field containing specified component of the vector field
    !
    class(array_vector_field), intent(in) :: this
    integer, intent(in) :: comp
    class(scalar_field), pointer :: res
    real(r8), allocatable, dimension(:,:) :: tmp
    call this%guard_temp()
    call this%allocate_scalar_field(res)
    select type(res)
    class is(array_scalar_field)
      call res%assign_meta_data(this)
      call res%allocate_deriv(this)
      if (comp <= this%vector_dims) then
        res%field_data = this%field_data(:,comp)
        if (res%has_deriv) res%deriv_data = this%deriv_data(:,comp)
      else
        res%field_data = 0.0_r8
        if (res%has_deriv) res%deriv_data = 0._r8
      end if
    class default
      error stop ('Non-array_scalar_field type allocated by '//&
                 '`allocate_scalar_field` routine.')
    end select
    call this%clean_temp()
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
    class(vector_field), pointer :: res
    integer :: i
    call this%guard_temp()
    call this%allocate_vector_field(res)
    select type(res)
    class is(array_vector_field)
      call res%assign_meta_data(this)
      call res%allocate_deriv(this)
      if (res%has_deriv) then
        do i = 1, this%vector_dims
          res%field_data(:,i) = this%array_dx(this%field_data(:,i), &
                                                         dir, order)
          res%deriv_data(:,i) = this%array_dx(this%deriv_data(:,i), &
                                                         dir, order)
        end do
      else
        do i = 1, this%vector_dims
          res%field_data(:,i) = this%array_dx(this%field_data(:,i), &
                                                         dir, order)
        end do
      end if
    class default
      error stop ('Non-array_scalar_field type allocated by '//&
                 '`allocate_scalar_field` routine.')
    end select
    call this%clean_temp()
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
    class(scalar_field), pointer :: res
    call this%guard_temp()
    call this%allocate_scalar_field(res)
    select type(res)
    class is(array_scalar_field)
      call res%assign_meta_data(this)
      call res%allocate_deriv(this)
      if (component <= this%vector_dims) then
        res%field_data = this%array_dx(this%field_data(:,component), &
                                       dir, order)
        if (res%has_deriv) res%deriv_data = this%array_dx(this%deriv_data(:,component), &
                                                         dir, order)
      else
        res%field_data = 0._r8
        if (res%has_deriv) res%deriv_data = 0._r8
      end if
    class default
      error stop ('Non-array_scalar_field type allocated by '//&
                 '`allocate_scalar_field` routine.')
    end select
    call this%clean_temp()
  end function array_vector_component_d_dx

  function array_vector_laplacian(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\nabla^2 \vec{\rm field}\)
    !
    class(array_vector_field), intent(in) :: this
    class(vector_field), pointer :: res !! The result of this operation
    integer :: i, j
    call this%guard_temp()
    call this%allocate_vector_field(res)
    select type(res)
    class is(array_vector_field)
      call res%assign_meta_data(this, .false.)
      call allocate_like(res%field_data, this%field_data)
      call res%allocate_deriv(this)
      do i = 1, this%vector_dims
        res%field_data(:,i) = this%array_dx(this%field_data(:,i),1,2)
        if (res%has_deriv) then
          res%deriv_data(:,i) = this%array_dx(this%deriv_data(:,i),1,2) 
          do j = 2, this%dimensions()
            res%field_data(:,i) = res%field_data(:,i) + &
                                  this%array_dx(this%field_data(:,i),j,2)
            res%deriv_data(:,i) = res%deriv_data(:,i) + &
                                  this%array_dx(this%deriv_data(:,i),j,2)
          end do
        else
          do j = 2, this%dimensions()
            res%field_data(:,i) = res%field_data(:,i) + &
                                  this%array_dx(this%field_data(:,i),j,2)
          end do
        end if
      end do
    class default
      error stop ('Non-array_vector_field type allocated by '//&
                 '`allocate_vector_field` routine.')
    end select
    call this%clean_temp()
  end function array_vector_laplacian
  
  function array_vector_divergence(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\nabla\cdot \vec{\rm field}\)
    !
    class(array_vector_field), intent(in) :: this
    class(scalar_field), pointer :: res !! The result of this operation
    integer :: i
    call this%guard_temp()
    call this%allocate_scalar_field(res)
    select type(res)
    class is(array_scalar_field)
      call res%assign_meta_data(this)
      call res%allocate_deriv(this)
      res%field_data = this%array_dx(this%field_data(:,1),1,1)
      if (res%has_deriv) then
        res%deriv_data = this%array_dx(this%deriv_data(:,1),1,1)
        do i = 2, this%dimensions()
          res%field_data = res%field_data + &
                           this%array_dx(this%field_data(:,i),i,1)
          res%deriv_data = res%deriv_data + &
                           this%array_dx(this%deriv_data(:,i),i,1)
        end do
      else
        do i = 2, this%dimensions()
          res%field_data = res%field_data + &
                           this%array_dx(this%field_data(:,i),i,1)
        end do
      end if
    class default
      error stop ('Non-array_scalar_field type allocated by '//&
                 '`allocate_scalar_field` routine.')
    end select
    call this%clean_temp()
  end function array_vector_divergence
  
  function array_vector_curl(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\nabla\times \vec{\rm field}\)
    !
    class(array_vector_field), intent(in) :: this
    class(vector_field), pointer :: res !! The result of this operation
    logical, dimension(3) :: been_set
    integer :: i
    call this%guard_temp()
    been_set = .false.
    call this%allocate_vector_field(res)
    select type(res)
    class is(array_vector_field)
      call res%assign_meta_data(this,.false.)
      call res%allocate_deriv(this)
      res%vector_dims = 3
      if (allocated(res%field_data)) then
        if (size(res%field_data, 1) /= this%numpoints .or. &
            size(res%field_data, 2) /= 3) then
          deallocate(res%field_data)
          allocate(res%field_data(res%numpoints,3))
        end if
      else
        allocate(res%field_data(res%numpoints,3))
      end if
      if (this%dimensions() >= 3) then
        res%field_data(:,2) = this%array_dx(this%field_data(:,1),3)
        if (res%has_deriv) res%deriv_data(:,2) = this%array_dx(this%deriv_data(:,1),3)
        been_set(2) = .true.
      end if
      if (this%dimensions() >= 2) then
        res%field_data(:,3) = -this%array_dx(this%field_data(:,1),2)
        if (res%has_deriv) res%deriv_data(:,3) = this%array_dx(this%deriv_data(:,1),1)
        been_set(3) = .true.
      end if
      if (this%vector_dims >= 2) then
        if (this%dimensions() >= 3) then
          res%field_data(:,1) = -this%array_dx(this%field_data(:,2),3)
          if (res%has_deriv) res%deriv_data(:,1) = -this%array_dx(this%deriv_data(:,2),3)
          been_set(1) = .true.
        end if
        if (been_set(3)) then
          res%field_data(:,3) = res%field_data(:,3) &
                             + this%array_dx(this%field_data(:,2),1)
          if (res%has_deriv) res%deriv_data(:,3) = res%deriv_data(:,3) &
                             + this%array_dx(this%deriv_data(:,2),1)
        else
          res%field_data(:,3) = this%array_dx(this%field_data(:,2),1)
          if (res%has_deriv) res%deriv_data(:,3) = this%array_dx(this%deriv_data(:,2),1)
          been_set(3) = .true.
        end if
      end if
      if (this%vector_dims >= 3) then
        if (this%dimensions() >= 2) then
          if (been_set(1)) then
            res%field_data(:,1) = res%field_data(:,1) &
                               + this%array_dx(this%field_data(:,3),2)
            if (res%has_deriv) res%deriv_data(:,1) = res%deriv_data(:,1) &
                               + this%array_dx(this%deriv_data(:,3),2)
          else
            res%field_data(:,1) = this%array_dx(this%field_data(:,3),2)
            if (res%has_deriv) res%deriv_data(:,1) = this%array_dx(this%deriv_data(:,3),2)
            been_set(1) = .true.
          end if
        end if
        if (been_set(2)) then
          res%field_data(:,2) = res%field_data(:,2) &
                             - this%array_dx(this%field_data(:,3),1)
          if (res%has_deriv) res%deriv_data(:,2) = res%deriv_data(:,2) &
                             + this%array_dx(this%deriv_data(:,3),1)
        else
          res%field_data(:,2) = -this%array_dx(this%field_data(:,3),1)
          if (res%has_deriv) res%deriv_data(:,2) = this%array_dx(this%deriv_data(:,3),1)
          been_set(2) = .true.
        end if
      end if
      do i = 1, 3
        if (.not. been_set(i)) then
          res%field_data(:,i) = 0.0_r8
          if (res%has_deriv) res%deriv_data(:,i) = 0.0_r8
        end if
      end do
    class default
      error stop ('Non-array_scalar_field type allocated by '//&
                  'allocate_scalar_field routine.')
    end select
    call this%clean_temp()
  end function array_vector_curl

  function array_vector_vf_cross_vf(this,rhs) result(res)
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
    class(vector_field), pointer :: res !! The restult of this operation
    real(r8), dimension(3) :: vec1, vec2, dvec1, dvec2
    real(r8), dimension(:), allocatable :: tmp
    integer :: i, dims1, dims2
    call this%guard_temp(); call rhs%guard_temp()
#:if defined('DEBUG')
    call this%check_compatible(rhs)
#:endif
    call this%allocate_vector_field(res)
    call res%assign_meta_data(this,.false.)
    select type(res)
    class is(array_vector_field)
      res%vector_dims = 3
      if (allocated(res%field_data)) then
        if (size(res%field_data, 1) /= this%numpoints .or. &
            size(res%field_data, 2) /= 3) then
          deallocate(res%field_data)
          allocate(res%field_data(this%numpoints,3))
        end if
      else
        allocate(res%field_data(this%numpoints,3))
      end if
      call res%allocate_deriv(this, rhs)
      vec1 = 0
      vec2 = 0
      dvec1 = 0
      dvec2 = 0
      dims1 = min(3,this%vector_dims)
      select type(rhs)
      class is(array_vector_field)
        dims2 = min(3,rhs%vector_dimensions())
        if (res%has_deriv) then
          do concurrent(i=1:this%numpoints)
            vec1(:dims1) = this%field_data(i,:dims1)
            vec2(:dims2) = rhs%field_data(i,:dims2)
            dvec1(:dims1) = this%deriv_data(i,:dims1)
            dvec2(:dims2) = rhs%deriv_data(i,:dims2)
            res%field_data(i,:) = [vec1(2)*vec2(3) - vec2(2)*vec1(3), &
                                   vec1(3)*vec2(1) - vec2(3)*vec1(1), &
                                   vec1(1)*vec2(2) - vec2(1)*vec1(2)]
            res%deriv_data(i,:) = [vec1(2)*dvec2(3) + dvec1(2)*vec2(3) - &
                                   vec2(2)*dvec1(3) - dvec2(2)*vec1(3), &
                                   vec1(3)*dvec2(1) + dvec1(3)*vec2(1) - &
                                   vec2(3)*dvec1(1) - dvec2(3)*vec1(1), &
                                   vec1(1)*dvec2(2) + dvec1(1)*vec2(2) - &
                                   vec2(1)*dvec1(2) - dvec2(2)*vec1(2)]
          end do
        else
          do concurrent(i=1:this%numpoints)
            vec1(:dims1) = this%field_data(i,:dims1)
            vec2(:dims2) = rhs%field_data(i,:dims2)
            res%field_data(i,:) = [vec1(2)*vec2(3) - vec2(2)*vec1(3), &
                                     vec1(3)*vec2(1) - vec2(3)*vec1(1), &
                                     vec1(1)*vec2(2) - vec2(1)*vec1(2)]
          end do
        end if
      class is(uniform_vector_field)
        dims2 = min(3,rhs%vector_dimensions())
        tmp = rhs%get_value()
        vec2(:dims2) = tmp(:dims2)
        if (res%has_deriv) then
          do concurrent(i=1:this%numpoints)
            vec1(:dims1) = this%field_data(i,:dims1)
            dvec1(:dims1) = this%deriv_data(i,:dims1)
            res%field_data(i,:) = [vec1(2)*vec2(3) - vec2(2)*vec1(3), &
                                   vec1(3)*vec2(1) - vec2(3)*vec1(1), &
                                   vec1(1)*vec2(2) - vec2(1)*vec1(2)]
            res%deriv_data(i,:) = [dvec1(2)*vec2(3) - vec2(2)*dvec1(3), &
                                   dvec1(3)*vec2(1) - vec2(3)*dvec1(1), &
                                   dvec1(1)*vec2(2) - vec2(1)*dvec1(2)]
          end do
        else
          do concurrent(i=1:this%numpoints)
            vec1(:dims1) = this%field_data(i,:dims1)
            res%field_data(i,:) = [vec1(2)*vec2(3) - vec2(2)*vec1(3), &
                                   vec1(3)*vec2(1) - vec2(3)*vec1(1), &
                                   vec1(1)*vec2(2) - vec2(1)*vec1(2)]
          end do
        end if
      end select
    class default
      error stop ('Non-array_vector_field type allocated by '//&
                 'allocate_vector_field routine.')
   end select
   call this%clean_temp(); call rhs%clean_temp()
  end function array_vector_vf_cross_vf

  function array_vector_vf_cross_vr(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \(\vec{\rm field} \times \vec{\rm real}\)
    !
    ! The returned vector will always consist of three components,
    ! regardless of the number of components in the arguments.
    !
    class(array_vector_field), intent(in) :: this
    real(r8), dimension(:), intent(in) :: rhs
    class(vector_field), pointer :: res !! The result of this operation
    logical :: res2d
    real(r8), dimension(3) :: vec1, vec2, dvec1
    integer :: i, dims1, dims2
    call this%guard_temp()
    call this%allocate_vector_field(res)
    call res%assign_meta_data(this,.false.)
    select type(res)
    class is(array_vector_field)
      res2d = size(rhs) >= 3
      if (res2d) res2d = rhs(1) == 0._r8 .and. rhs(2) == 0._r8
      if (res2d) then
        res%vector_dims = 2
        if (.not. allocated(res%field_data)) then
          allocate(res%field_data(this%numpoints,2))
        else if (size(res%field_data,1) /= this%numpoints .or. &
                 size(res%field_data,2) /= 2) then
          deallocate(res%field_data)
          allocate(res%field_data(this%numpoints,2))
        end if
        call res%allocate_deriv(this)
        if (this%vector_dims >= 2) then
          res%field_data(:,1) = this%field_data(:,2)*rhs(3)
          if (res%has_deriv) res%deriv_data(:,1) = this%deriv_data(:,2)*rhs(3)
        else
          res%field_data(:,1) = 0._r8
          if (res%has_deriv) res%deriv_data(:,1) = this%deriv_data(:,2)*rhs(3)
        end if
        if (this%vector_dims >= 1) then
          res%field_data(:,2) = -this%field_data(:,1)*rhs(3)
          if (res%has_deriv) res%deriv_data(:,2) = -this%deriv_data(:,1)*rhs(3)
        else
          res%field_data(:,2) = 0._r8
          if (res%has_deriv) res%deriv_data(:,2) = 0._r8
        end if
      else
        res%vector_dims = 3
        if (.not. allocated(res%field_data)) then
          allocate(res%field_data(this%numpoints,3))
        else if (size(res%field_data,1) /= this%numpoints .or. &
                 size(res%field_data,2) /= 3) then
          deallocate(res%field_data)
          allocate(res%field_data(this%numpoints,3))
        end if
        call res%allocate_deriv(this)
        dims1 = min(3,this%vector_dims)
        dims2 = min(3,size(rhs))
        vec1 = 0.0_r8
        vec2 = 0.0_r8
        dvec1 = 0.0_r8
        vec2(:dims2) = rhs(:dims2)
        if (res%has_deriv) then
          do concurrent(i=1:this%numpoints)
            vec1(:dims1) = this%field_data(i,:dims1)
            dvec1(:dims1) = this%deriv_data(i,:dims1)
            res%field_data(i,:) = [vec1(2)*vec2(3) - vec2(2)*vec1(3), &
                                   vec1(3)*vec2(1) - vec2(3)*vec1(1), &
                                   vec1(1)*vec2(2) - vec2(1)*vec1(2)]
            res%deriv_data(i,:) = [dvec1(2)*vec2(3) - vec2(2)*dvec1(3), &
                                   dvec1(3)*vec2(1) - vec2(3)*dvec1(1), &
                                   dvec1(1)*vec2(2) - vec2(1)*dvec1(2)]
          end do
        else
          do concurrent(i=1:this%numpoints)
            vec1(:dims1) = this%field_data(i,:dims1)
            res%field_data(i,:) = [vec1(2)*vec2(3) - vec2(2)*vec1(3), &
                                   vec1(3)*vec2(1) - vec2(3)*vec1(1), &
                                   vec1(1)*vec2(2) - vec2(1)*vec1(2)]
          end do
        end if
      end if
    class default
      error stop ('Non-array_vector_field type allocated by '//&
                 'allocate_vector_field routine.')
    end select
    call this%clean_temp()
  end function array_vector_vf_cross_vr

  function array_vector_vr_cross_vf(lhs,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \(\vec{\rm real} \times \vec{\rm field}\)
    !
    ! The returned vector will always consist of three components,
    ! regardless of the number of components in the arguments.
    !
    real(r8), dimension(:), intent(in) :: lhs
    class(array_vector_field), intent(in) :: rhs
    class(vector_field), pointer :: res !! The result of this operation
    logical :: res2d
    real(r8), dimension(3) :: vec1, vec2, dvec2
    integer :: i, dims1, dims2
    call rhs%guard_temp()
    call rhs%allocate_vector_field(res)
    call res%assign_meta_data(rhs,.false.)
    select type(res)
    class is(array_vector_field)
      res2d = size(lhs) >= 3
      if (res2d) res2d = lhs(1) == 0._r8 .and. lhs(2) == 0._r8
      if (res2d) then
        res%vector_dims = 2
        if (.not. allocated(res%field_data)) then
          allocate(res%field_data(rhs%numpoints,2))
        else if (size(res%field_data,1) /= rhs%numpoints .or. &
                 size(res%field_data,2) /= 2) then
          deallocate(res%field_data)
          allocate(res%field_data(rhs%numpoints,2))
        end if
        call res%allocate_deriv(rhs)
        if (rhs%vector_dims >= 2) then
          res%field_data(:,1) = -lhs(3)*rhs%field_data(:,2)
          if (res%has_deriv) res%deriv_data(:,1) = -lhs(3)*rhs%deriv_data(:,2)
        else
          res%field_data(:,1) = 0._r8
          if (res%has_deriv) res%deriv_data(:,1) = 0._r8
        end if
        if (rhs%vector_dims >= 1) then
          res%field_data(:,2) = lhs(3)*rhs%field_data(:,1)
          if (res%has_deriv) res%deriv_data(:,2) = lhs(3)*rhs%deriv_data(:,1)
        else
          res%field_data(:,2) = 0._r8
          if (res%has_deriv) res%deriv_data(:,2) = 0._r8
        end if
      else
        res%vector_dims = 3
        if (.not. allocated(res%field_data)) then
          allocate(res%field_data(rhs%numpoints,3))
        else if (size(res%field_data,1) /= rhs%numpoints .or. &
                 size(res%field_data,2) /= 3) then
          deallocate(res%field_data)
          allocate(res%field_data(rhs%numpoints,3))
        end if
        call res%allocate_deriv(rhs)
        dims1 = min(3,size(lhs))
        dims2 = min(3,rhs%vector_dims)
        vec1 = 0.0_r8
        vec2 = 0.0_r8
        dvec2 = 0.0_r8
        vec1(:dims1) = lhs(:dims1)
        if (res%has_deriv) then
          do concurrent(i=1:rhs%numpoints)
            vec2(:dims2) = rhs%field_data(i,:dims2)
            dvec2(:dims2) = rhs%deriv_data(i,:dims2)
            res%field_data(i,:) = [vec1(2)*vec2(3) - vec2(2)*vec1(3), &
                                   vec1(3)*vec2(1) - vec2(3)*vec1(1), &
                                   vec1(1)*vec2(2) - vec2(1)*vec1(2)]
            res%deriv_data(i,:) = [vec1(2)*dvec2(3) - dvec2(2)*vec1(3), &
                                   vec1(3)*dvec2(1) - dvec2(3)*vec1(1), &
                                   vec1(1)*dvec2(2) - dvec2(1)*vec1(2)]
          end do
        else
          do concurrent(i=1:rhs%numpoints)
            vec2(:dims2) = rhs%field_data(i,:dims2)
            res%field_data(i,:) = [vec1(2)*vec2(3) - vec2(2)*vec1(3), &
                                   vec1(3)*vec2(1) - vec2(3)*vec1(1), &
                                   vec1(1)*vec2(2) - vec2(1)*vec1(2)]
          end do
        end if
      end if
    class default
      error stop ('Non-array_vector_field type allocated by '//&
                 'allocate_vector_field routine.')
    end select
    call rhs%clean_temp()
  end function array_vector_vr_cross_vf

  function array_vector_vf_dot_vf(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\vec{\rm field} \cdot \vec{\rm field}\)
    !
    class(array_vector_field), intent(in) :: this
    class(vector_field), intent(in) :: rhs
    class(scalar_field), pointer :: res !! The restult of this operation
    integer :: i, min_dims
    real(r8), dimension(:), allocatable :: tmp
    call this%guard_temp(); call rhs%guard_temp()
#:if defined('DEBUG')
    call this%check_compatible(rhs)
#:endif
    call this%allocate_scalar_field(res)
    select type(res)
    class is(array_scalar_field)
      call res%assign_meta_data(this)
      call res%allocate_deriv(this, rhs)
      select type(rhs)
      class is(array_vector_field)
        min_dims = min(this%vector_dims,rhs%vector_dimensions())
        res%field_data = &
          sum(this%field_data(:,1:min_dims)*rhs%field_data(:,1:min_dims),2)
        if (this%has_deriv .and. rhs%has_deriv) then
          res%deriv_data = &
            sum(this%field_data(:,1:min_dims)*rhs%deriv_data(:,1:min_dims) + &
                this%deriv_data(:,1:min_dims)*rhs%field_data(:,1:min_dims), 2)
        else if (this%has_deriv) then
          res%deriv_data = &
            sum(this%deriv_data(:,1:min_dims)*rhs%field_data(:,1:min_dims), 2)
        else if (rhs%has_deriv) then
          res%deriv_data = &
            sum(this%field_data(:,1:min_dims)*rhs%deriv_data(:,1:min_dims), 2)
        end if
      class is(uniform_vector_field)
        min_dims = min(this%vector_dims,rhs%vector_dimensions())
        tmp = rhs%get_value()
        if (res%has_deriv) then
          do concurrent(i=1:this%numpoints)
            res%field_data(i) = dot_product(this%field_data(i,1:min_dims),tmp(1:min_dims))
            res%deriv_data(i) = dot_product(this%deriv_data(i,1:min_dims),tmp(1:min_dims))
          end do
        else
          do concurrent(i=1:this%numpoints)
            res%field_data(i) = dot_product(this%field_data(i,1:min_dims),tmp(1:min_dims))
          end do
        end if
      end select
    class default
      error stop ('Non-array_scalar_field type allocated by '//&
                 'allocate_scalar_field routine.')
    end select
    call this%clean_temp(); call rhs%clean_temp()
  end function array_vector_vf_dot_vf

  function array_vector_vf_dot_vr(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \(\vec{\rm field} \cdot \vec{\rm real}\)
    !
    class(array_vector_field), intent(in) :: this
    real(r8), dimension(:), intent(in) :: rhs
    class(scalar_field), pointer :: res !! The restult of this operation
    integer :: i, min_dims
    call this%guard_temp()
    call this%allocate_scalar_field(res)
    select type(res)
    class is(array_scalar_field)
      call res%assign_meta_data(this)
      call res%allocate_deriv(this)
      min_dims = min(this%vector_dims,size(rhs))
      if (res%has_deriv) then
        do concurrent(i=1:this%numpoints)
          res%field_data(i) = sum(this%field_data(i,1:min_dims)*rhs(1:min_dims))
          res%deriv_data(i) = sum(this%deriv_data(i,1:min_dims)*rhs(1:min_dims))
        end do
      else
        do concurrent(i=1:this%numpoints)
          res%field_data(i) = sum(this%field_data(i,1:min_dims)*rhs(1:min_dims))
        end do
      end if
    class default
      error stop ('Non-array_scalar_field type allocated by '//&
                 'allocate_scalar_field routine.')
    end select
    call this%clean_temp()
  end function array_vector_vf_dot_vr

  function array_vector_vr_dot_vf(lhs,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \(\vec{\rm real} \cdot \vec{\rm field}\)
    !
    real(r8), dimension(:), intent(in) :: lhs
    class(array_vector_field), intent(in) :: rhs
    class(scalar_field), pointer :: res !! The restult of this operation
    integer :: i, min_dims
    call rhs%guard_temp()
    call rhs%allocate_scalar_field(res)
    select type(res)
    class is(array_scalar_field)
      call res%assign_meta_data(rhs)
      call res%allocate_deriv(rhs)
      min_dims = min(size(lhs),rhs%vector_dims)
      if (res%has_deriv) then
        do concurrent(i=1:rhs%numpoints)
          res%field_data(i) = sum(lhs(1:min_dims)*rhs%field_data(i,1:min_dims))
          res%deriv_data(i) = sum(lhs(1:min_dims)*rhs%deriv_data(i,1:min_dims))
        end do
      else
        do concurrent(i=1:rhs%numpoints)
          res%field_data(i) = sum(lhs(1:min_dims)*rhs%field_data(i,1:min_dims))
        end do
      end if
    class default
      error stop ('Non-array_scalar_field type allocated by '//&
                  'allocate_scalar_field routine.')
    end select
    call rhs%clean_temp()
  end function array_vector_vr_dot_vf

  logical function array_vector_is_equal(this,rhs) result(iseq)
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
    real(r8), dimension(:), allocatable :: tmp
    call this%guard_temp(); call rhs%guard_temp()
    iseq = .true.
    select type(rhs)
    class is(array_vector_field)
      if (this%numpoints/=rhs%numpoints) then
        iseq = .false.
        call this%clean_temp(); call rhs%clean_temp()
        return
      end if
      if (this%vector_dims > rhs%vector_dims) then
        dims = rhs%vector_dims
        iseq = all(abs(this%field_data(:,dims:)) < get_tol())
      else if (this%vector_dims < rhs%vector_dims) then
        dims = this%vector_dims
        iseq = all(abs(rhs%field_data(:,dims:)) < get_tol())
      else
        dims = this%vector_dims
      end if
      if (.not. iseq) then
         call this%clean_temp(); call rhs%clean_temp()
         return
       end if
      do i=1,this%numpoints
        normalization = norm2(this%field_data(i,:dims))
        if (normalization < tiny(normalization)) normalization = 1.0_r8
        iseq = (norm2(this%field_data(i,1:dims) - rhs%field_data(i,1:dims))/ &
                normalization < get_tol()) .or. &
               all(is_nan(this%field_data(i,1:dims)) .eqv. &
                   is_nan(rhs%field_data(i,1:dims)))
        ! FIXME: This handling of NaNs will claim that things are
        ! equal when they aren't, just because they have a NAN in the
        ! same location.
        if (.not. iseq) then
          call this%clean_temp(); call rhs%clean_temp()
          return
        end if
      end do
    class is(uniform_vector_field)
      tmp = rhs%get_value()
      if (this%vector_dims > rhs%vector_dimensions()) then
        dims = rhs%vector_dimensions()
        iseq = all(abs(this%field_data(:,dims:)) < get_tol())
      else if (this%vector_dims < rhs%vector_dimensions()) then
        dims = this%vector_dimensions()
        iseq = all(abs(tmp(dims:)) < get_tol())
      else
        dims = this%vector_dims
      end if
      if (.not. iseq) then
        call this%clean_temp(); call rhs%clean_temp()
        return
      end if
      do i=1,this%numpoints
        normalization = norm2(this%field_data(i,:dims))
        if (normalization < get_tol()) normalization = 1.0_r8
        iseq = (norm2(this%field_data(i,1:dims) - tmp(1:dims))/normalization &
                < get_tol()) .or. &
               all(is_nan(this%field_data(i,1:dims)) .eqv. is_nan(tmp(1:dims))) 
        ! FIXME: This handling of NaNs will claim that things are
        ! equal when they aren't, just because they have a NAN in the
        ! same location.
        if (.not. iseq) then
          return
          call this%clean_temp(); call rhs%clean_temp()
        end if
      end do
    class default
      iseq = .false.
    end select
    call this%clean_temp(); call rhs%clean_temp()
  end function array_vector_is_equal

  subroutine array_vector_assign_meta_data(this, rhs, alloc)
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
    logical :: al
    call rhs%guard_temp()
    if (present(alloc)) then
      al = alloc
    else
      al = .true.
    end if
    call this%assign_subtype_meta_data(rhs)
    select type(rhs)
    class is(array_scalar_field)
      this%numpoints = rhs%numpoints
      this%vector_dims = rhs%dimensions()
      call allocate_like(this%field_data, rhs%field_data)
      this%has_deriv = rhs%has_deriv
    class is(array_vector_field)
      this%numpoints = rhs%numpoints
      this%vector_dims = rhs%vector_dims
      call allocate_like(this%field_data, rhs%field_data)
      this%has_deriv = rhs%has_deriv
    class default
      this%has_deriv = .false.
    end select
    call rhs%clean_temp()
  end subroutine array_vector_assign_meta_data

  subroutine array_vector_read_hdf(this, hdf_id, dataset_name, dims, &
                                   error)
    !* Author: Chris MacMackin
    !  Date: April 2017
    !
    ! Reads the contents of the field from a dataset in an HDF file.
    !
    class(array_vector_field), intent(inout)   :: this
    integer(hid_t), intent(in)                 :: hdf_id
      !! The identifier for the HDF file/group from which the field
      !! data is to be read
    character(len=*), intent(in)               :: dataset_name
      !! The name of the dataset to be read from the HDF file
    integer(hsize_t), dimension(:), intent(in) :: dims
      !! The dimensions of the data set to be read in. The last of
      !! these will be the vector dimension.
    integer, intent(out)                       :: error
      !! An error code which, upon succesful completion of the
      !! routine, is 0. Otherwise, contains the error code returned
      !! by the HDF library.
    call this%guard_temp()
    error = 0
    this%numpoints = int(product(dims(1:size(dims)-1)))
    this%vector_dims = int(dims(size(dims)))
    if (.not. allocated(this%field_data)) then
      allocate(this%field_data(this%numpoints, this%vector_dims))
    else if (size(this%field_data,1) /= this%numpoints .and. &
             size(this%field_data,2) /= this%vector_dims) then
      deallocate(this%field_data)
      allocate(this%field_data(this%numpoints,this%vector_dims))
    end if
    call h5ltread_dataset_double_f(hdf_id, dataset_name,        &
                                   this%field_data, dims, &
                                   error)
    call this%clean_temp()
  end subroutine array_vector_read_hdf

  subroutine array_vector_write_hdf(this, hdf_id, dataset_name, dims, &
                                    error)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Writes the contents of the field to a dataset in an HDF file.
    !
    class(array_vector_field), intent(in)      :: this
    integer(hid_t), intent(in)                 :: hdf_id
      !! The identifier for the HDF file/group in which the field
      !! data is to be written.
    character(len=*), intent(in)               :: dataset_name
      !! The name of the dataset to be created in the HDF file
      !! containing this field's data.
    integer(hsize_t), dimension(:), intent(in) :: dims
      !! The number of elements in each dimension of the fields
    integer, intent(out)                       :: error
      !! An error code which, upon succesful completion of the
      !! routine, is 0. Otherwise, contains the error code returned
      !! by the HDF library.
    call this%guard_temp()
#:if defined('DEBUG')
    if (product(dims) /= this%numpoints) &
      error stop ('Dimensions in requested output array do not match '// &
                  'number of data-points in the field.')
#:endif
    error = 0
    call h5ltmake_dataset_double_f(hdf_id, dataset_name, size(dims)+1,   &
                                   [dims,int(this%vector_dims,hsize_t)], &
                                   this%field_data, error)
    call this%clean_temp()
  end subroutine array_vector_write_hdf
  
  subroutine array_vector_compatible(this,other)
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
    call this%guard_temp(); call other%guard_temp()
    call this%check_subtype_compatible(other)
    select type(other)
    class is(array_scalar_field)
      if (this%numpoints /= other%numpoints) &
           error stop (err_message//'    different resolutions.')
      if (.not.(allocated(this%field_data).and.allocated(other%field_data))) &
           error stop (err_message//'    uninitialised fields.')
    class is(array_vector_field)
      if (this%numpoints /= other%numpoints) &
           error stop (err_message//'    different resolutions.')
      if (.not.(allocated(this%field_data).and.allocated(other%field_data))) &
           error stop (err_message//'    uninitialised fields.')
    class is(uniform_scalar_field)
      if (.not.(allocated(this%field_data))) &
           error stop (err_message//'    uninitialised fields.')
    class is(uniform_vector_field)
      if (.not.(allocated(this%field_data))) &
           error stop (err_message//'    uninitialised fields.')
    class default
      error stop (err_message//'    incompatible types.')
    end select
    call this%clean_temp(); call other%clean_temp()
  end subroutine array_vector_compatible

  function array_vector_get_element_vec(this,element) result(val)
    !* Author: Chris MacMackin
    !  Date: January 2017
    !
    ! Returns a vector of the field corresponding to the provided ID 
    ! number.
    !
    class(array_vector_field), intent(in) :: this
    integer, intent(in) :: element
      !! The ID number of the field element to be returned
    real(r8), allocatable, dimension(:) :: val
      !! The vector in the field corresponding to the specified ID
    call this%guard_temp()
    val = this%field_data(element,:)
    call this%clean_temp()
  end function array_vector_get_element_vec
  
  function array_vector_get_element_comp(this,element,component) result(val)
    !* Author: Chris MacMackin
    !  Date: January 2017
    !
    ! Returns a component of the vector of the field corresponding to 
    ! the provided ID number.
    !
    class(array_vector_field), intent(in) :: this
    integer, intent(in) :: element
      !! The ID number of the field element to be returned
    integer, intent(in) :: component
      !! The number of the vector component to be returned
    real(r8) :: val
      !! The vector component in the field corresponding to the 
      !! specified ID
    call this%guard_temp()
    val = this%field_data(element,component)
    call this%clean_temp()
  end function array_vector_get_element_comp

  subroutine array_vector_set_element_vec(this,element,val)
    !* Author: Chris MacMackin
    !  Date: January 2017
    !
    ! Sets the element of the field corresponding to the given ID to
    ! the given vector value.
    !
    class(array_vector_field), intent(inout) :: this
    integer, intent(in) :: element
      !! The ID number of the field element to be set
    real(r8), dimension(:), intent(in) :: val
      !! The new vector value the field element is to be set to
    call this%guard_temp()
    this%field_data(element,:) = val
    call this%clean_temp()
  end subroutine array_vector_set_element_vec

  subroutine array_vector_set_element_comp(this,element,component,val)
    !* Author: Chris MacMackin
    !  Date: January 2017
    !
    ! Sets the element of the field corresponding to the given ID to
    ! the given vector value.
    !
    class(array_vector_field), intent(inout) :: this
    integer, intent(in) :: element
      !! The ID number of the field element to be set
    integer, intent(in) :: component
      !! The number of the vector component to be returned
    real(r8), intent(in) :: val
      !! The new value of the vector component in the field element
    call this%guard_temp()
    this%field_data(element,component) = val
    call this%clean_temp()
  end subroutine array_vector_set_element_comp

  subroutine array_vector_set_deriv(this, deriv)
    !* Author: Chris MacMackin
    !  Date: March 2018
    !
    ! Sets a differential/derivative value for the field, which will
    ! be propagated through operations using automatic differentiation
    ! (tangent mode).
    !
    class(array_vector_field), intent(inout) :: this
    class(vector_field), intent(in)          :: deriv
      !! The derivative (i.e. differential) value of this field
    integer :: i
    call this%guard_temp(); call deriv%guard_temp()
#:if defined('DEBUG')
    call this%check_compatible(deriv)
#:endif
    select type(deriv)
    class is(array_vector_field)
      this%deriv_data = deriv%field_data
    class is(uniform_vector_field)
      if (allocated(this%field_data)) then
        call allocate_like(this%deriv_data, this%field_data)
        do i=1, this%numpoints
          this%deriv_data(i,:) = deriv%get_value()
        end do
      end if
    class default
      error stop ('Non-array_scalar_field and non-uniform_scalar_field '//&
                  'type passed as derivative value.')
    end select
    this%has_deriv = .true.
    call this%clean_temp(); call deriv%clean_temp()
  end subroutine array_vector_set_deriv
  
  function array_vector_get_deriv(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2018
    !
    ! Gets the differential/derivative value for this field. If no
    ! such information is available, returns a field with zero values.
    !
    class(array_vector_field), intent(in) :: this
    class(vector_field), pointer    :: res
    call this%guard_temp()
    call this%allocate_vector_field(res)
    call res%assign_meta_data(this)
    select type (res)
    class is(array_vector_field)
      if (allocated(this%deriv_data)) then
        res%field_data = this%deriv_data
      else
        res%field_data = 0._r8 * this%field_data
      end if   
      res%has_deriv = .false.
    class default
      error stop ('Non-array_scalar_field type allocated by '//&
                  '`allocate_scalar_field` routine.')
    end select
    call this%clean_temp()
  end function array_vector_get_deriv

  subroutine array_vector_alloc_deriv(this, field1, field2)
    !* Author: Chris MacMackin
    !  Date: March 2018
    !
    ! Helper function which decides whether a derivative is going to
    ! be calculated and allocates the array if necessary.
    !
    class(array_vector_field), intent(inout)  :: this
    class(abstract_field), intent(in)           :: field1
      !! The first field which is checked to see if it has a
      !! derivative.
    class(abstract_field), optional, intent(in) :: field2
      !! The optional second field which is checked to see if it has a
      !! derivative.
    logical :: alloc
    select type(field1)
    class is(array_scalar_field)
      alloc = field1%has_deriv
      if (present(field2) .and. .not. alloc) then
        select type(field2)
        class is(array_scalar_field)
          alloc = field2%has_deriv
        class is(array_vector_field)
          alloc = field2%has_deriv
        end select
      end if
    class is(array_vector_field)
      alloc = field1%has_deriv
      if (present(field2) .and. .not. alloc) then
        select type(field2)
        class is(array_scalar_field)
          alloc = field2%has_deriv
        class is(array_vector_field)
          alloc = field2%has_deriv
        end select
      end if
    class default
      alloc = .false.
    end select
    if (alloc) then
      call allocate_like(this%deriv_data, this%field_data)
      this%has_deriv = allocated(this%deriv_data)
    else
      this%has_deriv = .false.
    end if
  end subroutine array_vector_alloc_deriv

  function array_vector_get_bound(this,boundary,depth) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Returns a field containing the specified boundary
    ! information. For a array field this just means that a copy is
    ! returned.
    !
    class(array_vector_field), intent(in) :: this
    integer, intent(in) :: boundary
      !! Specifies which boundary is to be returned. The boundary
      !! will be the one normal to dimension of number
      !! `abs(boundary)`. If the argument is negative, then the
      !! lower boundary is returned. If positive, then the upper
      !! boundary is returned.
    integer, intent(in) :: depth
      !! The number of layers of data-points to return at the
      !! specified boundary.
    class(vector_field), pointer :: res
      !! A field, of the same type as `this` and with the same
      !! resolution, number of dimensions etc., but containing only
      !! the points within the specified number of layers of cells
      !! adjecent to the specified boundary.
    integer, dimension(:,:), allocatable :: slices
    integer :: i, j, k
    call this%guard_temp()
    if (boundary == 0) then
      allocate(res, mold=this)
      res = this
        call this%clean_temp()
      return
    end if
    call this%allocate_vector_field(res)
    select type(res)
    class is(array_vector_field)
      call res%subtype_boundary(this,boundary,depth,slices)
      res%vector_dims = this%vector_dims
      res%numpoints = sum(elements_in_slice(slices(1,:), slices(2,:), slices(3,:)))
      if (.not. allocated(res%field_data)) then
        allocate(res%field_data(res%numpoints,res%vector_dims))
      else if (size(res%field_data,1) /= res%numpoints .or. &
               size(res%field_data,2) /= res%vector_dims) then
        deallocate(res%field_data)
        allocate(res%field_data(res%numpoints,res%vector_dims))
      end if
      res%has_deriv = this%has_deriv
      if (res%has_deriv) then
        if (.not. allocated(res%deriv_data)) then
          allocate(res%deriv_data(res%numpoints,res%vector_dims))
        else if (size(res%deriv_data,1) /= res%numpoints .or. &
                 size(res%deriv_data,2) /= res%vector_dims) then
          deallocate(res%deriv_data)
          allocate(res%deriv_data(res%numpoints,res%vector_dims))
        end if
      else
        if (allocated(res%deriv_data)) deallocate(res%deriv_data)
      end if
      do concurrent (j=1:res%vector_dims)
        res%field_data(:,j) = &
             [(this%field_data(slices(1,i):slices(2,i):slices(3,i),j), &
             i=1, size(slices,2))]
        if (res%has_deriv) then
          res%deriv_data(:,j) = &
               [(this%deriv_data(slices(1,i):slices(2,i):slices(3,i),j), &
               i=1, size(slices,2))]
        end if
      end do
    class default
      error stop ('Non-array_scalar_field type allocated by '//&
                 'allocate_scalar_field routine.')
    end select
    call this%clean_temp()
  end function array_vector_get_bound

  function array_vector_is_allocated(this) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Returns a logical indicating whether both the array pointer
    ! _and_ the array which it holds are allocated.
    !
    class(array_vector_field), intent(in) :: this
    logical :: res
    call this%guard_temp()
    res = allocated(this%field_data)
    call this%clean_temp()
  end function array_vector_is_allocated

  elemental subroutine array_vector_finalise(this)
    !* Author: Chris MacMackin
    !  Date: January 2017
    !
    ! Deallocates the field data for this object.
    !
    class(array_vector_field), intent(inout) :: this
    if (allocated(this%field_data)) deallocate(this%field_data)
    call this%unset_temp()
  end subroutine array_vector_finalise

  function array_vector_interp(this, location) result(val)
    !* Author: Chris MacMackin
    !  Date: July 2017
    !
    ! Interpolates the value of the field at the specified location.
    !
    class(array_vector_field), intent(in) :: this
    real(r8), dimension(:), intent(in)    :: location
      !! The location at which to calculate the interpolated value.
    real(r8), dimension(:), allocatable   :: val
    call this%guard_temp()
    val = this%array_interpolate(this%field_data, location)
    call this%clean_temp()
  end function array_vector_interp

end module array_fields_mod
