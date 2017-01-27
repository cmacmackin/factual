!
!  uniform_fields.f90
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

module uniform_fields_mod
  !* Author: Chris MacMackin
  !  Date: November 2016
  !  License: GPLv3
  !
  ! Provides extensions of the fields in [[abstract_fields_mod]]
  ! represents a field which is spatially uniform. This can be used to
  ! optimise memory management and differentiation in simulations
  ! without having to change interfaces to use scalar real values.
  !
  use iso_fortran_env, only: r8 => real64
  use abstract_fields_mod
  use utils_mod, only: is_nan, check_set_from_raw, elements_in_slice
  use h5lt, only: hid_t, size_t, hsize_t, h5ltmake_dataset_double_f, &
                  h5ltset_attribute_string_f, h5ltset_attribute_int_f
  implicit none
  private

  character(len=20), parameter, public :: hdf_scalar_name = 'uniform_scalar_field'
  character(len=20), parameter, public :: hdf_vector_name = 'uniform_vector_field'

  integer, parameter :: grid_spacing_size = 10
  
  $:public_unary()
  public :: minval
  public :: maxval

  type, extends(scalar_field), public :: uniform_scalar_field
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! An abstract data type representing a mathematical field of
    ! scalar values, in which the field is uniform across all of
    ! space. This allows for reduced memory usage and avoiding having
    ! to apply a differentiation algorithm, compared to using another
    ! field type with all values set to be the same. Additionally,
    ! there is no need to change data types or interfaces in client
    ! code, as would be required if one were to switch from using a
    ! field to a real scalar. Assuming other fields have been properly
    ! implemented, this field type should be interoperable with all of
    ! them.
    !
    private
    real(r8) :: field_data
      !! The value of the scalar field.
  contains
    private
    procedure, public :: domain => uniform_scalar_domain
      !! Provides array with upper and lower limits of field's domain
    procedure, non_overridable, public :: elements => uniform_scalar_elements
      !! Specifies the number of individual data points present in this field.
    procedure, public :: dimensions => uniform_scalar_dimensions
      !! Returns number of dimensions of the field
    procedure, public :: raw_size => uniform_scalar_raw_size
      !! Provides the number of pieces of data needed to represent the
      !! field, i.e. the size of the uniform returned by `raw`.
    procedure, public :: raw => uniform_scalar_raw
      !! Returns uniform of data representing state of field. Can be
      !! useful for passing to nonlinear solvers
    procedure, non_overridable, public :: get_value => uniform_scalar_get_value
      !! Returns the scalar value held in the uniform field
    procedure, public :: set_from_raw => uniform_scalar_set_from_raw
      !! Assigns raw data, such as that produced by 
      !! [[uniform_scalar_field:raw]], to the field
    procedure, public :: resolution => uniform_scalar_resolution
      !! Returns array containing number of datapoints in each dimension.
    procedure :: field_multiply_field => uniform_scalar_sf_m_sf
      !! \({\rm field} \times {\rm field}\)
    procedure :: field_multiply_vecfield => uniform_scalar_sf_m_vf
      !! \({\rm field} \times {\rm \vec{field}}\)
    procedure, pass(rhs) :: real_multiply_field => uniform_scalar_r_m_sf
      !! \({\rm real} \times {\rm field}\)
    procedure, pass(rhs) :: real_array_multiply_field => uniform_scalar_vr_m_sf
      !! \(\vec{\rm real} \times {\rm field}\)
    procedure :: field_multiply_real => uniform_scalar_sf_m_r
      !! \({\rm field} \times {\rm real}\)
    procedure :: field_multiply_real_array => uniform_scalar_sf_m_vr
      !! \({\rm field} \times \vec{\rm real}\)
    procedure :: field_divide_field => uniform_scalar_sf_d_sf
      !! \(\frac{\rm field}{\rm field}\)
    procedure, pass(rhs) :: real_divide_field => uniform_scalar_r_d_sf
      !! \(\frac{\rm real}{\rm field}\)
    procedure, pass(rhs) :: real_array_divide_field => uniform_scalar_vr_d_sf
      !! \(\frac{\vec{\rm real}}{\rm field}\)
    procedure :: field_divide_real => uniform_scalar_sf_d_r
      !! \(\frac{\rm field}{\rm real}\)
    procedure :: field_add_field => uniform_scalar_sf_a_sf
      !! \({\rm field} + {\rm field}\)
    procedure, pass(rhs) :: real_add_field => uniform_scalar_r_a_sf
      !! \({\rm real} + {\rm field}\)
    procedure :: field_add_real => uniform_scalar_sf_a_r
      !! \({\rm field} + {\rm real}\)
    procedure :: field_sub_field => uniform_scalar_sf_s_sf
      !! \({\rm field} - {\rm field}\)
    procedure, pass(rhs) :: real_sub_field => uniform_scalar_r_s_sf
      !! \({\rm real} - {\rm field}\)
    procedure :: field_sub_real => uniform_scalar_sf_s_r
      !! \({\rm field} - {\rm real}\)
    procedure :: field_pow_real => uniform_scalar_sf_p_r
      !! \({\rm field}^{\rm real}\)
    procedure :: field_pow_real4 => uniform_scalar_sf_p_r4
      !! \({\rm field}^{\rm real}\)
    procedure :: field_pow_int => uniform_scalar_sf_p_i
      !! \({\rm field}^{\rm int}\)
#:for FUNC, TEX in UNARY_FUNCTIONS
    $:unary_binding(FUNC, TEX, 'uniform_scalar')
#:endfor
    procedure :: minval => uniform_scalar_minval
      !! \(\min({\rm field})\)
    procedure :: maxval => uniform_scalar_maxval
      !! \(\max({\rm field})\)
    procedure, public :: d_dx => uniform_scalar_d_dx
      !! \(\frac{\partial^n}{\partial x_i^n}({\rm field})\)
    procedure :: gradient => uniform_scalar_gradient
      !! \(\nabla {\rm field}\)
    procedure :: laplacian => uniform_scalar_laplacian
      !! \(\nabla^2 {\rm field}\)
    procedure :: is_equal => uniform_scalar_is_equal
      !! Checks fields are equal within a tolerance
    procedure :: assign_field => uniform_scalar_assign
      !! \({\rm field} = {\rm field}\)
    procedure, public :: assign_meta_data => uniform_scalar_assign_meta_data
      !! Copies all data other than values stored in field from another
      !! field object to this one.
    procedure, public :: allocate_scalar_field => uniform_scalar_allocate_scalar
      !! Allocates a scalar field with concrete type [[uniform_scalar_field]]
    procedure, public :: allocate_vector_field => uniform_scalar_allocate_vector
      !! Allocates a vector field with concrete type [[uniform_vector_field]]
    procedure, public :: id_to_position => uniform_scalar_id_to_pos
      !! Given the ID number of a location in the field, returns the
      !! coordinates of that position
    procedure, public :: write_hdf => uniform_scalar_write_hdf
      !! Write field data to a new dataset in an HDF5 file.
    procedure, public :: grid_spacing => uniform_scalar_grid_spacing
      !! Returns a vector field where the values are the size of cells
      !! in the grid of the field at each point, in each
      !! direction. This can be useful, e.g., when calculating time
      !! steps for integration.
    procedure, public :: get_element => uniform_scalar_get_element
      !! Returns one of the constituent values of the field, i.e. the 
      !! field's value at a particular location.
    procedure, public :: set_element => uniform_scalar_set_element
      !! Sets one of the constituent values of the field, i.e. the 
      !! field's value at a particular location.
    procedure, public :: get_boundary => uniform_scalar_get_bound
      !! Returns a field of the same type, containing only the
      !! specified ammount of data at the specified boundary.
  end type uniform_scalar_field

  interface uniform_scalar_field
    module procedure uniform_scalar_constructor
  end interface uniform_scalar_field


  type, extends(vector_field), public :: uniform_vector_field
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! An abstract data type representing a mathematical field of
    ! vector values, in which the field is uniform across all of
    ! space. This allows for reduced memory usage and avoiding having
    ! to apply a differentiation algorithm, compared to using another
    ! field type with all values set to be the same. Additionally,
    ! there is no need to change data types or interfaces in client
    ! code, as would be required if one were to switch from using a
    ! field to a real scalar. Assuming other fields have been properly
    ! implemented, this field type should be interoperable with all of
    ! them.
    !
    private
    real(r8), dimension(:), allocatable :: field_data
      !! The value of the vector field across all of space.
    integer :: vector_dims = 0
      !! The number of vector components
  contains
    private
    procedure, public :: domain => uniform_vector_domain
      !! Provides array with upper and lower limits of field's domain
    procedure, non_overridable, public :: &
         vector_dimensions => uniform_vector_vector_dimensions
      !! Returns dimension of the vectors in the field
    procedure, public :: elements => uniform_vector_elements
      !! Specifies the number of individual data points present in this field. 
    procedure, public :: dimensions => uniform_vector_dimensions
      !! Returns number of dimensions of the field
    procedure, public :: raw_size => uniform_vector_raw_size
      !! Provides the number of pieces of data needed to represent the
      !! field, i.e. the size of the uniform returned by `raw`.
    procedure, public :: raw => uniform_vector_raw
      !! Returns uniform of data representing state of field. Can be
      !! useful for passing to nonlinear solvers.
    procedure, non_overridable, public :: get_value => uniform_vector_get_value
      !! Returns the vector value held in the uniform field
    procedure, public :: set_from_raw => uniform_vector_set_from_raw
      !! Assigns raw data, such as that produced by 
      !! [[uniform_vector_field:raw]], to the field
    procedure, public :: resolution => uniform_vector_resolution
      !! Returns array containing number of datapoints in each dimension.
    procedure :: field_multiply_field => uniform_vector_vf_m_sf
      !! \({\rm field} \times {\rm field}\)
    procedure, pass(rhs) :: real_multiply_field => uniform_vector_r_m_vf
      !! \({\rm real}  \times {\rm field}\)
    procedure :: field_multiply_real => uniform_vector_vf_m_r
      !! \({\rm field} \times {\rm real}\)
    procedure :: field_divide_field => uniform_vector_vf_d_sf
      !! \(\frac{\rm field}{\rm field}\)
    procedure :: field_divide_real => uniform_vector_vf_d_r
      !! \(\frac{\rm field}{\rm real}\)
    procedure :: field_add_field => uniform_vector_vf_a_vf
      !! \({\rm field} + {\rm field}\)
    procedure, pass(rhs) :: real_add_field => uniform_vector_r_a_vf
      !! \({\rm real} + {\rm field}\)
    procedure :: field_add_real => uniform_vector_vf_a_r
      !! \({\rm field} + {\rm real}\)
    procedure :: field_sub_field => uniform_vector_vf_s_vf
      !! \({\rm field} - {\rm field}\)
    procedure, pass(rhs) :: real_sub_field => uniform_vector_r_s_vf
      !! \({\rm real} - {\rm field}\)
    procedure :: field_sub_real => uniform_vector_vf_s_r
      !! \({\rm field} - {\rm real}\)
    procedure, public :: norm => uniform_vector_norm
      !! \(\lVert {\rm \vec{field}} \rVert\)
    procedure, public :: component => uniform_vector_component
      !! Returns a scalar field containing the specified component of 
      !! the vector field
    procedure, public :: d_dx => uniform_vector_d_dx
      !! \(\frac{\partial^n}{\partial x_i^n}({\rm \vec{field}})\)
    procedure, public :: component_d_dx => uniform_vector_component_d_dx
      !! \(\frac{\partial^n}{\partial x_i^n}({\rm field_j})\)
    procedure :: divergence => uniform_vector_divergence
      !! \(\nabla\cdot {\rm field}\)
    procedure :: curl => uniform_vector_curl
      !! \(\nabla\times {\rm field}\)
    procedure :: laplacian => uniform_vector_laplacian
      !! \(\nabla^2 {\rm field}\)
    procedure :: field_dot_field => uniform_vector_vf_dot_vf
      !! \({\rm \vec{field}} \cdot {\rm \vec{field}}\)
    procedure :: field_dot_real => uniform_vector_vf_dot_vr
      !! \({\rm \vec{field}} \cdot {\rm \vec{real}}\)
    procedure, pass(rhs) :: real_dot_field => uniform_vector_vr_dot_vf
      !! \({\rm \vec{real}} \cdot {\rm \vec{field}}\)
    procedure :: field_cross_field => uniform_vector_vf_cross_vf
      !! \({\rm\vec{field}} \times {\rm\vec{field}}\)
    procedure :: field_cross_real => uniform_vector_vf_cross_vr
      !! \({\rm\vec{field}} \times {\rm\vec{real}}\)
    procedure, pass(rhs) :: real_cross_field => uniform_vector_vr_cross_vf
      !! \({\rm\vec{real}} \times {\rm\vec{field}}\)
    procedure :: assign_field => uniform_vector_assign
      !! \({\rm field} = {\rm field}\)
    procedure :: is_equal => uniform_vector_is_equal
      !! Checks fields are equal within a tolerance
    procedure, public :: assign_meta_data => uniform_vector_assign_meta_data
      !! Copies all data other than values stored in field from another
      !! field object to this one.
    procedure, public :: allocate_scalar_field => uniform_vector_allocate_scalar
      !! Allocates a scalar field with concrete type [[uniform_scalar_field]]
    procedure, public :: allocate_vector_field => uniform_vector_allocate_vector
      !! Allocates a vector field with concrete type [[uniform_vector_field]]
    procedure, public :: id_to_position => uniform_vector_id_to_pos
      !! Given the ID number of a location in the field, returns the
      !! coordinates of that position
    procedure, public :: write_hdf => uniform_vector_write_hdf
      !! Write field data to a new dataset in an HDF5 file.
    procedure, public :: grid_spacing => uniform_vector_grid_spacing
      !! Returns a vector field where the values are the size of cells
      !! in the grid of the field at each point, in each
      !! direction. This can be useful, e.g., when calculating time
      !! steps for integration.
    procedure :: get_element_vector => uniform_vector_get_element_vec
      !! Returns ones of the constituent vectors of the field, i.e. the 
      !! field's value at a particular location.
    procedure :: get_element_component => uniform_vector_get_element_comp
      !! Returns one of the components of a constituent vector of the 
      !! field, i.e. the component of the field's value at a particular 
      !! location.
    procedure :: set_element_vector => uniform_vector_set_element_vec
      !! Sets ones of the constituent vectors of the field, i.e. the 
      !! field's value at a particular location.
    procedure :: set_element_component => uniform_vector_set_element_comp
      !! Sets one of the components of a constituent vector of the 
      !! field, i.e. the component of the field's value at a particular 
      !! location.
    procedure, public :: get_boundary => uniform_vector_get_bound
      !! Returns a field of the same type, containing only the
      !! specified ammount of data at the specified boundary.
    final :: uniform_vector_finalize
      !! Deallocates all field contents for this object. It is a
      !! workaround for the fact that `gfortran` fails to deallocate
      !! polymorphic allocatable function results. While some memory
      !! will still be leaked, this will minimize the ammount.
  end type uniform_vector_field

  interface uniform_vector_field
    module procedure uniform_vector_constructor
  end interface uniform_vector_field

  
contains


  !=====================================================================
  ! Scalar Field Methods
  !=====================================================================

  function uniform_scalar_constructor(val) result(this)
    !* Author: Chris MacMackin
    !  Date: October 2016
    !
    ! Creates a new scalar field with a uniform value across all of space.
    !
    real(r8), intent(in) :: val
      !! The value of the field
    type(uniform_scalar_field) :: this
      !! A scalar field initiated based on the arguments to this function.
    this%field_data = val
  end function uniform_scalar_constructor
  
  pure function uniform_scalar_elements(this) result(elements)
    !* Author: Chris MacMackin
    !  Date: October 2016
    !
    ! Gives the number of individual data points present in the field.
    ! For a uniform field this is 1.
    !
    class(uniform_scalar_field), intent(in) :: this
    integer :: elements
    elements = 1
  end function uniform_scalar_elements

  pure function uniform_scalar_domain(this) result(domain)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Returns an array specifying the upper and lower bounds of the
    ! data in this field. As this field is uniform across all of
    ! space, the concept of limits is not meaningful.
    !
    class(uniform_scalar_field), intent(in) :: this
    real(r8), dimension(:,:), allocatable :: domain
    allocate(domain(1,2))
    domain(1,1) = 0.0_r8
    domain(1,2) = 0.0_r8
  end function uniform_scalar_domain

  elemental function uniform_scalar_dimensions(this) result(dims)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Gives the number of dimensions of the field. For a uniform
    ! field, this is effectively 0.
    !
    class(uniform_scalar_field), intent(in) :: this
    integer :: dims
    dims = 0
  end function uniform_scalar_dimensions

  pure function uniform_scalar_raw_size(this,exclude_lower_bound, &
                                        exclude_upper_bound) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Compute how many elements are in the raw representation of this
    ! field. For a uniform field this is 1.
    !
    class(uniform_scalar_field), intent(in)      :: this
    integer, dimension(:), optional, intent(in) :: exclude_lower_bound
      !! Specifies how many layers of data points should be excluded
      !! from the result at the lower boundary for each dimension. The
      !! number in element `n` of the array indicates how many layers
      !! of cells at the lower boundary normal to dimension `n` will
      !! be ignored. For a uniform field, this argument has no effect.
    integer, dimension(:), optional, intent(in) :: exclude_upper_bound
      !! Specifies how many layers of data points should be excluded
      !! from the result at the upper boundary for each dimension. The
      !! number in element `n` of the array indicates how many layers
      !! of cells at the upper boundary normal to dimension `n` will
      !! be ignored. For a uniform field, this argument has no effect.
    integer :: res
    res = 1
  end function uniform_scalar_raw_size
  
  pure function uniform_scalar_raw(this,exclude_lower_bound, &
                                   exclude_upper_bound) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Returns a representation of this field in the form of a 1D uniform
    ! of real numbers this allows, e.g. for manipulation by solver
    ! routines.
    !
    ! @BUG The returned value has length `this%raw_size()`, but
    ! a bug in gfortran 4.8 (fixed by version 5) caused the compiler
    ! to segfault if it was declared as such. As a workaround, it is
    ! allocatable isntead.
    !
    class(uniform_scalar_field), intent(in) :: this
    integer, dimension(:), optional, intent(in) :: exclude_lower_bound
      !! Specifies how many layers of data points should be excluded
      !! from the result at the lower boundary for each dimension. The
      !! number in element `n` of the array indicates how many layers
      !! of cells at the lower boundary normal to dimension `n` will
      !! be ignored. For a uniform field, this argument has no effect.
    integer, dimension(:), optional, intent(in) :: exclude_upper_bound
      !! Specifies how many layers of data points should be excluded
      !! from the result at the upper boundary for each dimension. The
      !! number in element `n` of the array indicates how many layers
      !! of cells at the upper boundary normal to dimension `n` will
      !! be ignored. For a uniform field, this argument has no effect.
    real(r8), dimension(:), allocatable :: res
      !! Uniform containing data needed to describe field
    allocate(res(1))
    res(1) = this%field_data
  end function uniform_scalar_raw

#:if defined('DEBUG')
#undef pure 
#undef elemental 
#:endif

  pure function uniform_scalar_get_value(this) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Returns the scalar value of the field. This is useful for
    ! implementing the interactions of other field types with the
    ! uniform field.
    !
    class(uniform_scalar_field), intent(in) :: this
    real(r8) :: res
    res = this%field_data
  end function uniform_scalar_get_value

#:if defined('DEBUG')
#define pure 
#define elemental 
#:endif

  pure subroutine uniform_scalar_set_from_raw(this,raw,provide_lower_bound, &
                                            provide_upper_bound)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Assigns raw data, such as that produced by 
    ! [[uniform_scalar_field:raw]], to the field. The routine will
    ! stop with an error if the uniform is the wrong length for a field
    ! of this resolution and with these boundary conditions. 
    !
    class(uniform_scalar_field), intent(inout) :: this
    real(r8), dimension(:), intent(in) :: raw
      !! The raw data to be stored in this uniform.
    integer, dimension(:), optional, intent(in) :: provide_lower_bound
      !! Specifies how many layers of data points are excluded from
      !! the raw data at the lower boundary for each dimension. The
      !! number in element `n` of the array indicates how many layers
      !! of cells at the lower boundary normal to dimension `n` are
      !! missed. For a uniform field, this argument has no effect.
    integer, dimension(:), optional, intent(in) :: provide_upper_bound
      !! Specifies how many layers of data points are excluded
      !! from the raw data at the upper boundary for each
      !! dimension. The number in element `n` of the array indicates
      !! how many layers of cells at the upper boundary normal to
      !! dimension `n` are missed. Defaults to 0 for all.
    call check_set_from_raw(this,raw,provide_lower_bound,provide_upper_bound)
    this%field_data = raw(1)
  end subroutine uniform_scalar_set_from_raw

  pure function uniform_scalar_resolution(this) result(resolution)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Returns an array specifying the number of data-points in each
    ! dimension. This is not a particularly meaningful quantity for a
    ! uniform array, as only a single data-point is needed.
    !
    class(uniform_scalar_field), intent(in)  :: this
    integer, dimension(:), allocatable :: resolution
      !! Array specifying the number of data points in each dimension.
    allocate(resolution(0))
  end function uniform_scalar_resolution
  
  pure function uniform_scalar_sf_m_sf(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \({\rm field} \times {\rm field}\)
    !
    class(uniform_scalar_field), intent(in) :: this
    class(scalar_field), intent(in) :: rhs
    class(scalar_field), allocatable :: res !! The restult of this operation
    allocate(res, mold=rhs)
    res = this%field_data * rhs
  end function uniform_scalar_sf_m_sf

  pure function uniform_scalar_sf_m_vf(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \({\rm field} \times {\rm \vec{field}}\)
    !
    class(uniform_scalar_field), intent(in) :: this
    class(vector_field), intent(in) :: rhs
    class(vector_field), allocatable :: res !! The result of this operation
    allocate(res, mold=rhs)
    res = this%field_data * rhs
  end function uniform_scalar_sf_m_vf

  pure function uniform_scalar_r_m_sf(lhs,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \({\rm real} \times {\rm field}\)
    !
    real(r8), intent(in) :: lhs
    class(uniform_scalar_field), intent(in) :: rhs
    class(scalar_field), allocatable :: res !! The result of this operation
    type(uniform_scalar_field), allocatable :: local
    allocate(local)
    local%field_data = lhs * rhs%field_data
    call move_alloc(local,res)
  end function uniform_scalar_r_m_sf

  pure function uniform_scalar_vr_m_sf(lhs,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \(\vec{\rm real} \times {\rm field}\)
    !
    real(r8), dimension(:), intent(in) :: lhs
    class(uniform_scalar_field), intent(in) :: rhs
    class(vector_field), allocatable :: res !! The result of this operation
    type(uniform_vector_field), allocatable :: local
    allocate(local)
    allocate(local%field_data(size(lhs)))
    local%vector_dims = size(lhs)
    local%field_data = lhs * rhs%field_data
    call move_alloc(local,res)
  end function uniform_scalar_vr_m_sf

  pure function uniform_scalar_sf_m_r(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \({\rm field} \times {\rm real}\)
    !
    class(uniform_scalar_field), intent(in) :: this
    real(r8), intent(in) :: rhs
    class(scalar_field), allocatable :: res !! The result of this operation
    type(uniform_scalar_field), allocatable :: local
    allocate(local)
    local%field_data = this%field_data * rhs
    call move_alloc(local,res)
  end function uniform_scalar_sf_m_r

  pure function uniform_scalar_sf_m_vr(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \({\rm field} \times \vec{\rm real}\)
    !
    class(uniform_scalar_field), intent(in) :: this
    real(r8), dimension(:), intent(in) :: rhs
    class(vector_field), allocatable :: res !! The result of this operation
    type(uniform_vector_field), allocatable :: local
    allocate(local)
    allocate(local%field_data(size(rhs)))
    local%vector_dims = size(rhs)
    local%field_data = this%field_data * rhs
    call move_alloc(local,res)
  end function uniform_scalar_sf_m_vr
  
  pure function uniform_scalar_sf_d_sf(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \({\rm field} / {\rm field}\)
    !
    class(uniform_scalar_field), intent(in) :: this
    class(scalar_field), intent(in) :: rhs
    class(scalar_field), allocatable :: res !! The restult of this operation
    allocate(res, mold=rhs)
    res = this%field_data / rhs    
  end function uniform_scalar_sf_d_sf

  pure function uniform_scalar_r_d_sf(lhs,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \({\rm real} / {\rm field}\)
    !
    real(r8), intent(in) :: lhs
    class(uniform_scalar_field), intent(in) :: rhs
    class(scalar_field), allocatable :: res !! The result of this operation
    type(uniform_scalar_field), allocatable :: local
    allocate(local)
    local%field_data = lhs / rhs%field_data
    call move_alloc(local,res)
  end function uniform_scalar_r_d_sf

  pure function uniform_scalar_vr_d_sf(lhs,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \(\vec{\rm real} / {\rm field}\)
    !
    real(r8), dimension(:), intent(in) :: lhs
    class(uniform_scalar_field), intent(in) :: rhs
    class(vector_field), allocatable :: res !! The result of this operation
    type(uniform_vector_field), allocatable :: local
    allocate(local)
    allocate(local%field_data(size(lhs)))
    local%vector_dims = size(lhs)
    local%field_data = lhs / rhs%field_data
    call move_alloc(local,res)
  end function uniform_scalar_vr_d_sf

  pure function uniform_scalar_sf_d_r(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \({\rm field} / {\rm real}\)
    !
    class(uniform_scalar_field), intent(in) :: this
    real(r8), intent(in) :: rhs
    class(scalar_field), allocatable :: res !! The result of this operation
    type(uniform_scalar_field), allocatable :: local
    allocate(local)
    local%field_data = this%field_data / rhs
    call move_alloc(local,res)
  end function uniform_scalar_sf_d_r
  
  pure function uniform_scalar_sf_s_sf(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \({\rm field} - {\rm field}\)
    !
    class(uniform_scalar_field), intent(in) :: this
    class(scalar_field), intent(in) :: rhs
    class(scalar_field), allocatable :: res !! The restult of this operation
    allocate(res, mold=rhs)
    res = this%field_data - rhs
  end function uniform_scalar_sf_s_sf

  pure function uniform_scalar_r_s_sf(lhs,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \({\rm real} - {\rm field}\)
    !
    real(r8), intent(in) :: lhs
    class(uniform_scalar_field), intent(in) :: rhs
    class(scalar_field), allocatable :: res !! The result of this operation
    type(uniform_scalar_field), allocatable :: local
    allocate(local)
    local%field_data = lhs - rhs%field_data
    call move_alloc(local,res)
  end function uniform_scalar_r_s_sf

  pure function uniform_scalar_sf_s_r(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \({\rm field} - {\rm real}\)
    !
    class(uniform_scalar_field), intent(in) :: this
    real(r8), intent(in) :: rhs
    class(scalar_field), allocatable :: res !! The result of this operation
    type(uniform_scalar_field), allocatable :: local
    allocate(local)
    local%field_data = this%field_data - rhs
    call move_alloc(local,res)
  end function uniform_scalar_sf_s_r
  
  pure function uniform_scalar_sf_a_sf(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \({\rm field} + {\rm field}\)
    !
    class(uniform_scalar_field), intent(in) :: this
    class(scalar_field), intent(in) :: rhs
    class(scalar_field), allocatable :: res !! The restult of this operation
    allocate(res, mold=this)
    res = this%field_data + rhs
  end function uniform_scalar_sf_a_sf

  pure function uniform_scalar_r_a_sf(lhs,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \({\rm real} + {\rm field}\)
    !
    real(r8), intent(in) :: lhs
    class(uniform_scalar_field), intent(in) :: rhs
    class(scalar_field), allocatable :: res !! The result of this operation
    type(uniform_scalar_field), allocatable :: local
    allocate(local)
    local%field_data = lhs + rhs%field_data
    call move_alloc(local,res)
  end function uniform_scalar_r_a_sf

  pure function uniform_scalar_sf_a_r(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \({\rm field} + {\rm real}\)
    !
    class(uniform_scalar_field), intent(in) :: this
    real(r8), intent(in) :: rhs
    class(scalar_field), allocatable :: res !! The result of this operation
    type(uniform_scalar_field), allocatable :: local
    allocate(local)
    local%field_data = this%field_data + rhs
    call move_alloc(local,res)
  end function uniform_scalar_sf_a_r

  pure function uniform_scalar_sf_p_r(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \({\rm field}^{\rm real}\)
    !
    class(uniform_scalar_field), intent(in) :: this
    real(r8), intent(in) :: rhs
    class(scalar_field), allocatable :: res !! The result of this operation
    type(uniform_scalar_field), allocatable :: local
    allocate(local)
    local%field_data = this%field_data ** rhs
    call move_alloc(local,res)
  end function uniform_scalar_sf_p_r

  pure function uniform_scalar_sf_p_r4(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \({\rm field}^{\rm real}\)
    !
    class(uniform_scalar_field), intent(in) :: this
    real, intent(in) :: rhs
    class(scalar_field), allocatable :: res !! The result of this operation
    type(uniform_scalar_field), allocatable :: local
    allocate(local)
    local%field_data = this%field_data ** rhs
    call move_alloc(local,res)
  end function uniform_scalar_sf_p_r4

  pure function uniform_scalar_sf_p_i(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \({\rm field}^{\rm integer}\)
    !
    class(uniform_scalar_field), intent(in) :: this
    integer, intent(in) :: rhs
    class(scalar_field), allocatable :: res !! The result of this operation
    type(uniform_scalar_field), allocatable :: local
    allocate(local)
    local%field_data = this%field_data ** rhs
    call move_alloc(local,res)
  end function uniform_scalar_sf_p_i

#:for FUNC, TEX in UNARY_FUNCTIONS
  pure function uniform_scalar_${FUNC}$(this) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \(${TEX}$({\rm field})\)
    !
    class(uniform_scalar_field), intent(in)  :: this
    class(scalar_field), allocatable        :: res !! The result of this operation
    type(uniform_scalar_field), allocatable :: local
    allocate(local)
    local%field_data = ${FUNC}$(this%field_data)
    call move_alloc(local, res)
  end function uniform_scalar_${FUNC}$

#:endfor

  pure function uniform_scalar_minval(this) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \(\min({\rm field})\)
    !
    class(uniform_scalar_field), intent(in) :: this
    real(r8) :: res !! The result of this operation
    res = this%field_data
  end function uniform_scalar_minval

  pure function uniform_scalar_maxval(this) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \(\max({\rm field})\)
    !
    class(uniform_scalar_field), intent(in) :: this
    real(r8) :: res !! The result of this operation
    res = this%field_data
  end function uniform_scalar_maxval

  function uniform_scalar_d_dx(this, dir, order) result(res)
    !* Author: Chris MacMackin
    !  Date: October 2016
    !
    ! \(\frac{\partial^{\rm order}}{\partial x_{\rm dir}^{\rm order}}{\rm field}\)
    !
    class(uniform_scalar_field), intent(in) :: this
    integer, intent(in) :: dir !! Direction in which to differentiate
    integer, optional, intent(in) :: order !! Order of the derivative, default = 1
    class(scalar_field), allocatable :: res
    type(uniform_scalar_field), allocatable :: local
    allocate(local)
    local%field_data = 0.0_r8
    call move_alloc(local, res)
  end function uniform_scalar_d_dx

  function uniform_scalar_laplacian(this) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \(\nabla^2 {\rm field}\)
    !
    class(uniform_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(uniform_scalar_field), allocatable :: local
    allocate(local)
    local%field_data = 0.0_r8
    call move_alloc(local, res)
  end function uniform_scalar_laplacian
  
  function uniform_scalar_gradient(this) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \(\nabla{\rm field}\)
    !
    class(uniform_scalar_field), intent(in) :: this
    class(vector_field), allocatable :: res !! The result of this operation
    type(uniform_vector_field), allocatable :: local
    allocate(local)
    allocate(local%field_data(1))
    local%field_data = 0.0_r8
    call move_alloc(local, res)
  end function uniform_scalar_gradient
  
  elemental subroutine uniform_scalar_assign(this,rhs)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \({\rm field} = {\rm field}\)
    !
    class(uniform_scalar_field), intent(inout) :: this
    class(scalar_field), intent(in) :: rhs
    select type(rhs)
    class is(uniform_scalar_field)
      this%field_data = rhs%field_data
    class default
      error stop ('Assigning incompatible type to uniform_scalar_field')
    end select
  end subroutine uniform_scalar_assign

  pure logical function uniform_scalar_is_equal(this,rhs) result(iseq)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Evaluates whether two scalar fields are equal within a tolerance,
    ! specified by [[set_tol]].
    !
    class(uniform_scalar_field), intent(in) :: this
    class(scalar_field), intent(in) :: rhs
    real(r8) :: normalization
    iseq = .true.
    select type(rhs)
    class is(uniform_scalar_field)
      normalization = abs(this%field_data)
      if (normalization < get_tol()) normalization = 1.0_r8
      iseq = iseq .and.( ((this%field_data-rhs%field_data)/normalization < &
                           get_tol()) .or. (is_nan(this%field_data).and. &
                                            is_nan(rhs%field_data)) )
    class default
      iseq = (rhs == this)
    end select
  end function uniform_scalar_is_equal

  pure subroutine uniform_scalar_assign_meta_data(this, rhs, alloc)
    !* Author: Chris MacMackin
    !  Date: December 2016
    !
    ! Assigns all data in `rhs` to `this` except for the actual values
    ! of the field. For a uniform scalar field, this does not actually
    ! do anything.
    !
    class(uniform_scalar_field), intent(inout) :: this
    class(abstract_field), intent(in) :: rhs
      !! The field whose metadata (domain, resolution, etc) is to be
      !! copied
    logical, optional, intent(in) :: alloc
      !! If present and false, do not allocate the array of `this`.
  end subroutine uniform_scalar_assign_meta_data

  pure subroutine uniform_scalar_allocate_scalar(this, new_field)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Allocates an abstract [[scalar_field(type)]] to have a type
    ! compatible for operations on this field and to be the type
    ! returned by this field's methods which produce scalars.
    !
    class(uniform_scalar_field), intent(in)         :: this
    class(scalar_field), allocatable, intent(inout) :: new_field
      !! A field which, upon return, is allocated to be of the same
      !! concrete type as scalar fields produced by `this`.
    if (allocated(new_field)) deallocate(new_field)
    allocate(uniform_scalar_field :: new_field)
  end subroutine uniform_scalar_allocate_scalar

  pure subroutine uniform_scalar_allocate_vector(this, new_field)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Allocates an abstract [[vector_field(type)]] to have a type
    ! compatible for operations on this field and to be the type
    ! returned by this field's methods which produce vectors.
    !
    class(uniform_scalar_field), intent(in)         :: this
    class(vector_field), allocatable, intent(inout) :: new_field
      !! A field which, upon return, is allocated to be of the same
      !! concrete type as vector fields produced by `this`.
    if (allocated(new_field)) deallocate(new_field)
    allocate(uniform_vector_field :: new_field)
  end subroutine uniform_scalar_allocate_vector

  pure function uniform_scalar_id_to_pos(this, id) result(pos)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Returns the coordinates for a particular data-value held in the
    ! field, based on its ID number. As this field is uniform
    ! throughout space, the result has little meaning.
    !
    class(uniform_scalar_field), intent(in) :: this
    integer, intent(in)                     :: id
      !! The ID number for some location in the field
    real(r8), dimension(:), allocatable :: pos
      !! The coordinates for this location in the field
    pos = [0.0_r8]
  end function uniform_scalar_id_to_pos
  
  subroutine uniform_scalar_write_hdf(this, hdf_id, dataset_name, &
                                            error)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Writes the contents of the uniform field to a dataset in an HDF
    ! file. The dataset will have an attribute specifying the name of
    ! the field type and an attribute with a zero value, indicating
    ! that this is not a vector field.
    !
    class(uniform_scalar_field), intent(in) :: this
    integer(hid_t), intent(in)              :: hdf_id
      !! The identifier for the HDF file/group in which the field
      !! data is to be written.
    character(len=*), intent(in)            :: dataset_name
      !! The name of the dataset to be created in the HDF file
      !! containing this field's data.
    integer, intent(out)                    :: error
      !! An error code which, upon succesful completion of the
      !! routine, is 0. Otherwise, contains the error code returned
      !! by the HDF library.
    error = 0
    call h5ltmake_dataset_double_f(hdf_id, dataset_name, 1, [1_hsize_t], &
                                   [this%field_data], error)
    if (error /= 0) return
    call h5ltset_attribute_string_f(hdf_id, dataset_name, hdf_field_type_attr, &
                                    hdf_scalar_name, error)
    if (error /= 0) return
    call h5ltset_attribute_int_f(hdf_id, dataset_name, hdf_vector_attr, [0], &
                                 1_size_t, error)
  end subroutine uniform_scalar_write_hdf

  pure function uniform_scalar_grid_spacing(this) result(grid)
    !* Author: Chris MacMackin
    !  Date: December 2016
    !
    ! Returns a field containing the grid spacing at each point. For a
    ! uniform field, this is not a meaningful concept
    !
    class(uniform_scalar_field), intent(in) :: this
    class(vector_field), allocatable        :: grid
      !! A field where the values indicate the grid spacing that
      !! point. Each vector dimension representes the spacing of the
      !! grid in that direction.
    type(uniform_vector_field), allocatable :: local
    allocate(local)
    allocate(local%field_data(grid_spacing_size))
    local%vector_dims = grid_spacing_size
    local%field_data = huge(local%field_data)
  end function uniform_scalar_grid_spacing

  pure function uniform_scalar_get_element(this,element) result(val)
    !* Author: Chris MacMackin
    !  Date: January 2017
    !
    ! Returns an element of the field corresponding to the provided ID
    ! number. For a uniform field, this will always be the same
    ! regardless of the ID number.
    !
    class(uniform_scalar_field), intent(in) :: this
    integer, intent(in) :: element
      !! The ID number of the field element to be returned
    real(r8) :: val
      !! The value of the field corresponding to the specified ID
    val = this%field_data
  end function uniform_scalar_get_element

  pure subroutine uniform_scalar_set_element(this,element,val)
    !* Author: Chris MacMackin
    !  Date: January 2017
    !
    ! Sets the element of the field corresponding to the given ID to
    ! the given value. For a uniform field, this ammounts to setting
    ! the field value everywhere.
    !
    class(uniform_scalar_field), intent(inout) :: this
    integer, intent(in) :: element
      !! The ID number of the field element to be set
    real(r8), intent(in) :: val
      !! The new value the field element is to be set to
    this%field_data = val
  end subroutine uniform_scalar_set_element
  
  pure function uniform_scalar_get_bound(this,boundary,depth) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Returns a field containing the specified boundary
    ! information. For a uniform field this just means that a copy is
    ! returned.
    !
    class(uniform_scalar_field), intent(in) :: this
    integer, intent(in) :: boundary
      !! Specifies which boundary is to be returned. The boundary
      !! will be the one normal to dimension of number
      !! `abs(boundary)`. If the argument is negative, then the
      !! lower boundary is returned. If positive, then the upper
      !! boundary is returned.
    integer, intent(in) :: depth
      !! The number of layers of data-points to return at the
      !! specified boundary.
    class(scalar_field), allocatable :: res
      !! A field, of the same type as `this` and with the same
      !! resolution, number of dimensions etc., but containing only
      !! the points within the specified number of layers of cells
      !! adjecent to the specified boundary.
    allocate(res, source=this)
  end function uniform_scalar_get_bound


  !=====================================================================
  ! Vector Field Methods
  !=====================================================================

  
  function uniform_vector_constructor(val) result(this)
    !* Author: Chris MacMackin
    !  Date: October 2016
    !
    ! Creates a new vector field with a uniform value across all of space.
    !
    real(r8), dimension(:) :: val
      !! The value of the field
    type(uniform_vector_field), allocatable :: this
      !! A scalar field initiated based on the arguments to this function.
    integer :: i
    this%vector_dims = size(val)
    this%field_data = val
  end function uniform_vector_constructor

  elemental function uniform_vector_vector_dimensions(this) result(dims)
    !* Author: Chris MacMackin
    !  Date: October 2016
    !
    ! Returns the number of dimensions/components are in the vectors
    ! of this field.
    !
    class(uniform_vector_field), intent(in) :: this
    integer :: dims !! Number of vector components
    dims = this%vector_dims
  end function uniform_vector_vector_dimensions

  pure function uniform_vector_elements(this) result(elements)
    !* Author: Chris MacMackin
    !  Date: October 2016
    !
    ! Gives the number of individual data points present in the
    ! field. For a uniform field this is 1.
    !
    class(uniform_vector_field), intent(in) :: this
    integer :: elements
    elements = 1
  end function uniform_vector_elements

  elemental function uniform_vector_dimensions(this) result(dims)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Gives the number of dimensions of the field. For a uniform
    ! field, this is effectively 0.
    !
    class(uniform_vector_field), intent(in) :: this
    integer :: dims
    dims = 0
  end function uniform_vector_dimensions
  
  pure function uniform_vector_domain(this) result(domain)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Returns an array specifying the upper and lower bounds of the
    ! data in this field. As this field is uniform across all of
    ! space, the concept of limits is not meaningful.
    !
    class(uniform_vector_field), intent(in) :: this
    real(r8), dimension(:,:), allocatable :: domain
    allocate(domain(1,2))
    domain(1,1) = 0.0_r8
    domain(1,2) = 0.0_r8
  end function uniform_vector_domain

  pure function uniform_vector_raw_size(this,exclude_lower_bound, &
                                        exclude_upper_bound) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Compute how many elements are in the raw representation of this
    ! field. This would be the number of data points, adjusted based on
    ! how boundary conditions are accounted for.
    !
    class(uniform_vector_field), intent(in)       :: this
    integer, dimension(:), optional, intent(in) :: exclude_lower_bound
      !! Specifies how many layers of data points should be excluded
      !! from the result at the lower boundary for each dimension. The
      !! number in element `n` of the array indicates how many layers
      !! of cells at the lower boundary normal to dimension `n` will
      !! be ignored. For a uniform field, this argument has no effect.
    integer, dimension(:), optional, intent(in) :: exclude_upper_bound
      !! Specifies how many layers of data points should be excluded
      !! from the result at the upper boundary for each dimension. The
      !! number in element `n` of the array indicates how many layers
      !! of cells at the upper boundary normal to dimension `n` will
      !! be ignored. For a uniform field, this argument has no effect.
    integer :: res
    integer, dimension(:,:), allocatable :: slices
    integer :: i
    if (allocated(this%field_data)) then
      res = this%vector_dims
    else
      res = 0
    end if
  end function uniform_vector_raw_size
  
  pure function uniform_vector_raw(this,exclude_lower_bound, &
                                   exclude_upper_bound) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Returns a representation of this field in the form of a 1D uniform
    ! of real numbers this allows, e.g. for manipulation by solver
    ! routines.
    !
    ! @BUG The returned value has length `this%raw_size()`, but
    ! a bug in gfortran 4.8 (fixed by version 5) caused the compiler
    ! to segfault if it was declared as such. As a workaround, it is
    ! allocatable isntead.
    !
    class(uniform_vector_field), intent(in) :: this
    integer, dimension(:), optional, intent(in) :: exclude_lower_bound
      !! Specifies how many layers of data points should be excluded
      !! from the result at the lower boundary for each dimension. The
      !! number in element `n` of the array indicates how many layers
      !! of cells at the lower boundary normal to dimension `n` will
      !! be ignored. For a uniform field, this argument has no effect.
    integer, dimension(:), optional, intent(in) :: exclude_upper_bound
      !! Specifies how many layers of data points should be excluded
      !! from the result at the upper boundary for each dimension. The
      !! number in element `n` of the array indicates how many layers
      !! of cells at the upper boundary normal to dimension `n` will
      !! be ignored. For a uniform field, this argument has no effect.
    real(r8), dimension(:), allocatable :: res
      !! Uniform containing data needed to describe field
    if (allocated(this%field_data)) then
      res = this%field_data
    else
      allocate(res(0))
    end if
  end function uniform_vector_raw

#:if defined('DEBUG')
#undef pure 
#undef elemental 
#:endif

  pure function uniform_vector_get_value(this) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Returns the vector value of the field. This is useful for
    ! implementing the interactions of other field types with the
    ! uniform field.
    !
    class(uniform_vector_field), intent(in) :: this
    real(r8), dimension(:), allocatable :: res
    res = this%field_data
  end function uniform_vector_get_value

#:if defined('DEBUG')
#define pure 
#define elemental 
#:endif
  
  pure subroutine uniform_vector_set_from_raw(this,raw,provide_lower_bound, &
                                            provide_upper_bound)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Assigns raw data, such as that produced by 
    ! [[uniform_scalar_field:raw]], to the field. The routine will
    ! stop with an error if the uniform is the wrong length for a field
    ! of this resolution and with these boundary conditions. 
    !
    class(uniform_vector_field), intent(inout) :: this
    real(r8), dimension(:), intent(in) :: raw
      !! The raw data to be stored in this uniform.
    integer, dimension(:), optional, intent(in) :: provide_lower_bound
      !! Specifies how many layers of data points are excluded from
      !! the raw data at the lower boundary for each dimension. The
      !! number in element `n` of the array indicates how many layers
      !! of cells at the lower boundary normal to dimension `n` are
      !! missed. For a uniform field, this argument has no effect.
    integer, dimension(:), optional, intent(in) :: provide_upper_bound
      !! Specifies how many layers of data points are excluded
      !! from the raw data at the upper boundary for each
      !! dimension. The number in element `n` of the array indicates
      !! how many layers of cells at the upper boundary normal to
      !! dimension `n` are missed. Defaults to 0 for all.
    call check_set_from_raw(this,raw,provide_lower_bound,provide_upper_bound)
    this%field_data = raw
  end subroutine uniform_vector_set_from_raw

  pure function uniform_vector_resolution(this) result(resolution)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Returns an array specifying the number of data-points in each
    ! dimension. This is not a particularly meaningful quantity for a
    ! uniform array, as only a single data-point is needed.
    !
    class(uniform_vector_field), intent(in)  :: this
    integer, dimension(:), allocatable :: resolution
      !! Array specifying the number of data points in each dimension.
    allocate(resolution(0))
  end function uniform_vector_resolution

  pure function uniform_vector_vf_m_sf(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \(\vec{\rm field} \times {\rm field}\)
    !
    class(uniform_vector_field), intent(in) :: this
    class(scalar_field), intent(in) :: rhs
    class(vector_field), allocatable :: res !! The restult of this operation
    call rhs%allocate_vector_field(res)
    res = this%field_data * rhs
  end function uniform_vector_vf_m_sf

  pure function uniform_vector_r_m_vf(lhs,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \({\rm real} \times \vec{\rm field}\)
    !
    real(r8), intent(in) :: lhs
    class(uniform_vector_field), intent(in) :: rhs
    class(vector_field), allocatable :: res !! The result of this operation
    type(uniform_vector_field), allocatable :: local
    allocate(local)
    local%vector_dims = rhs%vector_dims
    local%field_data = lhs * rhs%field_data
    call move_alloc(local,res)
  end function uniform_vector_r_m_vf

  pure function uniform_vector_vf_m_r(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \(\vec{\rm field} \times {\rm real}\)
    !
    class(uniform_vector_field), intent(in) :: this
    real(r8), intent(in) :: rhs
    class(vector_field), allocatable :: res !! The result of this operation
    type(uniform_vector_field), allocatable :: local
    allocate(local)
    local%vector_dims = this%vector_dims
    local%field_data = this%field_data * rhs
    call move_alloc(local,res)
  end function uniform_vector_vf_m_r
  
  pure function uniform_vector_vf_d_sf(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \(\vec{\rm field} / {\rm field}\)
    !
    class(uniform_vector_field), intent(in) :: this
    class(scalar_field), intent(in) :: rhs
    class(vector_field), allocatable :: res !! The restult of this operation
    call rhs%allocate_vector_field(res)
    res = this%field_data / rhs
  end function uniform_vector_vf_d_sf

  pure function uniform_vector_vf_d_r(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \(\vec{\rm field} / {\rm real}\)
    !
    class(uniform_vector_field), intent(in) :: this
    real(r8), intent(in) :: rhs
    class(vector_field), allocatable :: res !! The result of this operation
    type(uniform_vector_field), allocatable :: local
    allocate(local)
    local%vector_dims = this%vector_dims
    local%field_data = this%field_data / rhs
    call move_alloc(local,res)
  end function uniform_vector_vf_d_r
  
  pure function uniform_vector_vf_s_vf(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \(\vec{\rm field} - \vec{\rm field}\)
    !
    class(uniform_vector_field), intent(in) :: this
    class(vector_field), intent(in) :: rhs
    class(vector_field), allocatable :: res !! The restult of this operation
    allocate(res, mold=rhs)
    res = this%field_data - rhs
  end function uniform_vector_vf_s_vf

  pure function uniform_vector_r_s_vf(lhs,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \(\vec{\rm real} - \vec{\rm field}\)
    !
    real(r8), dimension(:), intent(in) :: lhs
    class(uniform_vector_field), intent(in) :: rhs
    class(vector_field), allocatable :: res !! The result of this operation
    type(uniform_vector_field), allocatable :: local
    integer :: min_dims, max_dims
    allocate(local)
    min_dims = min(rhs%vector_dims, size(lhs))
    max_dims = max(rhs%vector_dims, size(lhs))
    allocate(local%field_data(max_dims))
    local%vector_dims = max_dims
    local%field_data(:min_dims) = lhs(:min_dims) - rhs%field_data(:min_dims)
    if (rhs%vector_dims > size(lhs)) then
      local%field_data(min_dims+1:) = -rhs%field_data(min_dims+1:)
    else
      local%field_data(min_dims+1:) = lhs(min_dims+1:)
    end if
    call move_alloc(local,res)
  end function uniform_vector_r_s_vf

  pure function uniform_vector_vf_s_r(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \(\vec{\rm field} - \vec{\rm real}\)
    !
    class(uniform_vector_field), intent(in) :: this
    real(r8), dimension(:), intent(in) :: rhs
    class(vector_field), allocatable :: res !! The result of this operation
    type(uniform_vector_field), allocatable :: local
    integer :: min_dims, max_dims
    allocate(local)
    min_dims = min(this%vector_dims, size(rhs))
    max_dims = max(this%vector_dims, size(rhs))
    allocate(local%field_data(max_dims))
    local%vector_dims = max_dims
    local%field_data(:min_dims) = this%field_data(:min_dims) - rhs(:min_dims)
    if (this%vector_dims > size(rhs)) then
      local%field_data(min_dims+1:) = this%field_data(min_dims+1:)
    else
      local%field_data(min_dims+1:) = -rhs(min_dims+1:)
    end if
    call move_alloc(local,res)
  end function uniform_vector_vf_s_r
  
  pure function uniform_vector_vf_a_vf(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \(\vec{\rm field} - \vec{\rm field}\)
    !
    class(uniform_vector_field), intent(in) :: this
    class(vector_field), intent(in) :: rhs
    class(vector_field), allocatable :: res !! The restult of this operation
    integer :: max_dims, min_dims
    allocate(res, mold=this)
    res = this%field_data + rhs
  end function uniform_vector_vf_a_vf

  pure function uniform_vector_r_a_vf(lhs,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \(\vec{\rm real} + \vec{\rm field}\)
    !
    real(r8), dimension(:), intent(in) :: lhs
    class(uniform_vector_field), intent(in) :: rhs
    class(vector_field), allocatable :: res !! The result of this operation
    type(uniform_vector_field), allocatable :: local
    integer :: min_dims, max_dims
    allocate(local)
    min_dims = min(rhs%vector_dims, size(lhs))
    max_dims = max(rhs%vector_dims, size(lhs))
    allocate(local%field_data(max_dims))
    local%vector_dims = max_dims
    local%field_data(:min_dims) = lhs(:min_dims) + rhs%field_data(:min_dims)
    if (rhs%vector_dims > size(lhs)) then
      local%field_data(min_dims+1:) = rhs%field_data(min_dims+1:)
    else
      local%field_data(min_dims+1:) = lhs(min_dims+1:)
    end if
    call move_alloc(local,res)
  end function uniform_vector_r_a_vf

  pure function uniform_vector_vf_a_r(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \(\vec{\rm field} + \vec{\rm real}\)
    !
    class(uniform_vector_field), intent(in) :: this
    real(r8), dimension(:), intent(in) :: rhs
    class(vector_field), allocatable :: res !! The result of this operation
    type(uniform_vector_field), allocatable :: local
    integer :: min_dims, max_dims
    allocate(local)
    min_dims = min(this%vector_dims, size(rhs))
    max_dims = max(this%vector_dims, size(rhs))
    allocate(local%field_data(max_dims))
    local%vector_dims = max_dims
    local%field_data(:min_dims) = this%field_data(:min_dims) + rhs(:min_dims)
    if (this%vector_dims > size(rhs)) then
      local%field_data(min_dims+1:) = this%field_data(min_dims+1:)
    else
      local%field_data(min_dims+1:) = rhs(min_dims+1:)
    end if
    call move_alloc(local,res)
  end function uniform_vector_vf_a_r

  elemental subroutine uniform_vector_assign(this,rhs)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \(\vec{\rm field} = \vec{\rm field}\)
    !
    class(uniform_vector_field), intent(inout) :: this
    class(vector_field), intent(in) :: rhs
    select type(rhs)
    class is(uniform_vector_field)
      if (allocated(rhs%field_data)) then
        this%field_data = rhs%field_data
        this%vector_dims = rhs%vector_dims
      else
        this%vector_dims = 0
      end if
    class default
      error stop ('Assigning incompatible type to uniform_vector_field')
    end select
  end subroutine uniform_vector_assign

  pure function uniform_vector_norm(this) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \(\lVert \vec{\rm field} \rVert\)
    !
    class(uniform_vector_field), intent(in) :: this
    class(scalar_field), allocatable :: res
    type(uniform_scalar_field), allocatable :: local
    allocate(local)
    local%field_data = sqrt(sum(this%field_data**2))
    call move_alloc(local, res)
  end function uniform_vector_norm

  pure function uniform_vector_component(this,comp) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Returns a calar field containing specified component of the vector field
    !
    class(uniform_vector_field), intent(in) :: this
    integer, intent(in) :: comp
    class(scalar_field), allocatable :: res
    type(uniform_scalar_field), allocatable :: local
    allocate(local)
    if (comp <= this%vector_dims) then
      local%field_data = this%field_data(comp)
    else
      local%field_data = 0.0_r8
    end if
    call move_alloc(local, res)
  end function uniform_vector_component

  function uniform_vector_d_dx(this, dir, order) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \(\frac{\partial^{\rm order}}{\partial x_{\rm dir}^{\rm order}}\vec{\rm field}\)
    !
    class(uniform_vector_field), intent(in) :: this
    integer, intent(in) :: dir !! Direction in which to differentiate
    integer, optional, intent(in) :: order !! Order of the derivative, default = 1
    class(vector_field), allocatable :: res
    type(uniform_vector_field), allocatable :: local
    allocate(local)
    local%vector_dims = 1
    local%field_data = [0.0_r8]
    call move_alloc(local, res)
  end function uniform_vector_d_dx

  function uniform_vector_component_d_dx(this, dir, component, order) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \(\frac{\partial^{\rm order}}{\partial x_{\rm dir}^{\rm order}}{\rm field_{component}}\)
    !
    class(uniform_vector_field), intent(in) :: this
    integer, intent(in) :: dir !! Direction in which to differentiate
    integer, intent(in) :: component !! Which component of the vector is being differentiated
    integer, optional, intent(in) :: order !! Order of the derivative, default = 1
    class(scalar_field), allocatable :: res
    type(uniform_scalar_field), allocatable :: local
    allocate(local)
    local%field_data = 0.0_r8
    call move_alloc(local, res)
  end function uniform_vector_component_d_dx

  function uniform_vector_laplacian(this) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \(\nabla^2 \vec{\rm field}\)
    !
    class(uniform_vector_field), intent(in) :: this
    class(vector_field), allocatable :: res !! The result of this operation
    type(uniform_vector_field), allocatable :: local
    allocate(local)
    local%vector_dims = 1
    local%field_data = [0.0_r8]
    call move_alloc(local, res)
  end function uniform_vector_laplacian
  
  function uniform_vector_divergence(this) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \(\nabla\cdot \vec{\rm field}\)
    !
    class(uniform_vector_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(uniform_scalar_field), allocatable :: local
    allocate(local)
    local%field_data = 0.0_r8
  end function uniform_vector_divergence
  
  function uniform_vector_curl(this) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \(\nabla\times \vec{\rm field}\)
    !
    class(uniform_vector_field), intent(in) :: this
    class(vector_field), allocatable :: res !! The result of this operation
    type(uniform_vector_field), allocatable :: local
    allocate(local)
    local%vector_dims = 1
    local%field_data = [0.0_r8]
    call move_alloc(local, res)
  end function uniform_vector_curl

  pure function uniform_vector_vf_cross_vf(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \(\vec{\rm field} \times \vec{\rm field}\)
    !
    ! The returned vector will always consist of three components,
    ! regardless of the number of components in the arguments.
    !
    class(uniform_vector_field), intent(in) :: this
    class(vector_field), intent(in) :: rhs
    class(vector_field), allocatable :: res !! The restult of this operation
    allocate(res, mold=this)
    res = this%field_data .cross. rhs
  end function uniform_vector_vf_cross_vf

  pure function uniform_vector_vf_cross_vr(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \(\vec{\rm field} \times \vec{\rm real}\)
    !
    ! The returned vector will always consist of three components,
    ! regardless of the number of components in the arguments.
    !
    class(uniform_vector_field), intent(in) :: this
    real(r8), dimension(:), intent(in) :: rhs
    class(vector_field), allocatable :: res !! The result of this operation
    type(uniform_vector_field), allocatable :: local
    real(r8), dimension(3) :: vec1, vec2
    integer :: dims1, dims2
    allocate(local)
    allocate(local%field_data(3))
    local%vector_dims = 3
    dims1 = min(3,this%vector_dims)
    dims2 = min(3,size(rhs))
    vec1 = 0.0_r8
    vec2 = 0.0_r8
    vec1(:dims1) = this%field_data(:dims1)
    vec2(:dims2) = rhs(:dims2)
    local%field_data = [vec1(2)*vec2(3) - vec2(2)*vec1(3), &
                        vec1(3)*vec2(1) - vec2(3)*vec1(1), &
                        vec1(1)*vec2(2) - vec2(1)*vec1(2)]
    call move_alloc(local, res)
  end function uniform_vector_vf_cross_vr

  pure function uniform_vector_vr_cross_vf(lhs,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \(\vec{\rm real} \times \vec{\rm field}\)
    !
    ! The returned vector will always consist of three components,
    ! regardless of the number of components in the arguments.
    !
    real(r8), dimension(:), intent(in) :: lhs
    class(uniform_vector_field), intent(in) :: rhs
    class(vector_field), allocatable :: res !! The result of this operation
    type(uniform_vector_field), allocatable :: local
    real(r8), dimension(3) :: vec1, vec2
    integer :: dims1, dims2
    allocate(local)
    allocate(local%field_data(3))
    local%vector_dims = 3
    dims1 = min(3,size(lhs))
    dims2 = min(3,rhs%vector_dims)
    vec1 = 0.0_r8
    vec2 = 0.0_r8
    vec1(:dims1) = lhs(:dims1)
    vec2(:dims2) = rhs%field_data(:dims2)
    local%field_data = [vec1(2)*vec2(3) - vec2(2)*vec1(3), &
                        vec1(3)*vec2(1) - vec2(3)*vec1(1), &
                        vec1(1)*vec2(2) - vec2(1)*vec1(2)]
    call move_alloc(local, res)
  end function uniform_vector_vr_cross_vf
  
  pure function uniform_vector_vf_dot_vf(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \(\vec{\rm field} \cdot \vec{\rm field}\)
    !
    class(uniform_vector_field), intent(in) :: this
    class(vector_field), intent(in) :: rhs
    class(scalar_field), allocatable :: res !! The restult of this operation
    call rhs%allocate_scalar_field(res)
    res = this%field_data .dot. rhs
  end function uniform_vector_vf_dot_vf

  pure function uniform_vector_vf_dot_vr(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \(\vec{\rm field} \cdot \vec{\rm real}\)
    !
    class(uniform_vector_field), intent(in) :: this
    real(r8), dimension(:), intent(in) :: rhs
    class(scalar_field), allocatable :: res !! The restult of this operation
    type(uniform_scalar_field), allocatable :: local
    integer :: dims
    allocate(local)
    dims = min(this%vector_dims,size(rhs))
    local%field_data = sum(this%field_data(:dims)*rhs(:dims))
    call move_alloc(local, res)
  end function uniform_vector_vf_dot_vr

  pure function uniform_vector_vr_dot_vf(lhs,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! \(\vec{\rm real} \cdot \vec{\rm field}\)
    !
    real(r8), dimension(:), intent(in) :: lhs
    class(uniform_vector_field), intent(in) :: rhs
    class(scalar_field), allocatable :: res !! The restult of this operation
    type(uniform_scalar_field), allocatable :: local
    integer :: dims
    allocate(local)
    dims = min(size(lhs),rhs%vector_dims)
    local%field_data = sum(lhs(:dims) * rhs%field_data(:dims))
    call move_alloc(local, res)
  end function uniform_vector_vr_dot_vf

  pure logical function uniform_vector_is_equal(this,rhs) result(iseq)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Evaluates whether two vector fields are equal within a tolerance,
    ! specified by [[set_tol]].
    !
    class(uniform_vector_field), intent(in) :: this
    class(vector_field), intent(in) :: rhs
    integer :: dims
    real(r8) :: normalization
    iseq = .true.
    select type(rhs)
    class is(uniform_vector_field)
      if (this%vector_dims > rhs%vector_dims) then
        dims = rhs%vector_dims
        iseq = all(abs(this%field_data(dims:)) < get_tol())
      else if (this%vector_dims < rhs%vector_dims) then
        dims = this%vector_dims
        iseq = all(abs(rhs%field_data(dims:)) < get_tol())
      else
        dims = this%vector_dims
      end if
      if (.not. iseq) return
      iseq = (norm2(this%field_data(:dims) - rhs%field_data(:dims))/ &
              normalization < get_tol()) .or. &
             all(is_nan(this%field_data(1:dims)) .eqv. &
                 is_nan(rhs%field_data(1:dims)))
        ! FIXME: This handling of NaNs will claim that things are
        ! equal when they aren't, just because they have a NAN in the
        ! same location.
    class default
      iseq = (rhs == this)
    end select
  end function uniform_vector_is_equal

  pure subroutine uniform_vector_assign_meta_data(this, rhs, alloc)
    !* Author: Chris MacMackin
    !  Date: December 2016
    !
    ! Assigns all data in `rhs` to `this` except for the actual values
    ! of the field.
    !
    class(uniform_vector_field), intent(inout) :: this
    class(abstract_field), intent(in) :: rhs
      !! The field whose metadata (domain, resolution, etc) is to be
      !! copied
    logical, optional, intent(in) :: alloc
      !! If present and false, do not allocate the array of `this`.
    select type(rhs)
    class is(uniform_scalar_field)
      this%vector_dims = 1
    class is(uniform_vector_field)
      this%vector_dims = rhs%vector_dims
    end select
  end subroutine uniform_vector_assign_meta_data

  pure subroutine uniform_vector_allocate_scalar(this, new_field)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Allocates an abstract [[scalar_field(type)]] to have a type
    ! compatible for operations on this field and to be the type
    ! returned by this field's methods which produce scalars.
    !
    class(uniform_vector_field), intent(in)         :: this
    class(scalar_field), allocatable, intent(inout) :: new_field
      !! A field which, upon return, is allocated to be of the same
      !! concrete type as scalar fields produced by `this`.
    if (allocated(new_field)) deallocate(new_field)
    allocate(uniform_scalar_field :: new_field)
  end subroutine uniform_vector_allocate_scalar

  pure subroutine uniform_vector_allocate_vector(this, new_field)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Allocates an abstract [[vector_field(type)]] to have a type
    ! compatible for operations on this field and to be the type
    ! returned by this field's methods which produce vectors.
    !
    class(uniform_vector_field), intent(in)         :: this
    class(vector_field), allocatable, intent(inout) :: new_field
      !! A field which, upon return, is allocated to be of the same
      !! concrete type as vector fields produced by `this`.
    if (allocated(new_field)) deallocate(new_field)
    allocate(uniform_vector_field :: new_field)
  end subroutine uniform_vector_allocate_vector

  pure function uniform_vector_id_to_pos(this, id) result(pos)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Returns the coordinates for a particular data-value held in the
    ! field, based on its ID number. As this field is uniform
    ! throughout space, the result has little meaning.
    !
    class(uniform_vector_field), intent(in) :: this
    integer, intent(in)                     :: id
      !! The ID number for some location in the field
    real(r8), dimension(:), allocatable :: pos
      !! The coordinates for this location in the field
    pos = [0.0_r8]
  end function uniform_vector_id_to_pos
  
  subroutine uniform_vector_write_hdf(this, hdf_id, dataset_name, &
                                      error)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Writes the contents of the uniform field to a dataset in an HDF
    ! file. The dataset will have an attribute specifying the name of
    ! the field type and an attribute with a non-zero value,
    ! indicating that this is a vector field.
    !
    class(uniform_vector_field), intent(in) :: this
    integer(hid_t), intent(in)              :: hdf_id
      !! The identifier for the HDF file/group in which the field
      !! data is to be written.
    character(len=*), intent(in)            :: dataset_name
      !! The name of the dataset to be created in the HDF file
      !! containing this field's data.
    integer, intent(out)                    :: error
      !! An error code which, upon succesful completion of the
      !! routine, is 0. Otherwise, contains the error code returned
      !! by the HDF library.
    error = 0
    call h5ltmake_dataset_double_f(hdf_id, dataset_name, 1, &
                                   [int(this%vector_dims,hsize_t)], &
                                   [this%field_data], error)
    if (error /= 0) return
    call h5ltset_attribute_string_f(hdf_id, dataset_name, &
                                    hdf_field_type_attr, hdf_vector_name, error)
    if (error /= 0) return
    call h5ltset_attribute_int_f(hdf_id, dataset_name, hdf_vector_attr, [1], &
                                 1_size_t, error)
  end subroutine uniform_vector_write_hdf
  
  pure function uniform_vector_grid_spacing(this) result(grid)
    !* Author: Chris MacMackin
    !  Date: December 2016
    !
    ! Returns a field containing the grid spacing at each point. For a
    ! uniform field, this is not a meaningful concept
    !
    class(uniform_vector_field), intent(in) :: this
    class(vector_field), allocatable        :: grid
    !! A field where the values indicate the grid spacing that
    !! point. Each vector dimension representes the spacing of the
    !! grid in that direction.
    type(uniform_vector_field), allocatable :: local
    allocate(local)
    allocate(local%field_data(grid_spacing_size))
    local%vector_dims = grid_spacing_size
    local%field_data = huge(local%field_data)
  end function uniform_vector_grid_spacing

  pure function uniform_vector_get_element_vec(this,element) result(val)
    !* Author: Chris MacMackin
    !  Date: January 2017
    !
    ! Returns a vector of the field corresponding to the provided ID
    ! number. For a uniform field, this will always be the same
    ! regardless of the ID number.
    !
    class(uniform_vector_field), intent(in) :: this
    integer, intent(in) :: element
      !! The ID number of the field element to be returned
    real(r8), allocatable, dimension(:) :: val
      !! The vector in the field corresponding to the specified ID
    val = this%field_data
  end function uniform_vector_get_element_vec
  
  pure function uniform_vector_get_element_comp(this,element,component) result(val)
    !* Author: Chris MacMackin
    !  Date: January 2017
    !
    ! Returns a component of the vector of the field corresponding to
    ! the provided ID number. For a uniform field, this will always be
    ! the same regardless of the ID number.
    !
    class(uniform_vector_field), intent(in) :: this
    integer, intent(in) :: element
      !! The ID number of the field element to be returned
    integer, intent(in) :: component
      !! The number of the vector component to be returned
    real(r8) :: val
      !! The vector component in the field corresponding to the 
      !! specified ID
    val = this%field_data(component)
  end function uniform_vector_get_element_comp

  pure subroutine uniform_vector_set_element_vec(this,element,val)
    !* Author: Chris MacMackin
    !  Date: January 2017
    !
    ! Sets the element of the field corresponding to the given ID to
    ! the given vector value. For a uniform field, this ammounts to
    ! setting the value everywhere.
    !
    class(uniform_vector_field), intent(inout) :: this
    integer, intent(in) :: element
      !! The ID number of the field element to be set
    real(r8), dimension(:), intent(in) :: val
      !! The new vector value the field element is to be set to
    this%field_data(:) = val
  end subroutine uniform_vector_set_element_vec

  pure subroutine uniform_vector_set_element_comp(this,element,component,val)
    !* Author: Chris MacMackin
    !  Date: January 2017
    !
    ! Sets the element of the field corresponding to the given ID to
    ! the given vector value. For a uniform field, this ammounts to
    ! setting the component everywhere.
    !
    class(uniform_vector_field), intent(inout) :: this
    integer, intent(in) :: element
      !! The ID number of the field element to be set
    integer, intent(in) :: component
      !! The number of the vector component to be returned
    real(r8), intent(in) :: val
      !! The new value of the vector component in the field element
    this%field_data(component) = val
  end subroutine uniform_vector_set_element_comp

  pure function uniform_vector_get_bound(this,boundary,depth) result(res)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Returns a field containing the specified boundary
    ! information. For a uniform field this just means that a copy is
    ! returned.
    !
    class(uniform_vector_field), intent(in) :: this
    integer, intent(in) :: boundary
      !! Specifies which boundary is to be returned. The boundary
      !! will be the one normal to dimension of number
      !! `abs(boundary)`. If the argument is negative, then the
      !! lower boundary is returned. If positive, then the upper
      !! boundary is returned.
    integer, intent(in) :: depth
      !! The number of layers of data-points to return at the
      !! specified boundary.
    class(vector_field), allocatable :: res
      !! A field, of the same type as `this` and with the same
      !! resolution, number of dimensions etc., but containing only
      !! the points within the specified number of layers of cells
      !! adjecent to the specified boundary.
    allocate(res, source=this)
  end function uniform_vector_get_bound

  elemental subroutine uniform_vector_finalize(this)
    !* Author: Chris MacMackin
    !  Date: January 2017
    !
    ! Deallocates the field data for this object.
    !
    type(uniform_vector_field), intent(inout) :: this
    if (allocated(this%field_data)) deallocate(this%field_data)
  end subroutine uniform_vector_finalize

end module uniform_fields_mod
