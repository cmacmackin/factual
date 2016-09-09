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

module array_fields_mod
  !* Author: Chris MacMackin
  !  Date: September 2016
  !  License: GPLv3
  !
  ! Provides abstract extensions of the fileds in [[abstract_fields_mod]]
  ! which store field data in an array. These fields implement all of the
  ! non-calculus operations.
  !
  use iso_fortran_env, only: r8 => real64, stderr => error_unit
  use abstract_fields_mod, only: abstract_field, scalar_field, vector_field, &
                                 get_tol
  use utils_mod, only: is_nan, check_set_from_raw
  implicit none
  private

  type, extends(scalar_field), abstract, public :: array_scalar_field
    !* Author: Chris MacMackin
    !  Date: September 2016
    !
    ! An abstract data type representing a mathematical field of scalar
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
    real(r8), dimension(:,:), allocatable :: extent
      !! The start and end values of the domain
    real(r8), dimension(:), allocatable :: colloc_points
      !! The location of the data-points
    real(r8), dimension(:), allocatable :: field_data
      !! The value of the scalar field at the data-points
  contains
    private
    procedure, public :: raw_size => array_scalar_raw_size
      !! Provides the number of pieces of data needed to represent the
      !! field, i.e. the size of the array returned by `get_raw`.
    procedure, public :: raw => array_scalar_raw
      !! Returns array of data representing state of field. Can be
      !! useful for passing to nonlinear solvers
    procedure, public :: set_from_raw => array_scalar_set_from_raw
      !! Assigns raw data, such as that produced by 
      !! [[array_scalar_field:raw]], to the field
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
    procedure :: sin => array_scalar_sin
      !! \(\sin({\rm field})\)
    procedure :: cos => array_scalar_cos
      !! \(\cos({\rm field})\)
    procedure :: tan => array_scalar_tan
      !! \(\tan({\rm field})\)
    procedure :: asin => array_scalar_asin
      !! \(\sin^{-1}({\rm field})\)
    procedure :: acos => array_scalar_acos
      !! \(\cos^{-1}({\rm field})\)
    procedure :: atan => array_scalar_atan
      !! \(\tan^{-1}({\rm field})\)
    procedure :: sinh => array_scalar_sinh
      !! \(\sinh({\rm field})\)
    procedure :: cosh => array_scalar_cosh
      !! \(\cosh({\rm field})\)
    procedure :: tanh => array_scalar_tanh
      !! \(\tanh({\rm field})\)
    procedure :: asinh => array_scalar_asinh
      !! \(\sinh^{-1}({\rm field})\)
    procedure :: acosh => array_scalar_acosh
      !! \(\cosh^{-1}({\rm field})\)
    procedure :: atanh => array_scalar_atanh
      !! \(\tanh^{-1}({\rm field})\)
    procedure :: log => array_scalar_log
      !! \(\ln ({\rm field})\)
    procedure :: log10 => array_scalar_log10
      !! \(\log ({\rm field})\)
    procedure :: exp => array_scalar_exp
      !! \(\e^{\rm field})\)
    procedure :: abs => array_scalar_abs
      !! \(\|{\rm field}|\)
    procedure :: sqrt => array_scalar_sqrt
      !! \(\sqrt{\rm field}\)
    procedure :: minval => array_scalar_minval
      !! \(\min({\rm field})\)
    procedure :: maxval => array_scalar_maxval
      !! \(\max({\rm field})\)
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
    procedure(sf_compatible), deferred :: check_compatible
      !! Tests whether two fields are suitable for binary operations together
    procedure(sf_vf_result), public, deferred :: vector_result_type
      !! Returns an unallocated vector field object of the same type returned
      !! by the [array_scalar_field(type):gradient] method
  end type array_scalar_field

  abstract interface
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
    real(r8), dimension(:,:), allocatable :: extent
      !! The start and end values of the domain
    real(r8), dimension(:), allocatable :: colloc_points
      !! The location of the data-points
    real(r8), dimension(:,:), allocatable :: field_data
      !! The value of the vector field at the data-points
    integer :: numdims = 0
      !! The number of vector components
  contains
    private
    procedure, public :: raw_size => array_vector_raw_size
      !! Provides the number of pieces of data needed to represent the
      !! field, i.e. the size of the array returned by `get_raw`.
    procedure, public :: raw => array_vector_raw
      !! Returns array of data representing state of field. Can be
      !! useful for passing to nonlinear solvers.
    procedure, public :: set_from_raw => array_vector_set_from_raw
      !! Assigns raw data, such as that produced by 
      !! [[array_vector_field:raw]], to the field
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
    procedure(vf_compatible), deferred :: check_compatible
      !! Tests whether two fields are suitable for binary operations together
    procedure(vf_sf_result), public, deferred :: scalar_result_type
      !! Returns an unallocated scalar field object of the same type returned
      !! by the [array_vector_field(type):divergence], 
      !! [array_vector_field(type):dot_product], and other methods resulting
      !! in a scalar field
  end type array_vector_field

  abstract interface
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
  end interface

    
  abstract interface
    pure function vf_sf_result(this)
      import :: array_vector_field
      import :: array_scalar_field
      class(array_vector_field), intent(in) :: this
      class(array_scalar_field), allocatable :: vf_sf_result
        !! An unallocated field of the same type as returned by methods
        !! whose result is a scalar field.
    end function vf_sf_result

    pure function sf_vf_result(this)
      import :: array_vector_field
      import :: array_scalar_field
      class(array_scalar_field), intent(in) :: this
      class(array_vector_field), allocatable :: sf_vf_result
        !! An unallocated field of the same type as returned by methods
        !! whose result is a vector field.
    end function sf_vf_result
  end interface


contains


  !=====================================================================
  ! Scalar Field Methods
  !=====================================================================


  pure function array_scalar_raw_size(this,return_lower_bound, &
                                       return_upper_bound) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! Compute how many elements are in the raw representation of this
    ! field. This would be the number of data points, adjusted based on
    ! how boundary conditions are accounted for.
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
    integer :: res
    res = this%numpoints + 1
    if (present(return_lower_bound)) then
      if (.not. return_lower_bound(1)) res = res - 1
    end if
    if (present(return_upper_bound)) then
      if (.not. return_upper_bound(1)) res = res - 1
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
    logical, dimension(:), optional, intent(in) :: provide_lower_bound
      !! Specifies whether raw data contains values at the lower
      !! boundary, for each dimension, with the index of the element
      !! corresponding to the dimension. Defaults to all `.true.`.
    logical, dimension(:), optional, intent(in) :: provide_upper_bound
      !! Specifies whether raw data contains values at the upper
      !! boundary, for each dimension, with the index of the element
      !! corresponding to the dimension. Defaults to all `.true.`.
    integer :: lower, upper
    call check_set_from_raw(this,raw,provide_lower_bound,provide_upper_bound)
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
  end subroutine array_scalar_set_from_raw
  
  function array_scalar_sf_m_sf(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm field} \times {\rm field}\)
    !
    class(array_scalar_field), intent(in) :: this
    class(scalar_field), intent(in) :: rhs
    class(scalar_field), allocatable :: res !! The restult of this operation
    class(array_scalar_field), allocatable :: local
#ifdef DEBUG
    call this%check_compatible(rhs)
#endif
    allocate(local, mold=this)
    call local%assign_meta_data(this)
    select type(rhs)
    class is(array_scalar_field)
      local%field_data = this%field_data * rhs%field_data
    end select
    call move_alloc(local,res)
  end function array_scalar_sf_m_sf

  function array_scalar_sf_m_vf(this,rhs) result(res)
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
#ifdef DEBUG
    call this%check_compatible(rhs)
#endif
    select type(rhs)
    class is(array_vector_field)
      allocate(local, mold=rhs)
      call local%assign_meta_data(rhs)
      forall(i=1:this%numpoints+1) &
        local%field_data(i,:) = this%field_data(i) * rhs%field_data(i,:)
    end select
    call move_alloc(local,res)
  end function array_scalar_sf_m_vf

  function array_scalar_r_m_sf(lhs,rhs) result(res)
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

  function array_scalar_sf_m_r(this,rhs) result(res)
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
  
  function array_scalar_sf_d_sf(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm field} / {\rm field}\)
    !
    class(array_scalar_field), intent(in) :: this
    class(scalar_field), intent(in) :: rhs
    class(scalar_field), allocatable :: res !! The restult of this operation
    class(array_scalar_field), allocatable :: local
#ifdef DEBUG
    call this%check_compatible(rhs)
#endif
    allocate(local, mold=this)
    call local%assign_meta_data(this)
    select type(rhs)
    class is(array_scalar_field)
      local%field_data = this%field_data / rhs%field_data
    end select
    call move_alloc(local,res)
  end function array_scalar_sf_d_sf

  function array_scalar_r_d_sf(lhs,rhs) result(res)
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

  function array_scalar_sf_d_r(this,rhs) result(res)
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
  
  function array_scalar_sf_s_sf(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm field} - {\rm field}\)
    !
    class(array_scalar_field), intent(in) :: this
    class(scalar_field), intent(in) :: rhs
    class(scalar_field), allocatable :: res !! The restult of this operation
    class(array_scalar_field), allocatable :: local
#ifdef DEBUG
    call this%check_compatible(rhs)
#endif
    allocate(local, mold=this)
    call local%assign_meta_data(this)
    select type(rhs)
    class is(array_scalar_field)
      local%field_data = this%field_data - rhs%field_data
    end select
    call move_alloc(local,res)
  end function array_scalar_sf_s_sf

  function array_scalar_r_s_sf(lhs,rhs) result(res)
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

  function array_scalar_sf_s_r(this,rhs) result(res)
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
  
  function array_scalar_sf_a_sf(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm field} - {\rm field}\)
    !
    class(array_scalar_field), intent(in) :: this
    class(scalar_field), intent(in) :: rhs
    class(scalar_field), allocatable :: res !! The restult of this operation
    class(array_scalar_field), allocatable :: local
#ifdef DEBUG
    call this%check_compatible(rhs)
#endif
    allocate(local, mold=this)
    call local%assign_meta_data(this)
    select type(rhs)
    class is(array_scalar_field)
      local%field_data = this%field_data + rhs%field_data
    end select
    call move_alloc(local,res)
  end function array_scalar_sf_a_sf

  function array_scalar_r_a_sf(lhs,rhs) result(res)
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

  function array_scalar_sf_a_r(this,rhs) result(res)
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

  function array_scalar_sf_p_r(this,rhs) result(res)
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

  function array_scalar_sf_p_r4(this,rhs) result(res)
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

  function array_scalar_sf_p_i(this,rhs) result(res)
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

  pure function array_scalar_sin(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\sin({\rm field})\)
    !
    class(array_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    class(array_scalar_field), allocatable :: local
    allocate(local, mold=this)
    call local%assign_meta_data(this)
    if (allocated(this%field_data)) then
      local%field_data = sin(this%field_data)
    end if
    call move_alloc(local, res)
  end function array_scalar_sin

  pure function array_scalar_cos(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\cos({\rm field})\)
    !
    class(array_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    class(array_scalar_field), allocatable :: local
    allocate(local, mold=this)
    call local%assign_meta_data(this)
    if (allocated(this%field_data)) then
      local%field_data = cos(this%field_data)
    end if
    call move_alloc(local, res)
  end function array_scalar_cos

  pure function array_scalar_tan(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\tan({\rm field})\)
    !
    class(array_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    class(array_scalar_field), allocatable :: local
    allocate(local, mold=this)
    call local%assign_meta_data(this)
    if (allocated(this%field_data)) then
      local%field_data = tan(this%field_data)
    end if
    call move_alloc(local, res)
  end function array_scalar_tan

  pure function array_scalar_asin(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\sin^{-1}({\rm field})\)
    !
    class(array_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    class(array_scalar_field), allocatable :: local
    allocate(local, mold=this)
    call local%assign_meta_data(this)
    if (allocated(this%field_data)) then
      local%field_data = asin(this%field_data)
    end if
    call move_alloc(local, res)
  end function array_scalar_asin

  pure function array_scalar_acos(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\cos^{-1}({\rm field})\)
    !
    class(array_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    class(array_scalar_field), allocatable :: local
    allocate(local, mold=this)
    call local%assign_meta_data(this)
    if (allocated(this%field_data)) then
      local%field_data = acos(this%field_data)
    end if
    call move_alloc(local, res)
  end function array_scalar_acos

  pure function array_scalar_atan(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\tan^{-1}({\rm field})\)
    !
    class(array_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    class(array_scalar_field), allocatable :: local
    allocate(local, mold=this)
    call local%assign_meta_data(this)
    if (allocated(this%field_data)) then
      local%field_data = atan(this%field_data)
    end if
    call move_alloc(local, res)
  end function array_scalar_atan

  pure function array_scalar_sinh(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\sinh({\rm field})\)
    !
    class(array_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    class(array_scalar_field), allocatable :: local
    allocate(local, mold=this)
    call local%assign_meta_data(this)
    if (allocated(this%field_data)) then
      local%field_data = sinh(this%field_data)
    end if
    call move_alloc(local, res)
  end function array_scalar_sinh

  pure function array_scalar_cosh(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\cosh({\rm field})\)
    !
    class(array_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    class(array_scalar_field), allocatable :: local
    allocate(local, mold=this)
    call local%assign_meta_data(this)
    if (allocated(this%field_data)) then
      local%field_data = cosh(this%field_data)
    end if
    call move_alloc(local, res)
  end function array_scalar_cosh

  pure function array_scalar_tanh(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\tanh({\rm field})\)
    !
    class(array_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    class(array_scalar_field), allocatable :: local
    allocate(local, mold=this)
    call local%assign_meta_data(this)
    if (allocated(this%field_data)) then
      local%field_data = tanh(this%field_data)
    end if
    call move_alloc(local, res)
  end function array_scalar_tanh

  pure function array_scalar_asinh(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\sinh^{-1}({\rm field})\)
    !
    class(array_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    class(array_scalar_field), allocatable :: local
    allocate(local, mold=this)
    call local%assign_meta_data(this)
    if (allocated(this%field_data)) then
      local%field_data = asinh(this%field_data)
    end if
    call move_alloc(local, res)
  end function array_scalar_asinh

  pure function array_scalar_acosh(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\cosh^{-1}({\rm field})\)
    !
    class(array_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    class(array_scalar_field), allocatable :: local
    allocate(local, mold=this)
    call local%assign_meta_data(this)
    if (allocated(this%field_data)) then
      local%field_data = acosh(this%field_data)
    end if
    call move_alloc(local, res)
  end function array_scalar_acosh

  pure function array_scalar_atanh(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\tanh^{-1}({\rm field})\)
    !
    class(array_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    class(array_scalar_field), allocatable :: local
    allocate(local, mold=this)
    call local%assign_meta_data(this)
    if (allocated(this%field_data)) then
      local%field_data = atanh(this%field_data)
    end if
    call move_alloc(local, res)
  end function array_scalar_atanh

  pure function array_scalar_log(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\ln({\rm field})\)
    !
    class(array_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    class(array_scalar_field), allocatable :: local
    allocate(local, mold=this)
    call local%assign_meta_data(this)
    if (allocated(this%field_data)) then
      local%field_data = log(this%field_data)
    end if
    call move_alloc(local, res)
  end function array_scalar_log

  pure function array_scalar_log10(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\log({\rm field})\)
    !
    class(array_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    class(array_scalar_field), allocatable :: local
    allocate(local, mold=this)
    call local%assign_meta_data(this)
    if (allocated(this%field_data)) then
      local%field_data = log10(this%field_data)
    end if
    call move_alloc(local, res)
  end function array_scalar_log10

  pure function array_scalar_exp(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(e^{\rm field}\)
    !
    class(array_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    class(array_scalar_field), allocatable :: local
    allocate(local, mold=this)
    call local%assign_meta_data(this)
    if (allocated(this%field_data)) then
      local%field_data = exp(this%field_data)
    end if
    call move_alloc(local, res)
  end function array_scalar_exp

  pure function array_scalar_abs(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(|{\rm field}|\)
    !
    class(array_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    class(array_scalar_field), allocatable :: local
    allocate(local, mold=this)
    call local%assign_meta_data(this)
    if (allocated(this%field_data)) then
      local%field_data = abs(this%field_data)
    end if
    call move_alloc(local, res)
  end function array_scalar_abs

  pure function array_scalar_sqrt(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\sqrt{{\rm field}}\)
    !
    class(array_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    class(array_scalar_field), allocatable :: local
    allocate(local, mold=this)
    call local%assign_meta_data(this)
    if (allocated(this%field_data)) then
      local%field_data = sqrt(this%field_data)
    end if
    call move_alloc(local, res)
  end function array_scalar_sqrt

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
  
  function array_scalar_laplacian(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\nabla^2 {\rm field}\)
    !
    class(array_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    class(array_scalar_field), allocatable :: local
    allocate(local, mold=this)
    call local%assign_meta_data(this)
    local = this%d_dx(1,2)
    call move_alloc(local, res)
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
    allocate(local, mold=this%vector_result_type())
    call local%assign_meta_data(this)
    local%field_data(:,1) = this%field_data
    call differentiate_1d(local%field_data(:,1),local%colloc_points)
    call move_alloc(local,res)
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
    iseq = .true.
    select type(rhs)
    class is(array_scalar_field)
      if (this%numpoints /= rhs%numpoints) then
        iseq = .false.
        return
      end if
      do i=1,this%numpoints+1
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
    select type(rhs)
    class is(array_scalar_field)
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
    class is(array_vector_field)
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
  end subroutine array_scalar_assign_meta_data


  !=====================================================================
  ! Vector Field Methods
  !=====================================================================


  pure function array_vector_raw_size(this,return_lower_bound, &
                                       return_upper_bound) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! Compute how many elements are in the raw representation of this
    ! field. This would be the number of data points, adjusted based on
    ! how boundary conditions are accounted for.
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
    integer :: res
    res = (this%numpoints + 1)
    if (present(return_lower_bound)) then
      if (.not. return_lower_bound(1)) res = res - 1
    end if
    if (present(return_upper_bound)) then
      if (.not. return_upper_bound(1)) res = res - 1
    end if
    res = res*(this%numdims)
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
  end function array_vector_raw

  subroutine array_vector_set_from_raw(this,raw,provide_lower_bound, &
                                        provide_upper_bound)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! Assigns raw data, such as that produced by 
    ! [[array_vector_field:raw]], to the field. The routine will
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
    integer :: lower, upper
    call check_set_from_raw(this,raw,provide_lower_bound,provide_upper_bound)
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
      allocate(this%field_data(this%numpoints + 1,this%numdims))
    this%field_data(lower:upper,:) = reshape(raw,[upper-lower+1,this%numdims])
  end subroutine array_vector_set_from_raw
  
  function array_vector_vf_m_sf(this,rhs) result(res)
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
#ifdef DEBUG
    call this%check_compatible(rhs)
#endif
    allocate(local, mold=this)
    call local%assign_meta_data(this)
    select type(rhs)
    class is(array_scalar_field)
      forall (i=1:this%numdims) &
        local%field_data(:,i) = this%field_data(:,i) * rhs%field_data
    end select
    call move_alloc(local,res)
  end function array_vector_vf_m_sf

  function array_vector_r_m_vf(lhs,rhs) result(res)
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

  function array_vector_vf_m_r(this,rhs) result(res)
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
  
  function array_vector_vf_d_sf(this,rhs) result(res)
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
#ifdef DEBUG
    call this%check_compatible(rhs)
#endif
    allocate(local, mold=this)
    call local%assign_meta_data(this)
    select type(rhs)
    class is(array_scalar_field)
      forall (i=1:this%numdims) &
        local%field_data(:,i) = this%field_data(:,i) / rhs%field_data
    end select
    call move_alloc(local,res)
  end function array_vector_vf_d_sf

  function array_vector_vf_d_r(this,rhs) result(res)
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
  
  function array_vector_vf_s_vf(this,rhs) result(res)
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
#ifdef DEBUG
    call this%check_compatible(rhs)
#endif
    allocate(local, mold=this)
    select type(rhs)
    class is(array_vector_field)
      if (rhs%numdims > this%numdims) then
        call local%assign_meta_data(rhs)
      else
        call local%assign_meta_data(this)
      end if
      min_dims = min(this%numdims,rhs%numdims)
      max_dims = max(this%numdims,rhs%numdims)
      local%field_data(:,1:min_dims) = this%field_data(:,1:min_dims) &
                                     - rhs%field_data(:,1:min_dims)
      if (rhs%numdims > this%numdims) then
        local%field_data(:,min_dims+1:max_dims) = -rhs%field_data(:,min_dims+1:max_dims)
      else
        local%field_data(:,min_dims+1:max_dims) = this%field_data(:,min_dims+1:max_dims)
      end if
    end select
    call move_alloc(local,res)
  end function array_vector_vf_s_vf

  function array_vector_r_s_vf(lhs,rhs) result(res)
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
    min_dims = min(rhs%numdims, size(lhs))
    max_dims = max(rhs%numdims, size(lhs))
    allocate(local%field_data(local%numpoints + 1, max_dims))
    forall(i=1:local%numpoints+1) &
      local%field_data(i,:min_dims) = lhs(:min_dims) - rhs%field_data(i,:min_dims)
    if (rhs%numdims > size(lhs)) then
      local%field_data(:,min_dims+1:) = -rhs%field_data(:,min_dims+1:)
    else
      local%numdims = size(lhs)
      forall(i=1:local%numpoints+1) &
        local%field_data(i,min_dims+1:) = lhs(min_dims+1:)
    end if
    call move_alloc(local,res)
  end function array_vector_r_s_vf

  function array_vector_vf_s_r(this,rhs) result(res)
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
    min_dims = min(this%numdims, size(rhs))
    max_dims = max(this%numdims, size(rhs))
    allocate(local%field_data(local%numpoints + 1, max_dims))
    forall(i=1:local%numpoints+1) &
      local%field_data(i,:min_dims) = this%field_data(i,:min_dims) - rhs(:min_dims)
    if (this%numdims > size(rhs)) then
      local%field_data(:,min_dims+1:) = this%field_data(:,min_dims+1:)
    else
      local%numdims = size(rhs)
      forall(i=1:local%numpoints+1) &
        local%field_data(i,min_dims+1:) = -rhs(min_dims+1:)
    end if
    call move_alloc(local,res)
  end function array_vector_vf_s_r
  
  function array_vector_vf_a_vf(this,rhs) result(res)
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
#ifdef DEBUG
    call this%check_compatible(rhs)
#endif
    allocate(local, mold=this)
    select type(rhs)
    class is(array_vector_field)
      if (rhs%numdims > this%numdims) then
        call local%assign_meta_data(rhs)
      else
        call local%assign_meta_data(this)
      end if
      min_dims = min(this%numdims,rhs%numdims)
      max_dims = max(this%numdims,rhs%numdims)
      local%field_data(:,1:min_dims) = this%field_data(:,1:min_dims) &
                                     + rhs%field_data(:,1:min_dims)
      if (rhs%numdims > this%numdims) then
        local%field_data(:,min_dims+1:max_dims) = rhs%field_data(:,min_dims+1:max_dims)
      else
        local%field_data(:,min_dims+1:max_dims) = this%field_data(:,min_dims+1:max_dims)
      end if
    end select
    call move_alloc(local,res)
  end function array_vector_vf_a_vf

  function array_vector_r_a_vf(lhs,rhs) result(res)
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
    min_dims = min(rhs%numdims, size(lhs))
    max_dims = max(rhs%numdims, size(lhs))
    allocate(local%field_data(local%numpoints + 1, max_dims))
    forall(i=1:local%numpoints+1) &
      local%field_data(i,:min_dims) = lhs(:min_dims) + rhs%field_data(i,:min_dims)
    if (rhs%numdims > size(lhs)) then
      local%field_data(:,min_dims+1:) = rhs%field_data(:,min_dims+1:)
    else
      local%numdims = size(lhs)
      forall(i=1:local%numpoints+1) &
        local%field_data(i,min_dims+1:) = lhs(min_dims+1:)
    end if
    call move_alloc(local,res)
  end function array_vector_r_a_vf

  function array_vector_vf_a_r(this,rhs) result(res)
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
    min_dims = min(this%numdims, size(rhs))
    max_dims = max(this%numdims, size(rhs))
    allocate(local%field_data(local%numpoints + 1, max_dims))
    forall(i=1:local%numpoints+1) &
      local%field_data(i,:min_dims) = this%field_data(i,:min_dims) + rhs(:min_dims)
    if (this%numdims > size(rhs)) then
      local%field_data(:,min_dims+1:) = this%field_data(:,min_dims+1:)
    else
      local%numdims = size(rhs)
      forall(i=1:local%numpoints+1) &
        local%field_data(i,min_dims+1:) = rhs(min_dims+1:)
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
    allocate(local, mold=this%scalar_result_type())
    call local%assign_meta_data(this)
    select type(this)
    class is(array_vector_field)
      local%field_data = sqrt(sum(this%field_data**2,2))
    end select
    call move_alloc(local,res)
  end function array_vector_norm

  function array_vector_component(this,comp) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! Scalar field containing specified component of the vector field
    !
    class(array_vector_field), intent(in) :: this
    integer, intent(in) :: comp
    class(scalar_field), allocatable :: res
    class(array_scalar_field), allocatable :: local
    allocate(local, mold=this%scalar_result_type())
    call local%assign_meta_data(this)
    select type(this)
    class is(array_vector_field)
      if (comp <= this%numdims) then
        local%field_data = this%field_data(:,comp)
      else
        local%field_data = 0.0_r8
      end if
    end select
    call move_alloc(local,res)
  end function array_vector_component


  function array_vector_laplacian(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\nabla^2 \vec{\rm field}\)
    !
    class(array_vector_field), intent(in) :: this
    class(vector_field), allocatable :: res !! The result of this operation
    class(array_vector_field), allocatable :: local
    allocate(local, mold=this)
    call local%assign_meta_data(this)
    call move_alloc(local, res)
  end function array_vector_laplacian
  
  function array_vector_divergence(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\nabla\cdot \vec{\rm field}\)
    !
    class(array_vector_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    class(array_scalar_field), allocatable :: local
    allocate(local, mold=this%scalar_result_type())
    call local%assign_meta_data(this)
    local%field_data = this%field_data(:,1)
    call differentiate_1d(local%field_data,local%colloc_points)
    call move_alloc(local, res)
  end function array_vector_divergence
  
  function array_vector_curl(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\nabla\times \vec{\rm field}\)
    !
    class(array_vector_field), intent(in) :: this
    class(vector_field), allocatable :: res !! The result of this operation
    class(array_vector_field), allocatable :: local
    allocate(local, mold=this)
    call local%assign_meta_data(this,.false.)
    local%numdims = 3
    allocate(local%field_data(local%numpoints+1,3))
    local%field_data(:,1) = 0.0_r8
    if (this%numdims >= 3) then
      local%field_data(:,2) = -this%field_data(:,3)
      call differentiate_1d(local%field_data(:,2),local%colloc_points)
    else
      local%field_data(:,2) = 0.0_r8
    end if
    if (this%numdims >= 2) then
      local%field_data(:,3) = this%field_data(:,2)
      call differentiate_1d(local%field_data(:,3),local%colloc_points)
    else
      local%field_data(:,3) = 0.0_r8
    end if
    call move_alloc(local, res)
  end function array_vector_curl

  function array_vector_cross_prod(this,rhs) result(res)
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
#ifdef DEBUG
    call this%check_compatible(rhs)
#endif
    allocate(local, mold=this)
    call local%assign_meta_data(this,.false.)
    local%numdims = 3
    allocate(local%field_data(this%numpoints+1,3))
    vec1 = 0
    vec2 = 0
    select type(rhs)
    class is(array_vector_field)
      dims1 = min(3,this%numdims)
      dims2 = min(3,rhs%numdims)
      do i=1,this%numpoints+1
        vec1(:dims1) = this%field_data(i,:dims1)
        vec2(:dims2) = rhs%field_data(i,:dims2)
        local%field_data(i,:) = [vec1(2)*vec2(3) - vec2(2)*vec1(3), &
                                 vec1(3)*vec2(1) - vec2(3)*vec1(1), &
                                 vec1(1)*vec2(2) - vec2(1)*vec1(2)]
      end do
    end select
    call move_alloc(local,res)
  end function array_vector_cross_prod

  function array_vector_dot_prod(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\vec{\rm field} \cdot \vec{\rm field}\)
    !
    class(array_vector_field), intent(in) :: this
    class(vector_field), intent(in) :: rhs
    class(scalar_field), allocatable :: res !! The restult of this operation
    class(array_scalar_field), allocatable :: local
    integer :: min_dims
#ifdef DEBUG
    call this%check_compatible(rhs)
#endif
    allocate(local, mold=this%scalar_result_type())
    call local%assign_meta_data(this)
    select type(rhs)
    class is(array_vector_field)
      min_dims = min(this%numdims,rhs%numdims)
      local%field_data = &
        sum(this%field_data(:,1:min_dims)*rhs%field_data(:,1:min_dims),2)
    end select
    call move_alloc(local,res)
  end function array_vector_dot_prod

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
    iseq = .true.
    select type(rhs)
    class is(array_vector_field)
      if (this%numpoints/=rhs%numpoints) then
        iseq = .false.
        return
      end if
      if (this%numdims > rhs%numdims) then
        dims = rhs%numdims
        iseq = all(this%field_data(:,dims:)==0.0_r8)
        if (.not. iseq) return
      else if (this%numdims < rhs%numdims) then
        dims = this%numdims
        iseq = all(rhs%field_data(:,dims:)==0.0_r8)
        if (.not. iseq) return
      else
        dims = this%numdims
      end if
      do i=1,this%numpoints+1
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
   select type(rhs)
    class is(array_scalar_field)
      this%numpoints = rhs%numpoints
      this%extent = rhs%extent
      this%numdims = this%dimensions()
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
    class is(array_vector_field)
      this%numpoints = rhs%numpoints
      this%extent = rhs%extent
      this%numdims = rhs%numdims
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
  end subroutine array_vector_assign_meta_data

end module array_fields_mod
