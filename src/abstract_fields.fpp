!
!  abstract_fields.f90
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

module abstract_fields_mod
  !* Author: Chris MacMackin
  !  Date: March 2016
  !  License: GPLv3
  !
  ! Provides an abstract data types representing mathematical fields.
  ! Inspiration is taken from Rouson, Xia, and Xu (2011), in that this
  ! type encapsulates many of the mathematical operations which are
  ! performed on fields.
  !
  !####Bibliography
  ! *Scientific Software Design: The Object-Oriented Way*, Rouson, 
  ! Damian and Xia, Jim and Xu, Xiaofeng, 2011, ISBN 9780521888134, 
  ! Cambridge University Press, New York, NY, USA.
  !
  use iso_fortran_env, only: r8 => real64
  implicit none
  private
  
  real(r8) :: tolerance = 1e-10_r8
  public :: set_tol, get_tol
  
  type, abstract, public :: abstract_field
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! An abstract data type representing a mathematical field. 
    ! Inspiration is taken from the abstract calculus pattern described
    ! by Rouson, Xia, and Xu (2011). This particular type is meant only
    ! to provide whatever common properties are shared by fields of 
    ! scalar and vector quantities. These are type-bound procedures 
    ! which work with intrinsic types, either to provide information
    ! about the field or to provide the field data in a raw form.
    !
    !####Bibliography
    ! *Scientific Software Design: The Object-Oriented Way*, Rouson, 
    ! Damian and Xia, Jim and Xu, Xiaofeng, 2011, ISBN 9780521888134, 
    ! Cambridge University Press, New York, NY, USA.
    !
  contains
    procedure :: elements
      !! Specifies the number of individual data points present in this field.
    procedure(f_ret_r), deferred :: domain
      !! Provides array with upper and lower limits of field's domain
    procedure(f_ret_i), deferred :: dimensions
      !! Returns number of dimensions of the field
    procedure(f_rawsize), deferred :: raw_size
      !! Provides the number of pieces of data needed to represent the
      !! field, i.e. the size of the array returned by `get_raw`.
    procedure(f_raw), deferred :: raw
      !! Returns array of data representing state of field. Can be
      !! useful for passing to nonlinear solvers.
    procedure(f_eq_raw), deferred :: set_from_raw
      !! Set the field's values from raw data, such as that produced by 
      !! [[abstract_field:raw]]
    procedure(f_res), deferred :: resolution
      !! Returns array containing number of datapoints in each dimension.
    !procedure(f_eq_meta), deferred :: assign_meta_data
    !  !! Copies all data other than values stored in field from another
    !  !! field object to this one.
    procedure(scalar_factory), deferred :: allocate_scalar_field
      !! Allocates a scalar field to be of the same concrete type as
      !! those returned by type-bound procedures of this field which
      !! produce scalar fields.
    procedure(vector_factory), deferred :: allocate_vector_field
      !! Allocates a vector field to be of the same concrete type as
      !! those returned by type-bound procedures of this field which
      !! produce vector fields.
    procedure(id_pos), deferred :: id_to_position
      !! Given the ID number of a location in the field, returns the
      !! coordinates of that position
  end type abstract_field
  
  type, extends(abstract_field), abstract, public :: scalar_field
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! An abstract data type representing a mathematical field of scalar
    ! values. Inspiration is taken from the abstract calculus pattern
    ! described by Rouson, Xia, and Xu (2011), in that this type
    ! encapsulates all of Fortran's intrinsic mathematical operations,
    ! allowing for a very natural syntax to be used when manipulating
    ! fields. Support is also given for calculus operations.
    !
    ! Note that, when performing an operation on two fields, an error
    ! will occur if one of the fields is of the incorrect type or if
    ! the sizes of the fields do not match.
    !
    !####Bibliography
    ! *Scientific Software Design: The Object-Oriented Way*, Rouson, 
    ! Damian and Xia, Jim and Xu, Xiaofeng, 2011, ISBN 9780521888134, 
    ! Cambridge University Press, New York, NY, USA.
    !
   contains
    private
    procedure(sf_sf), deferred :: field_multiply_field
      !! \({\rm field} \times {\rm field}\)
    procedure(sf_vf), deferred :: field_multiply_vecfield
      !! \({\rm field} \times {\rm \vec{field}}\)
    procedure(r_sf), pass(rhs), deferred :: real_multiply_field
      !! \({\rm real}  \times {\rm field}\)
    procedure(r1d_sf), pass(rhs), deferred :: real_array_multiply_field
      !! \(\vec{\rm real}  \times {\rm field}\)
    procedure(sf_r), deferred :: field_multiply_real
      !! \({\rm field} \times {\rm real}\)
    procedure(sf_r1d), deferred :: field_multiply_real_array
      !! \( {\rm field}  \times \vec{\rm real} \)
    procedure(sf_sf), deferred :: field_divide_field
      !! \(\frac{\rm field}{\rm field}\)
    procedure(r_sf), pass(rhs), deferred :: real_divide_field
      !! \(\frac{\rm real}{\rm field}\)
    procedure(r1d_sf), pass(rhs), deferred :: real_array_divide_field
      !! \(\frac{\vec{\rm real}}{\rm field}\)
    procedure(sf_r), deferred :: field_divide_real
      !! \(\frac{\rm field}{\rm real}\)
    procedure(sf_sf), deferred :: field_add_field
      !! \({\rm field} + {\rm field}\)
    procedure(r_sf), pass(rhs), deferred :: real_add_field
      !! \({\rm real} + {\rm field}\)
    procedure(sf_r), deferred :: field_add_real
      !! \({\rm field} + {\rm real}\)
    procedure(sf_sf), deferred :: field_sub_field
      !! \({\rm field} - {\rm field}\)
    procedure(r_sf), pass(rhs), deferred :: real_sub_field
      !! \({\rm real} - {\rm field}\)
    procedure(sf_r), deferred :: field_sub_real
      !! \({\rm field} - {\rm real}\)
    procedure(sf_r), deferred :: field_pow_real
      !! \({\rm field}^{\rm real}\)
    procedure(sf_r4), deferred :: field_pow_real4
      !! \({\rm field}^{\rm real}\)
    procedure(sf_i), deferred :: field_pow_int
      !! \({\rm field}^{\rm int}\)
    !~ procedure(sf_elem), public, deferred :: get_element
      !! Returns one of the constituent values of the field, i.e. the 
      !! field's value at a particular location
#:for FUNC, TEX in UNARY_FUNCTIONS
    procedure(sf_ret_sf), deferred :: ${FUNC}$
      !! \(${TEX}$({\rm field})\)
#:endfor
    procedure(sf_ret_r), deferred :: minval
      !! \(\ ({\rm field})\)
    procedure(sf_ret_r), deferred :: maxval
      !! \(\ ({\rm field})\)
    procedure(sf_dx), public, deferred :: d_dx
      !! \(\frac{\partial^n}{\partial x_i^n}({\rm field})\)
    procedure(sf_grad), deferred :: gradient
      !! \(\nabla {\rm field}\)
    procedure(sf_laplace), deferred :: laplacian
      !! \(\nabla^2 {\rm field}\)
    procedure(sf_eq_sf), deferred :: assign_field
      !! \({\rm field} = {\rm field}\)
    generic, public :: operator(*) => field_multiply_field, &
        field_multiply_vecfield, real_multiply_field, &
        field_multiply_real, real_array_multiply_field, &
        field_multiply_real_array
    generic, public :: operator(/) => field_divide_field, &
        field_divide_real, real_divide_field, real_array_divide_field
    generic, public :: operator(+) => field_add_field, field_add_real, &
        real_add_field
    generic, public :: operator(-) => field_sub_field, field_sub_real, &
        real_sub_field
    generic, public :: operator(**) => field_pow_real, field_pow_int, &
        field_pow_real4
    generic, public :: operator(.grad.) => gradient
    generic, public :: operator(.laplacian.) => laplacian
    generic, public :: assignment(=) => assign_field
    procedure(sf_is_equal), deferred :: is_equal
      !! Checks fields are equal within a tolerance
    generic, public :: operator(==) => is_equal
    procedure :: negation => scalar_field_negation
    generic, public :: operator(-) => negation
  end type scalar_field
  
  type, extends(abstract_field), abstract, public :: vector_field
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! An abstract data type representing a mathematical field of vector
    ! values. Inspiration is taken from the abstract calculus pattern
    ! described by Rouson, Xia, and Xu (2011), in that this type
    ! encapsulates all of Fortran's intrinsic mathematical operations,
    ! allowing for a very natural syntax to be used when manipulating
    ! fields. Support is also given for calculus operations.
    !
    !####Bibliography
    ! *Scientific Software Design: The Object-Oriented Way*, Rouson, 
    ! Damian and Xia, Jim and Xu, Xiaofeng, 2011, ISBN 9780521888134, 
    ! Cambridge University Press, New York, NY, USA.
    !
  contains
    private
    procedure(vf_ret_i), public, deferred :: vector_dimensions
      !! Returns dimension of the vectors in the field
    procedure(vf_sf), deferred :: field_multiply_field
      !! \({\rm\vec{field}} \times {\rm field}\)
    procedure(r_vf), pass(rhs), deferred :: real_multiply_field
      !! \({\rm real}  \times {\rm \vec{field}}\)
    procedure(vf_r), deferred :: field_multiply_real
      !! \({\rm \vec{field}} \times {\rm real}\)
    procedure(vf_sf), deferred :: field_divide_field
      !! \(\frac{\rm \vec{field}}{\rm field}\)
    procedure(vf_r), deferred :: field_divide_real
      !! \(\frac{\rm \vec{field}}{\rm real}\)
    procedure(vf_vf), deferred :: field_add_field
      !! \({\rm \vec{field}} + {\rm \vec{field}}\)
    procedure(vr_vf), pass(rhs), deferred :: real_add_field
      !! \({\rm \vec{real}} + {\rm \vec{field}}\)
    procedure(vf_vr), deferred :: field_add_real
      !! \({\rm \vec{field}} + {\rm \vec{real}}\)
    procedure(vf_vf), deferred :: field_sub_field
      !! \({\rm \vec{field}} - {\rm \vec{field}}\)
    procedure(vr_vf), pass(rhs), deferred :: real_sub_field
      !! \({\rm \vec{real}} - {\rm \vec{field}}\)
    procedure(vf_vr), deferred :: field_sub_real
      !! \({\rm \vec{field}} - {\rm \vec{real}}\)
!~     procedure(vf_elem_vec), deferred :: get_element_vector
!~       !! Returns ones of the constituent vectors of the field, i.e. the 
!~       !! field's value at a particular location
!~     procedure(vf_elem_comp), deferred :: get_element_component
!~       !! Returns one of the components of a constituent vector of the 
!~       !! field, i.e. the component of the field's value at a particular 
!~       !! location
!~     generic, public :: get_element => get_element_vector, get_element_component
!~       !! Returns a constituent value of the field, i.e. the vector or 
!~       !! vector component giving the field's value at a particular 
!~       !! location
    procedure(vf_norm), public, deferred :: norm
      !! \(\lVert {\rm \vec{field}} \rVert\)
    procedure(vf_comp), public, deferred :: component
      !! Returns a scalar field containing the specified component of 
      !! the vector field
    procedure(vf_dx), public, deferred :: d_dx
      !! \(\frac{\partial^n}{\partial x_i^n}({\rm \vec{field}})\)
    procedure(vf_comp_dx), public, deferred :: component_d_dx
      !! \(\frac{\partial^n}{\partial x_i^n}({\rm field_j})\)
    procedure(impure_vf_ret_sf), deferred :: divergence
      !! \(\nabla\cdot {\rm \vec{field}}\)
    procedure(impure_vf_ret_vf), deferred :: curl
      !! \(\nabla\times {\rm \vec{field}}\)
    procedure(impure_vf_ret_vf), deferred :: laplacian
      !! \(\nabla^2 {\rm \vec{field}}\)
    procedure(vf_vf_ret_sf), deferred :: field_dot_field
      !! \({\rm \vec{field}} \cdot {\rm \vec{field}}\)
    procedure(vf_vr_ret_sf), deferred :: field_dot_real
      !! \({\rm \vec{field}} \cdot {\rm \vec{real}}\)
    procedure(vr_vf_ret_sf), pass(rhs), deferred :: real_dot_field
      !! \({\rm \vec{real}} \cdot {\rm \vec{field}}\)
    procedure(vf_vf), deferred :: field_cross_field
      !! \({\rm\vec{field}} \times {\rm\vec{field}}\)
    procedure(vf_vr), deferred :: field_cross_real
      !! \({\rm\vec{field}} \times {\rm\vec{real}}\)
    procedure(vr_vf), pass(rhs), deferred :: real_cross_field
      !! \({\rm\vec{real}} \times {\rm\vec{field}}\)
    procedure(vf_eq_vf), deferred :: assign_field
      !! \({\rm \vec{field}} = {\rm \vec{field}}\)
    generic, public :: operator(*) => field_multiply_field, &
        real_multiply_field, field_multiply_real
    generic, public :: operator(/) => field_divide_field, field_divide_real
    generic, public :: operator(+) => field_add_field, field_add_real, &
        real_add_field
    generic, public :: operator(-) => field_sub_field, field_sub_real, &
        real_sub_field
    generic, public :: operator(.div.) => divergence
    generic, public :: operator(.curl.) => curl
    generic, public :: operator(.laplacian.) => laplacian
    generic, public :: operator(.dot.) => field_dot_field, &
         real_dot_field, field_dot_real
    generic, public :: operator(.cross.) => field_cross_field, &
         real_cross_field, field_cross_real
    generic, public :: assignment(=) => assign_field
    procedure(vf_is_equal), deferred :: is_equal
      !! Checks fields are equal within a tolerance
    generic, public :: operator(==) => is_equal
    procedure :: negation => vector_field_negation
    generic, public :: operator(-) => negation
  end type vector_field


  abstract interface
    pure function f_ret_r(this)
      import :: abstract_field
      import :: r8
      class(abstract_field), intent(in) :: this
      real(r8), dimension(:,:), allocatable :: f_ret_r
        !* A 2D array of shape \(n \times 2\), where n is the number of
        !  dimensions of the field. Each row contains the lower and the
        !  upper extent of the fields domain in the corresponding 
        !  dimension
    end function f_ret_r

    elemental function f_ret_i(this)
      import :: abstract_field
      class(abstract_field), intent(in) :: this
      integer :: f_ret_i
    end function f_ret_i

    elemental function vf_ret_i(this)
      import :: vector_field
      class(vector_field), intent(in) :: this
      integer :: vf_ret_i
    end function vf_ret_i

    pure function f_rawsize(this,exclude_lower_bound,exclude_upper_bound)
      import :: abstract_field
      class(abstract_field), intent(in) :: this
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
      integer :: f_rawsize
        !! The number of elements in the array returned by `this%raw()`
        !! when given these values of `return_lower_bound` and 
        !! `return_upper_bound`.
    end function f_rawsize
    
    pure function f_raw(this,exclude_lower_bound,exclude_upper_bound)
      !* @BUG The returned value has length `this%raw_size()`, but
      !  a bug in gfortran 4.8 (fixed by version 5) caused the compiler
      !  to segfault if it was declared as such. As a workaround, it is
      !  allocatable instead.
      !
      import :: abstract_field
      import :: r8
      class(abstract_field), intent(in) :: this
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
      real(r8), dimension(:), allocatable :: f_raw
        !! Array containing data needed to describe field
    end function f_raw
    
!~    pure function f_jacob(this,order)
!~      !* A Jacobian matrix \(\frac{\partial}{\partial y_j}f_i(\vec{y})\)
!~      ! is returned, where \(\vec{y}\) is the field represented as a 1D
!~      ! array of values (i.e. the output of [[abstract_field:raw]])
!~      ! and \( f_i(\vec{y}) = \frac{\partial^{n}y_i}{\partial x^n} \).
!~      ! Here, \(n\) is the `order` of the derivative to be taken, with 0
!~      ! corresponding to not taking any derivative.
!~      !
!~      ! @BUG The returned value has shape `(this%raw_size(),
!~      ! this%raw_size())`, but a bug in gfortran 4.8 (fixed by version
!~      ! 5) caused the compiler to segfault if it was declared as such.
!~      ! As a workaround, it is allocatable instead.
!~      !
!~      import :: abstract_field
!~      import :: r8
!~      class(abstract_field), intent(in) :: this
!~      integer, intent(in), optional :: order
!~        !! The order of the derivative of the field whose Jacobian is
!~        !! to be returned. Default is 0 (no differentiation)
!~      real(r8), dimension(:,:), allocatable :: f_jacob
!~        !! The resulting Jacobian matrix
!~    end function f_jacob

    pure function f_res(this)
      import :: abstract_field
      class(abstract_field), intent(in) :: this
      integer, dimension(:), allocatable :: f_res
        !! Array specifying the number of data points in each dimension.
    end function f_res

    pure subroutine f_eq_raw(this,raw,provide_lower_bound,provide_upper_bound)
      !! Assigns raw data, such as that produced by 
      !! [[abstract_field:raw]], to the field
      import :: abstract_field
      import :: r8
      class(abstract_field), intent(inout) :: this
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
    end subroutine f_eq_raw
    
    pure subroutine f_eq_meta(this, rhs, alloc)
      import :: abstract_field
      class(abstract_field), intent(inout) :: this
      class(abstract_field), intent(in) :: rhs
        !! The field whose metadata (domain, resolution, etc) is to be
        !! copied
      logical, optional, intent(in) :: alloc
        !! If present and false, do not allocate the array of `this`.
    end subroutine f_eq_meta

    pure subroutine scalar_factory(this, new_field)
      import :: abstract_field
      import :: scalar_field
      class(abstract_field), intent(in)               :: this
      class(scalar_field), allocatable, intent(inout) :: new_field
        !! A field which, upon return, is allocated to be of the same
        !! concrete type as scalar fields produced by `this`.
    end subroutine scalar_factory

    pure subroutine vector_factory(this, new_field)
      import :: abstract_field
      import :: vector_field
      class(abstract_field), intent(in)               :: this
      class(vector_field), allocatable, intent(inout) :: new_field
        !! A field which, upon return, is allocated to be of the same
        !! concrete type as vector fields produced by `this`.
    end subroutine vector_factory

    pure function id_pos(this, id)
      import :: abstract_field
      import :: r8
      class(abstract_field), intent(in) :: this
      integer, intent(in)               :: id
        !! The ID number for some location in the field
      real(r8), dimension(:), allocatable :: id_pos
        !! The coordinates for this location in the field
    end function id_pos
    
    pure function sf_sf(this,rhs)
      !! \({\rm field} [{\rm operator}] {\rm field}\)
      import :: scalar_field
      class(scalar_field), intent(in) :: this, rhs
      class(scalar_field), allocatable :: sf_sf !! The restult of this operation
    end function sf_sf

    pure function sf_vf(this,rhs)
      !! \({\rm field} [{\rm operator}] {\rm \vec{field}}\)
      import :: scalar_field
      import :: vector_field
      class(scalar_field), intent(in) :: this
      class(vector_field), intent(in) :: rhs
      class(vector_field), allocatable :: sf_vf !! The result of this operation
    end function sf_vf

    pure function r_sf(lhs,rhs)
      !! \({\rm real} [{\rm operator}] {\rm field}\)
      import :: scalar_field
      import :: r8
      real(r8), intent(in) :: lhs
      class(scalar_field), intent(in) :: rhs
      class(scalar_field), allocatable :: r_sf !! The result of this operation
    end function r_sf

    pure function r1d_sf(lhs,rhs)
      !! \(\vec{\rm real} [{\rm operator}] {\rm field}\)
      import :: scalar_field
      import :: vector_field
      import :: r8
      real(r8), dimension(:), intent(in) :: lhs
      class(scalar_field), intent(in) :: rhs
      class(vector_field), allocatable :: r1d_sf !! The result of this operation
    end function r1d_sf
  
    pure function sf_r(this,rhs)
      !! \({\rm field} [{\rm operator}] {\rm real}\)
      import :: scalar_field
      import :: r8
      class(scalar_field), intent(in) :: this
      real(r8), intent(in) :: rhs
      class(scalar_field), allocatable :: sf_r !! The result of this operation
    end function sf_r
  
    pure function sf_r1d(this,rhs)
      !! \({\rm field} [{\rm operator}] \vec{\rm real}\)
      import :: scalar_field
      import :: vector_field
      import :: r8
      class(scalar_field), intent(in) :: this
      real(r8), dimension(:), intent(in) :: rhs
      class(vector_field), allocatable :: sf_r1d !! The result of this operation
    end function sf_r1d
  
    pure function sf_r4(this,rhs)
      !! \({\rm field} [{\rm operator}] {\rm real}\)
      import :: scalar_field
      class(scalar_field), intent(in) :: this
      real, intent(in) :: rhs
      class(scalar_field), allocatable :: sf_r4 !! The result of this operation
    end function sf_r4
  
    pure function sf_i(this,rhs)
      !! \({\rm field} [{\rm operator}] {\rm integer}\)
      import :: scalar_field
      class(scalar_field), intent(in) :: this
      integer, intent(in) :: rhs
      class(scalar_field), allocatable :: sf_i !! The result of this operation
    end function sf_i
  
    pure function sf_elem(this,element)
      !! Returns an element of the field corresponding to the provided ID 
      !! number
      import :: scalar_field
      import :: r8
      class(scalar_field), intent(in) :: this
      integer, intent(in) :: element
        !! The ID number of the field element to be returned
      real(r8) :: sf_elem
        !! The value of the field corresponding to the specified ID
    end function sf_elem

    pure subroutine sf_set_elem(this,element,val)
      !! Sets the element of the field corresponding to the given ID to
      !! the given value.
      import :: scalar_field
      import :: r8
      class(scalar_field), intent(inout) :: this
      integer, intent(in) :: element
        !! The ID number of the field element to be set
      real(r8), intent(in) :: val
        !! The new value the field element is to be set to
    end subroutine sf_set_elem

    pure function sf_ret_sf(this)
      !! \([{\rm operator}] {\rm field}\)
      import :: scalar_field
      class(scalar_field), intent(in) :: this
      class(scalar_field), allocatable :: sf_ret_sf !! The result of this operation
    end function sf_ret_sf

    pure function sf_ret_r(this)
      !! \([{\rm operator}] {\rm field}\)
      import :: scalar_field
      import :: r8
      class(scalar_field), intent(in) :: this
      real(r8) :: sf_ret_r !! The result of this operation
    end function sf_ret_r
    
    function sf_grad(this)
      !! \([{\rm operator}] {\rm field}\)
      import :: scalar_field
      import :: vector_field
      class(scalar_field), intent(in) :: this
      class(vector_field), allocatable :: sf_grad !! The result of this operation
    end function sf_grad
  
    function sf_laplace(this)
      !! \([{\rm operator}] {\rm field}\)
      import :: scalar_field
      class(scalar_field), intent(in) :: this
      class(scalar_field), allocatable :: sf_laplace !! The result of this operation
    end function sf_laplace
    
    elemental subroutine sf_eq_sf(this,rhs)
      !! \({\rm field} = {\rm field}\)
      import :: scalar_field
      class(scalar_field), intent(inout) :: this
      class(scalar_field), intent(in) :: rhs
    end subroutine sf_eq_sf
    
    pure function vf_sf(this,rhs)
      !! \({\rm \vec{field}} [{\rm operator}] {\rm field}\)
      import :: scalar_field
      import :: vector_field
      class(vector_field), intent(in) :: this
      class(scalar_field), intent(in) :: rhs
      class(vector_field), allocatable :: vf_sf !! The restult of this operation
    end function vf_sf

    pure function vf_vf(this,rhs)
      !! \({\rm \vec{field}} [{\rm operator}] {\rm \vec{field}}\)
      import :: vector_field
      class(vector_field), intent(in) :: this, rhs
      class(vector_field), allocatable :: vf_vf !! The result of this operation
    end function vf_vf

    pure function r_vf(lhs,rhs)
      !! \({\rm real} [{\rm operator}] {\rm field}\)
      import :: vector_field
      import :: r8
      real(r8), intent(in) :: lhs
      class(vector_field), intent(in) :: rhs
      class(vector_field), allocatable :: r_vf !! The result of this operation
    end function r_vf
  
    pure function vf_r(this,rhs)
      !! \({\rm field} [{\rm operator}] {\rm real}\)
      import :: vector_field
      import :: r8
      class(vector_field), intent(in) :: this
      real(r8), intent(in) :: rhs
      class(vector_field), allocatable :: vf_r !! The result of this operation
    end function vf_r
    
    pure function vf_elem_vec(this,element)
      !! Returns a vector of the field corresponding to the provided ID 
      !! number
      import :: vector_field
      import :: r8
      class(vector_field), intent(in) :: this
      integer, intent(in) :: element
        !! The ID number of the field element to be returned
      real(r8), allocatable, dimension(:) :: vf_elem_vec
        !! The vector in the field corresponding to the specified ID
    end function vf_elem_vec
    
    pure function vf_elem_comp(this,element,component)
      !! Returns a component of the vector of the field corresponding to 
      !! the provided ID number
      import :: vector_field
      import :: r8
      class(vector_field), intent(in) :: this
      integer, intent(in) :: element
        !! The ID number of the field element to be returned
      integer, intent(in) :: component
        !! The number of the vector component to be returned
      real(r8), allocatable, dimension(:) :: vf_elem_comp
        !! The vector component in the field corresponding to the 
        !! specified ID
    end function vf_elem_comp

    pure subroutine vf_set_elem(this,element,val)
      !! Sets the element of the field corresponding to the given ID to
      !! the given vector value.
      import :: vector_field
      import :: r8
      class(vector_field), intent(inout) :: this
      integer, intent(in) :: element
        !! The ID number of the field element to be set
      real(r8), dimension(:), intent(in) :: val
        !! The new vector value the field element is to be set to
    end subroutine vf_set_elem

    pure function vf_ret_sf(this)
      !! \([{\rm operator}] {\rm field}\)
      import :: scalar_field
      import :: vector_field
      class(vector_field), intent(in) :: this
      class(scalar_field), allocatable :: vf_ret_sf !! The result of this operation
    end function vf_ret_sf

    function impure_vf_ret_sf(this)
      !! \([{\rm operator}] {\rm field}\)
      import :: scalar_field
      import :: vector_field
      class(vector_field), intent(in) :: this
      class(scalar_field), allocatable :: impure_vf_ret_sf !! The result of this operation
    end function impure_vf_ret_sf
    
    pure function vf_ret_vf(this)
      !! \([{\rm operator}] {\rm field}\)
      import :: vector_field
      class(vector_field), intent(in) :: this
      class(vector_field), allocatable :: vf_ret_vf !! The result of this operation
    end function vf_ret_vf
    
    function impure_vf_ret_vf(this)
      !! \([{\rm operator}] {\rm field}\)
      import :: vector_field
      class(vector_field), intent(in) :: this
      class(vector_field), allocatable :: impure_vf_ret_vf !! The result of this operation
    end function impure_vf_ret_vf
    
    pure function vf_norm(this)
      !! \([{\rm operator}] {\rm field}\)
      import :: scalar_field
      import :: vector_field
      class(vector_field), intent(in) :: this
      class(scalar_field), allocatable :: vf_norm !! The result of this operation
    end function vf_norm
    
    pure function vf_vf_ret_sf(this,rhs)
      !! \({\rm \vec{field}} [{\rm operator}] {\rm \vec{field}}\)
      import :: scalar_field
      import :: vector_field
      class(vector_field), intent(in) :: this, rhs
      class(scalar_field), allocatable :: vf_vf_ret_sf !! The result of this operation
    end function vf_vf_ret_sf

    pure function vr_vf_ret_sf(lhs,rhs)
      !! \({\rm \vec{real}} [{\rm operator}] {\rm \vec{field}}\)
      import :: scalar_field
      import :: vector_field
      import :: r8
      real(r8), dimension(:), intent(in) :: lhs
      class(vector_field), intent(in) :: rhs
      class(scalar_field), allocatable :: vr_vf_ret_sf !! The result of this operation
    end function vr_vf_ret_sf

    pure function vf_vr_ret_sf(this,rhs)
      !! \({\rm \vec{field}} [{\rm operator}] {\rm \vec{real}}\)
      import :: scalar_field
      import :: vector_field
      import :: r8
      class(vector_field), intent(in) :: this
      real(r8), dimension(:), intent(in) :: rhs
      class(scalar_field), allocatable :: vf_vr_ret_sf !! The result of this operation
    end function vf_vr_ret_sf

    function sf_dx(this, dir, order)
      !! \(\frac{\partial^{\rm order}}{\partial x_{\rm dir}^{\rm order}}{\rm field}\)
      import :: scalar_field
      class(scalar_field), intent(in) :: this
      integer, intent(in) :: dir !! Direction in which to differentiation
      integer, optional, intent(in) :: order !! Order of the derivative, default = 1
      class(scalar_field), allocatable :: sf_dx
    end function sf_dx
    
    function vf_dx(this, dir, order)
      !! \(\frac{\partial^{\rm order}}{\partial x_{\rm dir}^{\rm order}}{\rm \vec{field}}\)
      import :: vector_field
      class(vector_field), intent(in) :: this
      integer, intent(in) :: dir !! Direction in which to differentiation
      integer, optional, intent(in) :: order !! Order of the derivative, default = 1
      class(vector_field), allocatable :: vf_dx !! The derivative
    end function vf_dx

    function vf_comp_dx(this, dir, component, order)
      !! \(\frac{\partial^{\rm order}}{\partial x_{\rm dir}^{\rm order}}{\rm field_{component}}\)
      import :: vector_field
      import :: scalar_field
      class(vector_field), intent(in) :: this
      integer, intent(in) :: dir !! Direction in which to differentiation
      integer, intent(in) :: component !! Which component of the vector is being differentiated
      integer, optional, intent(in) :: order !! Order of the derivative, default = 1
      class(scalar_field), allocatable :: vf_comp_dx !! The derivative
    end function vf_comp_dx

    elemental subroutine vf_eq_vf(this,rhs)
      !! \({\rm \vec{field}} = {\rm \vec{field}}\)
      import :: vector_field
      class(vector_field), intent(inout) :: this
      class(vector_field), intent(in) :: rhs
    end subroutine vf_eq_vf

    pure function vr_vf(lhs,rhs)
      !! \({\rm \vec{real}} [{\rm operator}] {\rm field}\)
      import :: vector_field
      import :: r8
      real(r8), dimension(:), intent(in) :: lhs
      class(vector_field), intent(in) :: rhs
      class(vector_field), allocatable :: vr_vf !! The result of this operation
    end function vr_vf

    pure function vf_vr(this,rhs)
      !! \({\rm \vec{field}} [{\rm operator}] {\rm \vec{real}}\)
      import :: vector_field
      import :: r8
      class(vector_field), intent(in) :: this
      real(r8), dimension(:), intent(in) :: rhs
      class(vector_field), allocatable :: vf_vr !! The result of this operation
    end function vf_vr

    pure function vf_comp(this,comp)
      !! Scalar field containing specified component of the vector field
      import :: scalar_field
      import :: vector_field
      class(vector_field), intent(in) :: this
      integer, intent(in) :: comp !! The index of the component
      class(scalar_field), allocatable :: vf_comp !! Component number `comp`
    end function vf_comp

!~     subroutine sf_bound(this,direction,lower,free_bound,bound_val, &
!~                         bound_deriv)
!~       !* Sets boundary conditions and values. If a boundary is not
!~       !  explicitly set by calling this subroutine then it defaults to
!~       !  being free.
!~       !
!~       import :: scalar_field
!~       import :: r8
!~       class(scalar_field), intent(inout) :: this
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
!~     end subroutine sf_bound
!~ 
!~     subroutine vf_bound(this,direction,lower,free_bound,bound_val, &
!~                         bound_deriv)
!~       !* Sets boundary conditions and values. If a boundary is not
!~       !  explicitly set by calling this subroutine then it defaults to
!~       !  being free.
!~       !
!~       import :: vector_field
!~       import :: r8
!~       class(vector_field), intent(inout) :: this
!~       integer, intent(in) :: direction
!~         !! The number corresponding to the direction/dimension whose
!~         !! boundary condition is to be set
!~       logical, optional, intent(in) :: lower
!~         !! Sets lower boundary if true (default), upper if false
!~       logical, optional, intent(in) :: free_bound
!~         !! If true, makes this a free boundary. Any boundary values 
!~         !! passed will be ignored. Default is `.false.`.
!~       real(r8), dimension(:), optional, intent(in) :: bound_val
!~         !! Value of the field at the boundary. Default is 0.
!~       real(r8), dimension(:), optional, intent(in) :: bound_deriv
!~         !! Value of the first derivative of the field at the boundary.
!~         !! Default is 0.
!~     end subroutine vf_bound
!~     
!~     function sf_same_bounds(this,other)
!~       import :: scalar_field
!~       class(scalar_field), intent(in) :: this
!~       class(scalar_field), intent(in) :: other
!~         !! The fields whose boundary conditions are being compared to
!~         !! those of this one
!~       logical :: sf_same_bounds
!~         !! `.true.` if boundary conditions agree within tolerance,
!~         !! `.false.` otherwise
!~     end function sf_same_bounds
!~     
!~     function vf_same_bounds(this,other)
!~       import :: vector_field
!~       class(vector_field), intent(in) :: this
!~       class(vector_field), intent(in) :: other
!~         !! The fields whose boundary conditions are being compared to
!~         !! those of this one
!~       logical :: vf_same_bounds
!~         !! `.true.` if boundary conditions agree within tolerance,
!~         !! `.false.` otherwise
!~     end function vf_same_bounds

    pure function sf_is_equal(this,rhs) result(iseq)
      import :: scalar_field
      class(scalar_field), intent(in) :: this
      class(scalar_field), intent(in) :: rhs
      logical :: iseq
        !! True if two fields are equal, false otherwise
    end function sf_is_equal

    pure function vf_is_equal(this,rhs) result(iseq)
      import :: vector_field
      class(vector_field), intent(in) :: this
      class(vector_field), intent(in) :: rhs
      logical :: iseq
        !! True if two fields are equal, false otherwise
    end function vf_is_equal
  end interface

#:for FUNC, TEX in UNARY_FUNCTIONS
  interface ${FUNC}$
    module procedure :: scalar_field_${FUNC}$
  end interface ${FUNC}$

#:endfor

  interface minval
    module procedure :: scalar_field_minval
  end interface

  interface maxval
    module procedure :: scalar_field_maxval
  end interface

$:public_unary()
  public :: minval
  public :: maxval
  
contains

  pure function elements(this)
    !* Author: Chris MacMackin
    !  Date: June 2016
    !
    ! Gives the number of individual data points present in the field.
    !
    class(abstract_field), intent(in) :: this
    integer :: elements
    elements = product(this%resolution())
  end function elements

#:for FUNC, TEX in UNARY_FUNCTIONS
  pure function scalar_field_${FUNC}$(field) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! Computes and returns \(${TEX}$({\rm field})\) for a scalar field.
    !
    class(scalar_field), intent(in) :: field
    class(scalar_field), allocatable :: res 
    class(scalar_field), allocatable :: tmp
    allocate(res, mold=field)
    res = field%${FUNC}$()
  end function scalar_field_${FUNC}$
  
#:endfor

  pure function scalar_field_minval(field) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! Computes and returns the minimum of the values in a scalar field.
    !
    class(scalar_field), intent(in) :: field
    real(r8) :: res
    res = field%minval()
  end function scalar_field_minval

  pure function scalar_field_maxval(field) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! Computes and returns the maximum of the values in a scalar field.
    !
    class(scalar_field), intent(in) :: field
    real(r8) :: res
    res = field%maxval()
  end function scalar_field_maxval

  function scalar_field_negation(field) result(res)
    !* Author: Chris MacMackin
    !  Date: April 2016
    !
    ! Computes the negative version of the scalar field.
    !
    class(scalar_field), intent(in) :: field
    class(scalar_field), allocatable :: res
    allocate(res,mold=field)
    res = 0.0_r8 - field
  end function scalar_field_negation

  function vector_field_negation(field) result(res)
    !* Author: Chris MacMackin
    !  Date: April 2016
    !
    ! Computes the negative version of the vector field.
    !
    class(vector_field), intent(in) :: field
    class(vector_field), allocatable :: res
    allocate(res,mold=field)
    res = [0.0_r8] - field
  end function vector_field_negation
  
  subroutine set_tol(tol)
    !* Author: Chris MacMackin
    !  Date: April 2016
    !
    ! Sets the tolerance within which two fields must agree to be 
    ! considered equal.
    !
    real(r8), intent(in) :: tol
    tolerance = tol
  end subroutine set_tol
  
  pure real(r8) function get_tol()
    !* Author: Chris MacMackin
    !  Date: April 2016
    !
    ! Returns the tolerance within which two fields must agree to be
    ! considered equal.
    !
    get_tol = tolerance
  end function get_tol

end module abstract_fields_mod
