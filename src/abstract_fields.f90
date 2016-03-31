!
!  abstract_fields.f90
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

module abstract_fields_mod
  !* Author: Chris MacMackin
  !  Date: March 2016
  !  License: LGPLv3
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
  
  character(len=9), dimension(4), parameter, public :: &
      boundary_conditions = ['Dirichlet', 'Neumann  ', 'Cauchy   ', 'Periodic ']
      !! Code for each type of boundary condition is its index
  
  type, abstract, public :: abstract_field
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! An abstract data type representing a mathematical field. 
    ! Inspiration is taken from the abstract calculus pattern described
    ! by Rouson, Xia, and Xu (2011). This particular type is meant only
    ! to provide whatever common, although largely non-mathematical,
    ! properties are shared by fields of scalar and vector quantities.
    !
    !####Bibliography
    ! *Scientific Software Design: The Object-Oriented Way*, Rouson, 
    ! Damian and Xia, Jim and Xu, Xiaofeng, 2011, ISBN 9780521888134, 
    ! Cambridge University Press, New York, NY, USA.
    !
  contains
    procedure(f_ret_r), deferred :: domain
      !! Provides array with upper and lower limits of field's domain
    procedure(f_ret_i), deferred :: dimensions
      !! Returns number of dimensions of the field
    procedure(f_ret_i), deferred :: raw_size
      !! Provides the number of pieces of data needed to represent the
      !! field, i.e. the size of the array returned by `get_raw`.
    procedure(f_raw), public, deferred :: get_raw
      !! Returns array of data representing state of field. Can be
      !! useful for passing to nonlinear solvers.
    procedure(f_jacob), public, deferred :: d_dself
      !! Returns Jacobian for desired x-derivative of each value in
      !! field relative to all other values in field.
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
    procedure(sf_r), deferred :: field_multiply_real
      !! \({\rm field} \times {\rm real}\)
    procedure(sf_sf), deferred :: field_divide_field
      !! \(\frac{\rm field}{\rm field}\)
    procedure(r_sf), pass(rhs), deferred :: real_divide_field
      !! \(\frac{\rm real}{\rm field}\)
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
    procedure(sf_r), deferred :: field_pow_real4
      !! \({\rm field}^{\rm real}\)
    procedure(sf_i), deferred :: field_pow_int
      !! \({\rm field}^{\rm int}\)
    procedure(sf_ret_sf), deferred :: sin
      !! \(\sin({\rm field})\)
    procedure(sf_ret_sf), deferred :: cos
      !! \(\cos({\rm field})\)
    procedure(sf_ret_sf), deferred :: tan
      !! \(\tan({\rm field})\)
    procedure(sf_ret_sf), deferred :: asin
      !! \(\sin^{-1}({\rm field})\)
    procedure(sf_ret_sf), deferred :: acos
      !! \(\cos^{-1}({\rm field})\)
    procedure(sf_ret_sf), deferred :: atan
      !! \(\tan^{-1}({\rm field})\)
    procedure(sf_ret_sf), deferred :: sinh
      !! \(\sinh({\rm field})\)
    procedure(sf_ret_sf), deferred :: cosh
      !! \(\cosh({\rm field})\)
    procedure(sf_ret_sf), deferred :: tanh
      !! \(\tanh({\rm field})\)
    procedure(sf_ret_sf), deferred :: asinh
      !! \(\sinh^{-1}({\rm field})\)
    procedure(sf_ret_sf), deferred :: acosh
      !! \(\cosh^{-1}({\rm field})\)
    procedure(sf_ret_sf), deferred :: atanh
      !! \(\tanh^{-1}({\rm field})\)
    procedure(sf_ret_sf), deferred :: log
      !! \(\ ({\rm field})\)
    procedure(sf_ret_sf), deferred :: log10
      !! \(\ ({\rm field})\)
    procedure(sf_ret_sf), deferred :: exp
      !! \(\ ({\rm field})\)
    procedure(sf_ret_sf), deferred :: abs
      !! \(\ ({\rm field})\)
    procedure(sf_ret_sf), deferred :: sqrt
      !! \(\ ({\rm field})\)
    procedure(sf_ret_sf), deferred :: minval
      !! \(\ ({\rm field})\)
    procedure(sf_ret_sf), deferred :: maxval
      !! \(\ ({\rm field})\)
    procedure(sf_dx), public, deferred :: d_dx
      !! \(\frac{\partial^n}{\partial x_i^n}({\rm field})\)
    procedure(sf_ret_vf), deferred :: gradient
      !! \(\nabla {\rm field}\)
    procedure(sf_ret_sf), deferred :: laplacian
      !! \(\nabla^2 {\rm field}\)
    procedure(sf_eq_sf), deferred :: assign_field
      !! \({\rm field} = {\rm field}\)
    procedure(sf_eq_r), deferred :: assign_real
      !! \({\rm field} = {\rm real}\)
    generic, public :: operator(*) => field_multiply_field, &
        field_multiply_vecfield, real_multiply_field, &
        field_multiply_real
    generic, public :: operator(/) => field_divide_field, &
        field_divide_real, real_divide_field
    generic, public :: operator(+) => field_add_field, field_add_real, &
        real_add_field
    generic, public :: operator(-) => field_sub_field, field_sub_real, &
        real_sub_field
    generic, public :: operator(**) => field_pow_real, field_pow_int, &
        field_pow_real4
    generic, public :: operator(.grad.) => gradient
    generic, public :: operator(.laplacian.) => laplacian
    generic, public :: assignment(=) => assign_field, assign_real
    procedure(sf_bound), public, deferred :: set_boundaries
      !! Set the type and value of boundary conditions
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
    ! ####Bibliography
    ! *Scientific Software Design: The Object-Oriented Way*, Rouson, 
    ! Damian and Xia, Jim and Xu, Xiaofeng, 2011, ISBN 9780521888134, 
    ! Cambridge University Press, New York, NY, USA.
    !
  contains
    private
    procedure(vf_sf), deferred :: field_multiply_field
      !! \({\rm\vec{field}} \times {\rm field}\)
    procedure(r_vf), pass(rhs), deferred :: real_multiply_field
      !! \({\rm real}  \times {\rm \vec{field}}\)
    procedure(vf_r), deferred :: field_multiply_real
      !! \({\rm \vec{field}} \times {\rm real}\)
    procedure(vf_vf), deferred :: field_divide_field
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
    procedure(vf_ret_vf), public, deferred :: norm
      !! \(\lVert {\rm \vec{field}} \rVert\)
    procedure(vf_comp), public, deferred :: component
      !! Returns a scalar field containing the specified component of 
      !! the vector field
    procedure(vf_dx), public, deferred :: d_dx
      !! \(\frac{\partial^n}{\partial x_i^n}({\rm \vec{field}})\)
    procedure(vf_ret_sf), deferred :: divergence
      !! \(\nabla\cdot {\rm \vec{field}}\)
    procedure(vf_ret_vf), deferred :: curl
      !! \(\nabla\times {\rm \vec{field}}\)
    procedure(vf_ret_vf), deferred :: laplacian
      !! \(\nabla^2 {\rm \vec{field}}\)
    procedure(vf_vf_ret_sf), deferred :: dot_prod
      !! \({\rm \vec{field}} \cdot {\rm \vec{field}}\)
    procedure(vf_vf), deferred :: cross_prod
      !! \({\rm\vec{field}} \times {\rm\vec{field}}\)
    procedure(vf_eq_vf), deferred :: assign_field
      !! \({\rm \vec{field}} = {\rm \vec{field}}\)
    procedure(vf_eq_vr), deferred :: assign_real
      !! \({\rm \vec{field}} = {\rm \vec{real}}\)
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
    generic, public :: operator(.dot.) => dot_prod
    generic, public :: operator(.cross.) => cross_prod
    generic, public :: assignment(=) => assign_field, assign_real
    procedure(vf_bound), public, deferred :: set_boundaries
      !! Set the type and value of boundary conditions
  end type vector_field


  abstract interface
    function f_ret_r(this)
      import :: abstract_field
      import :: r8
      class(abstract_field), intent(in) :: this
      real(r8), dimension(:,:), allocatable :: f_ret_r
        !* A 2D array of shape \(n \times 2\), where n is the number of
        !  dimensions of the field. Each row contains the lower and the
        !  upper extent of the fields domain in the corresponding 
        !  dimension
    end function f_ret_r

    pure function f_ret_i(this)
      import :: abstract_field
      class(abstract_field), intent(in) :: this
      integer :: f_ret_i
    end function f_ret_i

    function sf_sf(this,rhs)
      !! \({\rm field} [{\rm operator}] {\rm field}\)
      import :: scalar_field
      class(scalar_field), intent(in) :: this, rhs
      class(scalar_field), allocatable :: sf_sf !! The restult of this operation
    end function sf_sf

    function sf_vf(this,rhs)
      !! \({\rm field} [{\rm operator}] {\rm \vec{field}}\)
      import :: scalar_field
      import :: vector_field
      import :: r8
      class(scalar_field), intent(in) :: this
      class(vector_field), intent(in) :: rhs
      class(vector_field), allocatable :: sf_vf !! The result of this operation
    end function sf_vf

    function r_sf(lhs,rhs)
      !! \({\rm real} [{\rm operator}] {\rm field}\)
      import :: scalar_field
      import :: r8
      real(r8), intent(in) :: lhs
      class(scalar_field), intent(in) :: rhs
      class(scalar_field), allocatable :: r_sf !! The result of this operation
    end function r_sf
  
    function sf_r(this,rhs)
      !! \({\rm field} [{\rm operator}] {\rm real}\)
      import :: scalar_field
      import :: r8
      class(scalar_field), intent(in) :: this
      real(r8), intent(in) :: rhs
      class(scalar_field), allocatable :: sf_r !! The result of this operation
    end function sf_r
  
    function sf_r4(this,rhs)
      !! \({\rm field} [{\rm operator}] {\rm real}\)
      import :: scalar_field
      class(scalar_field), intent(in) :: this
      real, intent(in) :: rhs
      class(scalar_field), allocatable :: sf_r4 !! The result of this operation
    end function sf_r4
  
    function sf_i(this,rhs)
      !! \({\rm field} [{\rm operator}] {\rm real}\)
      import :: scalar_field
      class(scalar_field), intent(in) :: this
      integer, intent(in) :: rhs
      class(scalar_field), allocatable :: sf_i !! The result of this operation
    end function sf_i
  
    function sf_ret_sf(this)
      !! \([{\rm operator}] {\rm field}\)
      import :: scalar_field
      class(scalar_field), intent(in) :: this
      class(scalar_field), allocatable :: sf_ret_sf !! The result of this operation
    end function sf_ret_sf
    
    function sf_ret_vf(this)
      !! \([{\rm operator}] {\rm field}\)
      import :: scalar_field
      import :: vector_field
      class(scalar_field), intent(in) :: this
      class(vector_field), allocatable :: sf_ret_vf !! The result of this operation
    end function sf_ret_vf
    
    pure subroutine sf_eq_sf(this,rhs)
      !! \({\rm field} = {\rm field}\)
      import :: scalar_field
      class(scalar_field), intent(inout) :: this
      class(scalar_field), intent(in) :: rhs
    end subroutine sf_eq_sf
    
    pure subroutine sf_eq_r(this,rhs)
      !! \({\rm field} = {\rm real}\)
      import :: scalar_field
      import :: r8      
      class(scalar_field), intent(inout) :: this
      real(r8), intent(in) :: rhs
    end subroutine sf_eq_r
    
    function vf_sf(this,rhs)
      !! \({\rm \vec{field} [{\rm operator}] {\rm field}\)
      import :: scalar_field
      import :: vector_field
      class(vector_field), intent(in) :: this
      class(scalar_field), intent(in) :: rhs
      class(vector_field), allocatable :: vf_sf !! The restult of this operation
    end function vf_sf

    function vf_vf(this,rhs)
      !! \({\rm \vec{field}} [{\rm operator}] {\rm \vec{field}}\)
      import :: vector_field
      class(vector_field), intent(in) :: this, rhs
      class(vector_field), allocatable :: vf_vf !! The result of this operation
    end function vf_vf

    function r_vf(lhs,rhs)
      !! \({\rm real} [{\rm operator}] {\rm field}\)
      import :: vector_field
      import :: r8
      real(r8), intent(in) :: lhs
      class(vector_field), intent(in) :: rhs
      class(vector_field), allocatable :: r_vf !! The result of this operation
    end function r_vf
  
    function vf_r(this,rhs)
      !! \({\rm field} [{\rm operator}] {\rm real}\)
      import :: vector_field
      import :: r8
      class(vector_field), intent(in) :: this
      real(r8), intent(in) :: rhs
      class(vector_field), allocatable :: vf_r !! The result of this operation
    end function vf_r
    
    function vf_ret_sf(this)
      !! \([{\rm operator}] {\rm field}\)
      import :: scalar_field
      import :: vector_field
      class(vector_field), intent(in) :: this
      class(scalar_field), allocatable :: vf_ret_sf !! The result of this operation
    end function vf_ret_sf
    
    function vf_ret_vf(this)
      !! \([{\rm operator}] {\rm field}\)
      import :: vector_field
      class(vector_field), intent(in) :: this
      class(vector_field), allocatable :: vf_ret_vf !! The result of this operation
    end function vf_ret_vf
    
    function vf_vf_ret_sf(this,rhs)
      !! \({\rm \vec{field}} [{\rm operator}] {\rm \vec{field}}\)
      import :: scalar_field
      import :: vector_field
      class(vector_field), intent(in) :: this, rhs
      class(scalar_field), allocatable :: vf_vf_ret_sf !! The result of this operation
    end function vf_vf_ret_sf

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

    pure subroutine vf_eq_vf(this,rhs)
      !! \({\rm \vec{field}} = {\rm \vec{field}}\)
      import :: vector_field
      class(vector_field), intent(inout) :: this
      class(vector_field), intent(in) :: rhs
    end subroutine vf_eq_vf

    pure subroutine vf_eq_vr(this,rhs)
      !! \({\rm \vec{field}} = {\rm \vec{real}}\)
      import :: vector_field
      import :: r8
      class(vector_field), intent(inout) :: this
      real(r8), dimension(:), intent(in) :: rhs
    end subroutine vf_eq_vr

    function vr_vf(lhs,rhs)
      !! \({\rm \vec{real}} [{\rm operator}] {\rm field}\)
      import :: vector_field
      import :: r8
      real(r8), dimension(:), intent(in) :: lhs
      class(vector_field), intent(in) :: rhs
      class(vector_field), allocatable :: vr_vf !! The result of this operation
    end function vr_vf

    function vf_vr(this,rhs)
      !! \({\rm \vec{field}} [{\rm operator}] {\rm \vec{real}}\)
      import :: vector_field
      import :: r8
      class(vector_field), intent(in) :: this
      real(r8), dimension(:), intent(in) :: rhs
      class(vector_field), allocatable :: vf_vr !! The result of this operation
    end function vf_vr

    function vf_comp(this,comp)
      !! Scalar field containing specified component of the vector field
      import :: scalar_field
      import :: vector_field
      class(vector_field), intent(in) :: this
      integer, intent(in) :: comp !! The index of the component
      class(scalar_field), allocatable :: vf_comp !! Component number `comp`
    end function vf_comp

    pure subroutine sf_bound(this,bound_type,boundary,lower,value1,value2)
      !* Sets boundary conditions and values. Boundary conditions 
      !  (numbered with their corresponding integer code) are:
      !
      ! 1. Dirichlet
      ! 2. Neumann
      ! 3. Cauchy
      ! 4. Periodic
      !
      ! If set to 0 then is a free boundary.
      !
      import :: scalar_field
      import :: r8
      class(scalar_field), intent(inout) :: this
      integer, intent(in) :: bound_type
        !! Integer code for the type of boundary
      integer, intent(in) :: boundary
        !! The number corresponding to the dimension whose boundary
        !! condition is to be set
      logical, optional, intent(in) :: lower
        !! Sets lower boundary if true (default), upper if false
      real(r8), optional, intent(in) :: value1
        !! Value of the boundary condition for Dirichlet, Neumann, and
        !! Cauchy conditions. In the lattermost case, it is the value of
        !! the field itself at the boundary (not the derivative). Default
        !! is 0.
      real(r8), optional, intent(in) :: value2
        !! Value of the first derivative of the field, when specifying
        !! Cauchy boundary conditions. Default is 0.
    end subroutine sf_bound

    pure subroutine vf_bound(this,bound_type,boundary,lower,value1,value2)
      !* Sets boundary conditions and values. Boundary conditions 
      !  (numbered with their corresponding integer code) are:
      !
      ! 1. Dirichlet
      ! 2. Neumann
      ! 3. Cauchy
      ! 4. Periodic
      !
      ! If set to 0 then is a free boundary.
      !
      import :: vector_field
      import :: r8
      class(vector_field), intent(inout) :: this
      integer, intent(in) :: bound_type
        !! Integer code for the type of boundary
      integer, intent(in) :: boundary
        !! The number corresponding to the dimension whose boundary
        !! condition is to be set
      logical, optional, intent(in) :: lower
        !! Sets lower boundary if true (default), upper if false
      real(r8), dimension(:), optional, intent(in) :: value1
        !! Value of the boundary condition for Dirichlet, Neumann, and
        !! Cauchy conditions. In the lattermost case, it is the value of
        !! the field itself at the boundary (not the derivative). Default
        !! is 0.
      real(r8), dimension(:), optional, intent(in) :: value2
        !! Value of the first derivative of the field, when specifying
        !! Cauchy boundary conditions. Default is 0.
    end subroutine vf_bound
    
    function f_raw(this)
      import :: abstract_field
      import :: r8
      class(abstract_field), intent(in) :: this
      real(r8), dimension(this%raw_size()) :: f_raw
        !! Array containing data needed to describe field
    end function f_raw
    
    function f_jacob(this,order)
      !* A Jacobian matrix \(\frac{\partial}{\partial y_j}f_i(\vec{y})\)
      ! where \(\vec{y}\) is the field represented as a 1D array of 
      ! values (i.e. the output of [[abstract_field:get_raw]]) and
      ! \( f_i(\vec{y}) = \frac{\partial^{n}y_i}{\partial x^n} \). Here,
      ! \(n\) is the `order` of the derivative to be taken, with 0
      ! corresponding to not taking any derivative.
      import :: abstract_field
      import :: r8
      class(abstract_field), intent(in) :: this
      integer, intent(in), optional :: order
        !! The order of the derivative of the field whose Jacobian is
        !! to be returned. Default is 0 (no differentiation)
      real(r8), dimension(this%raw_size(),this%raw_size()) :: f_jacob
        !! The resulting Jacobian matrix
    end function f_jacob
  end interface

  interface sin
    module procedure :: scalar_field_sin
  end interface

  interface cos
    module procedure :: scalar_field_cos
  end interface

  interface tan
    module procedure :: scalar_field_tan
  end interface

  interface asin
    module procedure :: scalar_field_asin
  end interface

  interface acos
    module procedure :: scalar_field_acos
  end interface

  interface atan
    module procedure :: scalar_field_atan
  end interface

  interface sinh
    module procedure :: scalar_field_sinh
  end interface

  interface cosh
    module procedure :: scalar_field_cosh
  end interface

  interface tanh
    module procedure :: scalar_field_tanh
  end interface

  interface asinh
    module procedure :: scalar_field_asinh
  end interface

  interface acosh
    module procedure :: scalar_field_acosh
  end interface

  interface atanh
    module procedure :: scalar_field_atanh
  end interface

  interface log
    module procedure :: scalar_field_log
  end interface

  interface log10
    module procedure :: scalar_field_log10
  end interface

  interface exp
    module procedure :: scalar_field_exp
  end interface

  interface abs
    module procedure :: scalar_field_abs
  end interface

  interface sqrt
    module procedure :: scalar_field_sqrt
  end interface

  interface minval
    module procedure :: scalar_field_minval
  end interface

  interface maxval
    module procedure :: scalar_field_maxval
  end interface
  
  public :: sin, cos, tan, asin, acos, atan, sinh, cosh, tanh, asinh, &
            acosh, atanh, log, log10, exp, abs, sqrt, minval, maxval

contains
  
  function scalar_field_sin(field) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! Computes and returns the sine of the values in a scalar field.
    !
    class(scalar_field), intent(in) :: field
    class(scalar_field), allocatable :: res 
    class(scalar_field), allocatable :: tmp
    allocate(tmp, mold=field)
    tmp = field%sin()
    call move_alloc(tmp, res)
  end function scalar_field_sin

  function scalar_field_cos(field) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! Computes and returns the cosine of the values in a scalar field.
    !
    class(scalar_field), intent(in) :: field
    class(scalar_field), allocatable :: res 
    class(scalar_field), allocatable :: tmp
    allocate(tmp, mold=field)
    tmp = field%cos()
    call move_alloc(tmp, res)
  end function scalar_field_cos

  function scalar_field_tan(field) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! Computes and returns the tangent of the values in a scalar field.
    !
    class(scalar_field), intent(in) :: field
    class(scalar_field), allocatable :: res 
    class(scalar_field), allocatable :: tmp
    allocate(tmp, mold=field)
    tmp = field%tan()
    call move_alloc(tmp, res)
  end function scalar_field_tan

  function scalar_field_asin(field) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! Computes and returns the inverse sine of the values in a scalar
    ! field.
    !
    class(scalar_field), intent(in) :: field
    class(scalar_field), allocatable :: res 
    class(scalar_field), allocatable :: tmp
    allocate(tmp, mold=field)
    tmp = field%asin()
    call move_alloc(tmp, res)
  end function scalar_field_asin

  function scalar_field_acos(field) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! Computes and returns the inverse cosine of the values in a scalar
    ! field.
    !
    class(scalar_field), intent(in) :: field
    class(scalar_field), allocatable :: res 
    class(scalar_field), allocatable :: tmp
    allocate(tmp, mold=field)
    tmp = field%acos()
    call move_alloc(tmp, res)
  end function scalar_field_acos

  function scalar_field_atan(field) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! Computes and returns the inverse tangent of the values in a scalar
    ! field.
    !
    class(scalar_field), intent(in) :: field
    class(scalar_field), allocatable :: res 
    class(scalar_field), allocatable :: tmp
    allocate(tmp, mold=field)
    tmp = field%atan()
    call move_alloc(tmp, res)
  end function scalar_field_atan

  function scalar_field_sinh(field) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! Computes and returns the hyperbolic sine of the values in a scalar
    ! field.
    !
    class(scalar_field), intent(in) :: field
    class(scalar_field), allocatable :: res 
    class(scalar_field), allocatable :: tmp
    allocate(tmp, mold=field)
    tmp = field%sinh()
    call move_alloc(tmp, res)
  end function scalar_field_sinh

  function scalar_field_cosh(field) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! Computes and returns the hyperbolic cosine of the values in a
    ! scalar field.
    !
    class(scalar_field), intent(in) :: field
    class(scalar_field), allocatable :: res 
    class(scalar_field), allocatable :: tmp
    allocate(tmp, mold=field)
    tmp = field%cosh()
    call move_alloc(tmp, res)
  end function scalar_field_cosh

  function scalar_field_tanh(field) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! Computes and returns the hyperbolic tangent of the values in a
    ! scalar field.
    !
    class(scalar_field), intent(in) :: field
    class(scalar_field), allocatable :: res 
    class(scalar_field), allocatable :: tmp
    allocate(tmp, mold=field)
    tmp = field%tanh()
    call move_alloc(tmp, res)
  end function scalar_field_tanh

  function scalar_field_asinh(field) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! Computes and returns the inverse byperbolic sine of the values in
    ! a scalar field.
    !
    class(scalar_field), intent(in) :: field
    class(scalar_field), allocatable :: res 
    class(scalar_field), allocatable :: tmp
    allocate(tmp, mold=field)
    tmp = field%asinh()
    call move_alloc(tmp, res)
  end function scalar_field_asinh

  function scalar_field_acosh(field) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! Computes and returns the inverse hyperbolic cosine of the values
    ! in a scalar field.
    !
    class(scalar_field), intent(in) :: field
    class(scalar_field), allocatable :: res 
    class(scalar_field), allocatable :: tmp
    allocate(tmp, mold=field)
    tmp = field%acosh()
    call move_alloc(tmp, res)
  end function scalar_field_acosh

  function scalar_field_atanh(field) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! Computes and returns the inverse hyperbolic tangent of the values
    ! in a scalar field.
    !
    class(scalar_field), intent(in) :: field
    class(scalar_field), allocatable :: res 
    class(scalar_field), allocatable :: tmp
    allocate(tmp, mold=field)
    tmp = field%atanh()
    call move_alloc(tmp, res)
  end function scalar_field_atanh

  function scalar_field_log(field) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! Computes and returns the natural logarithm of the values in a
    ! scalar field.
    !
    class(scalar_field), intent(in) :: field
    class(scalar_field), allocatable :: res 
    class(scalar_field), allocatable :: tmp
    allocate(tmp, mold=field)
    tmp = field%log()
    call move_alloc(tmp, res)
  end function scalar_field_log

  function scalar_field_log10(field) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! Computes and returns the base-10 logarithm of the values in a
    ! scalar field.
    !
    class(scalar_field), intent(in) :: field
    class(scalar_field), allocatable :: res 
    class(scalar_field), allocatable :: tmp
    allocate(tmp, mold=field)
    tmp = field%log10()
    call move_alloc(tmp, res)
  end function scalar_field_log10

  function scalar_field_exp(field) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! Computes and returns the exponent (base e) of the values in a
    ! scalar field.
    !
    class(scalar_field), intent(in) :: field
    class(scalar_field), allocatable :: res 
    class(scalar_field), allocatable :: tmp
    allocate(tmp, mold=field)
    tmp = field%exp()
    call move_alloc(tmp, res)
  end function scalar_field_exp

  function scalar_field_abs(field) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! Computes and returns the absolute value of the values in a scalar
    ! field.
    !
    class(scalar_field), intent(in) :: field
    class(scalar_field), allocatable :: res 
    class(scalar_field), allocatable :: tmp
    allocate(tmp, mold=field)
    tmp = field%abs()
    call move_alloc(tmp, res)
  end function scalar_field_abs

  function scalar_field_sqrt(field) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! Computes and returns the squar root of the values in a scalar
    ! field.
    !
    class(scalar_field), intent(in) :: field
    class(scalar_field), allocatable :: res 
    class(scalar_field), allocatable :: tmp
    allocate(tmp, mold=field)
    tmp = field%sqrt()
    call move_alloc(tmp, res)
  end function scalar_field_sqrt

  function scalar_field_minval(field) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! Computes and returns the minimum of the values in a scalar field.
    !
    class(scalar_field), intent(in) :: field
    class(scalar_field), allocatable :: res 
    class(scalar_field), allocatable :: tmp
    allocate(tmp, mold=field)
    tmp = field%minval()
    call move_alloc(tmp, res)
  end function scalar_field_minval

  function scalar_field_maxval(field) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! Computes and returns the maximum of the values in a scalar field.
    !
    class(scalar_field), intent(in) :: field
    class(scalar_field), allocatable :: res 
    class(scalar_field), allocatable :: tmp
    allocate(tmp, mold=field)
    tmp = field%maxval()
    call move_alloc(tmp, res)
  end function scalar_field_maxval

end module abstract_fields_mod
