!
!  cheb2d_fields.f90
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

module cheb1d_fields_mod
  !* Author: Chris MacMackin
  !  Date: March 2016
  !  License: LGPLv3
  !
  ! Provides concrete implementations of the [[scalar_field]] and
  ! [[vector_field]] types. This implementation is for 1D fields (with
  ! 2D vectors) only. It tracks the values of the fields at Chebyshev
  ! collocation nodes and uses a pseudospectral approach to 
  ! differentiate. It may be useful for vertically-integrated, symmetric
  ! models.
  !
  use iso_fortran_env, only: r8 => real64
  use abstract_fields_mod
  implicit none
  private

  type, extends(scalar_field), public :: cheb1d_scalar_field
    private
    integer :: numpoints
      !! The number of datapoints used
    real(r8), dimension(1,2) :: extent = reshape([0.0_r8, 0.0_r8],[1,2])
      !! The start and end values of the domain
    real(r8), dimension(:), allocatable :: field_data
    integer :: lower_boundary
      !! The boundary-value code for the boundary at the start of the domain
    real(r8) :: lower_boundar_value1
      !! Value of the lower boundary condition for Dirichlet, Neumann,
      !! and Cauchy conditions. In the lattermost case, it is the value
      !! of the field itself at the boundary (not the derivative).
    real(r8) :: lower_boundar_value2
      !! Lower boundary value of the first derivative of the field, when
      !! specifying Cauchy boundary conditions. Default is 0.
    integer :: upper_boundary
      !! The boundary-value code for the boundary at the end of the domain
    real(r8) :: upper_boundar_value1
      !! Value of the upper boundary condition for Dirichlet, Neumann,
      !! and Cauchy conditions. In the lattermost case, it is the value
      !! of the field itself at the boundary (not the derivative).
    real(r8) :: upper_boundar_value2
      !! Upper boundary value of the first derivative of the field, when
      !! specifying Cauchy boundary conditions. Default is 0.
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
    procedure, public :: d_dself => cheb1d_scalar_d_dshelf
      !! Returns Jacobian for desired x-derivative of each value in
      !! field relative to all other values in field.
    procedure, public :: resolution => cheb1d_scalar_resolution
      !! Returns array containing number of datapoints in each dimension.
    procedure :: assign_raw => cheb1d_scalar_assign_raw
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
      !! \(\ ({\rm field})\)
    procedure :: log10 => cheb1d_scalar_log10
      !! \(\ ({\rm field})\)
    procedure :: exp => cheb1d_scalar_exp
      !! \(\ ({\rm field})\)
    procedure :: abs => cheb1d_scalar_abs
      !! \(\ ({\rm field})\)
    procedure :: sqrt => cheb1d_scalar_sqrt
      !! \(\ ({\rm field})\)
    procedure :: minval => cheb1d_scalar_minval
      !! \(\ ({\rm field})\)
    procedure :: maxval => cheb1d_scalar_maxval
      !! \(\ ({\rm field})\)
    procedure, public :: d_dx => cheb1d_scalar_d_dx
      !! \(\frac{\partial^n}{\partial x_i^n}({\rm field})\)
    procedure :: gradient => cheb1d_scalar_gradient
      !! \(\nabla {\rm field}\)
    procedure :: laplacian => cheb1d_scalar_laplacian
      !! \(\nabla^2 {\rm field}\)
    procedure :: assign_field => cheb1d_scalar_assign
      !! \({\rm field} = {\rm field}\)
    procedure, public :: set_boundaries => cheb1d_scalar_set_boundaries
      !! Set the type and value of boundary conditions    
  end type

contains

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

  pure function cheb1d_scalar_dimensions(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! Accessor for the number of dimensions.
    !
    class(cheb1d_scalar_field), intent(in) :: this
    integer :: res
    res = 1
  end function cheb1d_scalar_dimensions

  pure function cheb1d_scalar_raw_size(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! Compute how many elements are in the raw representation of this
    ! field. This would be the number of data points, adjusted based on
    ! how boundary conditions are accounted for.
    !
    class(cheb1d_scalar_field), intent(in) :: this
    integer :: res
    res = this%numpoints ! Will need to adjust this based on boundary conditions
  end function cheb1d_scalar_raw_size
  
  pure function cheb1d_scalar_raw(this) result(res)
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
    integer :: i
    i = max(1,this%raw_size())
    allocate(res(i))
    res(1:this%numpoints) = this%field_data ! Will need to adjust this based on boundary conditions
  end function cheb1d_scalar_raw
  
  pure function cheb1d_scalar_d_dshelf(this,order) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! A Jacobian matrix \(\frac{\partial}{\partial y_j}f_i(\vec{y})\)
    ! is returned, where \(\vec{y}\) is the field represented as a 1D
    ! array of values (i.e. the output of [[cheb1d_scalar_field:raw]])
    ! and \( f_i(\vec{y}) = \frac{\partial^{n}y_i}{\partial x^n} \).
    ! Here, \(n\) is the `order` of the derivative to be taken, with 0
    ! corresponding to not taking any derivative.
    !
    ! @BUG The returned value has shape `(this%raw_size(),
    ! this%raw_size())`, but a bug in gfortran 4.8 (fixed by version
    ! 5) caused the compiler to segfault if it was declared as such.
    ! As a workaround, it is allocatable isntead.
    !
    class(cheb1d_scalar_field), intent(in) :: this
    integer, intent(in), optional :: order
      !! The order of the derivative of the field whose Jacobian is
      !! to be returned. Default is 0 (no differentiation)
    real(r8), dimension(:,:), allocatable :: res
      !! The resulting Jacobian matrix
    integer :: i, j
    i = max(1,this%raw_size())
    allocate(res(i,i))
    res = 0.0_r8
    forall (j = 1:i) res(j,j) = 1.0_r8 ! This only works for no derivative and if nothing going on at boundaries. Will need to adjust it
  end function cheb1d_scalar_d_dshelf

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

  subroutine cheb1d_scalar_assign_raw(this,rhs)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! Assigns raw data, such as that produced by 
    ! [[cheb1d_scalar_field:raw]], to the field. The routine will
    ! stop with an error if the array is the wrong length for a field
    ! of this resolution and with these boundary conditions. 
    !
    class(cheb1d_scalar_field), intent(inout) :: this
    real(r8), dimension(:), intent(in) :: rhs
    ! Will have to do something here, once I've figured out how to 
    ! structure the raw data.
  end subroutine cheb1d_scalar_assign_raw
  
  function cheb1d_scalar_sf_m_sf(this,rhs) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm field} \times {\rm field}\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), intent(in) :: rhs
    class(scalar_field), allocatable :: res !! The restult of this operation
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
  end function cheb1d_scalar_sf_p_i

  pure function cheb1d_scalar_sin(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\sin({\rm field})\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: tmp
    allocate(tmp)
    if (allocated(this%field_data)) tmp%field_data = sin(this%field_data)
    call move_alloc(tmp, res)
  end function cheb1d_scalar_sin

  pure function cheb1d_scalar_cos(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\cos({\rm field})\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: tmp
    allocate(tmp)
    if (allocated(this%field_data)) tmp%field_data = cos(this%field_data)
    call move_alloc(tmp, res)
  end function cheb1d_scalar_cos

  pure function cheb1d_scalar_tan(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\tan({\rm field})\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: tmp
    allocate(tmp)
    if (allocated(this%field_data)) tmp%field_data = tan(this%field_data)
    call move_alloc(tmp, res)
  end function cheb1d_scalar_tan

  pure function cheb1d_scalar_asin(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\sin^{-1}({\rm field})\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: tmp
    allocate(tmp)
    if (allocated(this%field_data)) tmp%field_data = asin(this%field_data)
    call move_alloc(tmp, res)
  end function cheb1d_scalar_asin

  pure function cheb1d_scalar_acos(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\cos^{-1}({\rm field})\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: tmp
    allocate(tmp)
    if (allocated(this%field_data)) tmp%field_data = acos(this%field_data)
    call move_alloc(tmp, res)
  end function cheb1d_scalar_acos

  pure function cheb1d_scalar_atan(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\tan^{-1}({\rm field})\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: tmp
    allocate(tmp)
    if (allocated(this%field_data)) tmp%field_data = atan(this%field_data)
    call move_alloc(tmp, res)
  end function cheb1d_scalar_atan

  pure function cheb1d_scalar_sinh(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\sinh({\rm field})\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: tmp
    allocate(tmp)
    if (allocated(this%field_data)) tmp%field_data = sinh(this%field_data)
    call move_alloc(tmp, res)
  end function cheb1d_scalar_sinh

  pure function cheb1d_scalar_cosh(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\cosh({\rm field})\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: tmp
    allocate(tmp)
    if (allocated(this%field_data)) tmp%field_data = cosh(this%field_data)
    call move_alloc(tmp, res)
  end function cheb1d_scalar_cosh

  pure function cheb1d_scalar_tanh(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\tanh({\rm field})\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: tmp
    allocate(tmp)
    if (allocated(this%field_data)) tmp%field_data = tanh(this%field_data)
    call move_alloc(tmp, res)
  end function cheb1d_scalar_tanh

  pure function cheb1d_scalar_asinh(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\sinh^{-1}({\rm field})\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: tmp
    allocate(tmp)
    if (allocated(this%field_data)) tmp%field_data = asinh(this%field_data)
    call move_alloc(tmp, res)
  end function cheb1d_scalar_asinh

  pure function cheb1d_scalar_acosh(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\cosh^{-1}({\rm field})\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: tmp
    allocate(tmp)
    if (allocated(this%field_data)) tmp%field_data = acosh(this%field_data)
    call move_alloc(tmp, res)
  end function cheb1d_scalar_acosh

  pure function cheb1d_scalar_atanh(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\tanh^{-1}({\rm field})\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: tmp
    allocate(tmp)
    if (allocated(this%field_data)) tmp%field_data = atanh(this%field_data)
    call move_alloc(tmp, res)
  end function cheb1d_scalar_atanh

  pure function cheb1d_scalar_log(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\log({\rm field})\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: tmp
    allocate(tmp)
    if (allocated(this%field_data)) tmp%field_data = log(this%field_data)
    call move_alloc(tmp, res)
  end function cheb1d_scalar_log

  pure function cheb1d_scalar_log10(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\log10({\rm field})\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: tmp
    allocate(tmp)
    if (allocated(this%field_data)) tmp%field_data = log10(this%field_data)
    call move_alloc(tmp, res)
  end function cheb1d_scalar_log10

  pure function cheb1d_scalar_exp(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(e^({\rm field})\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: tmp
    allocate(tmp)
    if (allocated(this%field_data)) tmp%field_data = exp(this%field_data)
    call move_alloc(tmp, res)
  end function cheb1d_scalar_exp

  pure function cheb1d_scalar_abs(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\abs({\rm field})\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: tmp
    allocate(tmp)
    if (allocated(this%field_data)) tmp%field_data = abs(this%field_data)
    call move_alloc(tmp, res)
  end function cheb1d_scalar_abs

  pure function cheb1d_scalar_sqrt(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\sqrt{{\rm field}}\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: tmp
    allocate(tmp)
    if (allocated(this%field_data)) tmp%field_data = sqrt(this%field_data)
    call move_alloc(tmp, res)
  end function cheb1d_scalar_sqrt

  pure function cheb1d_scalar_minval(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\min({\rm field})\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: tmp
    allocate(tmp)
    if (allocated(this%field_data)) tmp%field_data = minval(this%field_data)
    call move_alloc(tmp, res)
  end function cheb1d_scalar_minval

  pure function cheb1d_scalar_maxval(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\max({\rm field})\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    type(cheb1d_scalar_field), allocatable :: tmp
    allocate(tmp)
    if (allocated(this%field_data)) tmp%field_data = maxval(this%field_data)
    call move_alloc(tmp, res)
  end function cheb1d_scalar_maxval
  
  pure function cheb1d_scalar_laplacian(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\nabla^2 {\rm field}\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
  end function cheb1d_scalar_laplacian
  
  pure function cheb1d_scalar_gradient(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(\nabla{\rm field}\)
    !
    class(cheb1d_scalar_field), intent(in) :: this
    class(vector_field), allocatable :: res !! The result of this operation
  end function cheb1d_scalar_gradient
  
  pure subroutine cheb1d_scalar_assign(this,rhs)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \({\rm field} = {\rm field}\)
    !
    class(cheb1d_scalar_field), intent(inout) :: this
    class(scalar_field), intent(in) :: rhs
  end subroutine cheb1d_scalar_assign

  subroutine cheb1d_scalar_set_boundaries(this,bound_type,boundary, &
                                          lower,value1,value2)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! Sets boundary conditions and values. Boundary conditions 
    ! (numbered with their corresponding integer code) are:
    !
    ! 1. Dirichlet
    ! 2. Neumann
    ! 3. Cauchy
    ! 4. Periodic
    !
    ! If set to 0 then is a free boundary.
    !
    class(cheb1d_scalar_field), intent(inout) :: this
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
  end subroutine cheb1d_scalar_set_boundaries

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
  end function cheb1d_scalar_d_dx

end module cheb1d_fields_mod
