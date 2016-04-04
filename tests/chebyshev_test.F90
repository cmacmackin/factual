!
!  chebyshev_test.f90
!  This file is part of FIAT.
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

module chebyshev_test
  use iso_fortran_env, only: r8 => real64
  use pfunit_mod
  use chebyshev_mod
  implicit none
  
!@TestParameter
  type, extends(abstracttestparameter) :: cheb_case
    integer :: nodes
    real(r8) :: lower, upper
    real(r8), dimension(:), allocatable :: expected
  contains
    procedure :: tostring
  end type cheb_case

!@TestCase(testParameters={get_params()}, constructor=new_test)
  type, extends(parameterizedtestcase) :: test_cheb
    integer :: nodes
    real(r8) :: lower, upper
    real(r8), dimension(:), allocatable :: expected
  end type test_cheb

contains

  type(cheb_case) function constructor(nodes,lower,upper,expected)
    integer, intent(in) :: nodes
    real(r8), intent(in), optional :: lower, upper
    real(r8), dimension(:), intent(in) :: expected
    constructor%nodes = nodes
    if (present(lower)) then
      constructor%lower = lower
    else
      constructor%lower = -1.0
    end if
    if (present(upper)) then
      constructor%upper = upper
    else
      constructor%upper = -1.0
    end if
    constructor%expected = expected
  end function constructor
  
  type(test_cheb) function new_test(param)
    type(cheb_case), intent(in) :: param
    new_test%nodes = param%nodes
    new_test%lower = param%lower
    new_test%upper = param%upper
    new_test%expected = param%expected
  end function new_test
  
  function get_params() result(params)
    type(cheb_case), allocatable, dimension(:) :: params
    params = [ constructor(3,-1._r8,1._r8,[1._r8,0.5_r8,-0.5_r8,-1._r8]), &
               constructor(3,-1._r8,2._r8,[2._r8,1.25_r8,-0.25_r8,-1._r8]), &
               constructor(3,0._r8,1._r8,[1._r8,0.75_r8,0.25_r8,0._r8]), &
               constructor(3,-2._r8,2._r8,[2._r8,1._r8,-1._r8,-2._r8]) ]
  end function get_params
  
  function tostring(this) result(string)
    class(cheb_case), intent(in) :: this
    character(:), allocatable :: string
    character(len=54), parameter :: format_str = &
        "(g0,' chebyshev collocation nodes from ',g0,' to ',g0)"
    allocate(character(len=120) :: string)
    write(string,format_str) this%nodes, this%lower, this%upper
  end function tostring
  
!@Test
  subroutine test_collocation(this)
    class(test_cheb), intent(inout) :: this
    real(r8), dimension(:), allocatable :: points
    if (this%lower==-1.0) then
      if (this%upper==1.0) then
        points = collocation_points(this%nodes)
      else
        points = collocation_points(this%nodes,upper_bound=this%upper)
      end if
    else
      if (this%upper==1.0) then
        points = collocation_points(this%nodes,this%lower)
      else
        points = collocation_points(this%nodes,this%lower,this%upper)
      end if
    end if
#line 107 "chebyshev_test.pf"
  call assertEqual(this%nodes+1,size(points),message="Collocation points array wrong length.", &
 & location=SourceLocation( &
 & 'chebyshev_test.pf', &
 & 107) )
  if (anyExceptions()) return
#line 108 "chebyshev_test.pf"
#line 108 "chebyshev_test.pf"
  call assertEqual(this%expected,points,message="Calculated Chebyshev collocation points incorrect.",tolerance=1.e-9_r8, &
 & location=SourceLocation( &
 & 'chebyshev_test.pf', &
 & 108) )
  if (anyExceptions()) return
#line 109 "chebyshev_test.pf"
  end subroutine test_collocation
  
end module chebyshev_test

module Wrapchebyshev_test
   use pFUnit_mod
   use chebyshev_test
   implicit none
   private

   public :: WrapUserTestCase
   public :: makeCustomTest
   type, extends(test_cheb) :: WrapUserTestCase
      procedure(userTestMethod), nopass, pointer :: testMethodPtr
   contains
      procedure :: runMethod
   end type WrapUserTestCase

   abstract interface
     subroutine userTestMethod(this)
        use chebyshev_test
        class (test_cheb), intent(inout) :: this
     end subroutine userTestMethod
   end interface

contains

   subroutine runMethod(this)
      class (WrapUserTestCase), intent(inout) :: this

      call this%testMethodPtr(this)
   end subroutine runMethod

   function makeCustomTest(methodName, testMethod, testParameter) result(aTest)
#ifdef INTEL_13
      use pfunit_mod, only: testCase
#endif
      type (WrapUserTestCase) :: aTest
#ifdef INTEL_13
      target :: aTest
      class (WrapUserTestCase), pointer :: p
#endif
      character(len=*), intent(in) :: methodName
      procedure(userTestMethod) :: testMethod
      type (cheb_case), intent(in) :: testParameter
      aTest%test_cheb = new_test(testParameter)

      aTest%testMethodPtr => testMethod
#ifdef INTEL_13
      p => aTest
      call p%setName(methodName)
#else
      call aTest%setName(methodName)
#endif
      call aTest%setTestParameter(testParameter)
   end function makeCustomTest

end module Wrapchebyshev_test

function chebyshev_test_suite() result(suite)
   use pFUnit_mod
   use chebyshev_test
   use Wrapchebyshev_test
   type (TestSuite) :: suite

   integer, allocatable :: npes(:)

   type (cheb_case), allocatable :: testParameters(:)
   type (cheb_case) :: testParameter
   integer :: iParam 
   integer, allocatable :: cases(:) 
 
   suite = newTestSuite('chebyshev_test_suite')

   testParameters = get_params()

   do iParam = 1, size(testParameters)
      testParameter = testParameters(iParam)
   call suite%addTest(makeCustomTest('test_collocation', test_collocation, testParameter))
   end do


end function chebyshev_test_suite

