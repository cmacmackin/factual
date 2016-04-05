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

module cheb1d_scalar_test
  use iso_fortran_env, only: r8 => real64
  use pfunit_mod
  use chebyshev_mod
  use cheb1d_fields_mod
  implicit none
  
  integer, parameter :: length = 20
  
!@TestCase
  type, extends(testcase) :: test_field
    type(cheb1d_scalar_field) :: field1, field2
    integer :: nodes = length
    real(r8) :: lower1 = -1.0_r8, upper1 = 1.0_r8
    real(r8) :: lower2 = -1.0_r8, upper2 = 1.0_r8
  contains
    procedure :: setup
  end type test_field

contains

  subroutine setup(this)
    class(test_field), intent(inout) :: this
    this%field1 = cheb1d_scalar_field(this%nodes,f1)
    this%field2 = cheb1d_scalar_field(this%nodes,f2,this%lower2,this%upper2)
  end subroutine setup

  pure real(r8) function f1(x)
    real(r8), intent(in) :: x
    f1 = sin(x)
  end function f1
  
  pure real(r8) function df1(x)
    real(r8), intent(in) :: x
    df1 = cos(x)
  end function df1
  
  pure real(r8) function d2f1(x)
    real(r8), intent(in) :: x
    d2f1 = -sin(x)
  end function d2f1
  
  pure real(r8) function f2(x)
    real(r8), intent(in) :: x
    f2 = 1.0_r8 + x*exp(x)
  end function f2
  
  pure real(r8) function df2(x)
    real(r8), intent(in) :: x
    df2 = exp(x) + x*exp(x)
  end function df2

  pure real(r8) function d2f2(x)
    real(r8), intent(in) :: x
    d2f2 = 2._r8*exp(x) + x*exp(x)
  end function d2f2
  
!@Test
  subroutine test_domain(this)
    class(test_field), intent(inout) :: this
    real(r8), dimension(2) :: dom
    dom = pack(this%field1%domain(),.true.)
#line 85 "cheb1d_scalar_test.pf"
  call assertEqual([this%lower1,this%upper1],dom,message='Incorrect domain returned.', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 85) )
  if (anyExceptions()) return
#line 86 "cheb1d_scalar_test.pf"
    dom = pack(this%field2%domain(),.true.)
#line 87 "cheb1d_scalar_test.pf"
  call assertEqual([this%lower2,this%upper2],dom,message='Incorrect domain returned.', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 87) )
  if (anyExceptions()) return
#line 88 "cheb1d_scalar_test.pf"
    this%field1 = cheb1d_scalar_field(this%nodes,f1,-20._r8,1.5_r8)
    dom = pack(this%field1%domain(),.true.)
#line 90 "cheb1d_scalar_test.pf"
  call assertEqual([-20._r8,1.5_r8],dom,message='Incorrect domain returned.', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 90) )
  if (anyExceptions()) return
#line 91 "cheb1d_scalar_test.pf"
  end subroutine test_domain
  
!@Test
  subroutine test_dimensions(this)
    class(test_field), intent(inout) :: this
    integer, dimension(1) :: dims
    dims = this%field1%dimensions()
#line 98 "cheb1d_scalar_test.pf"
  call assertEqual(1,dims(1),message='Incorrect number of dimensions returned.', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 98) )
  if (anyExceptions()) return
#line 99 "cheb1d_scalar_test.pf"
    dims = this%field2%dimensions()
#line 100 "cheb1d_scalar_test.pf"
  call assertEqual(1,dims(1),message='Incorrect number of dimensions returned.', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 100) )
  if (anyExceptions()) return
#line 101 "cheb1d_scalar_test.pf"
  end subroutine test_dimensions

!@Test
  subroutine test_raw_size(this)
    class(test_field), intent(inout) :: this
#line 106 "cheb1d_scalar_test.pf"
  call assertEqual(this%nodes+1,this%field1%raw_size(),message='Wrong raw size returned.', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 106) )
  if (anyExceptions()) return
#line 107 "cheb1d_scalar_test.pf"
#line 107 "cheb1d_scalar_test.pf"
  call assertEqual(this%nodes,this%field1%raw_size([.true.],[.false.]),message='Wrong raw size returned.', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 107) )
  if (anyExceptions()) return
#line 108 "cheb1d_scalar_test.pf"
#line 108 "cheb1d_scalar_test.pf"
  call assertEqual(this%nodes,this%field1%raw_size([.false.],[.true.]),message='Wrong raw size returned.', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 108) )
  if (anyExceptions()) return
#line 109 "cheb1d_scalar_test.pf"
#line 109 "cheb1d_scalar_test.pf"
  call assertEqual(this%nodes-1,this%field1%raw_size([.false.],[.false.]),message='Wrong raw size returned.', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 109) )
  if (anyExceptions()) return
#line 110 "cheb1d_scalar_test.pf"
  end subroutine test_raw_size

!@Test
  subroutine test_raw(this)
    class(test_field), intent(inout) :: this
    real(r8), allocatable, dimension(:) :: raw, expected
    integer :: i
    raw = this%field1%raw()
    expected = collocation_points(this%nodes,this%lower1,this%upper1)
    forall (i = 1:this%nodes+1) expected(i) = f1(expected(i))
#line 120 "cheb1d_scalar_test.pf"
  call assertEqual(expected,raw,message='Wrong raw contents returned', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 120) )
  if (anyExceptions()) return
#line 121 "cheb1d_scalar_test.pf"
    raw = this%field1%raw([.false.])
#line 122 "cheb1d_scalar_test.pf"
  call assertEqual(expected(:this%nodes),raw,message='Wrong raw contents returned', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 122) )
  if (anyExceptions()) return
#line 123 "cheb1d_scalar_test.pf"
    raw = this%field1%raw([.true.],[.false.])
#line 124 "cheb1d_scalar_test.pf"
  call assertEqual(expected(2:),raw,message='Wrong raw contents returned', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 124) )
  if (anyExceptions()) return
#line 125 "cheb1d_scalar_test.pf"
    raw = this%field1%raw([.false.],[.false.])
#line 126 "cheb1d_scalar_test.pf"
  call assertEqual(expected(2:this%nodes),raw,message='Wrong raw contents returned', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 126) )
  if (anyExceptions()) return
#line 127 "cheb1d_scalar_test.pf"
  end subroutine test_raw

!@Test
  subroutine test_assign_raw(this)
    class(test_field), intent(inout) :: this
    real(r8), dimension(this%nodes+1) :: raw, xvals
    integer :: i
    xvals = collocation_points(this%nodes,this%lower1,this%upper1)
    forall (i = 1:this%nodes+1) raw(i) = f2(xvals(i))
    call this%field1%set_from_raw(raw)
#line 137 "cheb1d_scalar_test.pf"
  call assertEqual(raw,this%field2%raw(),message='Error setting field from raw.', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 137) )
  if (anyExceptions()) return
#line 138 "cheb1d_scalar_test.pf"
    forall (i = 1:this%nodes+1) raw(i) = f1(xvals(i))
    call this%field2%set_from_raw(raw(2:),[.false.],[.true.])
#line 140 "cheb1d_scalar_test.pf"
  call assertEqual(raw(2:),this%field2%raw([.false.],[.true.]), &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 140) )
  if (anyExceptions()) return
#line 141 "cheb1d_scalar_test.pf"
    forall (i = 1:this%nodes+1) raw(i) = f2(xvals(i))
    call this%field2%set_from_raw(raw(:this%nodes),[.true.],[.false.])
#line 143 "cheb1d_scalar_test.pf"
  call assertEqual(raw(:this%nodes),this%field2%raw([.true.],[.false.]),message='Error setting field from raw.', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 143) )
  if (anyExceptions()) return
#line 144 "cheb1d_scalar_test.pf"
    forall (i = 1:this%nodes+1) raw(i) = f1(xvals(i))
    call this%field2%set_from_raw(raw(2:this%nodes),[.false.],[.false.])
#line 146 "cheb1d_scalar_test.pf"
  call assertEqual(raw(2:this%nodes),this%field2%raw([.false.],[.false.]),message='Error setting field from raw.', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 146) )
  if (anyExceptions()) return
#line 147 "cheb1d_scalar_test.pf"
    forall (i = 1:this%nodes+1) raw(i) = f1(xvals(i))
    call this%field2%set_from_raw(raw)
#line 149 "cheb1d_scalar_test.pf"
  call assertEqual(raw,this%field2%raw(),message='Error setting field from raw', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 149) )
  if (anyExceptions()) return
#line 150 "cheb1d_scalar_test.pf"
  end subroutine test_assign_raw

!@Test
  subroutine test_resolution(this)
    class(test_field), intent(inout) :: this
    integer, dimension(1) :: res
    res = this%field1%resolution()
#line 157 "cheb1d_scalar_test.pf"
  call assertEqual(1,size(res),message='Resolution vector has incorrect size', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 157) )
  if (anyExceptions()) return
#line 158 "cheb1d_scalar_test.pf"
#line 158 "cheb1d_scalar_test.pf"
  call assertEqual(this%nodes,res(1),message='Incorrect resolution returned', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 158) )
  if (anyExceptions()) return
#line 159 "cheb1d_scalar_test.pf"
  end subroutine test_resolution

!@Test
  subroutine test_fmf(this)
    class(test_field), intent(inout) :: this
    this%field1 = this%field1 * this%field2
    this%field2 = cheb1d_scalar_field(this%nodes,func)
#line 166 "cheb1d_scalar_test.pf"
  call assertTrue(this%field2==this%field1,message='Error calculating field*field', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 166) )
  if (anyExceptions()) return
#line 167 "cheb1d_scalar_test.pf"
  contains
    real(r8) pure function func(x)
      real(r8), intent(in) :: x
      func = f1(x) * f2(x)
    end function func
  end subroutine test_fmf
  
!~ @Test
!~   subroutine test_fmvf(this)
!~     class(test_field), intent(inout) :: this
!~     ! TODO: write implementation
!~   end subroutine test_fmvf
  
!@Test
  subroutine test_rmf(this)
    class(test_field), intent(inout) :: this
    this%field1 = 2.0_r8 * this%field1
    this%field2 = cheb1d_scalar_field(this%nodes,func)
#line 185 "cheb1d_scalar_test.pf"
  call assertTrue(this%field2==this%field1,message='Error calculating real*field', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 185) )
  if (anyExceptions()) return
#line 186 "cheb1d_scalar_test.pf"
  contains
    real(r8) pure function func(x)
      real(r8), intent(in) :: x
      func = f1(x) * 2.0_r8
    end function func
  end subroutine test_rmf
  
!@Test
  subroutine test_fmr(this)
    class(test_field), intent(inout) :: this
    this%field1 = this%field1 * 2._r8
    this%field2 = cheb1d_scalar_field(this%nodes,func)
#line 198 "cheb1d_scalar_test.pf"
  call assertTrue(this%field2==this%field1,message='Error calculating field*real', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 198) )
  if (anyExceptions()) return
#line 199 "cheb1d_scalar_test.pf"
  contains
    real(r8) pure function func(x)
      real(r8), intent(in) :: x
      func = f1(x) * 2.0_r8
    end function func
  end subroutine test_fmr
  
!@Test
  subroutine test_fdf(this)
    class(test_field), intent(inout) :: this
    this%field1 = this%field1 / this%field2
    this%field2 = cheb1d_scalar_field(this%nodes,func)
#line 211 "cheb1d_scalar_test.pf"
  call assertTrue(this%field2==this%field1,message='Error calculating field/field', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 211) )
  if (anyExceptions()) return
#line 212 "cheb1d_scalar_test.pf"
  contains
    real(r8) pure function func(x)
      real(r8), intent(in) :: x
      func = f1(x) / f2(x)
    end function func
  end subroutine test_fdf
  
!@Test
  subroutine test_rdf(this)
    class(test_field), intent(inout) :: this
    this%field1 = 5.0_r8 / this%field2
    this%field2 = cheb1d_scalar_field(this%nodes,func)
#line 224 "cheb1d_scalar_test.pf"
  call assertTrue(this%field2==this%field1,message='Error calculating real/field', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 224) )
  if (anyExceptions()) return
#line 225 "cheb1d_scalar_test.pf"
  contains
    real(r8) pure function func(x)
      real(r8), intent(in) :: x
      func = 5.0_r8 / f2(x)
    end function func
  end subroutine test_rdf
  
!@Test
  subroutine test_fdr(this)
    class(test_field), intent(inout) :: this
    this%field1 = this%field1 / 5.0_r8
    this%field2 = cheb1d_scalar_field(this%nodes,func)
#line 237 "cheb1d_scalar_test.pf"
  call assertTrue(this%field2==this%field1,message='Error calculating field/real', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 237) )
  if (anyExceptions()) return
#line 238 "cheb1d_scalar_test.pf"
  contains
    real(r8) pure function func(x)
      real(r8), intent(in) :: x
      func = f1(x)/5.0_r8
    end function func
  end subroutine test_fdr
  
!@Test
  subroutine test_faf(this)
    class(test_field), intent(inout) :: this
    this%field1 = this%field1 + this%field2
    this%field2 = cheb1d_scalar_field(this%nodes,func)
#line 250 "cheb1d_scalar_test.pf"
  call assertTrue(this%field2==this%field1,message='Error calculating field+field', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 250) )
  if (anyExceptions()) return
#line 251 "cheb1d_scalar_test.pf"
  contains
    real(r8) pure function func(x)
      real(r8), intent(in) :: x
      func = f1(x) + f2(x)
    end function func
  end subroutine test_faf
  
!@Test
  subroutine test_raf(this)
    class(test_field), intent(inout) :: this
    this%field1 = 1.0e3_r8 + this%field1
    this%field2 = cheb1d_scalar_field(this%nodes,func)
#line 263 "cheb1d_scalar_test.pf"
  call assertTrue(this%field2==this%field1,message='Error calculating real+field', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 263) )
  if (anyExceptions()) return
#line 264 "cheb1d_scalar_test.pf"
  contains
    real(r8) pure function func(x)
      real(r8), intent(in) :: x
      func = f1(x) + 1.0e3_r8
    end function func
  end subroutine test_raf
  
!@Test
  subroutine test_far(this)
    class(test_field), intent(inout) :: this
    this%field1 = this%field1 + 1.0e3_r8
    this%field2 = cheb1d_scalar_field(this%nodes,func)
#line 276 "cheb1d_scalar_test.pf"
  call assertTrue(this%field2==this%field1,message='Error calculating field+real', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 276) )
  if (anyExceptions()) return
#line 277 "cheb1d_scalar_test.pf"
  contains
    real(r8) pure function func(x)
      real(r8), intent(in) :: x
      func = f1(x) + 1.0e3_r8
    end function func
  end subroutine test_far
  
!@Test
  subroutine test_fsf(this)
    class(test_field), intent(inout) :: this
    this%field1 = this%field1 - this%field2
    this%field2 = cheb1d_scalar_field(this%nodes,func)
#line 289 "cheb1d_scalar_test.pf"
  call assertTrue(this%field2==this%field1,message='Error calculating field-field', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 289) )
  if (anyExceptions()) return
#line 290 "cheb1d_scalar_test.pf"
  contains
    real(r8) pure function func(x)
      real(r8), intent(in) :: x
      func = f1(x) - f2(x)
    end function func
  end subroutine test_fsf
  
!@Test
  subroutine test_fsr(this)
    class(test_field), intent(inout) :: this
    this%field1 = this%field1 - 4.323_r8
    this%field2 = cheb1d_scalar_field(this%nodes,func)
#line 302 "cheb1d_scalar_test.pf"
  call assertTrue(this%field2==this%field1,message='Error calculating field-real', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 302) )
  if (anyExceptions()) return
#line 303 "cheb1d_scalar_test.pf"
  contains
    real(r8) pure function func(x)
      real(r8), intent(in) :: x
      func = f1(x) - 4.323_r8
    end function func
  end subroutine test_fsr
  
!@Test
  subroutine test_rsf(this)
    class(test_field), intent(inout) :: this
    this%field1 = 10.2_r8 - this%field1
    this%field2 = cheb1d_scalar_field(this%nodes,func)
#line 315 "cheb1d_scalar_test.pf"
  call assertTrue(this%field2==this%field1,message='Error calculating real-field', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 315) )
  if (anyExceptions()) return
#line 316 "cheb1d_scalar_test.pf"
  contains
    real(r8) pure function func(x)
      real(r8), intent(in) :: x
      func = 10.2_r8 - f1(x)
    end function func
  end subroutine test_rsf
  
!@Test
  subroutine test_fpr(this)
    class(test_field), intent(inout) :: this
    this%field1 = this%field1**5.2_r8
    this%field2 = cheb1d_scalar_field(this%nodes,func)
#line 328 "cheb1d_scalar_test.pf"
  call assertTrue(this%field2==this%field1,message='Error calculating field**real', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 328) )
  if (anyExceptions()) return
#line 329 "cheb1d_scalar_test.pf"
  contains
    real(r8) pure function func(x)
      real(r8), intent(in) :: x
      func = f1(x) ** 5.2_r8
    end function func
  end subroutine test_fpr
  
!@Test
  subroutine test_fpr4(this)
    class(test_field), intent(inout) :: this
    this%field1 = this%field1**3.5
    this%field2 = cheb1d_scalar_field(this%nodes,func)
#line 341 "cheb1d_scalar_test.pf"
  call assertTrue(this%field2==this%field1,message='Error calculating field**real', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 341) )
  if (anyExceptions()) return
#line 342 "cheb1d_scalar_test.pf"
  contains
    real(r8) pure function func(x)
      real(r8), intent(in) :: x
      func = f1(x) ** 3.5
    end function func
  end subroutine test_fpr4
  
!@Test
  subroutine test_fpi(this)
    class(test_field), intent(inout) :: this
    this%field1 = this%field1**3
    this%field2 = cheb1d_scalar_field(this%nodes,func)
#line 354 "cheb1d_scalar_test.pf"
  call assertTrue(this%field2==this%field1,message='Error calculating field**int', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 354) )
  if (anyExceptions()) return
#line 355 "cheb1d_scalar_test.pf"
  contains
    real(r8) pure function func(x)
      real(r8), intent(in) :: x
      func = f1(x) ** 3
    end function func
  end subroutine test_fpi
  
!@Test
  subroutine test_sin(this)
    class(test_field), intent(inout) :: this
    this%field1 = sin(this%field1)
    this%field2 = cheb1d_scalar_field(this%nodes,func)
#line 367 "cheb1d_scalar_test.pf"
  call assertTrue(this%field2==this%field1,message='Error calculating sin(field)', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 367) )
  if (anyExceptions()) return
#line 368 "cheb1d_scalar_test.pf"
  contains
    real(r8) pure function func(x)
      real(r8), intent(in) :: x
      func = sin(f1(x))
    end function func
  end subroutine test_sin
  
!@Test
  subroutine test_cos(this)
    class(test_field), intent(inout) :: this
    this%field1 = cos(this%field1)
    this%field2 = cheb1d_scalar_field(this%nodes,func)
#line 380 "cheb1d_scalar_test.pf"
  call assertTrue(this%field2==this%field1,message='Error calculating cos(field)', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 380) )
  if (anyExceptions()) return
#line 381 "cheb1d_scalar_test.pf"
  contains
    real(r8) pure function func(x)
      real(r8), intent(in) :: x
      func = cos(f1(x))
    end function func
  end subroutine test_cos
  
!@Test
  subroutine test_tan(this)
    class(test_field), intent(inout) :: this
    this%field1 = tan(this%field1)
    this%field2 = cheb1d_scalar_field(this%nodes,func)
#line 393 "cheb1d_scalar_test.pf"
  call assertTrue(this%field2==this%field1,message='Error calculating tan(field)', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 393) )
  if (anyExceptions()) return
#line 394 "cheb1d_scalar_test.pf"
  contains
    real(r8) pure function func(x)
      real(r8), intent(in) :: x
      func = tan(f1(x))
    end function func
  end subroutine test_tan
  
!@Test
  subroutine test_asin(this)
    class(test_field), intent(inout) :: this
    this%field1 = asin(this%field1)
    this%field2 = cheb1d_scalar_field(this%nodes,func)
#line 406 "cheb1d_scalar_test.pf"
  call assertTrue(this%field2==this%field1,message='Error calculating asin(field)', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 406) )
  if (anyExceptions()) return
#line 407 "cheb1d_scalar_test.pf"
  contains
    real(r8) pure function func(x)
      real(r8), intent(in) :: x
      func = asin(f1(x))
    end function func
  end subroutine test_asin
  
!@Test
  subroutine test_acos(this)
    class(test_field), intent(inout) :: this
    this%field1 = acos(this%field1)
    this%field2 = cheb1d_scalar_field(this%nodes,func)
#line 419 "cheb1d_scalar_test.pf"
  call assertTrue(this%field2==this%field1,message='Error calculating acos(field)', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 419) )
  if (anyExceptions()) return
#line 420 "cheb1d_scalar_test.pf"
  contains
    real(r8) pure function func(x)
      real(r8), intent(in) :: x
      func = acos(f1(x))
    end function func
  end subroutine test_acos
  
!@Test
  subroutine test_atan(this)
    class(test_field), intent(inout) :: this
    this%field1 = atan(this%field1)
    this%field2 = cheb1d_scalar_field(this%nodes,func)
#line 432 "cheb1d_scalar_test.pf"
  call assertTrue(this%field2==this%field1,message='Error calculating atan(field)', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 432) )
  if (anyExceptions()) return
#line 433 "cheb1d_scalar_test.pf"
  contains
    real(r8) pure function func(x)
      real(r8), intent(in) :: x
      func = atan(f1(x))
    end function func
  end subroutine test_atan
  
!@Test
  subroutine test_sinh(this)
    class(test_field), intent(inout) :: this
    this%field1 = sinh(this%field1)
    this%field2 = cheb1d_scalar_field(this%nodes,func)
#line 445 "cheb1d_scalar_test.pf"
  call assertTrue(this%field2==this%field1,message='Error calculating sinh(field)', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 445) )
  if (anyExceptions()) return
#line 446 "cheb1d_scalar_test.pf"
  contains
    real(r8) pure function func(x)
      real(r8), intent(in) :: x
      func = sinh(f1(x))
    end function func
  end subroutine test_sinh
  
!@Test
  subroutine test_cosh(this)
    class(test_field), intent(inout) :: this
    this%field1 = cosh(this%field1)
    this%field2 = cheb1d_scalar_field(this%nodes,func)
#line 458 "cheb1d_scalar_test.pf"
  call assertTrue(this%field2==this%field1,message='Error calculating cosh(field)', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 458) )
  if (anyExceptions()) return
#line 459 "cheb1d_scalar_test.pf"
  contains
    real(r8) pure function func(x)
      real(r8), intent(in) :: x
      func = cosh(f1(x))
    end function func
  end subroutine test_cosh
  
!@Test
  subroutine test_tanh(this)
    class(test_field), intent(inout) :: this
    this%field1 = tanh(this%field1)
    this%field2 = cheb1d_scalar_field(this%nodes,func)
#line 471 "cheb1d_scalar_test.pf"
  call assertTrue(this%field2==this%field1,message='Error calculating tanh(field)', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 471) )
  if (anyExceptions()) return
#line 472 "cheb1d_scalar_test.pf"
  contains
    real(r8) pure function func(x)
      real(r8), intent(in) :: x
      func = tanh(f1(x))
    end function func
  end subroutine test_tanh
  
!@Test
  subroutine test_asinh(this)
    class(test_field), intent(inout) :: this
    this%field1 = asinh(this%field1)
    this%field2 = cheb1d_scalar_field(this%nodes,func)
#line 484 "cheb1d_scalar_test.pf"
  call assertTrue(this%field2==this%field1,message='Error calculating asinh(field)', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 484) )
  if (anyExceptions()) return
#line 485 "cheb1d_scalar_test.pf"
  contains
    real(r8) pure function func(x)
      real(r8), intent(in) :: x
      func = asinh(f1(x))
    end function func
  end subroutine test_asinh
  
!@Test
  subroutine test_acosh(this)
    class(test_field), intent(inout) :: this
    this%field1 = acosh(this%field2)
    this%field2 = cheb1d_scalar_field(this%nodes,func)
#line 497 "cheb1d_scalar_test.pf"
  call assertTrue(this%field2==this%field1,message='Error calculating acosh(field)', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 497) )
  if (anyExceptions()) return
#line 498 "cheb1d_scalar_test.pf"
  contains
    real(r8) pure function func(x)
      real(r8), intent(in) :: x
      func = acosh(f2(x))
    end function func
  end subroutine test_acosh
  
!@Test
  subroutine test_atanh(this)
    class(test_field), intent(inout) :: this
    this%field1 = atanh(this%field1)
    this%field2 = cheb1d_scalar_field(this%nodes,func)
#line 510 "cheb1d_scalar_test.pf"
  call assertTrue(this%field2==this%field1,message='Error calculating atanh(field)', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 510) )
  if (anyExceptions()) return
#line 511 "cheb1d_scalar_test.pf"
  contains
    real(r8) pure function func(x)
      real(r8), intent(in) :: x
      func = atanh(f1(x))
    end function func
  end subroutine test_atanh
  
!@Test
  subroutine test_log(this)
    class(test_field), intent(inout) :: this
    this%field1 = log(this%field1)
    this%field2 = cheb1d_scalar_field(this%nodes,func)
#line 523 "cheb1d_scalar_test.pf"
  call assertTrue(this%field2==this%field1,message='Error calculating log(field)', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 523) )
  if (anyExceptions()) return
#line 524 "cheb1d_scalar_test.pf"
  contains
    real(r8) pure function func(x)
      real(r8), intent(in) :: x
      func = log(f1(x))
    end function func
  end subroutine test_log
  
!@Test
  subroutine test_log10(this)
    class(test_field), intent(inout) :: this
    this%field1 = log10(this%field1)
    this%field2 = cheb1d_scalar_field(this%nodes,func)
#line 536 "cheb1d_scalar_test.pf"
  call assertTrue(this%field2==this%field1,message='Error calculating log10(field)', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 536) )
  if (anyExceptions()) return
#line 537 "cheb1d_scalar_test.pf"
  contains
    real(r8) pure function func(x)
      real(r8), intent(in) :: x
      func = log10(f1(x))
    end function func
  end subroutine test_log10
  
!@Test
  subroutine test_exp(this)
    class(test_field), intent(inout) :: this
    this%field1 = exp(this%field1)
    this%field2 = cheb1d_scalar_field(this%nodes,func)
#line 549 "cheb1d_scalar_test.pf"
  call assertTrue(this%field2==this%field1,message='Error calculating exp(field)', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 549) )
  if (anyExceptions()) return
#line 550 "cheb1d_scalar_test.pf"
  contains
    real(r8) pure function func(x)
      real(r8), intent(in) :: x
      func = exp(f1(x))
    end function func
  end subroutine test_exp
  
!@Test
  subroutine test_abs(this)
    class(test_field), intent(inout) :: this
    this%field1 = abs(this%field1)
    this%field2 = cheb1d_scalar_field(this%nodes,func)
#line 562 "cheb1d_scalar_test.pf"
  call assertTrue(this%field2==this%field1,message='Error calculating abs(field)', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 562) )
  if (anyExceptions()) return
#line 563 "cheb1d_scalar_test.pf"
  contains
    real(r8) pure function func(x)
      real(r8), intent(in) :: x
      func = abs(f1(x))
    end function func
  end subroutine test_abs
  
!@Test
  subroutine test_sqrt(this)
    class(test_field), intent(inout) :: this
    this%field1 = sqrt(this%field2)
    this%field2 = cheb1d_scalar_field(this%nodes,func)
#line 575 "cheb1d_scalar_test.pf"
  call assertTrue(this%field2==this%field1,message='Error calculating sqrt(field)', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 575) )
  if (anyExceptions()) return
#line 576 "cheb1d_scalar_test.pf"
  contains
    real(r8) pure function func(x)
      real(r8), intent(in) :: x
      func = sqrt(f2(x))
    end function func
  end subroutine test_sqrt
  
!@Test
  subroutine test_minval(this)
    class(test_field), intent(inout) :: this
#line 586 "cheb1d_scalar_test.pf"
  call assertEqual(f1(this%lower1),minval(this%field1),message='Error calculating minval(field)', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 586) )
  if (anyExceptions()) return
#line 587 "cheb1d_scalar_test.pf"
#line 587 "cheb1d_scalar_test.pf"
  call assertEqual(f2(this%lower2),minval(this%field2),message='Error calculating minval(field)', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 587) )
  if (anyExceptions()) return
#line 588 "cheb1d_scalar_test.pf"
  end subroutine test_minval
  
!@Test
  subroutine test_maxval(this)
    class(test_field), intent(inout) :: this
#line 593 "cheb1d_scalar_test.pf"
  call assertEqual(f1(this%upper1),maxval(this%field1),message='Error calculating maxval(field)', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 593) )
  if (anyExceptions()) return
#line 594 "cheb1d_scalar_test.pf"
#line 594 "cheb1d_scalar_test.pf"
  call assertEqual(f2(this%upper2),maxval(this%field2),message='Error calculating maxval(field)', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 594) )
  if (anyExceptions()) return
#line 595 "cheb1d_scalar_test.pf"
  end subroutine test_maxval
  
!@Test
  subroutine test_d_dx(this)
    class(test_field), intent(inout) :: this
    integer :: i
    this%field1 = this%field2%d_dx(1)
    this%field2 = cheb1d_scalar_field(this%nodes,func)
#line 603 "cheb1d_scalar_test.pf"
  call assertTrue(this%field2==this%field1,message='Error calculating d(field)/dx', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 603) )
  if (anyExceptions()) return
#line 604 "cheb1d_scalar_test.pf"
    this%field1 = this%field2%d_dx(1,2)
    this%field2 = cheb1d_scalar_field(this%nodes,func)
#line 606 "cheb1d_scalar_test.pf"
  call assertTrue(this%field2==this%field1,message='Error calculating d(field)/dx', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 606) )
  if (anyExceptions()) return
#line 607 "cheb1d_scalar_test.pf"
    do i = 2,4
      this%field1 = this%field2%d_dx(i)
      this%field2 = cheb1d_scalar_field(this%nodes,zeros)
#line 610 "cheb1d_scalar_test.pf"
  call assertTrue(this%field2==this%field1,message='Error calculating d(field)/dx', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 610) )
  if (anyExceptions()) return
#line 611 "cheb1d_scalar_test.pf"
    end do
  contains
    real(r8) pure function func(x)
      real(r8), intent(in) :: x
      func = df2(x)
    end function func
    real(r8) pure function func2(x)
      real(r8), intent(in) :: x
      func2 = d2f2(x)
    end function func2
    real(r8) pure function zeros(x)
      real(r8), intent(in) :: x
      zeros = 0.0_r8
    end function zeros
  end subroutine test_d_dx
  
!~ @Test
!~   subroutine test_grad(this)
!~     class(test_field), intent(inout) :: this
!~     ! TODO: write implementation
!~     @assertTrue(this%field2==this%field1,message='Error calculating gradient of field')
!~   contains
!~     real(r8) pure function func(x)
!~       real(r8), intent(in) :: x
!~       func = df2(x)
!~     end function func
!~   end subroutine test_grad
  
!@Test
  subroutine test_lap(this)
    class(test_field), intent(inout) :: this
    integer :: i
    this%field1 = .laplacian. this%field1
    this%field2 = cheb1d_scalar_field(this%nodes,func)
#line 645 "cheb1d_scalar_test.pf"
  call assertTrue(this%field2==this%field1,message='Error calculating laplacian of field', &
 & location=SourceLocation( &
 & 'cheb1d_scalar_test.pf', &
 & 645) )
  if (anyExceptions()) return
#line 646 "cheb1d_scalar_test.pf"
  contains
    real(r8) pure function func(x)
      real(r8), intent(in) :: x
      func = df1(x)
    end function func
  end subroutine test_lap

end module cheb1d_scalar_test

module Wrapcheb1d_scalar_test
   use pFUnit_mod
   use cheb1d_scalar_test
   implicit none
   private

   public :: WrapUserTestCase
   public :: makeCustomTest
   type, extends(test_field) :: WrapUserTestCase
      procedure(userTestMethod), nopass, pointer :: testMethodPtr
   contains
      procedure :: runMethod
   end type WrapUserTestCase

   abstract interface
     subroutine userTestMethod(this)
        use cheb1d_scalar_test
        class (test_field), intent(inout) :: this
     end subroutine userTestMethod
   end interface

contains

   subroutine runMethod(this)
      class (WrapUserTestCase), intent(inout) :: this

      call this%testMethodPtr(this)
   end subroutine runMethod

   function makeCustomTest(methodName, testMethod) result(aTest)
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
      aTest%testMethodPtr => testMethod
#ifdef INTEL_13
      p => aTest
      call p%setName(methodName)
#else
      call aTest%setName(methodName)
#endif
   end function makeCustomTest

end module Wrapcheb1d_scalar_test

function cheb1d_scalar_test_suite() result(suite)
   use pFUnit_mod
   use cheb1d_scalar_test
   use Wrapcheb1d_scalar_test
   type (TestSuite) :: suite

   integer, allocatable :: npes(:)

   suite = newTestSuite('cheb1d_scalar_test_suite')

   call suite%addTest(makeCustomTest('test_domain', test_domain))

   call suite%addTest(makeCustomTest('test_dimensions', test_dimensions))

   call suite%addTest(makeCustomTest('test_raw_size', test_raw_size))

   call suite%addTest(makeCustomTest('test_raw', test_raw))

   call suite%addTest(makeCustomTest('test_assign_raw', test_assign_raw))

   call suite%addTest(makeCustomTest('test_resolution', test_resolution))

   call suite%addTest(makeCustomTest('test_fmf', test_fmf))

   call suite%addTest(makeCustomTest('test_rmf', test_rmf))

   call suite%addTest(makeCustomTest('test_fmr', test_fmr))

   call suite%addTest(makeCustomTest('test_fdf', test_fdf))

   call suite%addTest(makeCustomTest('test_rdf', test_rdf))

   call suite%addTest(makeCustomTest('test_fdr', test_fdr))

   call suite%addTest(makeCustomTest('test_faf', test_faf))

   call suite%addTest(makeCustomTest('test_raf', test_raf))

   call suite%addTest(makeCustomTest('test_far', test_far))

   call suite%addTest(makeCustomTest('test_fsf', test_fsf))

   call suite%addTest(makeCustomTest('test_fsr', test_fsr))

   call suite%addTest(makeCustomTest('test_rsf', test_rsf))

   call suite%addTest(makeCustomTest('test_fpr', test_fpr))

   call suite%addTest(makeCustomTest('test_fpr4', test_fpr4))

   call suite%addTest(makeCustomTest('test_fpi', test_fpi))

   call suite%addTest(makeCustomTest('test_sin', test_sin))

   call suite%addTest(makeCustomTest('test_cos', test_cos))

   call suite%addTest(makeCustomTest('test_tan', test_tan))

   call suite%addTest(makeCustomTest('test_asin', test_asin))

   call suite%addTest(makeCustomTest('test_acos', test_acos))

   call suite%addTest(makeCustomTest('test_atan', test_atan))

   call suite%addTest(makeCustomTest('test_sinh', test_sinh))

   call suite%addTest(makeCustomTest('test_cosh', test_cosh))

   call suite%addTest(makeCustomTest('test_tanh', test_tanh))

   call suite%addTest(makeCustomTest('test_asinh', test_asinh))

   call suite%addTest(makeCustomTest('test_acosh', test_acosh))

   call suite%addTest(makeCustomTest('test_atanh', test_atanh))

   call suite%addTest(makeCustomTest('test_log', test_log))

   call suite%addTest(makeCustomTest('test_log10', test_log10))

   call suite%addTest(makeCustomTest('test_exp', test_exp))

   call suite%addTest(makeCustomTest('test_abs', test_abs))

   call suite%addTest(makeCustomTest('test_sqrt', test_sqrt))

   call suite%addTest(makeCustomTest('test_minval', test_minval))

   call suite%addTest(makeCustomTest('test_maxval', test_maxval))

   call suite%addTest(makeCustomTest('test_d_dx', test_d_dx))

   call suite%addTest(makeCustomTest('test_lap', test_lap))


end function cheb1d_scalar_test_suite

