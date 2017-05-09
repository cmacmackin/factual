!
!  vector_pool.f90
!  This file is part of FACTUAL.
!
!  Copyright 2017 Christopher MacMackin <cmacmackin@gmail.com>
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

module vector_pool_mod
  !* Author: Chris MacMackin
  !  Date: May 2017
  !  License: GPLv3
  !
  ! An [object
  ! pool](https://en.wikipedia.org/wiki/Object_pool_pattern) which can
  ! be used for "dynamically" allocating vector fields.
  !
  use abstract_fields_mod, only: vector_field
  implicit none

  private
  
  type, public :: vector_pool
    !* Author: Chris MacMackin
    !  Date: May 2017
    !
    ! Provides a number of vector fields which can be "dynamically"
    ! allocated as function results, without fear of memory
    ! leaks. This is an implementation of the [object
    ! pool](https://en.wikipedia.org/wiki/Object_pool_pattern)
    ! pattern.
    !
    integer                                    :: num_fields
      !! The maximum number of vector fields which can be used 
    class(vector_field), dimension(:), pointer :: pool => null()
      !! An array of vector fields which cane be used
    logical, dimension(:), allocatable         :: in_use
      !! Indicates whether a particular vector field is in use.
  contains
    procedure :: acquire
    procedure :: release
  end type vector_pool

  interface vector_pool
    module procedure :: constructor
  end interface vector_pool
  
contains

  function constructor(num_fields, prototype) result(this)
    !* Author: Chris MacMackin
    !  Date: May 2017
    !
    ! Creates a new vector pool of the specified size and type.
    !
    integer, intent(in)             :: num_fields
      !! Number of fields to be available in the new pool
    class(vector_field), intent(in) :: prototype
      !! Dynamic type of the fields in the new pool
    type(vector_pool)               :: this
    integer :: i
    this%num_fields = num_fields
    allocate(this%pool(num_fields), mold=prototype)
    allocate(this%in_use(num_fields))
    this%in_use = .false.
    do concurrent (i=1:num_fields)
      call this%pool(i)%set_pool_id(i)
    end do
  end function constructor

  function acquire(this) result(field)
    !* Author: Chris MacMackin
    !  Date: May 2017
    !
    ! Returns a reference to an available vector field.
    !
    class(vector_pool), intent(inout) :: this
    class(vector_field), pointer      :: field
    integer :: i
    do i = 1, this%num_fields
      if (.not. this%in_use(i)) then
        this%in_use(i) = .true.
        field => this%pool(i)
        call field%set_temp()
        return
      end if
    end do
    error stop ('Attempting to use more vector fields than '// &
                'allocated in this pool.')
  end function acquire

  pure subroutine release(this, id)
    class(vector_pool), intent(inout) :: this
    integer, intent(in)               :: id
      !! The ID number/array subscript of the field object to be
      !! released.
#:if defined('DEBUG')
    if (.not. this%in_use(id)) error stop ('Attempting to release '// &
                                           'unused field.')
#:endif
    this%in_use(id) = .false.
  end subroutine release
  
end module vector_pool_mod
