!
!  utils.f90
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
#:endif

module utils_mod
  !* Author: Chris MacMackin
  !  Date: September 2016
  !  License: GPLv3
  !
  ! Provides useful helper functions for internal use with FACTUAL.
  !
  use iso_fortran_env, only: r8 => real64, stderr => error_unit
  use abstract_fields_mod, only: abstract_field
!  use iso_c_binding, only: c_loc, c_ptr
  use hdf5
  use h5lt
  implicit none
  private

  public :: is_nan, check_set_from_raw, elements_in_slice, grid_to_spacing, &
            open_or_create_hdf_group, create_hdf_dset, lagrange_interp

contains

  elemental logical function is_nan(val)
    !* Author: Chris MacMackin
    !  Date: September 2016
    !
    ! Returns true if the argument is NaN, false otherwise. If using `gfortran`
    ! then the `isnan` extension is used. Otherwise, the value is tested for
    ! equality against itself.
    !
    real(r8), intent(in) :: val !! The value being checked to see if it is NaN
#ifdef __GFORTRAN__
    is_nan = isnan(val)
#else
    is_nan = val /= val
#endif
  end function is_nan

  pure subroutine check_set_from_raw(field,raw,provide_lower_bound,provide_upper_bound)
    !* Author: Chris MacMackin
    !  Date: September 2016
    !
    ! Checks to ensure that a raw array is the correct size to set the contents
    ! of a field object. If it is not then it will print an error message and
    ! stop the program.
    !
    class(abstract_field), intent(in) :: field
      !! The field being checked to see if it is compatible with the raw array.
    real(r8), dimension(:), intent(in) :: raw
      !! The raw data to be stored in this array.
    integer, dimension(:), optional, intent(in) :: provide_lower_bound
      !! Specifies how many layers of data points are excluded from
      !! the raw data at the lower boundary for each dimension. The
      !! number in element `n` of the array indicates how many layers
      !! of cells at the lower boundary normal to dimension `n` are
      !! missed.
    integer, dimension(:), optional, intent(in) :: provide_upper_bound
      !! Specifies how many layers of data points are excluded
      !! from the raw data at the upper boundary for each
      !! dimension. The number in element `n` of the array indicates
      !! how many layers of cells at the upper boundary normal to
      !! dimension `n` are missed.
#ifdef DEBUG
    integer :: expected
    expected = field%raw_size(provide_lower_bound,provide_upper_bound)
    if (expected /= size(raw)) then
      error stop('check_set_from_raw: Error, setting from raw array of '//&
                 'wrong size')
    end if
#endif
  end subroutine check_set_from_raw

  elemental function elements_in_slice(start, finish, stride) result(res)
    !* Author: Chris MacMackin
    !  Date: October 2016
    !
    ! Computes the number of array elements which would be returned
    ! from a slice of the form `start:finish:stride'.
    !
    integer, intent(in) :: start
      !! The index which the array slice starts counting from, inclusive
    integer, intent(in) :: finish
      !! The index at which the array slice stops reading elements, inclusive
    integer, intent(in) :: stride
      !! How frequently to read an element from the array. `stride =
      !! 1` would read every element, `stride = 2` would read every
      !! second element, etc.
    integer :: res
    res = (finish - start)/stride + 1
  end function elements_in_slice

  pure function grid_to_spacing(grid) result(space)
    !* Author: Chris MacMackin
    !  Date: December 2016
    !
    ! For a given specification of the location of grid-points,
    ! calculates the width of each grid cell. Grid cells are taken to
    ! be centred on the grid points.
    !
    real(r8), dimension(:), intent(in) :: grid
    real(r8), dimension(size(grid),1)  :: space
    integer :: i
    do concurrent (i=1:size(grid))
      if (i==1) then
        space(i,1) = grid(2) - grid(1)
      else if (i==size(grid)) then
        space(i,1) = grid(i) - grid(i-1)
      else
        space(i,1) = (grid(i+1) - grid(i-1))/2._r8
      end if
    end do
  end function grid_to_spacing

  subroutine open_or_create_hdf_group(loc_id, group_name, group_id, hdferr)
    !* Author: Paul Anton Letnes
    !  Date: April 2017
    !
    ! Checks if group with the given name exists. Create it if it doesn't, 
    ! open it if it does.
    !
    ! Copied fromp the [HDF
    ! forum](http://hdf-forum.184993.n3.nabble.com/Check-if-group-exists-td2601701.html)
    ! and slightly modified
    !
    integer(hid_t), intent(in)   :: loc_id
      !! File or group identifier
    character(len=*), intent(in) :: group_name 
      !! Name of the group
    integer(hid_t), intent(out)  :: group_id 
      !! Group identifier
    integer, intent(inout)       :: hdferr 
      !! Error code. 0 on success, -1 on failure

    logical :: group_exists 

    call h5lexists_f(loc_id, group_name, group_exists, hdferr) 
    if (group_exists) then 
      call h5gopen_f(loc_id, group_name, group_id, hdferr) 
    else 
      call h5gcreate_f(loc_id, group_name, group_id, hdferr) 
    end if 
  end subroutine open_or_create_hdf_group

  subroutine create_hdf_dset(loc_id, dset_name, rank, dims, buf, hdferr)
    !* Author: Chris MacMackin
    !  Date: April 2017
    !
    ! Creates a dataset if it doesn't already exist.
    !
    !
    integer(hid_t), intent(in)                 :: loc_id
      !! File or group identifier
    character(len=*), intent(in)               :: dset_name 
      !! Name of the data set
    integer, intent(in)                        :: rank
      !! Rank of the data set
    integer(hsize_t), dimension(*), intent(in) :: dims 
      !! Size of the data set  
    real(r8), intent(in), dimension(*)         :: buf  
      !! Buffer containing the data to be written to the data set.
    integer, intent(inout)                     :: hdferr 
      !! Error code. 0 on success, -1 on failure

    integer :: set_exists 

    set_exists = h5ltfind_dataset_f(loc_id, dset_name) 
    if (set_exists == 0) then 
      call h5ltmake_dataset_double_f(loc_id, dset_name, rank, dims, buf, hdferr) 
    end if 
  end subroutine create_hdf_dset

  subroutine get_hdf_dset_ref(loc_id, dset_name, rank, dims, buf, &
                              ref, hdferr)
    !* Author: Chris MacMackin
    !  Date: April 2017
    !
    ! Returns a reference to an HDF dataset of the given name. If such
    ! a dataset does not already exist, create it first.
    !
    integer(hid_t), intent(in)                 :: loc_id
      !! File or group identifier
    character(len=*), intent(in)               :: dset_name 
      !! Name of the data set
    integer, intent(in)                        :: rank
      !! Rank of the data set
    integer(hsize_t), dimension(*), intent(in) :: dims 
      !! Size of the data set  
    real(r8), intent(in), dimension(*)         :: buf  
      !! Buffer containing the data to be written to the data set.
    type(hobj_ref_t_f), intent(inout)          :: ref
      !! The reference to the data set
    integer, intent(inout)                     :: hdferr 
      !! Error code.
    
    integer :: exists
    
    exists = h5ltfind_dataset_f(loc_id, dset_name)
    if (exists /= 0) then
      call h5ltmake_dataset_double_f(loc_id, dset_name, rank, dims, &
                                     buf, hdferr)
    end if
    call h5rcreate_f(loc_id, dset_name, ref, hdferr)
  end subroutine get_hdf_dset_ref

  function lagrange_interp(x, f, xf) result(val)
    !* Author: Chris MacMackin
    !  Date: July 2017
    !
    ! Performs Lagrange interpolation. This is only suitable for
    ! nonuniformly space data-points, such as Chebyshev colocation
    ! nodes.
    !
    real(r8), intent(in)                     :: x
    real(r8), dimension(:), intent(in)       :: f
      !! The data values to interpolate between.
    real(r8), dimension(size(f)), intent(in) :: xf
      !! The coordinates of the data values. These must be either
      !! monotonically increasing or decreasing.
    real(r8)                                 :: val

    integer :: i, j, n
    real(r8) :: wx, denom
    logical :: increasing, recalc
    real(r8), dimension(:), allocatable, save :: w
    real(r8), save :: x1, xn

    val = 0._r8
    denom = 0._r8

    n = size(f)
    
    ! If asked to interpolate outside of domain, return nearest value
    if (xf(1) < xf(n)) then
      if (x < xf(1)) then
        val = f(1)
        return
      else if (x > xf(n)) then
        val = f(n)
        return
      end if
    else
      if (x > xf(1)) then
        val = f(1)
        return
      else if (x < xf(n)) then
        val = f(n)
        return
      end if
    end if

    if (allocated(w)) then
      if (size(w) /= n) then
        recalc = .true.
        deallocate(w)
        allocate(w(n))
      else
        recalc = (xf(1) == x1 .and. xf(n) == xn)        
      end if
    else
      recalc = .true.
      allocate(w(n))
    end if
    if (recalc) then
      x1 = xf(1)
      xn = xf(n)
      denom = 0._r8
      do i = 1, n
        w(i) = 1._r8
        do j = 1, n
          if (i /= j) w(i) = w(i) * (xf(i) - xf(j))
        end do
      end do
    end if

    do i = 1, n
      wx = 1._r8/(w(i)*(x - xf(i)))
      val = val + wx*f(i)
      denom = denom + wx
    end do
    val = val/denom
  end function lagrange_interp

!  subroutine set_attribute_ref(loc_id, dset_name, attr_name, ref, &
!                               errcode)
!    !* Author: Chris MacMackin
!    !  Date: April 2017
!    !
!    ! Creates an attribute on a dataset which will contain a reference
!    ! to another object.
!    !
!    integer(hid_t), intent(in)             :: loc_id
!      !! File or group identifier
!    character(len=*), intent(in)           :: dset_name 
!      !! Name of the data set
!    character(len=*), intent(in)           :: attr_name
!      !! Name of the attribute
!    type(hobj_ref_t_f), intent(in), target :: ref
!      !! The reference to be contained in the attribute
!    integer, intent(out)                   :: errcode 
!      !! Error code.
!    
!    integer(hid_t) :: dset_id, type_id, space_id, attr_id
!    type(c_ptr) :: reference
!
!    call h5dopen_f(loc_id, dset_name, dset_id, errcode)
!    call h5tcopy_f(H5T_REFERENCE_F, type_id, errcode)
!    call h5screate_f(H5S_SCALAR_F, space_id, errcode)
!    call h5acreate_f(dset_id, attr_name, type_id, space_id, &
!                     attr_id, errcode)
!    reference = c_loc(ref)
!    call h5awrite_f(attr_id, type_id, reference, errcode)
!    call h5aclose_f(attr_id, errcode)
!    call h5sclose_f(space_id, errcode)
!    call h5dclose_f(dset_id, errcode)
!  end subroutine set_attribute_ref

end module utils_mod
