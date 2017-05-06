#:setvar UNARY_FUNCTIONS [('sin','\sin'), ('cos','\cos'),              &
  & ('tan','\\tan'), ('asin','\sin^{-1}'), ('acos','\cos^{-1}'),       &
  & ('atan','\\tan^{-1}'), ('sinh','\sinh'), ('cosh','\cosh'),         &
  & ('tanh','\\tanh'), ('asinh','\sinh^{-1}'), ('acosh','\cosh^{-1}'), &
  & ('atanh','\\tanh^{-1}'), ('log','\ln'), ('log10','\log'),          &
  & ('exp','e^'), ('abs','\\abs'), ('sqrt','\sqrt')]

#:def unary_func(func,tex,field_name)
  function ${field_name}$_${func}$(this) result(res)
    !* Author: Chris MacMackin
    !  Date: March 2016
    !
    ! \(${tex}$({\rm field})\)
    !
    class(${field_name}$_field), intent(in) :: this
    class(scalar_field), allocatable :: res !! The result of this operation
    class(${field_name}$_field), allocatable :: local
    call this%guard_temp()
    allocate(local, mold=this)
    call local%assign_meta_data(this)
    if (allocated(this%field_data)) then
      local%field_data = ${func}$(this%field_data)
    end if
    call move_alloc(local, res)
    call res%set_temp()
    call this%clean_temp()
  end function ${field_name}$_${func}$
#:enddef

#:def unary_binding(func,tex,field_name)
    procedure :: ${func}$ => ${field_name}$_${func}$
      !! \(${tex}$({\rm field})\)
      
#:enddef

#:def public_unary()
#:for FUNC, TEX in UNARY_FUNCTIONS
  public :: ${FUNC}$
#:endfor
#:enddef

