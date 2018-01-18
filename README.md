# PaleoCalendarAdjust #

#### Work-in-progress versions of programs for calculating orbitally determined paleo month lengths, and using these to adjust CMIP5/PMIP3-formatted paleo simulation files that were created by summarizing data using a fixed modern calendar. ####

**Overview**

Code example:

```fortran	
subroutine hdaily(nm,nd,xm,monlen,no_negatives,xdh)

    implicit none

    integer(4), parameter   :: nh=6
    integer(4), intent(in)  :: nd,nm
    real(8), intent(in)     :: xm(nm)
    integer(4), intent(in)  :: monlen(nm)
    logical, intent(in)     :: no_negatives
    real(8), intent(out)    :: xdh(nd)

    real(8)             :: a(0:nh),b(0:nh)
    integer(4)          :: m

    if (debug) then
        do m=1,nm
            write (debug_unit,*) m,xm(m),monlen(m)
        end do
    end if

    ! interpolate daily values
    call harmonic_coeffs(nm,xm,a,b)
    call xdhat(nm,nd,monlen,a,b,xdh)
    if (no_negatives) call dzero(nm,nd,monlen,xm,xdh)

end subroutine hdaily
```

Some more text.
