    subroutine varset(y)
!   ********************
    implicit none

    include 'combine_foam.inc'
!   ****************************************
!   variable order:
!     all dissolved variables,
!     followed by all particulate variables,
!     followed by all isotopic variables
!   ****************************************

    k=0

    !========
    do i=1,5 ! dissolved variables
        do j=1,nbasin
            k=k+1
            var_diss(i,j) = y(k)/vol(k)
        end do
    end do
    !--------
    i = 6 ! Specific case for Oxygen => calculation done only in no-surface boxes...
    do j0=1,nnosurface
        j = jbox_nosurface(j0)
        var_diss(i,j) = y(k+j)/vol(k+j)
    end do
    var_diss(i,nbasin) = y(k+nbasin)/vol(k+nbasin) ! ... and in atmosphere box
    k = k+nbasin
    !--------
    do i=7,nvar_diss ! rest of dissolved variables
        do j=1,nbasin
            k=k+1
            var_diss(i,j) = y(k)/vol(k)
        end do
    end do

    !========
    do i=1,nvar_part ! particulate variables
        do j=1,nbasin
            k=k+1
            var_part(i,j) = y(k)/vol(k)
        end do
    end do

    !========
    do i=1,nvar_isot ! isotopic variables
        do j=1,nbasin
            k=k+1
            var_isot(i,j) = y(k)
        end do
    end do


    end subroutine
