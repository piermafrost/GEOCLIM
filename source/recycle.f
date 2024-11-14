    subroutine recycle(t)
!   *********************
    use constante, only: k_oxyd, KO2
    implicit none
    include 'combine_foam.inc'

    ranoxyd = 0

    ! Surface basins (O2 in equilibrium with atmosphere)
    do j0=1,nsurface
        j = jbox_surface(j0)
        roxyd(j) = 0
    enddo

    ! Deep and intermediate basins
    do j0=1,nnosurface
        j = jbox_nosurface(j0)
        roxyd(j) = k_oxyd*(1 - dexp(-var_diss(6,j)/KO2))
    end do

    ! Atmosphere
    roxyd(nbasin) = 0

    return
    end
