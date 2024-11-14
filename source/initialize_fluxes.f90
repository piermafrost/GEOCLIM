subroutine initialize_fluxes()

    ! Define initial value of fluxes
    ! > some of fluxes are held constant, other needs 0 value on certain basins

    use constante, only: xhydLiin
    implicit none
    include 'combine_foam.inc'

    fCO2_crust(:)      = 0
    fSO4_crust(:)      = 0
    fSO4_basin(:)      = 0
    fmor(:)            = 0
    FhydLi(:)          = 0
    ws(:)              = 0
    fin_sed(:)         = 0
    fodc(:)            = 0
    fodp(:)            = 0
    fO2_odc(:)         = 0
    fSulfRed(:)        = 0
    Corg_BE(:)         = 0
    fCO2atm_ocean(:)   = 0
    fC13atm_ocean(:)   = 0
    fC13ocean_atm(:)   = 0
    frain_sed(:)       = 0
    frain_POC(:)       = 0
    fin_sed(:)         = 0
    fin_POC(:)         = 0
    fexch_sed(:,:)     = 0

    fvol     = volin*clo
    fanthros = 0

    do j0=1,ndeep
        j = jbox_deep(j0)
        fmor(j)   = xMORin*clo/ndeep
        FhydLi(j) = xhydLiin*clo/ndeep
    end do

end subroutine
