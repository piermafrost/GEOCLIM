    subroutine biological_pump(t)
!   *****************************
    use constante, only: cpred
    implicit none
    include 'combine_foam.inc'


    ! Compute fluxes
    ! --------------

    ! initialization
    fbioP(:)      = 0
    fbioC(:)      = 0
    reff(:)       = 0
    rC_Corg(:)    = 0
    carb_ratio(:) = 0
    finorgC(:)    = 0
    finorgP(:)    = 0

    ! -> add continental inputs to biological P uptake flux (only in "app_cont = 1" basins)
    do j0=1,nappcont
        j = jbox_appcont(j0)
        fbioP(j) = fbioP(j) + fpw(j)
    end do

    ! -> rest of fluxes computation (in surface basins only)
    do j0=1,nsurface
        j = jbox_surface(j0)

        ! add ocean basins exchanges to biological P uptake flux
        do k = 1,nbasin-1
            fbioP(j) = fbioP(j) + F(k,j)*var_diss(3,k)
        end do

        ! Efficiency of biological P uptake (i.e., fraction of maximum potential P flux that is actually used)
        !                vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv   => CO2 limitation at very low pCO2
        reff(j)=max(0d0, (pco2_dissous(j)-0.2)/(pco2_dissous(j)-0.1))

        ! Fraction of calcifying primary producter (imposed)
        carb_ratio(j) = 0.15

        ! Reduce bioproductivity in polar basins (=> light limitation) 
        reff(j)       = (1 - 0.75*indice_polar(j))*reff(j)

        ! Reduce carbonate:org C productivity ratio in epicontinental basins
        carb_ratio(j) = (1 - 0.9*indice_epicont(j))*carb_ratio(j)

        ! -- Biological fluxes -- !
        fbioP(j) = reff(j)*fbioP(j) ! => only the "reff" fraction of P can be used by organisms
        fbioC(j) = fbioP(j)*cpred
        ! ----------------------- !

        ! -------- Biological inorganic fluxes (i.e., shells) -------- !
        if (omega_0(j).lt.1.0) then  !checking the saturation state of the ocean with respect to carbonates
            rC_Corg(j)=0.
        else
            rC_Corg(j)=carb_ratio(j)*(omega_0(j)-1)/(0.4+(omega_0(j)-1))
        endif

        ! Shelf-flag: calcifying organisms exist or not
        !   -> shelfal = 0.01 in non-epicontinental boxes if ishelfal==1 (i.e., calcifyers on shelves only)
        !      shelfal = 1    if ishelfal==0, and in epicontinental box regardless of ishelfal
        shelfal = 1 - 0.99*(1-indice_epicont(j))*ishelfal

        finorgC(j)=fbioC(j)*rC_Corg(j)*shells*shelfal
        finorgP(j)=0 ! => pas de P ds les coquilles ! !finorgC(j)/1000.
        ! ------------------------------------------------------------ !

    end do

    return
    end
