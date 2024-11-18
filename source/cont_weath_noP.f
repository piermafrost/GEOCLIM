subroutine cont_weath_noP(t)
!************************
use dynsoil, only: regolith_time_integration
use dynsoil_steady_state, only: steady_state_surface_regolith
use constante, only: PI_n_CO2_atm, Rgas, TSS_density, rsw_litho, CaMg_rock, P_rock, P2C_ker, P2C_carb, OC_rock, cp_cont, &
                     Sulf_rock, BASALT_LITHO_NUM, k_carb_wth
!use dynsoil_lithium, only: lithium
implicit none
double precision :: loc_granwth, sum_tss, loc_flx
include 'combine_foam.inc'
! WARNING: carbonate must be the last lithological class, whatever the chosen
! number of classes



! current atmospheric CO2 level (in PAL)
p=var_diss(7,nbasin)/PI_n_CO2_atm



! Total DISCHARGE of water to the ocean
!--------------------------------------

discharge(:) = 0.

do j0=1,ncontpxl
    j = list_cont_pixel(j0)
    kl = cont_basin_map(j)
    discharge(kl) = discharge(kl) + 1e10*runclim(j)*areaclimber(j)
end do !                            ^^^^^
! convert area [1e6km2 -> m2] and runoff [cm/y -> m/y]






!========================================================!
!========================================================!
!==               WEATHERING COMPUTATION               ==!
!========================================================!
!========================================================!



call FACCO2(t,pco2,fco2,fe)





!=========================================================================!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! VARIABLES DEPENDENT OF DYNSOIL USAGE !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



if (.not. coupling_dynsoil) then
! ######################################
! # if not coupled with dynsoil module #
! ######################################


    ! initialising:
    fsilw(:) = 0.
    fbasw(:) = 0.
    FrivLi(:) = 0
    dLiriv(:) = 0.
    fkerw(:) = 0.
    total_cont_POC_export(:) = 0.
    weighted_rsw(:) = 0.


    ! Riverine sediment discharge => parameteric power law of water discharge
    tss(:) = 2.9d6 * discharge(:)**0.5
    sum_tss = sum(tss)
     

    do j0=1,ncontpxl
        j = list_cont_pixel(j0)
        kl = cont_basin_map(j)

        ! =====================
        ! Basaltic weathering
        ! =====================
        ! "basalt" lithology => number BASALT_LITHO_NUM
        facbas=(-42300./Rgas)*(1./(Tclim(j)+273.15)-1./288.15)
        wth_litho(BASALT_LITHO_NUM,j) = (0.0483d+10)*dexp(facbas)*runclim(j) !0.044086200331d+10  !483685972.093689  !172421904
        wth_litho_wgh(BASALT_LITHO_NUM,j) = wth_litho(BASALT_LITHO_NUM,j) * litho_frac(BASALT_LITHO_NUM,j)
        fbasj = wth_litho_wgh(BASALT_LITHO_NUM,j) * areaclimber(j)
        fbasw(kl) =fbasw(kl) + fbasj*clo

        ! =====================
        ! Granitic weathering
        ! =====================
        ! skip basalts (BASALT_LITHO_NUM) and carbonates (last litho)
        facgra=(-48200./Rgas)*(1./(Tclim(j)+273.15)-1./288.15)
        loc_granwth = (0.11379d+10)*dexp(facgra)*runclim(j) * fco2*fe !0.101357275d+10
        fsilj = 0
        do k = 1,BASALT_LITHO_NUM-1 ! skip basalts
            wth_litho(k,j) = loc_granwth
            wth_litho_wgh(k,j) = wth_litho(k,j) * litho_frac(k,j)
            fsilj = fsilj + wth_litho_wgh(k,j)
            weighted_rsw(kl) = weighted_rsw(kl) + rsw_litho(k)*wth_litho_wgh(k,j)*areaclimber(j)
        end do
        do k = BASALT_LITHO_NUM+1,nlitho-1 ! skip carbonates
            wth_litho(k,j) = loc_granwth
            wth_litho_wgh(k,j) = wth_litho(k,j) * litho_frac(k,j)
            fsilj = fsilj + wth_litho_wgh(k,j)
            weighted_rsw(kl) = weighted_rsw(kl) + rsw_litho(k)*wth_litho_wgh(k,j)*areaclimber(j)
        end do
        fsilj = fsilj * areaclimber(j)
        fsilw(kl) = fsilw(kl) + fsilj*clo

        ! conversion /1e6km2 => /m2 for ALL silicate weathering rates
        wth_litho(1:nlitho-1,j)     = 1e-12*wth_litho(1:nlitho-1,j)
        wth_litho_wgh(1:nlitho-1,j) = 1e-12*wth_litho_wgh(1:nlitho-1,j)

        wth_allsil(j) = sum(wth_litho_wgh(1:nlitho-1,j))

        ! =====================
        ! Kerogen weathering
        ! =====================
        ! Old formulation:
        !! facker=(-43000./Rgas)*(1./(Tclim(j)+273.15)-1./288.15)
        !fker(j)=4.4184d+8*runclim(j)*areaclimber(j)   !*var_diss(6,nbasin)/38.d+18 !& !2.4 !1.432
        !                                              ! *dexp(facker)  !0.672d+10
        !fkerw(kl) = fkerw(kl) + fker(j)*clo
        fker(j) = 0.5 * sum(litho_frac(:,j)*OC_rock)*(sum_tss/TSS_density/areatot)
        !         ^^^ : oxidation efficiency
        if (.not. lock_oxygen_cycle)  fkerw(kl) = fkerw(kl) + fker(j)*areaclimber(j)*clo
        fker(j) = 1e-12*fker(j) ! conversion /1e6km2 => /m2

        ! Pyrite weathering (indicative value)
        pyr_wth(j) = 0.025*sum_tss/areatot

    end do


    ! ======================
    ! Continental POC export
    ! ======================

    ! erosion export (t/km2/an)
    eros_galy_unit=(sum_tss/areatot)*1e-6*1e-3  !moving from kg/1e6km2/yr to t/km2/yr
    do j0=1,ncontpxl
        j = list_cont_pixel(j0)
        kl = cont_basin_map(j)
        !scaling:
        POC_export_rate(j)=0.081*(eros_galy_unit**0.56)  !in t/km2/yr
        POC_export_rate(j)=POC_export_rate(j)/12.         !mol/m2/yr
        total_cont_POC_export(kl) = total_cont_POC_export(kl) + POC_export_rate(j)*areaclimber(j)*1.e+12  !mol/yr
    end do


    ! ===============================================================================
    ! Pyrite weathering and distribution in carb dissol, sil dissol and H2SO4 release
    ! ===============================================================================

    ! => molar flux of SO4^2- from pyrite oxidative weathering is "0.025*TSS":
    ! molar flux of Ca released to the ocean from carbonate dissolution by sulfuric acid
    fcarbsulfw(:)  =  0.724 * 0.025*tss(:)
    ! molar flux of Ca released to the ocean from silicate dissolution by sulfuric acid
    fsilsulfw(:)   =  0.276 * 0.025*tss(:)
    ! molar flux of H2SO4 from pyrite oxidation directly reaching the ocean
    fH2SO4sulfw(:) =  0! * 0.025*tss(:)


    ! lithology-averaged Strontium isotopic ratio
    where (fsilw /= 0d0)  weighted_rsw(:) = weighted_rsw(:)/fsilw(:)



else
! ######################################
! #   if coupled with dynsoil module   #
! ######################################


    ! increase time counter of dynsoil (=> asynchronous coupling) !
    icount_DS_int = icount_DS_int + 1
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!


    if ( icount_DS_int >= ijump_DS_integration ) then ! if integration time has been reached

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! integrate dynsoil module !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! counter reinitialization:
        icount_DS_int = 0

        ! =============================
        ! DYNSOIL (silicate weathering)
        ! =============================
        !
        if (use_dynsoil_steady_state) then
            call steady_state_surface_regolith(reg_thick, reg_x_surf, reg_tau_surf, reg_prod, reg_eros, reg_P_diss, &
                                               Tclim, runclim, slope, DS_timestep, veget_factor, veget_eros_factor, &
                                               list_cont_pixel, ncontpxl                                            )
        else
            call regolith_time_integration(xlevs, reg_thick, reg_x_surf, reg_tau_surf, reg_P_vol, reg_ktop,                 &
                                           reg_z_prof, reg_tau_prof,                                                        &
                                           reg_prod,reg_eros,reg_P_diss,reg_P_eros, reg_x_surf_eros, Tclim, runclim, slope, &
                                           DS_timestep,                                                                     &
                                           veget_factor, veget_eros_factor, list_cont_pixel, ncontpxl                       )
        end if
        !call lithium(reg_P_diss, Tclim, runclim, reg_thick, reg_Li_Friv, reg_Li_Fsp, reg_Li_driv, list_cont_pixel, ncontpxl)


        ! Continental sum:
        ! ----------------

        ! initialising:
        fsilw(:) = 0.
        fbasw(:) = 0.
        FrivLi(:) = 0.
        dLiriv(:) = 0.
        tss(:)   = 0.
        fkerw(:) = 0.
        total_cont_POC_export(:) = 0.
        weighted_rsw(:) = 0.
        fcarbsulfw(:)  = 0.
        fsilsulfw(:)   = 0.
        fH2SO4sulfw(:) = 0.

        do j0=1,ncontpxl
            j = list_cont_pixel(j0)
            kl = cont_basin_map(j)

            wth_allsil(j) = 0
            fsilj = 0
            do k = 1,nlitho-1 ! as k=nlitho is carbonate
                !
                wth_litho(k,j) = CaMg_rock(k) * reg_P_diss(k,j)
                wth_litho_wgh(k,j) = wth_litho(k,j) * litho_frac(k,j)
                wth_allsil(j) = wth_allsil(j) + wth_litho_wgh(k,j)
                fsilj = fsilj  +  wth_litho_wgh(k,j) * 1e12*areaclimber(j)
                !   WARNING: conversion: 1e6km2 -> m2: ^^^^
                weighted_rsw(kl) = weighted_rsw(kl) + rsw_litho(k)*wth_litho_wgh(k,j)*1e12*areaclimber(j)
                !FrivLi(kl) = FrivLi(kl) + reg_Li_Friv(k,j)*litho_frac(k,j)*1e12*areaclimber(j)   !mol/yr
                !dLiriv(kl) = dLiriv(kl) + reg_Li_driv(k,j)*(reg_Li_Friv(k,j)*litho_frac(k,j)*1e12*areaclimber(j))
            end do
            !
            ! Impose erosion of carbonates = erosion of sediments
            reg_eros(nlitho,j) = reg_eros(nlitho-1,j)
            !
            ! riverine sediment discharge
            tss(kl) = tss(kl) + TSS_density*sum(reg_eros(:,j)*litho_frac(:,j))*1e12*areaclimber(j)
            !
            ! Split silicates weathering in 2 fluxes: "basalt" weathering (mafic rocks) and "granite" weathering (rest of silicates)
            fbasj = wth_litho_wgh(BASALT_LITHO_NUM,j) * 1e12*areaclimber(j)
            fsilj = fsilj - fbasj
            fsilw(kl) = fsilw(kl) + fsilj*clo
            fbasw(kl) = fbasw(kl) + fbasj*clo
            !
            ! same for strontium isotopic ratio
            weighted_rsw(kl) = weighted_rsw(kl) - rsw_litho(BASALT_LITHO_NUM)*wth_litho_wgh(BASALT_LITHO_NUM,j)*1e12*areaclimber(j)


            ! =====================
            ! Kerogen weathering
            ! =====================
            ! Old formulation:
            ! facker=(-43000./Rgas)*(1./(Tclim(j)+273.15)-1./288.15)
            ! fker(j)=4.4184d+8*runclim(j)*areaclimber(j)   !*var_diss(6,nbasin)/38.d+18 !& !2.4 !1.432
                                                            ! *dexp(facker)  !0.672d+10
            fker(j) = 0.5 * sum(litho_frac(:,j)*OC_rock*reg_eros(:,j))
            !         ^^^ : oxidation efficiency
            if (.not. lock_oxygen_cycle)  fkerw(kl) = fkerw(kl) + fker(j)*1e12*areaclimber(j)*clo


            ! ======================
            ! Continental POC export
            ! ======================

            ! erosion export (t/km2/an)
            reg_eros_lithmean(j) = sum(litho_frac(:,j)*reg_eros(:,j)) / sum(litho_frac(:,j))
            !scaling:
            eros_galy_unit=reg_eros_lithmean(j)*TSS_density*1.e6*1.e-3  !moving from m/yr to t/km2/yr
            POC_export_rate(j)=0.081*(eros_galy_unit**0.56)  !in t/km2/yr
            POC_export_rate(j)=POC_export_rate(j)/12.         !mol/m2/yr
            total_cont_POC_export(kl) = total_cont_POC_export(kl) + POC_export_rate(j)*areaclimber(j)*1.e+12  !mol/yr


            ! ===============================================================================
            ! Pyrite weathering and distribution in carb dissol, sil dissol and H2SO4 release
            ! ===============================================================================

            ! molar flux of SO4^2- from pyrite oxidative weathering:
            pyr_wth(j) = sum(Sulf_rock(:) * reg_eros(:,j) * litho_frac(:,j))
            ! molar flux of Ca released to the ocean from carbonate dissolution by sulfuric acid
            fcarbsulfw(kl)  = fcarbsulfw(kl)    +  0.724*pyr_wth(j) * 1e12*areaclimber(j)
            ! molar flux of Ca released to the ocean from silicate dissolution by sulfuric acid
            fsilsulfw(kl)   = fsilsulfw(kl)     +  0.276*pyr_wth(j) * 1e12*areaclimber(j)
            ! molar flux of H2SO4 from pyrite oxidation directly reaching the ocean
            !fH2SO4sulfw(kl) = fH2SO4sulfw(kl)   +  0*pyr_wth(j)  *1e12*areaclimber(j) !=> All H2SO4 is neutralized

        end do


        ! Lithium riverine isotopic delta
        !dLiriv(:) = dLiriv(:)/FrivLi(:)

        ! lithology-averaged Strontium isotopic ratio
        where (fsilw /= 0d0)  weighted_rsw(:) = weighted_rsw(:)/fsilw(:)

    end if

end if




!=========================================================================!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! VARIABLES INDEPENDENT OF DYNSOIL USAGE !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! initialising:
fcarbw(:) = 0.
fpw(:)    = 0.

do j0=1,ncontpxl
    j = list_cont_pixel(j0)
    kl = cont_basin_map(j)

    ! =====================
    ! Carbonate weathering
    ! =====================
    !fcarbj=0.18577793d+11*dsqrt(runclim(j))*areaclimber(j)*fco2*fe !old
    call carbo(p,Tclim(j),runclim(j),carb_weath_conc,j)
    ! Carbonate is assumed to be the last lithological class (k = nlitho)
    wth_litho(nlitho,j) = k_carb_wth*runclim(j)*carb_weath_conc
    wth_litho_wgh(nlitho,j) = litho_frac(nlitho,j) * wth_litho(nlitho,j)
    fcarbj = 1d12*areaclimber(j) * wth_litho_wgh(nlitho,j)
    !        ^^^^^-- conversion because areaclimber is in 1e6km2 instead of m2:

    fcarbw(kl) = fcarbw(kl) + fcarbj*clo

    ! =====================
    ! Phosphorus weathering
    ! =====================
    ! => DO NOT UPDATE P WEATHERING RATES
    !fp(j) = sum((P_rock(1:nlitho-1)/CaMg_rock(1:nlitho-1)) * wth_litho_wgh(1:nlitho-1,j)) & ! silicate part
    !        +  P2C_carb * wth_litho_wgh(nlitho,j) & ! carbonate part
    !        +  P2C_ker * fker(j) ! kerogen part
    fpw(kl) = fpw(kl) + fp(j)*1d12*areaclimber(j)*clo*phosss

end do



! Substract Phosphorus that is in exported biospheric organic matter
! + reduce biospheric POC export if there is not enough Phosphorus (to avoid negative P fluxes)
where (total_cont_POC_export(:) > cp_cont*fpw(:))
    total_cont_POC_export(:) = cp_cont*fpw(:)
    fpw(:) = 0
else where
    fpw(:) = fpw(:) - total_cont_POC_export(:)/cp_cont
end where



! ****************************************************************** !
! Split global continental fluxes in multiple epicontinental basins  !
! proportionally to the basin areas in "unform routing" case         !
! ****************************************************************** !

if (uniform_routing) then
    ! note: for each variable, the entire flux was sent to the FIRST "appcont" box of the list (: kl)
    kl = jbox_appcont(1)
    do j0 = nappcont,1,-1 ! <-- reverse order so that the first "appcont" box is modified AT THE END OF THE LOOP
        j = jbox_appcont(j0)
        discharge(j)             = appcont_fract(j)*discharge(kl)
        tss(j)                   = appcont_fract(j)*tss(kl)
        fsilw(j)                 = appcont_fract(j)*fsilw(kl)
        fbasw(j)                 = appcont_fract(j)*fbasw(kl)
        fkerw(j)                 = appcont_fract(j)*fkerw(kl)
        total_cont_POC_export(j) = appcont_fract(j)*total_cont_POC_export(kl)
        fcarbsulfw(j)            = appcont_fract(j)*fcarbsulfw(kl)
        fsilsulfw(j)             = appcont_fract(j)*fsilsulfw(kl)
        fH2SO4sulfw(j)           = appcont_fract(j)*fH2SO4sulfw(kl)
        fpw(j)                   = appcont_fract(j)*fpw(kl)
        fcarbw(j)                = appcont_fract(j)*fcarbw(kl)
        FrivLi(j)                = appcont_fract(j)*FrivLi(kl)
    end do
end if




return
end

