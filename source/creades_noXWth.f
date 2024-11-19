  subroutine creades_noXWth(x)
! ****************************

    use omp_lib
    use constante
    use multidimensional_interpolation, only: climate_interpolation
    use seafloor_sedimentation, only: org_dep
    use water_chemistry, only: salinity, ocean_atm_flu, atmospheric_pco2, sea_omega, dc13_speciation 

    implicit none
    include 'combine_foam.inc'
    double precision:: loc_sum_flx


    ! asynchronous coupling:
    if (icount_climate==ijump_climate) then
        icount_climate = 0

        ! climate interpolation (from look-up table): oceanic temperature, oceanic circulation, land temperature and runoff
        ! -----------------------------------------------------------------------------------------------------------------

        p=var_diss(7,nbasin)/PI_n_CO2_atm ! atmospheric CO2 level (in PAL) 
        !
        !if (variable_watxch) then  !  ! DO NOT UPDATE WATER EXCHANGES ["Fxch_array=Fxch_clim, interp_Fxch=F"]
        !    call climate_interpolation(co2climber, clim_param_1, clim_param_2, clim_param_3, clim_param_4, clim_param_5,        &
        !                               p, cpvec, list_cont_pixel=list_cont_pixel, ncontpxl=ncontpxl,                            &
        !                               boxtemp_array=Toceclimber, interp_boxtemp=temp_box, Fxch_array=Fxch_clim, interp_Fxch=F, &
        !                               temp_array=Tairclimber, runf_array=Runclimber, interp_temp=Tclim, interp_runf=runclim    )
        !else
            call climate_interpolation(co2climber, clim_param_1, clim_param_2, clim_param_3, clim_param_4, clim_param_5,    &
                                       p, cpvec, list_cont_pixel=list_cont_pixel, ncontpxl=ncontpxl,                        &
                                       boxtemp_array=Toceclimber, interp_boxtemp=temp_box) ! continental climate fields are not needed
                                       !temp_array=Tairclimber, runf_array=Runclimber, interp_temp=Tclim, interp_runf=runclim)
        !end if

        !call cont_weath(x) !<-- do no update weathering variables

    else
        icount_climate = icount_climate + 1
    end if
    !
    call salinity(salin)
    call ocean_atm_flu()
    call atmospheric_pco2()
    call diss_oxygen()
    call sea_omega()
    call biological_pump(x)
    call carb_dep()
    call dc13_speciation()
    call bio_frac()
    call recycle(x)
    call anoxic()
    call org_dep()
    call degassing(x)
    call Phydrotherm(x)
    call phosphorite(x)
    call strontium_ratio(x)


    ! + + + + + + + + + + + !
    ! Lock elemental cycles !
    ! + + + + + + + + + + + !
    !
    if (lock_sulfur_cycle) then
        ! => keep sil. sulf. wth (because it eventually consumes CO2), set H2SO4 release to 0 and
        !    use carb. sulf. wth as adjustment variable to balance sulfate-reduction
        loc_sum_flx = sum(fSulfRed(1:nbasin-1))
        do j0=1,nappcont
            j = jbox_appcont(j0)
            fH2SO4sulfw(j) = 0
            fcarbsulfw(j) = appcont_fract(j)*loc_sum_flx - fsilsulfw(j)
        end do
    end if
    !
    if (lock_oxygen_cycle) then
        loc_sum_flx = sum(fO2_odc(1:nbasin-1))
        do j0=1,nappcont
            j = jbox_appcont(j0)
            fkerw(j) = appcont_fract(j)*loc_sum_flx - (15d0/8d0)*(fcarbsulfw(j) + fsilsulfw(j) + fH2SO4sulfw(j))
        end do
    end if



    ! Initialization of reaction rates
    ! --------------------------------
    !##########!
    R_diss = 0.
    R_part = 0.
    R_isot = 0.
    !##########!




!<><><><><><><><><><><><><><>!
!<>  ===================   <>!
!<>  DISSOLVED VARIABLES   <>!
!<>  ===================   <>!
!<><><><><><><><><><><><><><>!


!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    i=1 ! DIC
!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc

    do j=1,nbasin-1
        R_diss(i,j) = fCO2atm_ocean(j) - fbioC(j) &
                    - finorgC(j)+fdissol_carb(j)*var_part(5,j)*box_vol(j) &
                    + roxyd(j)*box_vol(j)*var_part(4,j) &
                    + (1-clo) * frain_PIC(j) &
                    + (fin_POC(j) - clo*fodc(j)) & ! reoxydated org. C from sediment
                    + fmor(j) + fCO2_crust(j) - freef(j)
    end do
    ! Add continental inputs
    do j0=1,nappcont
        j = jbox_appcont(j0)
        R_diss(i,j) = R_diss(i,j)  +  2*fsilw(j) + 2*fbasw(j) + 2*fcarbw(j) + fkerw(j) + fcarbsulfw(j)
    end do


!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    i=2 ! ALK
!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc

    do j=1,nbasin-1
        R_diss(i,j) = -2*finorgC(j) &
                      + 2*fdissol_carb(j)*var_part(5,j)*box_vol(j) &
                      + (1-clo) * 2*frain_PIC(j) &
                      - 2*freef(j) - 2*fSO4_basin(j) - 2*fSO4_crust(j) + 2*fSulfRed(j)
    end do
    ! Add continental inputs
    do j0=1,nappcont
        j = jbox_appcont(j0)
        R_diss(i,j) = R_diss(i,j)  +  2*fsilw(j) + 2*fbasw(j) + 2*fcarbw(j) - 2*fH2SO4sulfw(j)
    end do


!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    i=3 ! PO4
!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc

    do j=1,nbasin-1
        R_diss(i,j) = -fbioP(j) - finorgP(j) + fdissol_carb(j)*var_part(3,j)*box_vol(j) &
                      + roxyd(j)*var_part(2,j)*box_vol(j) &
                      + (fin_POP(j) - clo*fodp(j)) & ! P from reoxydation of org. matter in sediment
                      - fhydP(j) - fphos(j)
    end do
    ! Add continental inputs
    do j0=1,nappcont
        j = jbox_appcont(j0)
        R_diss(i,j) = R_diss(i,j)  +  fpw(j)
    end do


!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    i=4 ! Ca
!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc

    do j=1,nbasin-1
        R_diss(i,j) = -finorgC(j) + fdissol_carb(j)*var_part(5,j)*box_vol(j) &
                      + (1-clo) * frain_PIC(j) &
                      - freef(j)
    end do
    ! Add continental inputs
    do j0=1,nappcont
        j = jbox_appcont(j0)
        R_diss(i,j) = R_diss(i,j)  +  fsilw(j) + fbasw(j) + fcarbw(j) + fcarbsulfw(j) + fsilsulfw(j)
    end do


!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    i=5 ! Sr
!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc

    do j=1,nbasin-1
        R_diss(i,j) = -finorgC(j)*rSrdep(j) + fdissol_carb(j)*var_part(1,j)*box_vol(j) &
                      + (1-clo) * sink_veloc*var_part(1,j)*surf_sedi(j) &
                      - freef(j)*rSrdep(j)
    end do
    ! Add continental inputs
    do j0=1,nappcont
        j = jbox_appcont(j0)
        R_diss(i,j) = R_diss(i,j)  +  (fsilw(j)+fbasw(j)+fsilsulfw(j))*rSrsil + (fcarbw(j)+fcarbsulfw(j))*rSrCar
    end do


!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    i=6 ! Oxygen
!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc

!   non-superficial basins
    do j0=1,nnosurface
        j = jbox_nosurface(j0)
        R_diss(i,j) = -rO2C*roxyd(j)*box_vol(j)*var_part(4,j) &
                      - (fin_POC(j)*rO2C - clo*fO2_odc(j)) ! reoxydation of org. C in sediment 
    end do

!   atmosphere
    R_diss(i,nbasin) = sum(total_cont_POC_export)*rO2C - sum(fkerw) - (15d0/8d0)*sum(fcarbsulfw+fsilsulfw+fH2SO4sulfw)

!   !##>> Impose equilibrium between atmosphere and superficial basins <<##!
    !
    !  => keep R_diss=0 in superficial basins
    !
    !  => add to the atmospheric box the net fluxes that should be in superficial basins
    !
    do k0=1,nsurface ! -> add to atmosphere net "reaction" O2 flux of surface basins
        k = jbox_surface(k0)
        R_diss(i,nbasin) = R_diss(i,nbasin) + rO2C*fbioC(k) & ! note: roxyd=0 in surface basins
                           - (fin_POC(k)*rO2C - clo*fO2_odc(k)) ! reoxydation of org. C in sediment 
    end do
    do k=1,nbasin-1 ! -> add to atmosphere net advection O2 flux of surface basins
        do j0=1,nsurface
            j = jbox_surface(j0)
            R_diss(i,nbasin) = R_diss(i,nbasin) + F(k,j)*var_diss(i,k) - F(j,k)*var_diss(i,j)
        end do
    end do



!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    i=7 ! atmospheric PCO2
!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc

    if (nclimber/=1) then ! keep null "R_diss" for in case of fixed-CO2 run (nclimber==1)
        do k=1,nbasin-1
            R_diss(i,nbasin) = R_diss(i,nbasin) - fCO2atm_ocean(k)
        end do
        R_diss(i,nbasin) = R_diss(i,nbasin) + fvol + ftrap + fanthros - 2*sum(fsilw) - 2*sum(fbasw) - sum(fcarbw) &
                           - sum(total_cont_POC_export)
    end if


!ccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccc
    i=8 ! SO4^2-
!ccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccc        

    ! Epicontinental box(es):
    do j0=1,nappcont
        j = jbox_appcont(j0)
        R_diss(i,j) = fcarbsulfw(j) + fsilsulfw(j) + fH2SO4sulfw(j)
    end do

    ! Sedimentary boxes
    do j0=1,nsedi
        j = jbox_sedi(j0)
        R_diss(i,j) = R_diss(i,j) - fSulfRed(j)
    end do




!<><><><><><><><><><><><><><><>!
!<>  =====================   <>!
!<>  PARTICULATE VARIABLES   <>!
!<>  =====================   <>!
!<><><><><><><><><><><><><><><>!


!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    i=1 ! SrPIC
!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc

    do j=1,nbasin-1
        R_part(i,j) = finorgC(j)*rSrdep(j)  &
                      - fdissol_carb(j)*var_part(i,j)*box_vol(j) &
                      - sink_veloc*var_part(i,j)*oce_surf(j)
    end do
    ! add sinking particules from the upper box (if there is one)
    do j0=1,nnosurface
        j = jbox_nosurface(j0)
        R_part(i,j) = R_part(i,j) + sink_veloc*var_part(i,j-1)*(oce_surf(j-1)-surf_sedi(j-1))
    end do


!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    i=2 !POP
!cccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccc

    do j=1,nbasin-1
        R_part(i,j) = fbioP(j) - roxyd(j)*box_vol(j)*var_part(i,j) - sink_veloc*var_part(i,j)*oce_surf(j)
    end do
    ! add continental inputs (in affected boxes)
    do j0=1,nappcont
        j = jbox_appcont(j0)
        R_part(i,j) = R_part(i,j) + total_cont_POC_export(j)/cp_cont
    end do
    ! add sinking particules from the upper box (if there is one)
    do j0=1,nnosurface
        j = jbox_nosurface(j0)
        R_part(i,j) = R_part(i,j) + sink_veloc*var_part(i,j-1)*(oce_surf(j-1)-surf_sedi(j-1))
    end do


!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    i=3 !PIP
!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc

    do j=1,nbasin-1
        R_part(i,j) = finorgP(j) - fdissol_carb(j)*var_part(i,j)*box_vol(j) - sink_veloc*var_part(i,j)*oce_surf(j)
    end do
    ! add sinking particules from the upper box (if there is one)
    do j0=1,nnosurface
        j = jbox_nosurface(j0)
        R_part(i,j) = R_part(i,j) + sink_veloc*var_part(i,j-1)*(oce_surf(j-1)-surf_sedi(j-1))
    end do


!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    i=4 ! POC
!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc

    do j=1,nbasin-1
        R_part(i,j) = fbioC(j) - roxyd(j)*box_vol(j)*var_part(i,j) - sink_veloc*var_part(i,j)*oce_surf(j)
    end do
    ! add continental inputs (in affected boxes)
    do j0=1,nappcont
        j = jbox_appcont(j0)
        R_part(i,j) = R_part(i,j) + total_cont_POC_export(j)
    end do
    ! add sinking particules from the upper box (if there is one)
    do j0=1,nnosurface
        j = jbox_nosurface(j0)
        R_part(i,j) = R_part(i,j) + sink_veloc*var_part(i,j-1)*(oce_surf(j-1)-surf_sedi(j-1))
    end do


!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    i=5 ! PIC
!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc

    do j=1,nbasin-1
        R_part(i,j) = finorgC(j) - fdissol_carb(j)*var_part(i,j)*box_vol(j) - sink_veloc*var_part(i,j)*oce_surf(j)
    end do
    ! add sinking particules from the upper box (if there is one)
    do j0=1,nnosurface
        j = jbox_nosurface(j0)
        R_part(i,j) = R_part(i,j) + sink_veloc*var_part(i,j-1)*(oce_surf(j-1)-surf_sedi(j-1))
    end do




!<><><><><><><><><><><><><>!
!<>  ==================  <>!
!<>  ISOTOPIC VARIABLES  <>!
!<>  ==================  <>!
!<><><><><><><><><><><><><>!


!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    i=1 ! DIC dc13
!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc

    do j=1,nbasin-1
        R_isot(i,j) = fc13ocean_atm(j) &
                      - fbioC(j) * (dh2co3(j) - epsiC(j) - var_isot(i,j)) &
                      - finorgC(j) * (dco3(j) - var_isot(i,j)) &
                      + fdissol_carb(j)*var_part(5,j)*box_vol(j) * (var_isot(i+1,j) - var_isot(i,j)) &
                      + roxyd(j)*box_vol(j)*var_part(4,j) * (var_isot(i+2,j) - var_isot(i,j)) &
                      + (1-clo) * frain_PIC(j) * (var_isot(i+1,j) - var_isot(i,j)) &
                      + (fin_POC(j) - clo*fodc(j)) * (var_isot(i+2,j) - var_isot(i,j)) & ! d13C flux from reoxydated org. C from sediment
                      + fmor(j)*(dcmor - var_isot(i,j)) &
                      - freef(j)*(dco3(j) - var_isot(i,j))
    end do

    ! add continental inputs (in affected boxes)
    do j0=1,nappcont
        j = jbox_appcont(j0)
        R_isot(i,j) = R_isot(i,j) &
                      + (2*fsilw(j) + 2*fbasw(j) + fcarbw(j)) * (var_isot(4,nbasin) - var_isot(i,j)) &
                      + (fcarbw(j) + fcarbsulfw(j)) * (dccarbw - var_isot(i,j)) &
                      + fkerw(j) * (dckerw - var_isot(i,j))
    end do

    ! Put advection part in "R" term (because of specificity of isotopic equation) and divide by "main" variable
    do j=1,nbasin-1
        do k=1,nbasin-1 ! advection terms
            R_isot(i,j) = R_isot(i,j)  +  F(k,j)*var_diss(1,k)*(var_isot(i,k) - var_isot(i,j))
        end do
        R_isot(i,j) = R_isot(i,j) / (var_diss(1,j)*box_vol(j)) ! divide by "main" variable
    end do


!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    i=2 ! PIC dc13
!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc

    do j=1,nbasin-1
        R_isot(i,j) = finorgC(j)*(dco3(j) - var_isot(i,j))
    end do
    ! add sinking particules from the upper box (if there is one)
    do j0=1,nnosurface
        j = jbox_nosurface(j0)
        R_isot(i,j) = R_isot(i,j) &
                      + sink_veloc*var_part(5,j-1)*(oce_surf(j-1)-surf_sedi(j-1)) * (var_isot(i,j-1) - var_isot(i,j))
    end do
    ! divide by "main" variable
    do j=1,nbasin-1
        !if (var_part(5,j).gt.1.e-6) then
            R_isot(i,j) = R_isot(i,j) / (var_part(5,j)*box_vol(j))
        !else
        !    R_isot(i,j) = 0
        !end if
    end do

    ! Note: no advection of particules


!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    i=3 ! POC dc13
!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc

    do j=1,nbasin-1
        R_isot(i,j) = fbioC(j)*(dh2co3(j) - epsiC(j) - var_isot(i,j)) 
    end do
    ! add continental inputs (in affected boxes)
    do j0=1,nappcont
        j = jbox_appcont(j0)
        R_isot(i,j) = R_isot(i,j) &
                      + total_cont_POC_export(j) * (var_isot(4,nbasin) - epsiCont - var_isot(i,j))
    end do
    ! add sinking particules from the upper box (if there is one)
    do j0=1,nnosurface
        j = jbox_nosurface(j0)
        R_isot(i,j) = R_isot(i,j) &
                      + sink_veloc*var_part(4,j-1)*(oce_surf(j-1)-surf_sedi(j-1)) * (var_isot(i,j-1) - var_isot(i,j))
    end do
    ! divide by "main" variable
    do j=1,nbasin-1
        !if (var_part(4,j).gt.1.e-6) then
            R_isot(i,j) = R_isot(i,j) / (var_part(4,j)*box_vol(j))
        !else
        !    R_isot(i,j) = 0
        !end if
    end do

    ! Note: no advection of particules


!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    i=4 ! PCO2 dc13
!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc

    R_isot(i,nbasin) = ( -sum(fC13atm_ocean(1:nbasin-1)) &
                         + fvol*(dcvol - var_isot(i,nbasin)) &
                         + ftrap*(dctrap - var_isot(i,nbasin)) &
                         + epsiCont*sum(total_cont_POC_export) &
                       ) /var_diss(7,nbasin)


!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    i=5 ! 87Sr/86Sr
!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc

    do j=1,nbasin-1
        R_isot(i,j) = fdissol_carb(j)*var_part(1,j)*box_vol(j) * (var_isot(i+1,j) - var_isot(i,j)) / (9.43 + var_isot(i+1,j)) &
                      + (1-clo) * sink_veloc*var_part(1,j)*surf_sedi(j) * (var_isot(i+1,j) - var_isot(i,j)) &
                                                                                                   / (9.43 + var_isot(i+1,j)) &
                      + fmor(j)*rSrmor * (rmor - var_isot(i,j)) / (9.43 + rmor)
                    
    end do

    ! add continental inputs (in affected boxes)
    do j0=1,nappcont
        j = jbox_appcont(j0)
        R_isot(i,j) = R_isot(i,j) &
                      + (fsilw(j) + fsilsulfw(j))*rSrSil * (weighted_rsw(j) - var_isot(i,j)) / (9.43 + weighted_rsw(j)) &
                      + fbasw(j)*rSrSil * (rbas - var_isot(i,j)) / (9.43 + rbas) &
                      + (fcarbw(j) + fcarbsulfw(j))*rSrCar * (rcw - var_isot(i,j)) / (9.43 + rcw)
    end do

    ! Put advection part in "R" term (because of specificity of isotopic equation) and divide by "main" variable
    do j=1,nbasin-1
        do k=1,nbasin-1 ! advection terms
            R_isot(i,j) = R_isot(i,j)  +  F(k,j) * var_diss(5,k) * (var_isot(i,k) - var_isot(i,j))/(9.43 + var_isot(i,k))
        end do
        R_isot(i,j) = R_isot(i,j) * (9.43 + var_isot(i,j)) / (var_diss(5,j)*box_vol(j)) ! divide by "main" variable
    end do


!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
        i=6 ! PIC 87Sr/86Sr
!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc

    do j=1,nbasin-1
        R_isot(i,j) = finorgC(j)*rSrdep(j) * (var_isot(i-1,j) - var_isot(i,j)) / (9.43 + var_isot(i-1,j))
    end do
    ! add sinking particules from the upper box (if there is one)
    do j0=1,nnosurface
        j = jbox_nosurface(j0)
        R_isot(i,j) = R_isot(i,j) &
                     + sink_veloc*var_part(1,j-1)*(oce_surf(j-1)-surf_sedi(j-1)) * &
                                             (var_isot(i,j-1) - var_isot(i,j)) / (9.43 + var_isot(i,j-1))
    end do
    ! divide by "main" variable
    do j=1,nbasin-1
        !if (var_part(1,j).gt.1.e-6) then
            R_isot(i,j) = R_isot(i,j) * (9.43 + var_isot(i,j)) / (var_part(1,j)*box_vol(j))
        !else
        !    R_isot(i,j) = 0
        !end if
    end do

    ! Note: no advection of particules




  end subroutine
