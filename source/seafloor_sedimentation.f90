module seafloor_sedimentation
implicit none

contains


subroutine org_dep()

    use constante, only: rho_sed, hml, hsr, betahml, gammahsr, xmethan
    include 'combine_foam.inc'


    ! Compute sedimentation flux (m3/yr) and sedimenting POC flux (mol/yr), at seafloor, on each box
    ! ----------------------------------------------------------------------------------------------

    ! Initialization of sedim flux: 
    do j0=1,nsedi
        j = jbox_sedi(j0)
        fin_sed(j) = TSS(j)/rho_sed ! -> put terrestrial detritical sediments in epicontinental surface boxes
        fin_POC(j) = 0 !
        fin_POP(j) = 0 ! -> terrestrial POC|POP is already merged with oceanic POC|POP in water column
    end do

    ! ======================= !
    ! Sediment routing scheme !
    ! ======================= !

    ! 1. Shelf (epicontinental) boxes
    call drop_and_export(jbox_sediepicsurf, nsediepicsurf,                                                                 &
    !                    ____PIC_____   ____POC_____   ____POP_____
                         var_part(5,:), var_part(4,:), var_part(2,:), surf_sedi, sedbastype_fract, accumul_capacity, Xsed, &
                         frain_sed, frain_POC, frain_PIC, frain_POP, fin_sed, fin_POC, fin_POP, fexch_sed,                 &
                         joutbox=jbox_sediepicnosurf, noutbox=nsediepicnosurf, route_sedim=sed_routing_flag                )

    ! 2. Epicontinental deep boxes
    call drop_and_export(jbox_sediepicnosurf, nsediepicnosurf,                                                             &
    !                    ____PIC_____   ____POC_____   ____POP_____
                         var_part(5,:), var_part(4,:), var_part(2,:), surf_sedi, sedbastype_fract, accumul_capacity, Xsed, &
                         frain_sed, frain_POC, frain_PIC, frain_POP, fin_sed, fin_POC, fin_POP, fexch_sed,                 &
                         joutbox=jbox_sediintermed, noutbox=nsediintermed, route_sedim=sed_routing_flag                    )

    ! 3. Intermediate depth open ocean boxes
    call drop_and_export(jbox_sediintermed, nsediintermed,                                                                 &
    !                    ____PIC_____   ____POC_____   ____POP_____
                         var_part(5,:), var_part(4,:), var_part(2,:), surf_sedi, sedbastype_fract, accumul_capacity, Xsed, &
                         frain_sed, frain_POC, frain_PIC, frain_POP, fin_sed, fin_POC, fin_POP, fexch_sed,                 &
                         joutbox=jbox_sedideep, noutbox=nsedideep, route_sedim=sed_routing_flag                            )

    ! 4. Deep open ocean boxes
    call drop_and_export(jbox_sedideep, nsedideep,                                                                         &
    !                    ____PIC_____   ____POC_____   ____POP_____
                         var_part(5,:), var_part(4,:), var_part(2,:), surf_sedi, sedbastype_fract, accumul_capacity, Xsed, &
                         frain_sed, frain_POC, frain_PIC, frain_POP, fin_sed, fin_POC, fin_POP, fexch_sed                  )
                         !  do not export (infinite accumulation capacity)


!   Francois and Walker formalism, with TSS-dependent sedim rate:
!   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    do j0=1,nsedi
        j = jbox_sedi(j0)

!       sedimentation rate (m/yr) [from flux in m3/yr]
        ws(j) = fin_sed(j) / surf_sedi(j)
!       "fin_POC" is the org C input at the sediment top. Convert it: mol(C)/yr => mol(C)/m2/yr
!       Corg at the basis of the bioturbated layer (mixed layer):
        Corg_hml(j) = (fin_POC(j)/surf_sedi(j)) / (ws(j) + betahml*var_diss(6,j)*hml)
!       Corg at the basis of the sulfate reduction zone
        Corg_hsr(j) = (ws(j)*Corg_hml(j))/(ws(j) + gammahsr*var_diss(8,j)*hsr) ! var_diss(8,:) : [SO4^2-]

!       Convert in specific fluxes (per squared meters):
        fodc_noSR_m2 = ws(j)*Corg_hml(j)
        fodcm2       = ws(j)*Corg_hsr(j)
        ! Sulfate-Reduction flux:
        fSulfRed_m2  = (1d0/2d0)*(fodc_noSR_m2 - fodcm2)
        !              ^^^^^^^^^ 1 S reduced for 2 C oxidized

!       Oxygen flux: not affected by sulfate-reduction, except for the Fe^2+ part
        fO2_odc_m2   = fodc_noSR_m2 - (1d0/8d0)*fSulfRed_m2
!                                   ^^^^^^^^^^^^^^^^^^^^^^^   for 1 S reduced, 1/2 Fe^2+ released => 1/8 O2 eventually consumed

!       Substract C loss through methanogenesis 
        fO2_odc_m2   = fO2_odc_m2 - (1-xmethan)*fodcm2 ! leaking methane will eventually be oxidized by O2
        fodcm2       = xmethan*fodcm2

!       Total fluxes:
        fodc(j)     = fodcm2*surf_sedi(j)*clo
        fodp(j)     = fodc(j)/cp_burial(j)*phosss
        fSulfRed(j) = fSulfRed_m2*surf_sedi(j)*clo
        fO2_odc(j)  = fO2_odc_m2*surf_sedi(j)*clo

        ! Burial efficiency of C
        Corg_BE(j) = fodc(j)/fin_POC(j)

    end do


end subroutine



! ======================================== !



subroutine drop_and_export(jbox, nbox, PIC, POC, POP, surf_sedi, sedbastype_fract, accumul_capacity, sed_routing_matrix, &
                           frain_sed, frain_POC, frain_PIC, frain_POP, fin_sed, fin_POC, fin_POP, sed_exch_flux, &
                           joutbox, noutbox, route_sedim)

    use constante, only: sink_veloc, PIC_mmass, POC_mmass, rho_sed
    include 'shape.inc'
    integer, intent(in):: jbox(nbasin), nbox
    double precision, dimension(nbasin), intent(in):: PIC, POC, POP, surf_sedi, sedbastype_fract, accumul_capacity
    double precision, dimension(nbasin, nbasin), intent(in):: sed_routing_matrix
    double precision, dimension(nbasin), intent(inout):: frain_sed, frain_POC, frain_PIC, frain_POP, fin_sed, fin_POC, fin_POP
    double precision, dimension(nbasin, nbasin), intent(inout):: sed_exch_flux
    integer, intent(in), optional:: joutbox(nbasin), noutbox
    logical, intent(in), optional:: route_sedim
    double precision:: x, fout_sed, fout_POC, fout_POP, fPOC_jk, fPOP_jk
    integer:: j, j0, k, k0
    logical:: export, loc_route_sedim


    ! Options
    ! -------

    export = (present(joutbox) .and. present(noutbox))

    if (present(route_sedim)) then
        loc_route_sedim = route_sedim
    else
        loc_route_sedim = .false.
    end if


    do j0 = 1,nbox
        j = jbox(j0)

        ! Compute "raining" fluxes (i.e., particles from the water column) and add them to "upslope" fluxes
        ! -------------------------------------------------------------------------------------------------

        ! bulk sediment (m3)
        frain_sed(j) = sink_veloc*(PIC_mmass*PIC(j) + POC_mmass*POC(j))*surf_sedi(j) / rho_sed
        ! organic matter (mol)
        frain_POC(j) = sink_veloc*POC(j)*surf_sedi(j)
        ! particulate inorganic carbon  (mol) [All "raining" PIC is buried]
        frain_PIC(j) = sink_veloc*PIC(j)*surf_sedi(j)
        ! Phosphorus associated with organic matter (mol)
        frain_POP(j) = sink_veloc*POP(j)*surf_sedi(j)

        ! > add raining fluxes to fluxes from upslope (initial "fin" values)
        fin_sed(j) = fin_sed(j) + frain_sed(j)
        fin_POC(j) = fin_POC(j) + frain_POC(j)
        fin_POP(j) = fin_POP(j) + frain_POP(j)
        ! no computation for PIC (all is buried)


        if (export) then

            ! Compute staying and exported fractions
            ! --------------------------------------

            ! Michaelis-like saturation: export what exceeds accumulation capacity
            x = 1  -  1 / (1 + fin_sed(j)/accumul_capacity(j)) ! <-- exported fraction

            ! Exported sediment fluxes
            fout_sed = x*fin_sed(j)
            fout_POC = x*fin_POC(j)
            fout_POP = x*fin_POP(j)

            ! Sediment fluxes that stay in the box's sediment
            fin_sed(j) = fin_sed(j) - fout_sed
            fin_POC(j) = fin_POC(j) - fout_POC
            fin_POP(j) = fin_POP(j) - fout_POP


            ! Export surplus to downslope ocean basins
            ! ----------------------------------------

            do k0 = 1,noutbox
                k = joutbox(k0)

                ! Compute export flux from box #j to box #k
                ! - - - - - - - - - - - - - - - - - - - - -

                if (loc_route_sedim) then !--> use routing matrix to export sedim
                                          !    ("sed_routing_matrix" represents the FRACTION of exported fluxes)

                    sed_exch_flux(j,k) = sed_routing_matrix(j,k) * fout_sed
                    fPOC_jk            = sed_routing_matrix(j,k) * fout_POC
                    fPOP_jk            = sed_routing_matrix(j,k) * fout_POP

                else !--> "spread" exported sedim to all "immediate" downslope basins
                     !    (weighted by their area fraction "sedbastype_fract")

                    sed_exch_flux(j,k) = sedbastype_fract(k) * fout_sed
                    fPOC_jk            = sedbastype_fract(k) * fout_POC
                    fPOP_jk            = sedbastype_fract(k) * fout_POP

                end if

                ! Put summed exported fluxes as "next" boxes initial value of "fin"
                fin_sed(k) = fin_sed(k) + sed_exch_flux(j,k)
                fin_POC(k) = fin_POC(k) + fPOC_jk
                fin_POP(k) = fin_POP(k) + fPOP_jk

            end do

        end if

    end do

end subroutine



end module
