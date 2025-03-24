module dynsoil_lithium
implicit none

!##################################################################################################!
!#                                                                                                #!
!# LITHIUM MODEL FROM "Caves Rugentstein et al., Nature, 2019,                                    #!
!#                     'Neogene cooling driven by land surface reactivity ...'                    #!
!#                     https://doi.org/10.1038/s41586-019-1332-y"                                 #!
!#                                                                                                #!
!#   >> approximate "weathering intensity" as "1 - primary phases proportion at top of regolith"  #!
!#                                                                                                #!
!##################################################################################################!


! =================================== !
! Caves Rugenstein's model parameters !
! =================================== !

! Weathering Intensity coefficients for d7Li riv (from Caves Rugenstein et al., 2019)
double precision, parameter:: a_WI_CR19 = 4.5d-3 !(dimensionless, NOT permil!)
double precision, parameter:: b_WI_CR19 = 0.0575 !(note: is actually 10 times lower than indicated in CR2019)
! Fractionation factor from Caves Rugenstein et al. (2019)
double precision, parameter:: Dsec_CR19 = 17d-3 !(dimensionless, NOT permil!)


contains


! ================================== !
! Caves Rugenstein's model equations !
! ================================== !

function d_d7Li_riv(wi)
    ! Compute "riverine (dissolved) d7Li MINUS bedrock d7Li" based on weahtering intensity (wi)
    double precision, intent(in):: wi
    double precision:: d_d7Li_riv
    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
    d_d7Li_riv = ((1 - wi)/wi) * a_WI_CR19 * dexp(-1*b_WI_CR19/wi)
    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
end function

function f_Li_diss(dd7Li_r)
    ! Compute the fraction of weathered Li that is not retained in secondary phases (ie, dissolved),
    ! based on "dd7Li_r" (see above), with Rayleigh distillation assumption (Caves Rugenstein et al., 2019)
    double precision, intent(in):: dd7Li_r
    double precision:: f_Li_diss
    !<><><><><><><><><><><><><><><><><>!
    f_Li_diss = dexp(-dd7Li_r/Dsec_CR19)
    !<><><><><><><><><><><><><><><><><>!
end function


! ============ !
! Main routine !
! ============ !

subroutine lithium(reg_prod, reg_P_diss, list_cont_pxl, cont_basin_map, ncontpxl, lithfrac, area, Li_Friv, Li_driv)

  use dynsoil_physical_parameters, only: rho_UCC, Li_UCC, dLi_UCC
  use constante, only: CaMg_rock
  include 'shape.inc' ! => nlon, nlat, nlitho, nbasin
  integer, parameter:: npxl = nlon*nlat

  ! parameter to convert volumetric weathering into massic TDS weahtering
  double precision, parameter, dimension(nlitho):: TDS_sil = (550d9/2500d9) * CaMg_rock
  !  Ratio between TDS sil. (kg/yr) and Ca-Mg sil. (mol/yr) --^^^^^^^^^^^^
  !  (from Gaillardet et al., 1999)

  double precision, dimension(nlitho,npxl), intent(in)::  reg_prod, reg_P_diss, lithfrac
  double precision, dimension(npxl),        intent(in)::  area
  double precision, dimension(nbasin),      intent(out):: Li_Friv, Li_driv
  integer,          dimension(npxl),        intent(in)::  list_cont_pxl, cont_basin_map
  integer,                                  intent(in)::  ncontpxl
  !
  double precision, dimension(nbasin):: wgh_ave_WI, denud_flux, Li_wth_flux
  double precision:: chem_den
  integer:: i, j, j0, kl


  ! Compute, on each watershed,
  !   the integrated lithium "congruent" weathering flux (before incorporation in secondary phases),
  !   the erosion-weighted average weathering intensity (with the approximation: WI ~= 1 - xP_surf)
  ! ------------------------------------------------------------------------------------------------

  Li_wth_flux = 0
  denud_flux  = 0
  wgh_ave_WI  = 0

  do j0 = 1,ncontpxl
    j  = list_cont_pxl(j0)
    kl = cont_basin_map(j)

    do i = 1,nlitho-1 ! Skip last lithology class => carbonates

      chem_den  =  reg_P_diss(i,j) * lithfrac(i,j) * 1d12*area(j)
      !
      Li_wth_flux(kl)  =  Li_wth_flux(kl)  +  chem_den
      wgh_ave_WI(kl)   =  wgh_ave_WI(kl)   +  TDS_sil(i) * chem_den
      denud_flux(kl)   =  denud_flux(kl)   +  rho_UCC * reg_prod(i,j) * lithfrac(i,j) * 1d12*area(j)

    end do
  end do

  do kl = 1,nbasin
    if (Li_wth_flux(kl) /= 0d0) then

      ! conversion in "true" units
      Li_wth_flux(kl) = Li_UCC * rho_UCC * Li_wth_flux(kl)
      wgh_ave_WI(kl) = wgh_ave_WI(kl) / denud_flux(kl)

      ! "Actual" Li dissolved (riverine) flux and d7Li, on each watershed
      ! => Caves Rugenstein lithium model
      ! -----------------------------------------------------------------
      Li_driv(kl) = d_d7Li_riv(wgh_ave_WI(kl))
      Li_Friv(kl) = Li_wth_flux(kl) * f_Li_diss(Li_driv(kl))

    else

      Li_driv(kl) = 0
      Li_Friv(kl) = 0

    end if
  end do

end subroutine


end module
