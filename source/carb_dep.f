    subroutine carb_dep()
!   *********************
    use constante, only: akcr, akdiss, rarag
    implicit none
    include 'combine_foam.inc'


!   Reefal carbonate (aragonite) production on shelves [mol(CaCO3)/yr]
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    do j0=1,nappcont
        j = jbox_appcont(j0)
        if (omega_ara_0(j).ge.1.) then
            freef(j)  = akcr * surf_sedi(j) *(omega_ara_0(j)-1)**1.7 * clo * 1 !Toar=18.33  !Ceno=24.16
            freefP(j) = 0 !(1./1000.)*freef(j)*clo*phosss                      !Maas=16.66  !Berra=19.5
        else!           ^-- pas de P dans les coquilles                        !PT=22.0     !Carn=18.33
            freef(j)  = 0                                                      !Aptien=16.6 !Triasinf=13.33
            freefP(j) = 0                                                      !Rhetien=18.33
        endif
    end do



!   Open ocean Pelagic Carbonate dissolution rate [m^-3.yr^-1]
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!   (area_** expressed in percentage of total seafloor)
!   no CCD correction !!

    do j=1,nbasin-1

        ! obsolete => carbonate deposition computed as sink_veloc*[PIC]*surf_sedi
        !fdepc(j)=0.
        !fdepa(j)=0.

        if (dplysc(j) <= depth_top_box(j)) then
            f_diss_c(j) = akdiss*(1-rarag) ! dissolve on the full thickness
        elseif (dplysc(j) < depth_box(j)) then
            f_diss_c(j) = akdiss * (depth_box(j)-dplysc(j))/box_thick(j) * (1-rarag)
        else
            f_diss_c(j) = 0
        endif

        if (dplysa(j) <= depth_top_box(j)) then
            f_diss_a(j) = akdiss*rarag ! dissolve on the full thickness
        elseif (dplysa(j) < depth_box(j)) then
            f_diss_a(j) = akdiss * (depth_box(j)-dplysa(j))/box_thick(j) * rarag
        else
            f_diss_a(j) = 0
        endif

        fdissol_carb(j)  = f_diss_c(j) + f_diss_a(j)   !*shells ! fdownt
        fdissol_carbP(j) = 0 ! pas de P dans les coquilles !(1./1000.)*fdissol_carb(j)*shells

    end do

    end subroutine
