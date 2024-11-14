module water_chemistry

! ########################################################################### !
! This module contains the functions related to inorganic chemistry
! of seawater: thermodynamics constant, pH computation, carbonate speciation,
! gas dissolution, salinity...
! ########################################################################### !

implicit none

contains



    ! <><><><><><><><><> !
    ! Main 4 subroutines !
    ! <><><><><><><><><> !



    subroutine ocean_atm_flu()
!   ***************************
    use constante, only: dens
    include 'combine_foam.inc'

    do j=1,nbasin - 1 ! pour ne pas calculer ds l atmos

!       Thermodynamic constants:
        call eqcte(temp_box(j), salin(j), press_box(j), bco2(j), ak1(j), ak2(j), akb(j), akc_0(j))

!       carbonate speciation:
        call phbor(var_diss(1,j), var_diss(2,j), ch(j), salin(j), dens,ak1(j), ak2(j), akb(j))
        call chimie(var_diss(1,j), var_diss(2,j), ch(j), h2co3(j), hco3(j), co3(j), ak1(j), ak2(j))

        ph(j)           = -dlog10(ch(j)*1.d-3)
        pco2_dissous(j) = h2co3(j)/bco2(j)

    end do

    end subroutine



    ! ========================================= !



    subroutine atmospheric_pco2()
!   *****************************
    use constante, only: akk0, phias, phisa, PI_n_CO2_atm
    include 'combine_foam.inc'

    do j0=1,nsurface
        j = jbox_surface(j0)

        fCO2atm_ocean(j) = akk0*(var_diss(7,nbasin)/PI_n_CO2_atm  - pco2_dissous(j))*oce_surf(j)
        fC13atm_ocean(j) = akk0*(phias*var_diss(7,nbasin)/PI_n_CO2_atm - &
                                 (var_isot(1,j)-var_isot(4,nbasin) + phisa)*pco2_dissous(j)) &
                           *oce_surf(j)
        fC13ocean_atm(j) = akk0*((var_isot(4,nbasin) - var_isot(1,j) + phias)*var_diss(7,nbasin)/PI_n_CO2_atm - &
                                 phisa*pco2_dissous(j)) &
                           *oce_surf(j)

    enddo

    end subroutine



    ! ========================================= !



    subroutine dc13_speciation()
!   ****************************
    include 'combine_foam.inc'

    do j=1,nbasin-1
        edb(j) = (24.12-9866/temp_box(j))*1.d-3 !H2CO3 -> HCO3 
        ebc(j) = ((653.627/(temp_box(j)-233.45)**2) + 0.22)*1.d-3 ! HCO3 -> CO3
        dco3(j)   = ( var_diss(1,j)*var_isot(1,j)-ebc(j)*hco3(j) - ebc(j)*h2co3(j)-edb(j)*h2co3(j) ) &
                    / var_diss(1,j)
        dhco3(j)  = dco3(j)+ebc(j)
        dh2co3(j) = dhco3(j)+edb(j)
    end do

    end subroutine



    ! ========================================= !



    subroutine sea_omega()
!   **********************
    use constante, only: dens
    include 'combine_foam.inc'

    do j=1,nbasin-1

        ! Carbonate saturation at Pressure=0
        omega_0(j)     = var_diss(4,j)*co3(j)/akc_0(j)
        omega_ara_0(j) = omega_0(j)/1.5

        call lyso(omega_0(j), temp_box(j), salin(j), dens, dplysc(j), dplysa(j))

    enddo

    end subroutine



    ! ========================================= !
    ! ========================================= !



    ! <><><><><><><><><><><> !
    ! Auxiliary subroutines  !
    ! <><><><><><><><><><><> !



    subroutine lyso(omc_0, tem, sal, dens, dplysc, dplysa)
!   ------------------------------------------------------
!
!   Compute lysocline depth for calcite and aragonite given *current box* conditions
!   (temperature, salinity, ... and saturation at pressure=0)

    double precision, parameter:: grav = 9.806
    double precision, intent(in):: omc_0, tem, sal, dens
    double precision, intent(out):: dplysc,  dplysa
    double precision tcel, rt, cc, dvc, dkc, aa, bb, reali, pratmc, pratma

    tcel = tem - 273.15
    rt   = 82.056*tem

    ! Solve 2nd-order equation "aa*P^2 + bb*P + cc = 0"

    ! definition of 3 coeffs
    cc = -dlog(omc_0)
    !
    dvc = -(48.76 - 0.5304*tcel)
    dkc = -1.e-3*(11.76 - 0.3692*tcel)
    aa = 0.5*dkc/rt
    !
    bb = -dvc/rt

    ! Solving
    reali = bb*bb - 4.*aa*cc
    pratmc = (-bb + dsqrt(reali))/(2.*aa)

    ! Same for aragonite

    cc = -dlog(omc_0/1.5)
    !
    dvc = -(48.76 - 0.5304*tcel - 2.8)
    bb = -dvc/rt

    reali = bb*bb - 4.*aa*cc
    pratma = (-bb + dsqrt(reali))/(2.*aa)

    !,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    dplysc = 1.e-3*pratmc*1.013e5/(dens*grav)
    dplysa = 1.e-3*pratma*1.013e5/(dens*grav)
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    end subroutine



    ! ========================================= !



    function omega_with_pressure(omega_0, press, temp, aragonite)
!   -------------------------------------------------------------
!
!   Compute the calcite (or aragonite) saturation at a given pressure
!   given the saturation at pressure=0 (omega_0)

    double precision:: omega_with_pressure
    double precision, intent(in):: omega_0, press, temp
    logical, intent(in), optional:: aragonite
    double precision tcel, rt, dvc, dkc, aa, bb

    tcel = temp - 273.15
    rt   = 82.056*temp

    dvc = -(48.76 - 0.5304*tcel)
    if (present(aragonite)) then
        if (aragonite) dvc = dvc + 2.8
    end if
    dkc = -1.e-3*(11.76 - 0.3692*tcel)
    aa = 0.5*dkc/rt
    bb = -dvc/rt

    omega_with_pressure = omega_0*dexp(-aa*press**2 - bb*press)

    end function



    ! ========================================= !



    subroutine phbor(ct, alk, ah, sal, dens, ak1, ak2, akb)
!   -------------------------------------------------------
!
!   Compute pH by inversion of borate and carbonate acid-base equations

    double precision, intent(in):: ct, alk, sal, dens, ak1, ak2, akb
    double precision, intent(out):: ah
    double precision cl, tbor, ba, sa, a, b, c, fah, dfah, residu, err

    cl=(sal-0.030)/1.8050
    tbor=2.19e-5*cl*dens
    ba=tbor/alk
    sa=ct/alk
    a=(akb*(1.-ba)+ak1*(1.-sa))*1.e+5
    b=(akb*(1.-ba-sa)+ak2*(1.-sa-sa))*ak1*1.e+10
    c=(1.-ba-sa-sa)*ak1*ak2*akb*1.e+15

!   methode de newton-raphson amelioree

    ah=max(abs(c),1.+abs(b),1.+abs(a))
    fah=c+ah*(b+ah*(a+ah))
    dfah=b+ah*(a+a+3.*ah)

!   commencer par des pas doubles jusqua ce que fah<0

 1    continue
    ah=ah-2.*fah/dfah
    fah=c+ah*(b+ah*(a+ah))
    dfah=b+ah*(a+a+3.*ah)
    if(fah.gt.0.)go to 1

!   continuer par newton-raphson classique

 2    continue
    residu=-fah/dfah
    ah=ah+residu
    fah=c+ah*(b+ah*(a+ah))
    dfah=b+ah*(a+a+3.*ah)
    err=abs(residu/ah)
    if(err.gt.1.e-3)go to 2

    ah=ah*1.e-5

    end subroutine



    ! ========================================= !



    subroutine chimie(co2t, alk, ch, h2co3, hco3, co3, ak1, ak2)
!   ------------------------------------------------------------
      
    double precision, intent(in):: co2t, alk, ch, ak1, ak2
    double precision, intent(out):: h2co3, hco3, co3

    hco3  = co2t/(ch/ak1 + 1. + ak2/ch)
    co3   = ak2*hco3/ch
    h2co3 = hco3*ch/ak1

    end subroutine



    ! ========================================= !



    subroutine eqcte(tem, sal, pr, bc, ak1, ak2, akb, akc_0)
!   --------------------------------------------------------
!
!   this subroutine returns the equilibrium constants in mol/m3
!   for the dissociation of h2co3 (ak1), hco3 (ak2) and h3bo3 (akb)
!   as a function of temperature (tem) in kelvin, salinity (sal)
!   in per mil and pressure (pr) in bar. it also returns
!   henrys law constant bc for co2,(weiss (1974) for co2 solubility
!   in water, in zeebe & wolf-gladrow -co2 in seawater: equilibrium,
!   kinetics, isotopes-) in mol/m3/pal, and the saturation product
!   (akc_0) for calcite **at pr=0** in mol/m3*mol/m3.

    double precision, intent(in):: tem, sal, pr
    double precision, intent(out):: bc, ak1, ak2, akb, akc_0
    double precision:: ak0, pratm, tcel, rt, dv1, dk1, arg, dv2, dk2, dvb, dkb

    ak0 = -60.2409+9345.17/tem+23.3585*dlog(tem/100) + sal*(0.023517-0.00023656*tem+0.0047036*(tem*tem/10000))
    ak0 = dexp(ak0)*1.e+3
    bc = 280.e-6*ak0

    ak1 = 290.9097-14554.21/tem-45.0575*dlog(tem) + dsqrt(sal)*(0.0221+34.02/tem)
    ak1 = dexp(ak1)*1.e+3

    ak2 = 207.6548-11843.79/tem-33.6485*dlog(tem) + dsqrt(sal)*(0.9805-92.65/tem)-0.03294*sal
    ak2 = dexp(ak2)*1.e+3

    akb = 148.0248-8966.90/tem-24.4344*dlog(tem) + dsqrt(sal)*(0.5998-75.25/tem)-0.01767*sal
    akb = dexp(akb)*1.e+3

    akc_0 = 303.1308-13348.09/tem-48.7537*dlog(tem) + dsqrt(sal)*(1.6233-118.64/tem)-0.06999*sal
    akc_0 = dexp(akc_0)*1.e+6

    if(pr /= 0d0) then

        pratm = pr/1.013
        tcel  = tem - 273.15
        rt    = 82.056*tem

        dv1 = -(25.50+0.151*(sal-34.8)-0.1271*tcel)
        dk1 = -1.e-3*(3.08+0.578*(sal-34.8)-0.0877*tcel)
        arg = -dv1*pratm/rt + 0.5*dk1*pratm*pratm/rt
        ak1 = ak1*dexp(arg)

        dv2 = -(15.82-0.321*(sal-34.8)+0.0219*tcel)
        dk2 = -1.e-3*(-1.13+0.314*(sal-34.8)+0.1475*tcel)
        arg = -dv2*pratm/rt + 0.5*dk2*pratm*pratm/rt
        ak2 = ak2*dexp(arg)

        dvb = -(29.48-0.295*(sal-34.8)-0.1622*tcel+0.002608*tcel*tcel)
        dkb = -1.e-3*(2.84-0.354*(sal-34.8))
        arg = -dvb*pratm/rt + 0.5*dkb*pratm*pratm/rt
        akb = akb*dexp(arg)

    endif

    end subroutine



    ! ========================================= !



    subroutine salinity(salin)
!   **************************
    use constante, only: sal
    include 'shape.inc'
    double precision, dimension(nbasin), intent(out):: salin
    integer:: j

    do j=1,nbasin
        salin(j)=sal
    enddo

    end subroutine



end module
