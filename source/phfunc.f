    subroutine phfunc(xo,fsoil,dfsoil,temper)
!   -----------------------------------------
!   This routine calculate pH for a solution at equilibrium with
!   soil PCO2 and carbonate rock. Used to estimate carbonate weathering

    !use constante, only: deltah0, xloak25
    use water_chemistry, only: eqcte
    implicit none

    include 'combine_foam.inc'


    salsoil=0.
    prsoil=0.
    
    call eqcte(temper,salsoil,prsoil,betasoil,ak1soil &
               ,ak2soil,akbsoil,akcsoil)

!   ak25=10**xloak25

!   aksp=ak25*exp(-deltah0/8.314*(1/temper-1/298))
!   aksp=akcsoil
    aksp=-0.38031396+0.0037566543*temper &
         -1.2063922d-5*temper*temper &
         +1.2696055d-8*temper*temper*temper
    if (temper.lt.278.15) then
        aksp=0.004460795   !to maintain aksp constant below 278.15 K
    endif 
    
    Hsoil=xo

    fsoil=2.*aksp*xo*xo*xo*xo/(ak1soil*ak2soil*betasoil*pco2soil) &
          -ak1soil*betasoil*pco2soil*xo &
          -2*ak1soil*ak2soil*betasoil*pco2soil &
          -2.*SO4soil*xo*xo
    dfsoil=2.*4*aksp*xo*xo*xo/(ak1soil*ak2soil*betasoil*pco2soil)- &
           ak1soil*betasoil*pco2soil-4*SO4soil*xo

    return
    end
