    subroutine Phydrotherm(t)
!   *************************
!   scavenging of P on iron oxides (hydrothermal)
    use constante, only: akPhyd
    implicit none
    include 'combine_foam.inc'

    do i=1,nbasin
        fhydP(i)=akPhyd*var_diss(3,i)*clo*phosss*indice_deep(i)
    enddo
    return
    end
