    subroutine newton2(x,fx,dfx,ermaxa)
!    -----------------------------------
    implicit none
    integer num
    double precision x,fx,dfx,ermaxa,eror,fx0
    double precision fxfx0,residu,err


!   methode de newton-raphson amelioree

    eror=1.d-10

    call phfunc2(x,fx0,dfx)
    fx=fx0

    num=1

!   commencer par des pas doubles jusqu'a ce que fx*fx0<0

 1    continue
    x=x-2.*fx/dfx
    call phfunc2(x,fx,dfx)
    num=num+1
    fxfx0=fx*fx0
    if(num.ge.51)then
        write(74,*)'num=',num,' x=',x
        write(74,*)'fx=',fx,' dfx=',dfx
        write(74,*)'program stop'
        stop 77
    endif
    if(fxfx0.gt.0.)go to 1

!   continuer par newton-raphson classique

 2    continue
    residu=-fx/dfx
    x=x+residu
    call phfunc2(x,fx,dfx)
    num=num+1
    err=dabs(residu/x)
    if(num.ge.101)then
        write(74,*)'num=',num,' x=',x
        write(74,*)'fx=',fx,' dfx=',dfx
        write(74,*)'program stop'
        stop 77
    endif
    if(err.gt.eror)go to 2

    return
    end
