    subroutine rtbis(x1,x2,xacc,out,temper)
!   ---------------------------------------
    implicit none
    integer, parameter :: jmax=80
    double precision f1,f2,xnothing
    double precision dx,xmid,fmid,fxmid

    include 'combine_foam.inc'

    call phfunc(x2,f2,xnothing,temper)
    fmid=f2
    call phfunc(x1,f1,xnothing,temper)
    fsoil=f1
    if(fsoil*fmid.ge.0.) then
        write(6,*)x1,fsoil,x2,fmid
        write(6,*)'root must be bracketed:rtbis'
        stop 77
    endif
    if(fsoil.lt.0.) then
        out=x1
        dx=x2-x1
    else
        out=x2
        dx=x1-x2
    endif

    do j=1,jmax
        dx=dx*0.5
        xmid=out+dx
        call phfunc(xmid,fxmid,xnothing,temper)
        fmid=fxmid
        if(fmid.le.0.)out=xmid
        if(abs(dx).lt.xacc.or.fmid.eq.0.) return
    end do

    write(6,*) 'too many bisection in rtbis'
    return
    end
