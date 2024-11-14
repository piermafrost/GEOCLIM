    subroutine rtbis2(x1,x2,xacc,out)
!   ---------------------------------
    implicit none
    integer, parameter :: jmax=800
    double precision f2,xnothing,fmid,f1,dx,xmid,fxmid
    include 'combine_foam.inc'

    call phfunc2(x2,f2,xnothing)
    fmid=f2
    call phfunc2(x1,f1,xnothing)
    fph=f1
    if(fph*fmid.ge.0.) then
        write(6,*)x1,fph,x2,fmid
        write(6,*)'root must be bracketed:rtbis'
        stop 77
    endif
    if(fph.lt.0.) then
        out=x1
        dx=x2-x1
    else
        out=x2
        dx=x1-x2
    endif

    do j=1,jmax
        dx=dx*0.5
        xmid=out+dx
        call phfunc2(xmid,fxmid,xnothing)
        fmid=fxmid
        if(fmid.le.0.)out=xmid
        if(abs(dx).lt.xacc.or.fmid.eq.0.) return
    end do

!    pause 'too many bisection in rtbis2'
    end
