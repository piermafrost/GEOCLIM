    subroutine derivs(x,y,dydx)
!   ******************************
    implicit none

    include 'combine_foam.inc'


    call varset(y)

!   Initialise les derives
    diff_diss = 0.
    diff_part = 0.
    diff_isot = 0.
    

!    hard code
!    if (x.ge.1.d+4.and.x.le.3d+4) then
!        sluggish=0.5
!    else
        sluggish=1.
!    endif

!   Calcul

!   ccccccccccccccccccccccccccccccccccccccccccccccccc
!   Dissolved variables - except oxygen (diss var #6)
!   ccccccccccccccccccccccccccccccccccccccccccccccccc
!   --------
    do i=1,5
        do j=1,nbasin
            do k=1,nbasin
                diff_diss(i,j)=diff_diss(i,j)+F(k,j)*var_diss(i,k)*sluggish  &
                                             -F(j,k)*var_diss(i,j)*sluggish
            enddo
            diff_diss(i,j)=diff_diss(i,j) + R_diss(i,j)
        enddo
    end do
!   --------
    do i=7,nvar_diss
        do j=1,nbasin
            do k=1,nbasin
                diff_diss(i,j)=diff_diss(i,j)+F(k,j)*var_diss(i,k)*sluggish  &
                                             -F(j,k)*var_diss(i,j)*sluggish
            enddo
            diff_diss(i,j)=diff_diss(i,j) + R_diss(i,j)
        enddo
    end do

!   cccccccccccccccccccccccccccccccccccccccccccc
!   Oxygen (diss var #6)
!   cccccccccccccccccccccccccccccccccccccccccccc
!   --------
    i=6 ! Oxygen (diss var #6)
    do j0=1,nnosurface
        j = jbox_nosurface(j0)
        do k=1,nbasin-1
            diff_diss(i,j)=diff_diss(i,j)+F(k,j)*var_diss(i,k)*sluggish  &
                                         -F(j,k)*var_diss(i,j)*sluggish
        enddo
        diff_diss(i,j)=diff_diss(i,j) + R_diss(i,j)
    end do
    diff_diss(i,nbasin)=diff_diss(i,nbasin) + R_diss(i,nbasin) ! atmosphere box
    do j0=1,nsurface !!! oxygN de surf. equilibr (diff_diss(i,j)=0)
        j = jbox_surface(j0)
        diff_diss(i,j)=R_diss(i,j)
    end do

!   cccccccccccccccccccccccccccccccccccccccccccc
!   Particulate variables
!   cccccccccccccccccccccccccccccccccccccccccccc
!   --------
    do i=1,nvar_part
        do j=1,nbasin
            diff_part(i,j)=diff_part(i,j) + R_part(i,j)
        enddo
    end do

!   cccccccccccccccccccccccccccccccccccccccccccc
!   Isotopic variables
!   cccccccccccccccccccccccccccccccccccccccccccc
!   --------
    do i=1,nvar_isot 
        do j=1,nbasin ! dissolved isotope are mixed in creades
            diff_isot(i,j)=diff_isot(i,j) + R_isot(i,j)
        end do
    end do


!   Acceleration factors
!   --------------------
    diff_diss(6,nbasin) = diff_diss(6,nbasin)*oxy_acc_fact ! atmospheric O2
    diff_diss(8,:) = diff_diss(8,:)*sulf_acc_fact          ! oceanic sulfate (all basins)


!   Put derivatives in a vector
!   ---------------------------

    k = 0
    do i = 1,nvar_diss
        do j = 1,nbasin
            k = k+1
            dydx(k) = diff_diss(i,j)
        end do
    end do
    do i = 1,nvar_part
        do j = 1,nbasin
            k = k+1
            dydx(k) = diff_part(i,j)
        end do
    end do
    do i = 1,nvar_isot
        do j = 1,nbasin
            k = k+1
            dydx(k) = diff_isot(i,j)
        end do
    end do

    end subroutine
