module uplift_scenario
implicit none

! Shape of the temporal uplift pulse = > lognormal function of time
! =================================================================

! -------------------------------------------------------------- !
! Shape parameters of log-normal uplift pulse
double precision, parameter:: t_peak = 3d6
double precision, parameter:: sigma_pulse = 0.3
double precision, parameter:: pulse_amplitude = 3
double precision, parameter:: pulse_baseline = 1
double precision, parameter:: t_mu = t_peak * dexp(sigma_pulse**2)
! -------------------------------------------------------------- !

contains

function pulse_shape(t)
    ! lognormal function, with module's parameters
    double precision, intent(in):: t
    double precision:: pulse_shape
    if (t==0d0) then
       pulse_shape = pulse_baseline
    else
        pulse_shape = pulse_baseline + (pulse_amplitude-1)*pulse_baseline * &
                    dexp(-(sigma_pulse**2)/2) * (t_mu/t) * dexp(-dlog(t/t_mu)**2 / (2*(sigma_pulse**2)))
    end if
end function


! ================= !
! Uplift subroutine !
! ================= !

subroutine uplift(t, init_slope, mask, curr_slope, missingpoints)

    include 'shape.inc' ! --> nlon, nlat
    integer, parameter:: npxl = nlon*nlat
    !
    double precision, intent(in):: t, init_slope(npxl), mask(npxl)
    double precision, intent(out):: curr_slope(npxl)
    logical, intent(in):: missingpoints(npxl)

    where (.not. missingpoints)
        curr_slope = (pulse_shape(t)*mask + pulse_baseline*(1-mask)) * init_slope
    end where

end subroutine

end module

