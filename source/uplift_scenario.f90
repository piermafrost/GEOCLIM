module uplift_scenario
implicit none

! Shape of the temporal uplift pulse = > lognormal function of time
! =================================================================

! ----------------------------------------------------------------------------------------------------- !
! Parameters meant to replicate "Baby et al., Terra Nova, 2019"'s sediment pulse on the Orange river fan,
! for a starting time (t=0) at 105 Ma, and ending time at 65 Ma
double precision, parameter:: t_peak = 15d6 ! (i.e., peak of the pulse at 90 Ma)
double precision, parameter:: sigma_pulse = 0.27
double precision, parameter:: pulse_amplitude = 15
double precision, parameter:: pulse_baseline = 1
double precision, parameter:: t_mu = t_peak * dexp(sigma_pulse**2)
! ----------------------------------------------------------------------------------------------------- !

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

