    subroutine degassing(t)
!   ***********************
    use constante, only: Gt_to_mol, tstart_deg
    implicit none
    include 'combine_foam.inc'


    ftrap=0.
    !TRAP DEGASSING
    if (ipeak.eq.1) then
        do j=1,n_peaks    
            if (t.gt.tstart_deg+pulse_start(j).and.t.lt.tstart_deg+pulse_end(j)) then
                peak_duration(j)=pulse_end(j)-pulse_start(j)      !duration of one degassing peak
                ftrap=amount_peak(j)*Gt_to_mol/peak_duration(j)   !degassing in moles/yr
                dctrap=dc13_peak(j)
            endif
        end do
    endif


    return
    end
