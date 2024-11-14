module climatic_parameters
implicit none

contains


    !=========================================================================!


    subroutine aperiodic_param_len(alen_p1, alen_p2, alen_p3, alen_p4, alen_p5)
        ! Determine the length of parameter removing 1 value for periodic ranges of parameter values
        include 'shape.inc'
        integer, intent(out):: alen_p1, alen_p2, alen_p3, alen_p4, alen_p5
        alen_p1 = len_p1
        alen_p2 = len_p2
        alen_p3 = len_p3
        alen_p4 = len_p4
        alen_p5 = len_p5
        if (p1_period/=0d0 .and. len_p1>1) alen_p1 = alen_p1-1
        if (p2_period/=0d0 .and. len_p2>1) alen_p2 = alen_p2-1
        if (p3_period/=0d0 .and. len_p3>1) alen_p3 = alen_p3-1
        if (p4_period/=0d0 .and. len_p4>1) alen_p4 = alen_p4-1
        if (p5_period/=0d0 .and. len_p5>1) alen_p5 = alen_p5-1
    end subroutine

    ! ---------- !

    subroutine check_axis_length(n, length, param_name)
        integer, intent(in):: n, length
        character(len=*), intent(in):: param_name
        if (n > length) then
            print *, 'ERROR: too many different values read in GCM configuration file for climatic parameter '//param_name
            print *, 'Length of current climatic parameter axis is ', length
            print *, 'Found ', n, ' different values'
            stop 89
        end if
    end subroutine

    ! ---------- !

    subroutine sort_single_param(param_list, param_axis, order, param_name)

        double precision, dimension(:), intent(in):: param_list
        double precision, dimension(:), intent(out):: param_axis
        integer, dimension(:), intent(out):: order
        character(len=*), intent(in):: param_name
        double precision:: p
        integer:: list_len, axis_len, k, imax, i, j
        logical:: assigned

        list_len = size(param_list)
        axis_len = size(param_axis)
        if (size(order)/=list_len) then
            print *, 'INTERNAL ERROR in module "climatic_parameters", subroutine "sort_single_param":'
            print *, 'ouput argument "order" must be of same size than input argument "param_list"'
            stop 421
        end if

        if (axis_len==1 .and. param_name/='CO2') then

            ! parameter does not exist, "axis_list" was not attributed, !expect for CO2!
            order = 1

        else

            imax = 0
            do k = 1,list_len

                p = param_list(k)
                assigned = .false.

                do i = 1,imax
                    if (.not. assigned) then
                        if (p == param_axis(i)) then

                            order(k) = i
                            assigned = .true.

                        elseif (p < param_axis(i)) then

                            imax = imax + 1
                            call check_axis_length(imax, axis_len, param_name)
                            param_axis(i+1:imax) = param_axis(i:imax-1)
                            param_axis(i) = p
                            assigned = .true.
                            do j = 1,k-1
                                if (order(j) >= i) order(j) = order(j) + 1
                            end do
                            order(k) = i

                        end if
                    end if
                end do

                if (.not. assigned) then
                    imax = imax + 1
                    call check_axis_length(imax, axis_len, param_name)
                    order(k) = imax
                    param_axis(imax) = p
                end if

            end do

        end if

    end subroutine

    ! ---------- !

    subroutine retrieve_clim_param_space(climparam_array, &
                                         CO2_levels, clim_param_1, clim_param_2, clim_param_3, clim_param_4, clim_param_5, &
                                         lp1, lp2, lp3, lp4, lp5, paramspace_filling_order, is_paramspace_filled)
        include 'shape.inc'
        integer, parameter:: paramspace_size = nclimber*len_p1*len_p2*len_p3*len_p4*len_p5
        !
        double precision, intent(in), dimension(paramspace_size,6):: climparam_array
        double precision, intent(out):: CO2_levels(nclimber), clim_param_1(len_p1), clim_param_2(len_p2), clim_param_3(len_p3), &
                                        clim_param_4(len_p4), clim_param_5(len_p5)
        integer, intent(out):: lp1, lp2, lp3, lp4, lp5
        integer, dimension(paramspace_size, 6), intent(out):: paramspace_filling_order
        logical, dimension(nclimber,len_p1,len_p2,len_p3,len_p4,len_p5), intent(out), optional:: is_paramspace_filled
        integer:: k, i1, i2, i3, i4, i5, i6, psize

        call aperiodic_param_len(lp1, lp2, lp3, lp4, lp5)
        psize = nclimber*lp1*lp2*lp3*lp4*lp5

        call sort_single_param(climparam_array(1:psize,1), CO2_levels,          paramspace_filling_order(1:psize,1), 'CO2')
        call sort_single_param(climparam_array(1:psize,2), clim_param_1(1:lp1), paramspace_filling_order(1:psize,2), '#1')
        call sort_single_param(climparam_array(1:psize,3), clim_param_2(1:lp2), paramspace_filling_order(1:psize,3), '#2')
        call sort_single_param(climparam_array(1:psize,4), clim_param_3(1:lp3), paramspace_filling_order(1:psize,4), '#3')
        call sort_single_param(climparam_array(1:psize,5), clim_param_4(1:lp4), paramspace_filling_order(1:psize,5), '#4')
        call sort_single_param(climparam_array(1:psize,6), clim_param_5(1:lp5), paramspace_filling_order(1:psize,6), '#5')

        if (present(is_paramspace_filled)) then
            is_paramspace_filled = .false.
            do k = 1,psize
                i1 = paramspace_filling_order(k,1)
                i2 = paramspace_filling_order(k,2)
                i3 = paramspace_filling_order(k,3)
                i4 = paramspace_filling_order(k,4)
                i5 = paramspace_filling_order(k,5)
                i6 = paramspace_filling_order(k,6)
                is_paramspace_filled(i1,i2,i3,i4,i5,i6) = .true.
            end do
        end if

        ! Periodic climatic parameter arrays
        if (lp1==len_p1-1)  clim_param_1(len_p1) = clim_param_1(1) + p1_period
        if (lp2==len_p2-1)  clim_param_2(len_p2) = clim_param_2(1) + p2_period
        if (lp3==len_p3-1)  clim_param_3(len_p3) = clim_param_3(1) + p3_period
        if (lp4==len_p4-1)  clim_param_4(len_p4) = clim_param_4(1) + p4_period
        if (lp5==len_p5-1)  clim_param_5(len_p5) = clim_param_5(1) + p5_period

    end subroutine


    !=========================================================================!


    function clim_param_modulo(param, lower_bound, period)
       double precision, intent(in):: param, lower_bound, period
       double precision:: clim_param_modulo
       clim_param_modulo = lower_bound + mod(param-lower_bound, period)
    end function

    ! ---------- !

    subroutine get_clim_param(t, tend, ijump_climparam, icount_climparam, prev_cpvec, next_cpvec, cpvec, cp_lowest_value, &
                              initialization)

        use utils, only: read_comment
        integer, parameter:: FILEID = 49 ! ID of climatic parameter file, open by suroutine "open_ascii_files"
        include 'shape.inc' ! => get variable "nclimparam" and "p*_periods"
        include 'coupler.inc' ! => get variables "climparam_loop" and "climparam_kill"
        double precision, intent(in):: t
        double precision, intent(inout):: tend
        integer, intent(in):: ijump_climparam
        integer, intent(inout):: icount_climparam
        double precision, intent(inout), dimension(5):: prev_cpvec, next_cpvec
        double precision, intent(out), dimension(5):: cpvec
        double precision, intent(in), dimension(5):: cp_lowest_value
        logical, intent(in), optional:: initialization
        logical:: loc_init
        character(len=1000):: line
        integer:: i_error
        double precision:: lam, p

        if (nclimparam > 0) then

            if (present(initialization)) then
                loc_init = initialization
            else
                loc_init = .false. ! default behaviour: no initialization
            end if


            if (loc_init) then

                ! only read starting climatic parameters
                read(unit=FILEID, fmt='(A)') line
                read(line, fmt=*) prev_cpvec(1:nclimparam)
                next_cpvec(1:nclimparam) = prev_cpvec(1:nclimparam)
                cpvec(1:nclimparam) = prev_cpvec(1:nclimparam)

            else

                ! (potentially) read next timestep' climatic parameters and interpolate "current" climatic parameters
                ! between current and next time steps

                if (icount_climparam == ijump_climparam) then

                    prev_cpvec(1:nclimparam) = next_cpvec(1:nclimparam)
                    cpvec(1:nclimparam) = prev_cpvec(1:nclimparam)
                    read(unit=FILEID, fmt='(A)', iostat=i_error) line

                    ! regular condition
                    if (i_error==0) then ! (no error)

                        read(line, fmt=*) next_cpvec(1:nclimparam)
                        icount_climparam = 0

                    ! end of file reached
                    elseif (i_error<0) then ! (end-of-file error)

                        if (climparam_loop) then ! restart from beginning of file
                            rewind(unit=FILEID)
                            call read_comment(FILEID) ! only accept comments at beginning of file
                            read(unit=FILEID, fmt='(A)') line
                            read(line, fmt=*) next_cpvec(1:nclimparam)
                            icount_climparam = 0
                        elseif (climparam_kill) then
                            tend = t ! signal to end GEOCLIM run
                        else
                            print *, 'reached end of climatic parameters file'
                            icount_climparam = ijump_climparam + 1 ! => program will never enter the "reading" if-loop again
                            ! note: "prev" and "next" cpvec are now identical => the linear interpolation will keep the last value
                        end if

                    ! other error
                    else

                        print *, 'error in "get_clim_param" (module "climatic_parameters") while reading climatic parameters file'
                        print *, 'IOstatus:', i_error
                        stop 97

                    end if

                else ! between 2 reading time steps

                    ! Linear interpolation between "previous" and "next" climatic parameters
                    lam = dble(icount_climparam)/dble(ijump_climparam)
                    cpvec(1:nclimparam) = (1-lam)*prev_cpvec(1:nclimparam) + lam*next_cpvec(1:nclimparam)

                    icount_climparam = icount_climparam + 1

                end if

            end if

            ! Climatic parameters with periodic ranges:
            if (p1_period /= 0d0)  cpvec(1) = clim_param_modulo(cpvec(1), cp_lowest_value(1), p1_period)
            if (p2_period /= 0d0)  cpvec(2) = clim_param_modulo(cpvec(2), cp_lowest_value(2), p2_period)
            if (p3_period /= 0d0)  cpvec(3) = clim_param_modulo(cpvec(3), cp_lowest_value(3), p3_period)
            if (p4_period /= 0d0)  cpvec(4) = clim_param_modulo(cpvec(4), cp_lowest_value(4), p4_period)
            if (p5_period /= 0d0)  cpvec(5) = clim_param_modulo(cpvec(5), cp_lowest_value(5), p5_period)

        end if

    end subroutine


    !=========================================================================!


end module
