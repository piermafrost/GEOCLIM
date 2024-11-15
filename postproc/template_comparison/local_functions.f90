module local_functions
implicit none

! Maximum allowed relative error
real, parameter:: MAX_ALLOWED_ERROR = 1e-5


  contains



    function compare_var1D(varname, var, var_ref)

        logical:: compare_var1D
        character(len=*), intent(in):: varname
        real, dimension(:), intent(in):: var, var_ref
        logical:: passed
        real:: xref, delta
        integer:: nt

        nt = size(var)

        xref = abs(sum(var_ref)) / nt
        delta = maxval(abs(var - var_ref))
        if (xref > 1e-12) delta = delta/xref
        compare_var1D = (delta <= MAX_ALLOWED_ERROR)

        ! Print message
        if (compare_var1D) then
            write(*, fmt='(A4,A50,ES10.3,A1)') '  * ', trim(varname)//': PASSED   (delta = ', delta, ')'
        else
            write(*, fmt='(A4,A50,ES10.3)')     '  * ', trim(varname)//': FAILED:   delta = ', delta
        end if
        write(*, fmt='(A)') '                          ----------'

    end function


    !--------!


    function compare_var2D(varname, var, var_ref, sedim)

        include 'shape.inc' ! => 'nbasin'
        integer, dimension(5), parameter:: sedim_basins = (/2, 5, 6, 7, 9/)

        logical:: compare_var2D
        character(len=*), intent(in):: varname
        real, dimension(:,:), intent(in):: var, var_ref
        logical, intent(in), optional:: sedim
        logical:: loc_sedim, passed
        real:: xref, loc_delta
        real, dimension(nbasin):: delta, delta_fail
        integer, dimension(nbasin):: failed
        integer:: j, j0, jend, nt, kfail

        nt = size(var, 2)
        kfail = 0

        if (present(sedim)) then
            loc_sedim = sedim
        else
            loc_sedim = .false.
        end if

        if (loc_sedim) then
            jend = size(sedim_basins)
            do j0 = 1,jend
                j = sedim_basins(j0)
                xref = abs(sum(var_ref(j,:))) / nt
                loc_delta = maxval(abs(var(j,:) - var_ref(j,:)))
                if (xref > 1e-12) loc_delta = loc_delta/xref
                passed = (loc_delta <= MAX_ALLOWED_ERROR)
                if (.not. passed) then
                    kfail = kfail + 1
                    failed(kfail) = j
                    delta_fail(kfail) = loc_delta
                end if
                delta(j0) = loc_delta
            end do
        else
            jend = nbasin-1 ! skip atmospheric box
            do j = 1,jend
                xref = abs(sum(var_ref(j,:))) / nt
                loc_delta = maxval(abs(var(j,:) - var_ref(j,:)))
                if (xref > 1e-12) loc_delta = loc_delta/xref
                passed = (loc_delta <= MAX_ALLOWED_ERROR)
                if (.not. passed) then
                    kfail = kfail + 1
                    failed(kfail) = j
                    delta_fail(kfail) = loc_delta
                end if
                delta(j) = loc_delta
            end do
        end if

        ! Print message
        if (kfail==0) then
            write(*, fmt='(A4,A50,ES10.3,A1)')     '  * ', trim(varname)//': PASSED   (delta = ', maxval(delta(1:jend)), ')'
            compare_var2D = .true.
        else
            write(*, fmt='(A4,A41,I0,A8,I0,A1)')   '  * ', trim(varname)//': FAILED  (', kfail, ' out of ', jend, ')'
            do j = 1,kfail
                write(*, fmt='(A24,I0,A9,ES10.3)') '       basin #', failed(j), 'delta = ', delta_fail(j)
            end do
            compare_var2D = .false.
        end if
        write(*, fmt='(A)') '                          ----------'

    end function


end module

