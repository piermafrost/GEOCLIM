module multidimensional_interpolation
implicit none

contains


    !==============================================================================================!


    subroutine climate_interpolation(co2axis, param1, param2, param3, param4, param5, c, pvec, &
                                     boxtemp_array, interp_boxtemp, Fxch_array, interp_Fxch, &
                                     list_cont_pixel, ncontpxl, temp_array, runf_array, interp_temp, interp_runf)

        include 'shape.inc'
        !   => get data shape: nlon, nlat, nbasin, nclimber, nclimparam, len_p1, len_p2, len_p3, len_p4 and len_p5
        integer, parameter:: npxl=nlon*nlat
        include 'coupler.inc' ! => get variable 'CO2_interpolation'
        !
        double precision, intent(in):: co2axis(nclimber), param1(len_p1), param2(len_p2), param3(len_p3), param4(len_p4), &
                                       param5(len_p5)
        double precision, intent(in):: c, pvec(5)
        !
        double precision, intent(in), dimension(nbasin, nclimber, len_p1, len_p2, len_p3, len_p4, len_p5), optional:: boxtemp_array
        double precision, intent(out), dimension(nbasin), optional:: interp_boxtemp
        !
        double precision, intent(in), dimension(nbasin**2, nclimber, len_p1, len_p2, len_p3, len_p4, len_p5), optional:: Fxch_array
        double precision, intent(out), dimension(nbasin, nbasin), optional:: interp_Fxch
        !
        integer, intent(in), dimension(npxl), optional:: list_cont_pixel
        integer, intent(in), optional:: ncontpxl
        double precision, intent(in), dimension(npxl, nclimber, len_p1, len_p2, len_p3, len_p4, len_p5), optional:: temp_array, &
                                                                                                                    runf_array
        double precision, intent(out), dimension(npxl), optional:: interp_temp, interp_runf
        !
        integer, dimension(6):: subshp
        integer, dimension(6,2):: intidx
        double precision, dimension(2,2,2,2,2,2):: interp_coeff
        integer:: j0, j, i, k

        ! get interpolation coefficients and indices
        ! ------------------------------------------

        select case (CO2_interpolation)
            case ('linear')
                call get_multilininterp_coeff(co2axis, param1, param2, param3, param4, param5, &
                                              c, pvec(1), pvec(2), pvec(3), pvec(4), pvec(5), &
                                              subshp, intidx, interp_coeff)
            case ('log')
                call get_multilininterp_coeff(dlog(co2axis), param1, param2, param3, param4, param5, &
                                              dlog(c), pvec(1), pvec(2), pvec(3), pvec(4), pvec(5), &
                                              subshp, intidx, interp_coeff)
            case default
                print *
                print *, 'ERROR in "source/coupler.inc": illegal interpolation mode "'//CO2_interpolation//'"'
                print *, 'Expect "linear" or "log"'
                stop 31
        end select

        ! compute multinear interpolation
        ! -------------------------------

        if (present(list_cont_pixel) .and. present(ncontpxl)) then

            if (present(temp_array) .and. present(interp_temp)) then
                do j0 = 1,ncontpxl
                    j = list_cont_pixel(j0)
                    !
                    interp_temp(j) = sum( &
                        interp_coeff(1:subshp(1), 1:subshp(2), 1:subshp(3), 1:subshp(4), 1:subshp(5), 1:subshp(6)) &
                        * &
                        temp_array(j, intidx(1,1):intidx(1,2), intidx(2,1):intidx(2,2), intidx(3,1):intidx(3,2), &
                                      intidx(4,1):intidx(4,2), intidx(5,1):intidx(5,2), intidx(6,1):intidx(6,2)) &
                                        )
                end do
            end if

            if (present(runf_array) .and. present(interp_runf)) then
                do j0 = 1,ncontpxl
                    j = list_cont_pixel(j0)
                    !
                    interp_runf(j) = sum( &
                        interp_coeff(1:subshp(1), 1:subshp(2), 1:subshp(3), 1:subshp(4), 1:subshp(5), 1:subshp(6)) &
                        * &
                        runf_array(j, intidx(1,1):intidx(1,2), intidx(2,1):intidx(2,2), intidx(3,1):intidx(3,2), &
                                      intidx(4,1):intidx(4,2), intidx(5,1):intidx(5,2), intidx(6,1):intidx(6,2)) &
                                        )
                    ! avoid negative runoff in out-of range extrapolation
                    interp_runf(j) = max(0d0, interp_runf(j))
                end do
            end if

        end if

        if (present(boxtemp_array) .and. present(interp_boxtemp)) then
            do j = 1,nbasin
                interp_boxtemp(j) = sum( &
                    interp_coeff(1:subshp(1), 1:subshp(2), 1:subshp(3), 1:subshp(4), 1:subshp(5), 1:subshp(6)) &
                    * &
                    boxtemp_array(j, intidx(1,1):intidx(1,2), intidx(2,1):intidx(2,2), intidx(3,1):intidx(3,2), &
                                     intidx(4,1):intidx(4,2), intidx(5,1):intidx(5,2), intidx(6,1):intidx(6,2)) &
                                       )
            end do
        end if

        if (present(Fxch_array) .and. present(interp_Fxch)) then
            k = 0
            do j = 1,nbasin
                do i = 1,nbasin
                    k = k + 1
                    interp_Fxch(i, j) = sum( &
                        interp_coeff(1:subshp(1), 1:subshp(2), 1:subshp(3), 1:subshp(4), 1:subshp(5), 1:subshp(6)) &
                        * &
                        Fxch_array(k, intidx(1,1):intidx(1,2), intidx(2,1):intidx(2,2), intidx(3,1):intidx(3,2), &
                                      intidx(4,1):intidx(4,2), intidx(5,1):intidx(5,2), intidx(6,1):intidx(6,2)) &
                                           )
                end do
            end do
        end if



    end subroutine


    !==============================================================================================!


    subroutine get_multilininterp_coeff(vector1, vector2, vector3, vector4, vector5, vector6, val1, val2, val3, val4, val5, val6, &
                                        hypercube_subshape, interp_indices, hypercube_interp_coeff)

        double precision, intent(in), dimension(:):: vector1, vector2, vector3, vector4, vector5, vector6
        double precision, intent(in):: val1, val2, val3, val4, val5, val6
        integer, intent(out), dimension(6):: hypercube_subshape
        integer, dimension(6,2), intent(out):: interp_indices
        double precision, intent(out), dimension(2,2,2,2,2,2):: hypercube_interp_coeff
        double precision, dimension(6,2):: interp_coeff_1D
        integer:: i1, i2, i3, i4, i5, i6


        ! Find hybercube in the parameter spacing "framing" the current parameters value
        ! and compute 1-dimensional interpolation coefficients for each of the 6 dimensions 
        ! note: hypercube_subshape(i) = 1 if the i-th axis is length-1
        call find_interval(vector1, val1, interp_indices(1,1), interp_indices(1,2), interp_coeff_1D(1,2), hypercube_subshape(1))
        call find_interval(vector2, val2, interp_indices(2,1), interp_indices(2,2), interp_coeff_1D(2,2), hypercube_subshape(2))
        call find_interval(vector3, val3, interp_indices(3,1), interp_indices(3,2), interp_coeff_1D(3,2), hypercube_subshape(3))
        call find_interval(vector4, val4, interp_indices(4,1), interp_indices(4,2), interp_coeff_1D(4,2), hypercube_subshape(4))
        call find_interval(vector5, val5, interp_indices(5,1), interp_indices(5,2), interp_coeff_1D(5,2), hypercube_subshape(5))
        call find_interval(vector6, val6, interp_indices(6,1), interp_indices(6,2), interp_coeff_1D(6,2), hypercube_subshape(6))

        interp_coeff_1D(:,1) = 1 - interp_coeff_1D(:,2)
        !   => interpolation equation is (1-lambda)*vector(k) + lamda*vector(k+1)


        ! compute hypercube (6-dimensional) interpolation coefficient array
        ! ie, product of the 6 1-dimensional interpolation coefficients
        do i6 = 1,hypercube_subshape(6)
            do i5 = 1,hypercube_subshape(5)
                do i4 = 1,hypercube_subshape(4)
                    do i3 = 1,hypercube_subshape(3)
                        do i2 = 1,hypercube_subshape(2)
                            do i1 = 1,hypercube_subshape(1)
                                !
                                hypercube_interp_coeff(i1,i2,i3,i4,i5,i6) =  &
                                    interp_coeff_1D(1,i1)*interp_coeff_1D(2,i2)*interp_coeff_1D(3,i3)* &
                                    interp_coeff_1D(4,i4)*interp_coeff_1D(5,i5)*interp_coeff_1D(6,i6)
                                !
                            end do
                        end do
                    end do
                end do
            end do
        end do


    end subroutine
    

    !==============================================================================================!
    

    subroutine find_interval(vector, val, kinf, ksup, xi, interval_length)

        ! find the two consecutive elements of an *increasingly sorted* array 'vector' that frame
        ! 'val'. ie: find kinf and ksup such as vector(kinf) <= val < vector(ksup).
        ! If val is lower than the 1st element or higher than the last element, the first 2 (or
        ! last 2) indices are returned.
        ! Also returns 'xi', the linear interpolation coefficient

        double precision, intent(in), dimension(:):: vector
        double precision, intent(in)::               val
        integer, intent(out)::                       kinf, ksup
        double precision, intent(out)::              xi
        integer, intent(out), optional::             interval_length
        integer::                                    k

        if (size(vector)==1) then

            kinf = 1
            ksup = 1
            xi = 0d0
            if (present(interval_length)) interval_length = 1

        else

            kinf = 1
            ksup = size(vector)
            ! find by dichotomy:
            do while (ksup-kinf > 1)
                k = kinf + (ksup-kinf)/2
                if (val >= vector(k)) then
                    kinf = k
                else
                    ksup = k
                end if
            end do

            ! linear interpolation coefficient
            xi = (val - vector(kinf)) / (vector(ksup) - vector(kinf))

            if (present(interval_length)) interval_length = 2

        end if

    end subroutine


    !==============================================================================================!
    


end module
