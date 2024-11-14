module read_oceanic_temperature_mod
implicit none
contains

subroutine read_oceanic_temperature(IOunit, Toceclimber, co2_axis, order, ndata)

    include 'shape.inc'
    integer, parameter:: psize = nclimber*len_p1*len_p2*len_p3*len_p4*len_p5
    integer, intent(in) :: IOunit
    double precision, dimension(nbasin, nclimber, len_p1, len_p2, len_p3, len_p4, len_p5), intent(out) :: Toceclimber
    double precision, dimension(nclimber), intent(in), optional :: co2_axis
    integer, dimension(psize, 6), intent(in), optional:: order
    integer, intent(in), optional :: ndata
    double precision, dimension(nclimber) :: SST, deepT
    double precision :: co2, temp(nbasin-1)
    integer :: i1, i2, i3, i4, i5, i6, j, k, n
    character(len=10) :: line

    ! Try reading temperature generating mode (string) in input file
    read(unit=IOunit, fmt=*) line

    if (line == 'parametric') then

        if (.not. present(co2_axis)) then
            print *
            print *, 'INTERNAL ERROR: in subroutine "read_oceanic_temperature",'
            print *, 'argument "co2_axis", needed to compute oceanic temperature ("parametric" mode),'
            print *, 'was not given to the subroutine'
            stop 421
        end if

        print *
        print *
        print *, 'Parametric CO2-dependent oceanic temperature' 
    
        ! Parametric CO2-temperature relationship
        ! ---------------------------------------

        !SST = 15. !6.8161*dlog(co2_axis*280.) - 22.089

        ! Pre-industrial
        SST = 20 + 4.88504*dlog(co2_axis)
        deepT = -6.58807 + 6.89026*dlog(co2_axis)
        ! Fixed oceanic temperature
        !SST = 20
        !deepT = -6.58807

        ! 265Ma
        !SST = 15.6279 + 4.88504*dlog(co2_axis)
        !deepT = -6.58807 + 6.89026*dlog(co2_axis)

        ! 365Ma
        !SST = 12.5071 + 5.98331*dlog(co2_axis)
        !deepT = -30.8490 + 16.9424*dlog(co2_axis)

        do k = 1,nclimber

            ! -----------------------------------------------------
            ! SURFACE:
            Toceclimber(1,k,:,:,:,:,:) = deepT(k)  !SST(k)-18
            Toceclimber(3,k,:,:,:,:,:) = SST(k)
            Toceclimber(6,k,:,:,:,:,:) = SST(k)
            Toceclimber(8,k,:,:,:,:,:) = deepT(k)  !SST(k)-18

            ! -----------------------------------------------------
            ! DEEP
            Toceclimber(2,k,:,:,:,:,:) = deepT(k)
            Toceclimber(4,k,:,:,:,:,:) = deepT(k)
            Toceclimber(5,k,:,:,:,:,:) = deepT(k)
            Toceclimber(7,k,:,:,:,:,:) = Toceclimber(6,k,:,:,:,:,:)
            Toceclimber(9,k,:,:,:,:,:) = deepT(k)

        end do


    else

        ! if input file does not contain the word "parametric",
        ! expect direct ascii-file temperature reading, one climatic parameter combination per line

        rewind(IOunit)

        if (present(order)) then

            print *
            print *
            write(*, fmt='(A,I0,A)') 'Line ordering in oceanic temperature input file (#',IOunit, &
                                     ') assumed to be the same than for continental climate'

            if (present(ndata)) then
                n = ndata
            else
                n = psize
            end if

            do k = 1,n
                i1 = order(k,1)
                i2 = order(k,2)
                i3 = order(k,3)
                i4 = order(k,4)
                i5 = order(k,5)
                i6 = order(k,6)
                read(unit=IOunit, fmt=*) co2, Toceclimber(1:nbasin-1, i1, i2, i3, i4, i5, i6)
                !                   skip atmospheric box -^^^^^^^^^^
            end do

        else

            print *
            print *
            print *,               'Oceanic temperature assumed to be dependent on CO2 ONLY'
            write(*, fmt='(A,I0)') '  -> expect 1 line per CO2 level, in decreasing order, in input file #',IOunit

            do k = nclimber,1,-1
                read(unit=IOunit, fmt=*) co2, temp
                do j = 1,nbasin-1 ! skip atmospheric box
                    Toceclimber(j,k,:,:,:,:,:) = temp(j)
                end do
                if (present(co2_axis)) then
                    ! Check that CO2 axis match first input of each line
                    ! (old security test, only valid in absence of any climatic parameters)
                    if (abs(co2/280. - co2_axis(k)) > 1d-6) then
                        print *
                        print *, 'Error while reading oceanic temperature'
                        print *, 'CO2 mismatch for level #', k
                        write(*, fmt='(A,F10.2,A,F10.2)') 'Got: ', real(co2), ', expect: ', 280*co2_axis(k)
                        stop 13
                    end if
                end if
            end do

        end if


    end if


    ! Security: avoid too low oceanic temperature because it makes the chemistry routines crash
    where (Toceclimber(1:nbasin-1,:,:,:,:,:,:) <= 1)  Toceclimber(1:nbasin-1,:,:,:,:,:,:) = 1d0


    ! Temperature conversion: deg C => Kelvin
    Toceclimber(1:nbasin-1,:,:,:,:,:,:) = Toceclimber(1:nbasin-1,:,:,:,:,:,:) + 273.15


    close(unit=IOunit)

end subroutine

end module
