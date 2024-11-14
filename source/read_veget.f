subroutine read_veget
    implicit none
    include 'combine_foam.inc'
    integer:: i_error

    if (coupling_veget) then

        icount_veget = icount_veget + 1

        if (icount_veget >= ijump_veget) then

            read(unit=48,fmt=*,iostat=i_error) veget_occup

            if (i_error==0) then ! if no error

                where (veget_occup==1)
                    veget_factor = 1
                    veget_eros_factor = 1
                else where (veget_occup==2)
                    veget_factor = 1.0
                    veget_eros_factor = 1./2 
                end where

                icount_veget = 0

            elseif (i_error<0) then ! if end-of-file error

                print *, 'reached end of vegetation file'
                icount_veget = -10 * int( (tfi-tsta)/ts )
!                 => make sure icount_veget will never reach ijump_veget value again

            else ! other error

                print *, 'error in "read_veget.f" while reading vegetation file'
                print *, 'IOstatus:',i_error
                stop 97

            end if

        end if

    end if

end subroutine
