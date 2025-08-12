! WEIMER_GRID.F

    subroutine weimer_grid

    use parameter_mod

    real :: weimer_lat(nf+1),weimer_lon(nlt+1)
    real :: blatpt(nzp1,nfp1,nlt)
    real :: blonpt(nzp1,nfp1,nlt)

    open ( unit=76, file='blatpu.dat',form='unformatted' )
    open ( unit=77, file='blonpu.dat',form='unformatted' )

    open ( unit=78, file='weimer_grid.dat',form='formatted')
!    open ( unit=79, file='weimer_lat.dat',form='formatted')
!    open ( unit=80, file='weimer_lon.dat',form='formatted')

    read(76) blatpt
    read(77) blonpt

    close(76)
    close(77)

    do j = 1,nf+1
!        print *,'j,blat ',j,blatpt(nz-1,j,1)
        weimer_lat(j) = blatpt(nz-1,j,1)
    enddo

    do k = 1,nlt
!        print *,'k,blon ',k,blonpt(1,1,k)
        weimer_lon(k) = blonpt(1,1,k)
    enddo

    weimer_lon(nlt+1) = weimer_lon(1) + 360.

    write(78,*) weimer_lat,weimer_lon

!    write(79,*) weimer_lat
!    write(80,*) weimer_lon

    close(78)
!    close(79)
!    close(80)

    end subroutine weimer_grid

