
    integer,parameter :: linesuv = 105

    integer :: iday,nday(8),iyear
    character(len=4) :: cyear,cday
    character(len=3) :: ccday
    character(len=1024) :: dir,fism_file
    real :: flux(linesuv)

    iyear   = 2002
    write(cyear,10) iyear
10  format(i4)

    print *,'iyear,cyear',iyear,cyear

    nday(1) = 1
    nday(2) = 11
    nday(3) = 132
    nday(4) = 33
    nday(5) = 34
    nday(6) = 35
    nday(7) = 36
    nday(8) = 37

    print *,'input iday '
    read *,iday

    if ( nday(iday) < 10 ) then
      write(cday,11) nday(iday)
      ccday = '00'//cday
      print *,'1',nday(iday)
    endif
    if ( 10 <= nday(iday) .AND. nday(iday) < 100 ) then
      write(cday ,12) nday(iday)
      ccday = '0'//cday
      print *,'2',nday(iday)
    endif
    if ( 100 <= nday(iday) ) then
      write(cday ,13) nday(iday)
      ccday = cday
      print *,'3',nday(iday)
    endif

    print *,'iday,nday',iday,nday(iday),cday,ccday

    dir       = './fism/'
    fism_file = trim(dir)//'FISM_'//trim(cyear)//'_'//ccday//'.inp'
    print *,'fism_file = ',trim(fism_file)
    open ( unit=65, file=fism_file  )  ! fism
          do i = 1,linesuv
            read (65,*) fism
            flux(i) = fism
          enddo
          print *,'max flux', maxval(flux)


11  format(i1)
12  format(i2)
13  format(i3)

    end
