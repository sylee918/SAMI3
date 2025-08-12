c this must be linked to Geopack.o (Geopack.f from NASA Goddard)
c This verson takes quantites from the ACE OMNI dataset (level 2 combined 1-min data)

      PARAMETER (nwind0=525601)

      CHARACTER*32 pwind,pwine
      real bygse(nwind0),bzgse(nwind0),vxgse(nwind0),denin(nwind0)
      real sec(nwind0),timein(nwind0),tilt(nwind0)
      real hrutd(nwind0),by(nwind0),bz(nwind0),vx(nwind0),den(nwind0)
      integer ihr(nwind0),minu(nwind0),id(nwind0)
      integer imonth(nwind0),iyear(nwind0)
      integer :: nlines = 525601


      COMMON /GEOPACK1/AA(16)
c
      earthr = 6380.

C Output to be used by the Weimer model
!      pwine='imf_for_weimer.txt'
       pwine='weimer_data.inp'

c *** These next few lines need to be filled in by hand ***
c     ** Downloaded data **
!      pwind = 'OMNI_HRO_1MIN_2009_03_01.txt'
!      pwind = 'OMNI_HRO_1MIN_2012_06_22.txt'
!      pwind = 'OMNI_HRO_1MIN_2001_03_30.txt'
      pwind = 'OMNI.txt'
c     ** Spacraft distance from Earth in Re **
c     ** (240 for ACE, 0 for "OMNI")        **
      xcraft = 0.0

c     ** Input header lines to skip **
      npskip = 70
c     ** Line number of last line of data **
!      npwind = 4391
!      npwind = 1511  ! one day
!      npwind = 7271  ! five day
!      npwind = 10151  ! seven day

c     Open IMF file
!      open(33,FILE=pwind,STATUS='OLD',FORM='FORMATTED')
!      do i = 1,npskip 
!        read(33,*)
!      enddo

c     ** Year and day of year corresponding to first line of data **
!      iyr = 2001
!      iday0 = 89
!      iyr = 2013
!      iday0 = 1

c     ** Here we set the cadence and the window width (minutes) **
      cadence = 15.0
!      window = 20.0
      window = 30.

c     Open IMF file
      open(33,FILE=pwind,STATUS='OLD',FORM='FORMATTED')
      do i = 1,npskip 
        read(33,*)
      enddo

!      open(33,file='OMNI.txt',form='formatted')

      imonth_1 = 12
      iday_1   = 1
      imonth_2 = 12
      iday_2   = 6
      iyear0   = 2002

      call julian_day(iday_1,imonth_1,iyear0,iday0)

      do i = 1,nlines
        read(33,13) id(i),imonth(i),iyear(i),ihr(i),minu(i),sec(i),
     .        isb,isp,bygse(i),bzgse(i),vxgse(i),denin(i)
        if ( imonth(i) == imonth_1 .AND. id(i) == iday_1 
     .         .AND. iyear(i) == iyear0
     .         .AND. ihr(i) == 0 .AND. minu(i) == 0 ) then
          i_start = i
          print *,'found start',i_start,imonth(i),id(i)
        endif
        if ( imonth(i) == imonth_2 .AND. id(i) == iday_2 
     .         .AND. iyear(i) == iyear0
     .         .AND. ihr(i) == 0 .AND. minu(i) == 0 ) then
          i_end = i
          print *,'found end',i_end,imonth(i),id(i)
        endif
      enddo

      print *,'i_start,i_end',i_start,i_end

 13   format(i2,1x,i2,1x,i4,1x,i2,1x,i2,1x,f6.3,i12,i15,
     .       2f14.5,f17.3,f15.5,f14.3)

      close(33)

c     re-Open IMF file
      open(33,FILE=pwind,STATUS='OLD',FORM='FORMATTED')
      do i = 1,npskip 
        read(33,*)
      enddo

c     Read the actual data
c     (in this code we have bygse, etc, but the OMNI By,Bz data comes in GSM, so it is GSM)
!      npread = npwind - npskip
      iday   = 0
!      npread = i_end - i_start
!
      do ii=i_start,i_end
        i = ii - i_start + 1
!        print *,'i,ii',i,ii
        read(33,3) id(i),ihr(i),minu(i),sec(i),isb,isp,
     @             bygse(i),bzgse(i),vxgse(i),denin(i)

c       CDAWeb sometimes changes the formatting and
c       so this may need to be adjusted 

!    3   format(i2,9x,i2,1x,i2,1x,f6.3,i12,i15,2f14.5,f17.3,f15.5)
!    3   format(i2,1x,i2,1x,i2,1x,f6.3,i12,i15,2f14.5,f17.3,f15.5,f14.3
!     .         ,f14.5,f14.4)

    3   format(i2,9x,i2,1x,i2,1x,f6.3,i12,i15,2f14.5,f17.3,f15.5,f14.3
     .         ,f14.5,f14.4)
 
c EPOCH_TIME              IMF_S/C_ID# PLASMA_S/C_ID#       BY,_GSM       BZ,_GSM VX_VELOCITY,_GSE PROTON_DENSITY
c dd-mm-yyyy hh:mm:ss.ms                                        nT            nT             km/s           n/cc
c 30-01-2001 00:00:00.000          71             71       1.07000       4.13000         -455.300        6.01000

!        if(i .eq.1) write(*,*) ' h:m:s ',ihr(i),minu(i),sec(i)
!        if(i .eq.1) write(*,*) ' ID 1  ',id(i)
!        if(i .eq.1) write(*,*) ' By Bz ',bygse(i),bzgse(i)
!        if(i .eq.1) write(*,*) ' Vx n  ',vxgse(i),denin(i)

        if(i .ne. 1) then 
          if (id(i) .ne. id(i-1)) iday = iday + 1
        endif

c       Fill in timeline (decimal minutes)
        timein(i) = iday*24.0*60.0 + ihr(i)*60.0 + minu(i) + sec(i)/60.
c
      enddo
c

      write(*,*) ' Completed the LOOP '
      write(*,*) ' h:m:s ',ihr(i-1),minu(i-1),sec(i-1)
      write(*,*) ' By Bz ',bygse(i-1),bzgse(i-1)
      write(*,*) ' Vx n  ',vxgse(i-1),denin(i-1)
      write(*,*) ' ID of END ',id(i-1)
      npts = i-1
      write(*,*) ' Data points read: ',npts
      write(*,*) ' '

c     Set up output file
      time1 = timein(1) + window/2.0
      time2 = timein(npts) - window/2.0
      nout = int((time2 - time1)/cadence) + 1

c     loop over times for time, by, bz
      do j=1,nout

c       time, time1, twin1, twin2 all in minutes
c       hrutd in decimal hours UT relative to first day (monotonically increasing)
        time = time1 + (j-1)*cadence
        twin1 = time - window/2.0
        twin2 = time + window/2.0
        hrutd(j) = time/60.0

c       Average
        navg = 0
        bysum = 0.0
        bzsum = 0.0
        i = 0
        do 
          i = i+1
c         Just in case, keep track of last good By, Bz before time=time
          if (timein(i) .le. time .and. bygse(i) .lt. 9999.0 
     &      .and. bzgse(i) .lt. 9999.0) then
            bylg = bygse(i) 
            bzlg = bzgse(i) 
            ilg = i
          endif
          if (twin1 .le. timein(i) .and. timein(i) .le. twin2) then
            if (bygse(i) .lt. 9999.0 .and. bzgse(i) .lt. 9999.0) then
              navg = navg + 1
              bysum = bysum + bygse(i)
              bzsum = bzsum + bzgse(i)
            endif
          elseif (timein(i) .gt. twin2) then
            exit
          endif
        enddo    

c       By, Bz.  This craps out if the file begins or ends with a lot of bad data
c       (In those cases, we can modify the input file by adding a dummy "good" 
c       data point at the beginning or end of the file)
        if (navg .gt. 3) then
          by(j) = bysum/float(navg)
          bz(j) = bzsum/float(navg)
        else
c         If we have very sparse data, then simply interpolate.
c         Get first good data after time=time
          i = ilg
          do 
            i = i+1
            if (timein(i) .gt. time .and. bygse(i) .lt. 9999.0 
     &        .and. bzgse(i) .lt. 9999.0) then
              byfg = bygse(i) 
              bzfg = bzgse(i) 
              ifg  = i
              exit
            endif
          enddo
          fac = (time - timein(ilg))/(timein(ifg)- timein(ilg))
          by(j) = bylg - (bylg - byfg)*fac
          bz(j) = bzlg - (bzlg - bzfg)*fac
        endif
      enddo    

c     loop over times for density
      do j=1,nout

c       time, time1, twin1, twin2, timein all in minutes
        time = time1 + (j-1)*cadence
        twin1 = time - window/2.0
        twin2 = time + window/2.0

c       Average
        navg = 0
        dsum = 0.0
        i = 0
        do 
          i = i+1
c         Just in case, keep track of last good den before time=time
          if (timein(i) .le. time .and. denin(i) .lt. 999.0) then
            denlg = denin(i) 
            ilg = i
          endif
          if (twin1 .le. timein(i) .and. timein(i) .le. twin2) then
            if (denin(i) .lt. 999.0) then
              navg = navg + 1
              dsum = dsum + denin(i)
            endif
          elseif (timein(i) .gt. twin2) then
            exit
          endif
        enddo    

c       n.  This craps out if the file begins or ends with a lot of bad data
c       (In those cases, we can modify the input file by adding a dummy "good" 
c       data point at the beginning or end of the file)
        if (navg .gt. 3) then
          den(j) = dsum/float(navg)
        else
c         If we have very sparse data, then simply interpolate.
c         Get first good data after time=time
          i = ilg
          do 
            i = i+1
            if (timein(i) .gt. time .and. denin(i) .lt. 999.0) then
              denfg = denin(i) 
              ifg  = i
              exit
            endif
          enddo
          fac = (time - timein(ilg))/(timein(ifg)- timein(ilg))
          den(j) = denlg - (denlg - denfg)*fac
        endif
      enddo    

c     loop over times for Vx = -<Vx_gse>
      do j=1,nout

c       time, time1, twin1, twin2, timein all in minutes
        time = time1 + (j-1)*cadence
        twin1 = time - window/2.0
        twin2 = time + window/2.0

c       Average
        navg = 0
        vxsum = 0.0
        i = 0
        do 
          i = i+1
c         Just in case, keep track of last good den before time=time
          if (timein(i) .le. time .and. vxgse(i) .lt. 99999.0) then
            vxlg = -vxgse(i) 
            ilg = i
          endif
          if (twin1 .le. timein(i) .and. timein(i) .le. twin2) then
            if (vxgse(i) .lt. 99999.0) then
              navg = navg + 1
              vxsum = vxsum + vxgse(i)
            endif
          elseif (timein(i) .gt. twin2) then
            exit
          endif
        enddo    

c       Vx.  This craps out if the file begins with a lot of bad data
c        if (navg .gt. 0) then
c          vx(j) = -vxsum/float(navg)
c        else
c          vx(j) = vx(j-1)
c        endif

c       Vx.  This craps out if the file begins or ends with a lot of bad data
c       (In those cases, we can modify the input file by adding a dummy "good" 
c       data point at the beginning or end of the file)
        if (navg .gt. 3) then
          vx(j) = -vxsum/float(navg)
        else
c         If we have very sparse data, then simply interpolate.
c         Get first good data after time=time
          i = ilg
          do 
            i = i+1
            if (timein(i) .gt. time .and. vxgse(i) .lt. 99999.0) then
              vxfg = -vxgse(i) 
              ifg  = i
              exit
            endif
          enddo
          fac = (time - timein(ilg))/(timein(ifg)- timein(ilg))
          vx(j) = vxlg - (vxlg - vxfg)*fac
        endif
      enddo    

      vbigsum = 0.0
      do j=1,nout
        vbigsum = vbigsum + vx(j)
      enddo    
      vavg = vbigsum/nout

c     shift times based on average speed]
      tshift = (xcraft - 0.0)*earthr/(vavg*3600.0)
      hrutd = hrutd + tshift

      print *, ' Vavg (km/s), time shift (h) = ',vavg,tshift 

c     Dipole tilt angle

      do j=1,nout
        iday = int(hrutd(j)/24.0)
        ihour = int(hrutd(j) - iday*24.0)
        imin = int(hrutd(j)*60.0 - ihour*60.0 - iday*24.0*60.0)
        isec = 0.0
        iday_g = iday + iday0

c       recalc is Geopack routine , set the vgsex,... for computing dipole tilt angle
c       as we are using this, the value of vgsex doesn't matter so long as < 0
        vgsex = -400.
        vgsey = 0.
        vgsez = 0.
        iyr_i=iyear0
        if (iyr_i .gt. 2010) iyr_i = 2010
        call recalc_08(iyr_i,iday_g,ihour,imin,isec,vgsex,vgsey,vgsez)
        tilt(j) = aa(16)*180./3.14159

        if ( j .eq. 1 ) then
          print *,'stuff',iyr_i,iday_g,ihour,imin,isec,vgsex,
     .     vgsey,vgsez
        endif

      enddo

c     Write it all out
      close(33)
!      open(33,file=pwine,status='new',form='formatted')
      open(33,file=pwine,form='formatted')
!      write(33,*) ' Data comes from file(s): ',pwind       
!      write(33,*) ' Cadence, window (min): ',cadence,window
!      write(33,*) ' Npts, Year, DOY, <VSW>[km/s],',  
!     &            ' spacecraft-to-1AU time [h]'
!      write(33,*) ' hUT[decimal], By[nt], Bz[nt], -Vx[km/s],
!     & n[cm^-3], tilt[degrees]'

      write(33,44) nout,iyear0,iday0
   44 format(1x,3i8,2f15.5)
      do j=1,nout
        write(33,43) hrutd(j),by(j),bz(j),vx(j),den(j),tilt(j)
      enddo
   43 format(6f15.5)

c   43 format(11x,i2,1x,i2,1x,f6.3,3f14.1,3f14.2,f7.2)

c     End it all
      close(33)
      stop
      end


      subroutine julian_day(day,mon,year,jday)
    !============================    
    !filename: task.f95

      implicit none  

      integer day, jday, i, month(12), jd, year, mon  

      month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]    

      jday = 0  

      if (mod(year,4)==0) then
          month(2) = 29     
          do i = 1, mon-1
            jday = jday + month(i) 
          end do    
          jday = jday + day
      else    
          do i = 1, mon-1
            jday = jday + month(i) 
          end do    
          jday = jday + day
      end if  

      print *, jday   

      end 
