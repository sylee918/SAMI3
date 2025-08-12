!******************************************
!******************************************

!            photprod

!******************************************
!******************************************

! photoproduction rates

    subroutine photprod ( phprodr,nfl,nll,hrut)

    use parameter_mod
    use variable_mod
    use namelist_mod
    use photo_production_mod
    use atomic_mod
    use grid_mod
    use misc_mod
    use time_mod
    use message_passing_mod

    real :: phprodr(nz,nion),xmass(3),dtmask,tmod60
    integer :: idx(3),itmod60
    real(8) :: ch1
   ! real(8) :: atm_chapman
    real    :: xscale
    
    real(8) :: xscale8
    real(8) :: arg8
    
    logical :: event
    integer :: istep
    
    real  :: dflx(linesuv)
    real  :: nflx(linesnt)
    real  :: dphprodr
    real  :: nphprodr

  ! scale height of neutral atmosphere

    hcof = 1.e-5 * bolt / ( gzero * amu * re ** 2 )

    nll1 = nll

    phprodr = 0.0 ! phprodr(:,:) = 0.0
    do iz = 1,nz
        
        coschi = cx(iz,nfl,nll)
        
        ! only consider o, n2, o2 for absorption

        idx(1) = pto
        idx(2) = ptn2
        idx(3) = pto2

        rp    = alts(iz,nfl,nll1) + re
        rp2   = rp * rp
                 
        if ( coschi >= coschicrit(iz,nfl,nll1) ) then ! sun is up

        ! daytime deposition

            do i = 1,3
                hscale   = hcof * tn(iz,nfl,nll1) * rp2 / amn(idx(i))                
                xscale8   = rp / hscale
                arg8 = rtod*acos(coschi)
                ch1 = chapman(xscale8,arg8)                
                xmass(i) = denn(iz,nfl,nll1,idx(i)) * hscale * ch1 * 1.e5
            enddo

             do l = 1,linesuv
               exa =  xmass(1) * sigabsdt(l,1) + xmass(2) * sigabsdt(l,2) + xmass(3) * sigabsdt(l,3)
               exa = min(exa,35.0)
               dflx(l) = flux(l) * exp(-exa)
            enddo
            
            do j = nion1,nion2
                dphprodr = 0.0
                do l = 1,linesuv
                    dphprodr = dphprodr +  sigidt(l,j) * dflx(l)
                enddo
                phprodr(iz,j) = phprodr(iz,j) + dphprodr
            enddo
            
          !  do l=1,linesuv
          !      exa =  xmass(1) * sigabsdt(l,1) &
          !          + xmass(2) * sigabsdt(l,2) &
          !          + xmass(3) * sigabsdt(l,3)
          !      if(exa > 35.) exa = 35.
          !      ! add eclipse mask
          !      flx = emask(iz) * flux(l) * exp(-exa)
          !      u4(iz,nfl,nll)  = emask(iz)
          !      do j = nion1,nion2
          !          phprodr(iz,j) = phprodr(iz,j) + sigidt(l,j) * flx
          !      enddo
          !  enddo
            
            ! photoelectron ionization ????

            pei_rate = exp(-(alts(iz,nfl,nll) - 150.0)/40.0)
            if ( alts(iz,nfl,nll) .gt. 200.0 ) pei_rate = 0.2

            phprodr(iz,ptop)  = phprodr(iz,ptop)  * (1.0 + pei_rate)
            phprodr(iz,ptn2p) = phprodr(iz,ptn2p) * (1.0 + pei_rate)
            phprodr(iz,pto2p) = phprodr(iz,pto2p) * 1.2

        endif

    enddo
    
! ----------------------------------------------------------------------------
! High-energy ionization contributes both day and night.
! ----------------------------------------------------------------------------

    do iz = 1,nz
    
        coschi = cx(iz,nfl,nll)
        ang    = acos ( coschi )
        itheta0 = int ( ang / po180 ) - 90
        itheta  = int ( amax1 ( float(itheta0), 1. ) )
        ithetap1 = itheta + 1
        del     = ang/po180 - int(ang/po180)
        
        if (itheta0 .lt. 1) then
             do l = 1,linesnt
                nflx(l) = fluxnt(l,itheta,iz,nfl,nll)
             enddo
        else
            do l = 1,linesnt
                nflx(l) = fluxnt(l,itheta,iz,nfl,nll) * (1.0 - del) + &
                          fluxnt(l,ithetap1,iz,nfl,nll) * del
            enddo
        endif
        
        do j = nion1,nion2
            nphprodr = 0.0
            do l = 1,linesnt
                nphprodr = nphprodr + sigint(l,j) * nflx(l)
            enddo
            phprodr(iz,j) = phprodr(iz,j) + nphprodr
        enddo
        
        !ang    = acos ( coschi )
        !itheta0 = int ( ang / po180 ) - 90
        !itheta  = int ( amax1 ( float(itheta0), 1. ) )
        !del     = ang/po180 - int(ang/po180)
        !
        !do l = 1,linesnt
        !    do j=nion1,nion2
        !        if (itheta0 < 1) then
        !            fluxntt = fluxnt(iz,nfl,nll1,itheta,l)
        !        else
        !            fluxntt = fluxnt(iz,nfl,nll1,itheta,l) * (1.-del) + &
        !            fluxnt(iz,nfl,nll1,itheta+1,l) * del
        !        endif
        !        phprodr(iz,j) =   phprodr(iz,j) &
        !        + sigint(l,j) * fluxntt
        !                              
        !    enddo
        !enddo

        ! u4(iz,nfl,nll) = phprodr(iz,ptop)

    enddo

    return
    end subroutine photprod


!******************************************
!******************************************
!            f1026
!******************************************
!******************************************

! subroutine to calculate the nighttime flux of
! lyman beta (1026) (note: line = 1)

    subroutine sf1026 ( f,line,nfl,nll )

    use parameter_mod
    use variable_mod
    use photo_production_mod
    use grid_mod


    real :: f(nz,nf,nl,91)

    imax = 1

! determine f for the 4 known values of theta

    do i = 1,nz
        if ( alts(i,nfl,nll) < zaltnt(line,1) ) then
            do k = 1,4
                f( i,nfl,nll,int(thetant(line,k))+1-90 ) = 1.
            enddo
        elseif ( zaltnt(line,1) <= alts(i,nfl,nll) .AND. &
            alts(i,nfl,nll) <= zaltnt(line,2)       ) then
            f( i,nfl,nll,int(thetant(line,1))+1-90 ) = &
            &        1.4e8 * tanh ( (alts(i,nfl,nll) - 90.) / 50. )
            f( i,nfl,nll,int(thetant(line,2))+1-90 ) = &
            &        3.8e7 * tanh ( (alts(i,nfl,nll) - 90.) / 50. )
            f( i,nfl,nll,int(thetant(line,3))+1-90 ) = &
            &        1.4e7 * tanh ( (alts(i,nfl,nll) - 93.) / 55. )
            f( i,nfl,nll,int(thetant(line,4))+1-90 ) = &
            &        9.2e6 * tanh ( (alts(i,nfl,nll) - 94.) / 55. )
            imax = i
        else
            do k = 1,4
                f( i,nfl,nll,   int(thetant(line,k))+1-90 ) = &
                f( imax,nfl,nll,int(thetant(line,k))+1-90 )
            enddo
        endif
    enddo

    do k = 1,4
        do i = 1,nz
            f( i,nfl,nll,int(thetant(line,k))+1-90 ) = &
            amax1 ( 1., f( i,nfl,nll,int(thetant(line,k))+1-90 ) )
        enddo
    enddo

! now interpolate to all valuse of theta (90 - 180)
     
    do k = 1,91
        k90 = 90 + k - 1
        ji  = 1
        ki  = int(thetant(line,1))
        do j = 1,4
            if ( k90 > int(thetant(line,j)) ) then
                ji = j
                ki = int(thetant(line,ji))
            endif
        enddo
        jip1 = ji + 1
        kip1 = int(thetant(line,jip1))
        delk = float (   int(thetant(line,jip1)) &
        - int(thetant(line,ji  )) )
        do i = 1,nz
            flog =   alog10(f(i,nfl,nll,ki+1-90)) &
            + (k90 - ki) / delk &
            * (  alog10(f(i,nfl,nll,kip1+1-90)) &
            - alog10(f(i,nfl,nll,ki  +1-90)) )
            f(i,nfl,nll,k) = 10 ** flog
        enddo
    enddo

    return
    end subroutine sf1026

!******************************************
!******************************************
!            f584
!******************************************
!******************************************

! subroutine to calculate the nighttime flux of
! he i (584) (note: line = 2)

    subroutine sf584 ( f,line,nfl,nll )

    use parameter_mod
    use variable_mod
    use photo_production_mod
    use grid_mod


    real :: f(nz,nf,nl,91)

    imax = 1

! determine f for the 4 known values of theta

    do i = 1,nz
        if ( alts(i,nfl,nll) < zaltnt(line,1) ) then
            do k = 1,4
                f( i,nfl,nll,int(thetant(line,k))+1-90 ) = 1.
            enddo
        elseif ( zaltnt(line,1) <= alts(i,nfl,nll) .AND. &
            alts(i,nfl,nll) <= zaltnt(line,2)       ) then
            f( i,nfl,nll,int(thetant(line,1))+1-90 ) = &
            &        1.85e5 * ( alts(i,nfl,nll) - 170. ) ** 1.20
            f( i,nfl,nll,int(thetant(line,2))+1-90 ) = &
            &        2.60e4 * ( alts(i,nfl,nll) - 170. ) ** 1.25
            f( i,nfl,nll,int(thetant(line,3))+1-90 ) = &
            &        2.60e3 * ( alts(i,nfl,nll) - 170. ) ** 1.20
            f( i,nfl,nll,int(thetant(line,4))+1-90 ) = &
            &        2.60e2 * ( alts(i,nfl,nll) - 170. ) ** 1.20
            imax = i
        else
            do k = 1,4
                f( i   ,nfl,nll,int(thetant(line,k))+1-90 ) = &
                f( imax,nfl,nll,int(thetant(line,k))+1-90 )
            enddo
        endif
    enddo

    do k = 1,4
        do i = 1,nz
            f( i,nfl,nll,int(thetant(line,k))+1-90 ) = &
            amax1 ( 1., f( i,nfl,nll,int(thetant(line,k))+1-90 ) )
        enddo
    enddo

! now interpolate to all valuse of theta (90 - 180)
! set f(i,nfl,nll,theta=180) = 1.

    do k = 1,91
        k90 = 90 + k - 1
        ji  = 1
        ki  = int(thetant(line,1))
        do j = 1,4
            if ( k90 > int(thetant(line,j)) ) then
                ji = j
                ki = int(thetant(line,ji))
            endif
        enddo
        if ( ji /= 4 ) then
            jip1 = ji + 1
            kip1 = int(thetant(line,jip1))
            delk = float (   int(thetant(line,jip1)) &
            - int(thetant(line,ji  )) )
            do i = 1,nz
                flog =   alog10(f(i,nfl,nll,ki+1-90)) &
                + (k90 - ki) / delk &
                * (  alog10(f(i,nfl,nll,kip1+1-90)) &
                - alog10(f(i,nfl,nll,ki  +1-90)) )
                f(i,nfl,nll,k) = 10 ** flog
            enddo
        else
            delk = float (   180 &
            - int(thetant(line,ji  )) )
            do i = 1,nz
                flog =   alog10(f(i,nfl,nll,ki+1-90)) &
                + (k90 - ki) / delk &
                * (  alog10(1.) &
                - alog10(f(i,nfl,nll,ki  +1-90)) )
                f(i,nfl,nll,k) = 10 ** flog
            enddo
        endif
    enddo

    return
    end subroutine sf584

!******************************************
!******************************************

!            f304

!******************************************
!******************************************

! subroutine to calculate the nighttime flux of
! he ii (304) (note: line = 3)

    subroutine sf304 ( f,line,nfl,nll )

    use parameter_mod
    use variable_mod
    use photo_production_mod
    use grid_mod


    real :: f(nz,nf,nl,91)

    imax = 1

! determine f for the 4 known values of theta

    do i = 1,nz
        if ( alts(i,nfl,nll) < zaltnt(line,1) ) then
            do k = 1,4
                f( i,nfl,nll,int(thetant(line,k))+1-90 ) = 1.
            enddo
        elseif ( zaltnt(line,1) <= alts(i,nfl,nll) .AND. &
            alts(i,nfl,nll) <= zaltnt(line,2)       ) then
            f( i,nfl,nll,int(thetant(line,1))+1-90 ) = &
            &        3.8e6 * tanh ( (alts(i,nfl,nll) - 138.) / 80. )
            f( i,nfl,nll,int(thetant(line,2))+1-90 ) = &
            &        3.0e6 * tanh ( (alts(i,nfl,nll) - 138.) / 80. )
            f( i,nfl,nll,int(thetant(line,3))+1-90 ) = &
            &        2.5e6 * tanh ( (alts(i,nfl,nll) - 138.) / 80. )
            f( i,nfl,nll,int(thetant(line,4))+1-90 ) = &
            &        2.5e6 * tanh ( (alts(i,nfl,nll) - 138.) / 80. )
            imax = i
        else
            do k = 1,4
                f( i,   nfl,nll,int(thetant(line,k))+1-90 ) = &
                f( imax,nfl,nll,int(thetant(line,k))+1-90 )
            enddo
        endif
    enddo

    do k = 1,4
        do i = 1,nz
            f( i,nfl,nll,int(thetant(line,k))+1-90 ) = &
            amax1 ( 1., f( i,nfl,nll,int(thetant(line,k))+1-90 ) )
        enddo
    enddo

! now interpolate to all valuse of theta (90 - 180)
! set f(i,nfl,nll,theta=180) = 1.

    do k = 1,91
        k90 = 90 + k - 1
        ji  = 1
        ki  = int(thetant(line,1))
        do j = 1,4
            if ( k90 > int(thetant(line,j)) ) then
                ji = j
                ki = int(thetant(line,ji))
            endif
        enddo
        if ( ji /= 4 ) then
            jip1 = ji + 1
            kip1 = int(thetant(line,jip1))
            delk = float (   int(thetant(line,jip1)) &
            - int(thetant(line,ji  )) )
            do i = 1,nz
                flog =   alog10(f(i,nfl,nll,ki+1-90)) &
                + (k90 - ki) / delk &
                * (  alog10(f(i,nfl,nll,kip1+1-90)) &
                - alog10(f(i,nfl,nll,ki  +1-90)) )
                f(i,nfl,nll,k) = 10 ** flog
            enddo
        else
            delk = float (   180 &
            - int(thetant(line,ji  )) )
            do i = 1,nz
                flog =   alog10(f(i,nfl,nll,ki+1-90)) &
                + (k90 - ki) / delk &
                * (  alog10(1.) &
                - alog10(f(i,nfl,nll,ki  +1-90)) )
                f(i,nfl,nll,k) = 10 ** flog
            enddo
        endif
    enddo

    return
    end subroutine sf304

!******************************************
!******************************************

!            f1216

!******************************************
!******************************************

!     subroutine to calculate the nighttime flux of
!     lyman alpha (1216) (note: line = 4)

    subroutine sf1216 ( f,line,nfl,nll )

    use parameter_mod
    use variable_mod
    use photo_production_mod
    use grid_mod


    real :: f(nz,nf,nl,91)

    imax = 1

! determine f for the 4 known values of theta

    do i = 1,nz
        if ( alts(i,nfl,nll) < zaltnt(line,1) ) then
            do k = 1,4
                f( i,nfl,nll,int(thetant(line,k))+1-90 ) = 1.
            enddo
        elseif ( zaltnt(line,1) <= alts(i,nfl,nll) .AND. &
            alts(i,nfl,nll) <= zaltnt(line,2)       ) then
            f( i,nfl,nll,int(thetant(line,1))+1-90 ) = &
            &        1.2e10 * tanh ( (alts(i,nfl,nll) - 80.) / 50. ) + 3.e9
            f( i,nfl,nll,int(thetant(line,2))+1-90 ) = &
            &        4.0e9  * tanh ( (alts(i,nfl,nll) - 80.) / 50. ) + 1.e9
            f( i,nfl,nll,int(thetant(line,3))+1-90 ) = &
            &        2.0e9  * tanh ( (alts(i,nfl,nll) - 65.) / 50. ) + 1.e8
            f( i,nfl,nll,int(thetant(line,4))+1-90 ) = &
            &        1.5e9  * tanh ( (alts(i,nfl,nll) - 75.) / 50. ) + 1.e8
            imax = i
        else
            do k = 1,4
                f( i,   nfl,nll,int(thetant(line,k))+1-90 ) = &
                f( imax,nfl,nll,int(thetant(line,k))+1-90 )
            enddo
        endif
    enddo

    do k = 1,4
        do i = 1,nz
            f( i,nfl,nll,int(thetant(line,k))+1-90 ) = &
            amax1 ( 1., f( i,nfl,nll,int(thetant(line,k))+1-90 ) )
        enddo
    enddo

! now interpolate to all valuse of theta (90 - 180)
     
    do k = 1,91
        k90 = 90 + k - 1
        ji  = 1
        ki  = int(thetant(line,1))
        do j = 1,4
            if ( k90 > int(thetant(line,j)) ) then
                ji = j
                ki = int(thetant(line,ji))
            endif
        enddo
        jip1 = ji + 1
        kip1 = int(thetant(line,jip1))
        delk = float (   int(thetant(line,jip1)) &
        - int(thetant(line,ji  )) )
        do i = 1,nz
            flog =   alog10(f(i,nfl,nll,ki+1-90)) &
            + (k90 - ki) / delk &
            * (  alog10(f(i,nfl,nll,kip1+1-90)) &
            - alog10(f(i,nfl,nll,ki  +1-90)) )
            f(i,nfl,nll,k) = 10 ** flog
        enddo
    enddo

    return
    end subroutine sf1216

!******************************************
!******************************************

!            zenith

!******************************************
!******************************************

    subroutine zenith (hrut,nfl,nll)

    use parameter_mod
    use variable_mod
    use namelist_mod
    use photo_production_mod
    use grid_mod

! geometric variables


! bdec: magnetic declination angle
! sdec: solar zenith angle
! cx:  cos of the zenith angle

    do i = 1,nz
        hrl   = mod(hrut + glons(i,nfl,nll) / 15.,24.)
        sdec          = rtod * asin (  sin (2.*pie*(day-dayve)/sidyr) * sin (solinc/rtod)             )
        cossdec       = cos ( po180 * sdec )
        sinsdec       = sin ( po180 * sdec )
        clat          = cos ( po180 * glats(i,nfl,nll) )
        slat          = sin ( po180 * glats(i,nfl,nll) )
        cx(i,nfl,nll) =   slat * sinsdec - clat * cossdec * cos ( 15.0*po180*hrl )
        
    ! MS: Since we will be taking acos of this value in photprod, make
    ! sure that the absolute value does not minutely exceed 1 because of
    ! round-off error.
    !    if (abs(abs(cx(i,nfl,nll))-1.) < 1.e-6) &
    !    cx(i,nfl,nll) = sign(1.,cx(i,nfl,nll))
    
        enddo

    return
    end subroutine zenith


