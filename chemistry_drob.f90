!******************************************
!******************************************
!            chemrate
!******************************************
!******************************************

!     chemical producation and loss rates
!     bb: bailley and balan (red book, 1996)
!     etc....

    subroutine chemrate( chrate,nfl,nll )

    use parameter_mod
    use variable_mod

    integer,intent(in)  :: nfl,nll

    real,intent(out)    :: chrate(nz,nchem)
    real                :: ti300o
    real                :: x
    real                :: teff

    do iz = 1,nz

        ti300o = ti(iz,nfl,nll,ptop) / 300.

    ! h+ + o --> o+ + h (bb)

        ! Banks and Kockarts

        !   chrate (iz,1) = 2.2e-11 * sqrt( ti(iz,nfl,nll,pthp) )

        ! Stancil (1999)

        teff = (ti(iz,nfl,nll,pthp)*16.0 + tn(iz,nfl,nll)*1.0)/17.0

        x = teff/10000.0

        chrate(iz,1) = 1.26e-9 * x ** 0.517 * exp( -227.0 / teff) + 4.25e-10 * x ** 0.00669

    ! he+ + n2 --> n2+ + he (bb)

        chrate (iz,2) = 3.6e-10

    ! he+ + n2 --> n+ + n + he (schunk)

        chrate (iz,3) = 8.4e-10

    ! he+ + o2 --> o+ + o + he

        chrate (iz,4) = 6.4e-12

    ! he+ + o2 --> o2+ + he

        chrate (iz,5) = 1.0e-9

    ! n+ + o2 --> no+ + o  (schunk)

        chrate (iz,6) = 2.0e-10

    ! n+ + o2 --> o2+ + n(2d) (schunk)

        chrate (iz,7) = 4.0e-10

    ! n+ + O --> o+ + n

        chrate (iz,8) = 2.2e-12

    ! n+ + no --> no+ + o (schunk)

        teff = (ti(iz,nfl,nll,ptnp) * 30.0 + tn(iz,nfl,nll)* 14.0)/44.0
        chrate (iz,9) = 5.72e-9 / teff ** 0.44

    ! o+ + h --> h+ + o   (bb)

        ! Banks and Kockarts 

        ! chrate(iz,10) = 2.5e-11 * sqrt( tn(iz,nfl,nll) )

        ! Stancil (1999)

        teff = (ti(iz,nfl,nll,ptop)*1.0 + tn(iz,nfl,nll)*16.0)/17.0

        x = teff/10000.0

        chrate(iz,10) = 2.08e-9 * x ** 0.405 + 1.11e-11 / x ** 0.458

    ! o+ + n2 --> no+ + n (bb)

        teff = (ti(iz,nfl,nll,ptop)*28.0 + tn(iz,nfl,nll)*16.0)/44.0 
        ti300o = teff / 300.

        chrate(iz,11) = 1.533e-12 - &
        &               5.920e-13 * ti300o + &
        &               8.600e-14 * ti300o ** 2

        if ( teff > 1700 ) &
            chrate(iz,11) = 2.730e-12 - &
        &               1.155e-12 * ti300o + &
        &               1.483e-13 * ti300o ** 2

    ! o+ + o2 --> o2+ + o

        teff = (ti(iz,nfl,nll,ptop)*32.0 + tn(iz,nfl,nll)*16.0)/48.0 
        ti300o = teff / 300.

        chrate(iz,12) = 2.820e-11 - &
        &               7.740e-12 * ti300o + &
        &               1.073e-12 * ti300o ** 2 - &
        &               5.170e-14 * ti300o ** 3 + &
        &               9.650e-16 * ti300o ** 4

    ! o+ + no --> no+ + o

        teff = (ti(iz,nfl,nll,ptnop)*30.0 + tn(iz,nfl,nll)*16.0)/46.0 
        ti300o = teff / 300.

        chrate(iz,13) = 8.86e-13 - 2.02e-13 * ti300o + 6.95e-14 * ti300o ** 2

    ! n2+ + o --> no+ + n(2d) (bb)

        teff = (ti(iz,nfl,nll,ptn2p)*16.0 + tn(iz,nfl,nll)*28.0)/44.0 
        ti300o = teff / 300.

        chrate(iz,14) = 1.4e-10 / ti300o ** 0.44

        if ( teff > 1500.0 ) &

            chrate(iz,14) = 5.2e-11 * ti300o ** 0.2

   ! n2+ + o2 --> o2+ + n2 

        teff = (ti(iz,nfl,nll,ptn2p)*32.0 + tn(iz,nfl,nll)*28.0)/60.0 
        x = teff/300.0

        if (teff .lt. 1000.0) then
            chrate(iz,15) = 1.0e-10 * (0.1142 - 0.1171 * x + 5.047e-2 * x ** 2 - 9.551e-3 * x ** 3 + 6.344e-4 * x ** 4)
        else
            chrate(iz,15) = 1.0e-11 * (teff/1000.0) ** 0.912
        endif

 

    ! n2+ + o2 --> no+ + no

        chrate(iz,16) = 1.0e-14

    ! n2+ + no --> no+ + n2

        teff = (ti(iz,nfl,nll,ptn2p)*30.0 + tn(iz,nfl,nll)*28.0)/58.0 
        chrate(iz,17) = 7.5e-9 / teff ** 0.52

    ! o2+ + n --> no+ + o

        chrate(iz,18) = 1.0e-10

    ! o2+ + n(2d) --> n+ + o2

        chrate(iz,19) = 2.5e-10

    ! o2+ + no --> no+ + o2 (bb)

        chrate(iz,20) = 4.1e-10

    ! o2+ + n2 --> no+ + no (schunk)

        chrate(iz,21) = 5.0e-16

    enddo

    return

    end subroutine chemrate

!******************************************

!******************************************

!            recorate

!******************************************

!******************************************

! recombination rates
! bb: bailley and balan (red book, 1996)

    subroutine recorate ( relossr,nfl,nll )

    use parameter_mod
    use variable_mod
    use namelist_mod

    real :: relossr(nz,nion)

    do iz = 1,nz

        te300 = te(iz,nfl,nll) / 300.0 
        xe = 1.0 / te300

        relossr(iz,pthp) =   4.43e-12 / te300 ** 0.7
        relossr(iz,pthep) =  4.43e-12 / te300 ** 0.7 ! same as factor for relossr(iz,pthp)
        relossr(iz,ptnp)  =  4.43e-12 / te300 ** 0.7 ! same as factor for relossr(iz,pthp)
        relossr(iz,ptop)  =  4.43e-12 / te300 ** 0.7 ! same as factor for relossr(iz,pthp)

        if (te(iz,nfl,nll) .lt. 1200.0) then
            relossr(iz,ptn2p) =  2.2e-7 * xe ** 0.39       !  
            relossr(iz,ptnop) =  3.5e-7 * xe ** 0.69      !  
            relossr(iz,pto2p) =  1.95e-7 * xe ** 0.70     !  
        else
            relossr(iz,ptn2p) = 1.95e-7 * xe ** 0.57      !
            relossr(iz,ptnop) = 3.02e-7 * xe ** 0.56     !
            relossr(iz,pto2p) = 1.93e-7 * xe ** 0.61     !
        endif

!        relossr(iz,ptmgp) = relossfact(ptmgp)%now * 2.81e-12 * xe ** 0.855

 

    enddo

 

    return

    end subroutine recorate



!******************************************
!******************************************

!          chempl

!******************************************
!******************************************
            
! chemical loss (chloss) and production (chprod)

! chrate: chemical reaction rates calculated in chemrate
! ichem: input data file showing loss, neutral, production
!        species for each reaction

    subroutine chempl ( chrate,chloss,chprod,nfl,nll )

    use parameter_mod
    use variable_mod
    use namelist_mod
    use chemistry_mod

    real :: chrate(nz,nchem),chloss(nz,nion),chprod(nz,nion)

    do i = nion1,nion2
        do iz = 1,nz
            chloss(iz,i)   = 0.
            chprod(iz,i)   = 0.
        enddo
    enddo

    do k = 1,nchem
        il = ichem(k,1) ! ion species (reacting) loss
        in = ichem(k,2) ! neutral species reacting
        ip = ichem(k,3) ! ion species produced
        do iz = 1,nz
            chem  = denn(iz,nfl,nll,in) * chrate(iz,k)
            tdeni = deni(iz,nfl,nll,il) * chem
            chloss(iz,il) = tdeni + chloss(iz,il)
            chprod(iz,ip) = tdeni + chprod(iz,ip)
        enddo
    enddo

    return
    end subroutine chempl

