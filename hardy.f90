! hardy.f90

!      precipitation ionization rates  
 
!      Rees, Physics and chemistry of the upper atmosphere 
!            Cambridge Press, 1989 
!            pp. 39 - 41 
 
!      Rees, Auroral ionization and excitation by 
!            incident energetic electrons, 
!            Planet. Space Sci. 11, 1209, 1963 

!   particle flux and average electron energy from Hardy model 

!    subroutine hardy (nfl,nll,blats,alts,blons)

    subroutine hardy (hrut,nfl,nll)

      use hardy_mod
      use namelist_mod
      use grid_mod
      use parameter_mod

      integer :: nzh,nfl,nll,j,k,ii,ngrids,is,iit,nglam,ngridn,in
      real :: dele0,aven,value,part_flux
      real :: e0ps,rps,r0s
      real :: e0pn,rpn,r0n
      real :: thelat,themlt

!      real :: blats(nz,nf,nlt)
!      real ::  alts(nz,nf,nlt)
!      real :: blons(nz,nf,nlt)

      data iread_hardy / 1 /

      if ( iread_hardy == 1 ) then
         close(1001) 
         close(1002) 
         open(unit=1001,file='lambda_zor.inp') 
         open(unit=1002,file='atm_params_int.inp') 
         do ii = 1,no_lamzor 
           read(1001,*) zor(ii), lambda(ii) 
         enddo 
         do ii = 1,no_atm_ps 
           read(1002,*) height(ii), z(ii), rho(ii), nNn(ii), nM(ii) 
         enddo 
         close(1001) 
         close(1002) 
         iread_hardy = 0 
      endif

!      first initialize all ionization rates to 0 
 
       do i = 1,nz 
         do j = 1,nf 
           do k = 1,nl 
             preciprs(i,j,k) = 0. 
             preciprn(i,j,k) = 0. 
             iin(j,k)        = nz 
             iis(j,k)        = 1 
           enddo 
         enddo 
       enddo 

!      some parameters

!       kp    = 6

       dele0 = 35.e-3 
       nzh  = nz / 2 
       j    = nfl 
       k    = nll 
 
!      south pole 

       thelat = blats(1,j,k)
       themlt = mod(hrut + blons(nzh,nfl,nll) / 15.,24.)
       if ( lcr ) themlt = mod(blons(nzh,nfl,nll) / 15.,24.)

       call eavekp(kp,thelat,themlt,aven)
       call elekp(1,kp,thelat,themlt,value)
       part_flux = 10 ** value

       e0ps  = aven
       rps   = 4.57e-6 * e0ps ** 1.75 
       ii    = 1 
       do while ( z(ii) < rps ) 
         ii = ii + 1 
       enddo 
       ngrids = ii 
       r0s    = rps / rho(ngrids) 
       is     = 1 
!       print *,'is,nzh,alts(is,j,k)',is,nzh,alts(is,j,k)
       do while ( is < nzh .and. alts(is,j,k) <= 300. ) 
         iis(j,k) = is 
         if ( alts(is,j,k) < height(ngrids) ) then 
           preciprs(is,j,k) = 0. 
         else 
           iit      = ngrids  
           do while ( alts(is,j,k) > height(iit) ) 
             iit = iit - 1 
           enddo  
           if ( iit < ngrids .and. rps /= 0. ) then 
             nglam = ( z(iit) / rps ) * no_lamzor + 1 
             preciprs(is,j,k) = ( e0ps / r0s ) / dele0 * & 
                                lambda(nglam) *  part_flux * &
                                nM(iit)/nM(ngrids) 
           endif 
!           print *,'south',alts(is,j,k),preciprs(is,j,k)
         endif 
             is = is + 1 
       enddo 

!       print *,'max S',maxval(preciprs)

!      north pole 

       thelat = blats(nz,j,k)
       themlt = mod(hrut + blons(nzh,nfl,nll) / 15.,24.)
       if ( lcr ) themlt = mod(blons(nzh,nfl,nll) / 15.,24.)

       call eavekp(kp,thelat,themlt,aven)              
       call elekp(1,kp,thelat,themlt,value)
       part_flux = 10 ** value

!       print *,'aven,part_flux',aven,part_flux

       e0pn  = aven
       rpn   = 4.57e-6 * e0pn ** 1.75 
       ii    = 1 
       do while ( z(ii) < rpn ) 
         ii = ii + 1 
       enddo 
       ngridn = ii 
       r0n    = rpn / rho(ngridn) 
       in     = nz 
!       print *,'in,nzh,alts',in,nzh,alts(in,j,k)
       do while ( in > nzh .and. alts(in,j,k) <= 300. ) 
         iin(j,k) = in 
         if ( alts(in,j,k) < height(ngridn) ) then 
           preciprn(in,j,k) = 0. 
         else 
           iit      = ngridn 
           do while ( alts(in,j,k) > height(iit) ) 
             iit = iit - 1 
           enddo  
           if ( iit .lt. ngridn .and. rpn .ne. 0. ) then 
             nglam = ( z(iit) / rpn ) * no_lamzor + 1 
             preciprn(in,j,k) = ( e0pn / r0n ) / dele0 * & 
                                 lambda(nglam) * part_flux * &
                                 nM(iit)/nM(ngridn)
           endif 
!         print *,'north',alts(in,j,k),preciprn(in,j,k)
         endif 
         in       = in - 1 
       enddo 

!       print *,'max N',maxval(preciprn)

    end subroutine hardy


    SUBROUTINE EAVEKP(IKP,THELAT,THEMLT,AVEN)

!******************************************************************************

!   INPUTS:

!     IKP          INTEGRAL KP VALUE (0 THROUGH 6)
!     THELAT       CORRECTED GEOMAGNETIC LATITUDE
!     THEMLT       MAGNETIC LOCAL TIME

!   OUTPUTS:

!     AVEN         THE ELECTRON AVERAGE ENERGY FOR THIS KP MAP (KeV)

!   FILES :

!     EAVEKP.INP   ASCII FILE CONTAINING THE COEFFICIENTS FOR THE
!                  CHEBYSHEV-FOURIER EXPANSION

!   NOTES:

!     ROUTINE IS VALID ABOVE CORRECTED GEOMAGNETIC LATITUDE OF 50-DEG

!   REVISIONS:

!     16 SEP 98  ORIGINAL

!******************************************************************************

    use hardy_mod

    PARAMETER (LL=10)
    PARAMETER (MM=6)

    COMMON/SAVEKPE/A(7,LL),B(7,LL,MM),C(7,LL,MM)

    DIMENSION Y(20)

    DATA IFIRST/1/

    IF(IFIRST == 1)THEN
        IFIRST=0
        OPEN(10,FILE='eavekp.inp',FORM='FORMATTED', &
        STATUS='OLD',ERR=99)
        DO 10 I=1,7
            READ(10,*,ERR=98)JMOD,LL1,MM1
            IF(LL1 /= LL)THEN
                PRINT*,'NOT PROPER L VALUE....'
                STOP
            ENDIF
            IF(MM1 /= MM)THEN
                PRINT*,'NOT PROPER M VALUE....'
                STOP
            ENDIF
            READ(10,*,ERR=98)(A(JMOD,J),J=1,LL), &
            ((B(JMOD,J,K),J=1,LL),K=1,MM), &
            ((C(JMOD,J,K),J=1,LL),K=1,MM)
        10 END DO
        CLOSE(10)
    ENDIF

!   GET X VALUE FROM LATITUDE

    TTLAT=ABS(THELAT)
    IF(TTLAT < 50.0)THEN
        TTLAT=50.0
    ENDIF

!   MAX OUT THE LATITUDE AT THE 89-DEGREES

    IF(TTLAT > 89.0)THEN
        TTLAT=89.0
    ENDIF

    TWOTHIRDS=2.D0/3.D0
    XX=2.0*(TTLAT-50.0)/40.0-1.0
    IF(XX < 0.0)THEN
        SIG=-1.0
    ELSE
        SIG=1.0
    ENDIF
    XX=SIG*ABS(XX)**TWOTHIRDS

!   GET CHEVY POLYNOMIALS

    CALL GETCHEVY3(XX,LL,Y)

!   CHECK FOR VALID KP

    IF(IKP < 0) then 
      PRINT *,'KP MUST BE >= 0'
      stop
    endif

    IF(IKP > 6) then
      PRINT *,'KP MUST BE <= 6'
      stop
    endif

!   MODEL NUMBER IS IKP+1

    IMOD=IKP+1

    AV=0.0

!   GET COSINE TERMS

    TODEG=45.0/ATAN(1.0)
    DEGR=THEMLT*360.0/24.0/TODEG

    DO 20 I=1,LL
        AV=AV+Y(I)*A(IMOD,I)
        DO 20 J=1,MM
            FMM=FLOAT(J)
            AV=AV+Y(I)*B(IMOD,I,J)*SIN(FMM*DEGR)
            AV=AV+Y(I)*C(IMOD,I,J)*COS(FMM*DEGR)
    20 END DO

    AVEN=EXP(AV)

    RETURN

    99 PRINT*,'THE FILE iavekp.inp IS NOT HERE'
    STOP
    98 PRINT*,'THE FILE iavekp.inp IS NO GOOD'
    STOP

    END SUBROUTINE EAVEKP

    SUBROUTINE GETCHEVY3(X,IORD,CHEVY)

!   THIS ROUTINE EVALUATES THE DIRECT EXPANSION OF THE CHEVYS

    DIMENSION CHEVY(20)

    IF(IORD > 20) then
       PRINT*,'NOT MORE THAN 20 TERMS...'
       stop
    endif

    CHEVY(1)=1.
    CHEVY(2)=X
    IF(IORD <= 2) return
    DO 11 J=3,IORD
        CHEVY(J)=2.*X*CHEVY(J-1)-CHEVY(J-2)
    11 END DO

    end subroutine getchevy3


    SUBROUTINE ELEKP(IOPT,KP,THELAT,THEMLT,VALUE)

!   INPUTS:      IOPT    1=ELECTRON NUMBER FLUX
!                        2=ELECTRON ENERGY FLUX

!                KP     THE Kp MAP NUMBER
!                        0 0,0+
!                        1 1-,1,1+
!                        ...
!                        6 => 6-

!                THELAT  CORRECTED GEOMAGNETIC LATITUDE (DEG)
!                THEMLT  MAGNETIC LOCAL TIME (0->24)

!   OUTPUT:      VALUE   LOG10 OF THE FLUX (1/CM^2/SEC ?)

!   NOTES:       THE POLEWARD VALUES OF THE FLUXES ARE SET TO
!                A MINIMUM OF THE VALUES CONTAINED IN PRAIN()
!                WHICH ARE DERIVED FROM THE DATA.  EQUATORWARD
!                VALUES ARE SET AT UNITY.

!                QUESTIONS CONCERNING THE OPERATION OF THIS
!                MODEL SHOULD BE DIRECTED TO:

!                      BILL MCNEIL
!                      RADEX, INC.
!                      (617)275-6767
!                      MCNEIL@PLH.AF.MIL

!   REVISIONS:   11 SEP 98   ORIGINAL

    use hardy_mod

    COMMON/SAVKPE/CR0(2,7,17),CH0(2,7,17),CH1(2,7,17), &
    CS0(2,7,17),CR1(2,7,17),CS2(2,7,17), &
    PRAIN(2,7)

    CHARACTER(80) :: LINE

    DATA IFIRST/0/

    DATA PRAIN/7.07292,7.01250,7.02709, &
    &            7.18542,7.00417,6.90417,6.98750, &
    &            6.48542,6.68333,6.62083, &
    &            6.82708,6.54583,6.49167,6.46458/

!   READ THE DATA ON FIRST PASS

    IF(IFIRST == 0)THEN
        IFIRST=1
        OPEN(101,FILE='elekp.inp',FORM='FORMATTED', &
        STATUS='OLD',ERR=99)
        DO 10 JOPT=1,2
            DO 10 JKP=1,7
                READ(101,100,ERR=98)LINE
                100 FORMAT(A80)
                DO 10 JCOEF=1,17
                    READ(101,*,ERR=98)JNUM,CR0(JOPT,JKP,JCOEF), &
                    CH0(JOPT,JKP,JCOEF), &
                    CR1(JOPT,JKP,JCOEF),CH1(JOPT,JKP,JCOEF), &
                    CS0(JOPT,JKP,JCOEF),CS2(JOPT,JKP,JCOEF)
        10 END DO
        CLOSE(101)
    ENDIF

    IKP=KP+1

    IF(KP < 0) then 
      PRINT *,'KP MUST BE >= 0'
      stop
    endif

    IF(KP > 6) then
      PRINT *,'KP MUST BE <= 6'
      stop
    endif

    IKP=KP+1

    PI=3.14159265
    R0=CR0(IOPT,IKP,1)/2.0
    H0=CH0(IOPT,IKP,1)/2.0
    H1=CH1(IOPT,IKP,1)/2.0
    S0=CS0(IOPT,IKP,1)/2.0
    R1=CR1(IOPT,IKP,1)/2.0
    S2=CS2(IOPT,IKP,1)/2.0

    ARG=PI*THEMLT/12.0
    K=2
    DO 12 I=1,8
        AARG=ARG*I
        C=COS(AARG)
        S=SIN(AARG)
        R0=R0+C*CR0(IOPT,IKP,K)
        H0=H0+C*CH0(IOPT,IKP,K)
        H1=H1+C*CH1(IOPT,IKP,K)
        S0=S0+C*CS0(IOPT,IKP,K)
        R1=R1+C*CR1(IOPT,IKP,K)
        S2=S2+C*CS2(IOPT,IKP,K)
        R0=R0+S*CR0(IOPT,IKP,K+8)
        H0=H0+S*CH0(IOPT,IKP,K+8)
        H1=H1+S*CH1(IOPT,IKP,K+8)
        S0=S0+S*CS0(IOPT,IKP,K+8)
        R1=R1+S*CR1(IOPT,IKP,K+8)
        S2=S2+S*CS2(IOPT,IKP,K+8)
        K=K+1
    12 END DO

    B1=ALOG((1.0+EXP(H1-H0))/2.0)
    B2=ALOG(2.0/(1.+EXP(H0-H1)))
    S1=(R1-R0-S0*(H1-H0)+S0*B1-S2*B2)/(B1-B2)

    H=ABS(THELAT)
    IF(H < 50.0)H=50.0

    EH=R0+S0*(H-H0)
    EH=EH+(S1-S0)*ALOG((1.0+EXP(H-H0))/2.0)
    EH=EH+(S2-S1)*ALOG((1.0+EXP(H-H1))/ &
    (1.0+EXP(H0-H1)))

    IF(EH < 0.0)EH=0.0

    IF(H > H1 .AND. EH < PRAIN(IOPT,IKP))THEN
        EH=PRAIN(IOPT,IKP)
    ENDIF

    VALUE=EH

    RETURN

    99 PRINT*,'THE FILE elekp.inp MUST BE PRESENT'
    STOP
    98 PRINT*,'THE FILE elekp.inp IS CORRUPTED'
    STOP

    END SUBROUTINE ELEKP

