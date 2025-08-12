!******************************************
!******************************************

!       COM3_MPI-2.00.INC

!******************************************
!******************************************

module com3

    use param3

    real dt
!     s grid data

! Worker grid
    real :: alts(nz,nf,nl),grs(nz,nf,nl),glats(nz,nf,nl),glons(nz,nf,nl)
    real :: glons0(nz,nf,nl)
    real :: bms(nz,nf,nl),gs(nz,nl),ps(nz,nf,nl),gp(nz,nf,nl), &
    blats(nz,nf,nl),blons(nz,nf,nl),cfs(nz,nl)
    real :: gsthetax(nz,nf,nl),gsthetay(nz,nf,nl),gsthetaz(nz,nf,nl)
    real :: gsphix(nz,nf,nl),gsphiy(nz,nf,nl),gsphiz(nz,nf,nl)
    real :: gsrx(nz,nf,nl),gsry(nz,nf,nl),gsrz(nz,nf,nl)
    real :: ds(nz,nf,nl),d2s(nz,nf,nl),d22s(nz,nf,nl)
    real :: dels(nz,nf,nl)
    integer :: iz300s(nf,nl),iz300n(nf,nl)
    real :: grad_inp (nfp1)
    real :: xnorms(nzp1,nf,nl),ynorms(nzp1,nf,nl),znorms(nzp1,nf,nl)
    real :: xnormp(nz,nfp1,nl),ynormp(nz,nfp1,nl),znormp(nz,nfp1,nl)
    real :: xnormh(nz,nf,nlp1),ynormh(nz,nf,nlp1),znormh(nz,nf,nlp1)
    real :: xrg(nz,nf,nl),xthg(nz,nf,nl),xphig(nz,nf,nl)
    real :: qs(nz,nf,nl),brs(nz,nf,nl)

    real :: delpp(nz,nfp1,nl),delhp(nz,nfp1,nl)
    real :: delph(nz,nf,nlp1),delhh(nz,nf,nl)
    real :: delps(nzp1,nf,nl),delhs(nzp1,nf,nl)

    real :: phihp(nfp1,nlp1)
    real :: phish(nfp1,nl)
    real :: phisp(nf,nlp1)
    real :: eps(nzp1,nf,nl),eph(nz,nf,nlp1)
    real :: ehs(nzp1,nf,nl),ehp(nz,nfp1,nl)
    real :: jp(nz,nf,nl),jphi(nz,nf,nl)

    real :: bdirhx(nz,nf,nlp1), &
    bdirhy(nz,nf,nlp1), &
    bdirhz(nz,nf,nlp1)
    real :: bdirsx(nz,nf,nl), &
    bdirsy(nz,nf,nl), &
    bdirsz(nz,nf,nl)
    real :: bdirpx(nz,nfp1,nl), &
    bdirpy(nz,nfp1,nl), &
    bdirpz(nz,nfp1,nl)
    real :: blatss(nfp1,nlp1),blonss(nfp1,nlp1),bradss(nfp1,nlp1), &
    blatsp(nf,nlp1),blonsp(nf,nlp1),bradsp(nf,nlp1), &
    blatsh(nfp1,nl),blonsh(nfp1,nl),bradsh(nfp1,nl)



!       alts     altitude  (in km) on s mesh
!       grs      radial geographic distance to field line on s mesh
!       glats    geographic latitude on s mesh
!       glons    geographic longitude on s mesh
!       bms      normalized magnetic field on s mesh (b/b0)
!       ds,d2,
!       d22s     differential `distances' used in diff eqs
!       dels     actual arc length of grid in s direction
!       iz300s   index of grid at 300 km in south
!       iz300n   index of grid at 300 km in north
     
!     p grid data

! Worker grid

    real :: vol(nz,nf,nl)
    real :: areap(nz,nfp1,nl),areas(nzp1,nf,nl),areah(nz,nf,nlp1)
    real :: vpnx(nzp1,nfp1,nlp1),vpny(nzp1,nfp1,nlp1), &
    vpnz(nzp1,nfp1,nlp1)
    real :: vpsnx(nz,nf,nl),vpsny(nz,nf,nl),vpsnz(nz,nf,nl)
    real :: vhsnx(nz,nf,nl),vhsny(nz,nf,nl),vhsnz(nz,nf,nl)
    real :: vhnx(nzp1,nfp1,nlp1),vhny(nzp1,nfp1,nlp1), &
    vhnz(nzp1,nfp1,nlp1)
    real :: xdels(nz,nfp1,nlp1),xdelp(nzp1,nf,nlp1),xdelh(nzp1,nfp1,nl)
    real :: vexbs(nzp1,nf,nl),vexbp(nz,nfp1,nl),vexbh(nz,nf,nlp1)

    real :: vexbs_phi(nzp1,nf,nl), &
    vexbp_phi(nz,nfp1,nl), &
    vexbh_phi(nz,nf,nlp1)

    real :: vexb(nzp1,nfp1,nlp1)
    real :: hipc(nf,nl),hihc(nf,nl)
    real :: hipcp(nf,nl),hipcphi(nf,nl),hihcm(nf,nl)
    real :: hidv(nf,nl)
    real :: hidpv(nf,nl),hidphiv(nf,nl)
    real :: hidpg(nf,nl),hidphig(nf,nl)
    real :: sigmap(nz,nf,nl),sigmah(nz,nf,nl)
    real :: sigmapic(nz,nf,nl),sigmahic(nz,nf,nl)
    real :: blatp(nzp1,nfp1,nlp1),blatpt(nzp1,nfp1,nlt)
    real :: blonpt(nzp1,nfp1,nlt)
    real :: pp(nzp1,nfp1,nlp1),qp(nzp1,nfp1,nlp1),blonp(nzp1,nfp1,nlp1)
    real :: baltp(nzp1,nfp1,nlp1),brp(nzp1,nfp1,nlp1)

    real :: bmpf(nz,nfp1,nl),bmhf(nz,nf,nlp1),bmsf(nzp1,nf,nl)
    real :: epsx(nzp1,nf,nl),epsy(nzp1,nf,nl),epsz(nzp1,nf,nl)
    real :: ehsx(nzp1,nf,nl),ehsy(nzp1,nf,nl),ehsz(nzp1,nf,nl)
    real :: ehpx(nz,nfp1,nl),ehpy(nz,nfp1,nl),ehpz(nz,nfp1,nl)
    real :: ephx(nz,nf,nlp1),ephy(nz,nf,nlp1),ephz(nz,nf,nlp1)

    real :: deni_mnp(nz,nion),ti_mnp(nz,nion),te_mnp(nz)


!     vol        volume (i.e., area) of cell

!     chemical reaction data

    integer :: ichem(nchem,3)
!      real ireact(nion,nneut,nchem)

    common / chem / ichem

!     variables

! Total
! Total matrices are now allocated in main for the master
! well almost all

!      real u1pt(nz,nf,nlt),u2st(nz,nf,nlt),u3ht(nz,nf,nlt)
    real :: hipcpt(nf,nlt),hipcphit(nf,nlt),hihcmt(nf,nlt)
    real :: hipct(nf,nlt),hihct(nf,nlt),hidvt(nf,nlt)
    real :: hidpvt(nf,nlt),hidphivt(nf,nlt)
    real :: hidpgt(nf,nlt),hidphigt(nf,nlt)
!      real sigmapt(nz,nf,nlt),sigmaht(nz,nf,nlt)
!      real sigmapict(nz,nf,nlt),sigmahict(nz,nf,nlt)
!      real vnqt(nz,nf,nlt),vnpt(nz,nf,nlt),vnphit(nz,nf,nlt)
!      real xpt(nzp1,nfp1,nlt), ypt(nzp1,nfp1,nlt),
!     .                            zpt(nzp1,nfp1,nlt)
!      real d22st(nz,nf,nlt),pst(nz,nf,nlt),bmst(nz,nf,nlt)

!      real blonp0t(nlt+3),pst(nz,nf,nlt)
    real :: ppt(nzp1,nfp1,nlt)
    real*8 :: blonp0t(nlt+3)

!      real jpt(nz,nf,nlt),jphit(nz,nf,nlt)

! Worker
    real :: deni(nz,nf,nl,nion),ne(nz,nf,nl),denn(nz,nf,nl,nneut)
    real :: denni(nz,nf,nl,nneut),dennf(nz,nf,nl,nneut)
    real :: vsi(nz,nf,nl,nion),vsid(nz,nf,nl,nion)
    real :: sumvsi(nz,nf,nl,nion),vsic(nz,nf,nl,nion)
    real :: te(nz,nf,nl),ti(nz,nf,nl,nion),tn(nz,nf,nl)
    real :: tni(nz,nf,nl),tnf(nz,nf,nl)
    real :: u(nz,nf,nl),v(nz,nf,nl),vpi(nz,nf,nl),w(nz,nf,nl)
    real :: ui(nz,nf,nl),vi(nz,nf,nl),wi(nz,nf,nl)
    real :: uf(nz,nf,nl),vf(nz,nf,nl),wf(nz,nf,nl)
    real :: vnq(nz,nf,nl),vnp(nz,nf,nl),vnphi(nz,nf,nl)
    real :: vhi(nz,nf,nl),nuen(nz,nf,nl)

    real :: hruti,hrutf

!     velocity change positions

    real :: dxv(nz,nf,nl,nion),dyv(nz,nf,nl,nion),dzv(nz,nf,nl,nion)

!     atomic masses

    real :: ami(nion),amn(nneut),alpha0(nneut),aap(7)

!     zenith data

    real :: cosbdec(nz,nf,nl),sinbdec(nz,nf,nl),cx(nz,nf,nl), &
    coschicrit(nz,nf,nl)

    common / geom / cosbdec,sinbdec,cx,coschicrit

!     photodeposition rates
!     used 3 (absorption) and 7 (nion) explicitly
!     used 4 (number of angles in nighttime deposition)

    real :: sigabsdt(linesuv,3),flux(linesuv),sigidt(linesuv,7)
    real :: sigint(linesnt,7),fluxnt(nz,nf,nl,91,linesnt)
    real :: thetant(linesnt,4),zaltnt(linesnt,2)

    common / photdep / sigabsdt,flux,sigidt, &
    sigint,fluxnt,thetant,zaltnt

!     diagnostic variables

    real :: t1(nz,nf,nl,nion),t2(nz,nf,nl,nion),t3(nz,nf,nl,nion)
    real :: u1(nz,nf,nl),u2(nz,nf,nl),u3(nz,nf,nl), &
    u4(nz,nf,nl),u5(nz,nf,nl)
    real :: u1p(nz,nf,nl),u2s(nz,nf,nl),u3h(nz,nf,nl)

! Message passing stuff
           
! this is a buffer that holds deni(nz,nion), ti(nz,nion), and te(nz)

    real :: tl1s(nz,nf,nion+nion+1), tl1r(nz,nf,nion+nion+1)
    real :: tr1s(nz,nf,nion+nion+1), tr1r(nz,nf,nion+nion+1)

    integer,parameter :: MASTER = 0
    logical :: flagit, flagit1, flagit10
    integer :: taskid, source, dest, numworkers

end module com3
