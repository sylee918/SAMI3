
  use parameter_mod
  use grid_mod
  use misc_mod

  real,parameter :: nt = 337

  real :: blat(nz,nf,nlt)
  real :: blon(nz,nf,nlt)
  real :: balt(nz,nf,nlt)
  real :: zalt(nz,nf,nlt)
  real :: qst(nz,nf,nlt)
  real :: pst(nz,nf,nlt)
  real :: L_n(nz),Ln_inv(nz)
  real :: ne(nz,nf,nlt,nt),tne(nz,nf,nlt)

  real :: q0,qp1,qp2,qm1,qm2,slope_p
  real :: del_ne_p,ne_p0,slope_m,del_ne_m,ne_m0
  real :: del_p0,ne_0

  integer :: ip(nz),im(nz)

  open(101,file=trim(dir)//'blatu.dat',form='unformatted')
  open(102,file=trim(dir)//'blonu.dat',form='unformatted')
  open(103,file=trim(dir)//'baltu.dat',form='unformatted')
  open(104,file=trim(dir)//'deneu.dat',form='unformatted')

  open(201,file=trim(dir)//'L_n.dat',form='unformatted')

  read(101) blat
  read(102) blon
  read(103) balt

  zalt = balt - 6370.

  print *,'min balt',minval(zalt)

  do n = 1,nt 
    read(104) tne
    ne(:,:,:,n) = tne
  enddo

  print *,'max ne = ',maxval(ne)

  print *,'input nfl/nll'
  read  *,nfl,nll

  print *,'nfl,nf',nfl,nf

  qst = (re/(re+zalt))**2 * sin(blat*pie/180.)
  pst = (re+zalt)/re / cos(blat*pie/180.) / cos(blat*pie/180.)

!  do i = 1,nz
!    print *,'qst,p,m = ',i,qst(i,nfl,nll),qst(i,nfl+1,nll),qst(i,nfl-1,nll)
!  enddo
 
  do i0 = 1,nz
    q0 = qst(i0,nfl,nll)     
    do i = 1,nz
      if ( nfl /= nf .AND. nfl /= 1 ) then
        if ( qst(i,nfl+1,nll) < q0 ) ip(i0) = i
        if ( qst(i,nfl-1,nll) < q0 ) im(i0) = i
      endif
      if ( nfl == 1 ) then
        if ( qst(i,3,nll) < q0 ) ip(i0) = i
        if ( qst(i,1,nll) < q0 ) im(i0) = i
      endif
      if ( nfl == nf ) then
        if ( qst(i,nf,nll)   < q0 ) ip(i0) = i
        if ( qst(i,nf-2,nll) < q0 ) im(i0) = i
      endif
    enddo
  enddo

  do n = 100,100
    do i = 1,nz
      q0 = qst(i,nfl,nll)     
      if ( nfl /= nf .AND. nfl /= 1 ) then
        qp1      = qst(ip(i),nfl+1,nll)
        qp2      = qst(ip(i)+1,nfl+1,nll)
        qm1      = qst(im(i),nfl-1,nll)
        qm2      = qst(im(i)+1,nfl-1,nll)
        slope_p  = (q0 - qp1)/(qp2 - qp1)
        del_ne_p = ne(ip(i)+1,nfl+1,nll,n) - ne(ip(i),nfl+1,nll,n)
        ne_p0    = ne(ip(i),nfl+1,nll,n) + del_ne_p * slope_p
        slope_m  = (q0 - qm1)/(qm2 - qm1)
        del_ne_m = ne(im(i)+1,nfl-1,nll,n) - ne(im(i),nfl-1,nll,n)
        ne_m0    = ne(im(i),nfl-1,nll,n) + del_ne_m * slope_m
        del_p0   = pst(i,nfl+1,nll) - pst(i,nfl-1,nll)
        ne_0     = ne(i,nfl,nll,n)
        Ln_inv(i) = (1./ne_0)*(ne_p0 - ne_m0)/del_p0
        if ( zalt(i,nfl,nll) < 100. ) Ln_inv(i) = 0.
      endif
      print *,'i,Ln_inv',i,zalt(i,nfl,nll),Ln_inv(i)
    enddo        
    write(201) Ln_inv
  enddo

  end
