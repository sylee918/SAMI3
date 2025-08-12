!******************************************
!******************************************
!             courant
!******************************************
!******************************************

    subroutine courant

    use parameter_mod
    use variable_mod
    use namelist_mod
    use time_mod
    use exb_mod
    use grid_mod

! parallel motion

    dtnew = 1.e6
    do k = nion1,nion2
        do l = 1,nl
            do j = 1,nf
                do i = 1,nz
                    dt1 = dels(i,j,l) / amax1(1.,abs(vsi(i,j,l,k)))
                    if ( dt1 <= dtnew ) then
                        dtnew = dt1
                        i0    = i
                        j0    = j
                        k0    = k
                        l0    = l
                    endif
                enddo
            enddo
        enddo
    enddo

! perpendicular motion (traditionaly limited to O+)

    do k = 1,nl
        do j = 1,nf
            do i = 1,nz
                dts = xdels(i,j,k) / amax1(1.,abs(vexbs(i,j,k,ptop)))
                dtp = xdelp(i,j,k) / amax1(1.,abs(vexbp(i,j,k,ptop)))
                dth = xdelh(i,j,k) / amax1(1.,abs(vexbh(i,j,k,ptop)))
                dt1 = amin1 ( dts,dtp,dth )
                if ( dt1 <= dtnew ) then
                    dtnew = dt1
                    i0    = i
                    j0    = j
                    k0    = k
                    dtsnew = dts
                    dtpnew = dtp
                    dthnew = dth
                endif
            enddo
        enddo
    enddo

    ! Consider perpendicular motion of metals a limiting factor too (DPD)
    
    ! do k = 8,9
    !     do l = 1,nl
    !         do j = 1,nf
    !             do i = 1,nz
    !                 dts = xdels(i,j,l) / amax1(1.,abs(vexbs(i,j,l,k)))
    !                 dtp = xdelp(i,j,l) / amax1(1.,abs(vexbp(i,j,l,k)))
    !                 dth = xdelh(i,j,l) / amax1(1.,abs(vexbh(i,j,l,k)))
    !                 dt1 = amin1 ( dts,dtp,dth )
    !                 if ( dt1 <= dtnew ) then
    !                     dtnew = dt1
    !                     i0    = i
    !                     j0    = j
    !                     k0    = k
    !                     dtsnew = dts
    !                     dtpnew = dtp
    !                     dthnew = dth
    !                 endif
    !             enddo
    !         enddo
    !     enddo
    ! enddo
    
    if ( dtnew <= .01 ) then
        print *,' Time step too small'
        stop
    elseif ( dtnew >= 5e4 ) then
        print *,' Time step too big: dtnew',dtnew
        stop
    endif
    
    
    if ( dtnew <= dt ) then
        dt = dtnew
    else
        if ( dtnew > dt*1.1) then
            dt = dt * 1.1
        else
            dt = dtnew
        endif
    endif
        
    ! Added to mitgate problems with drift in
    ! split time step syncronization as the result
    ! round off error in single precision DPD
    
    !dt = anint(dt*100.0)*0.01 ! truncate to 0.01 second resolution

    return
    end subroutine courant


    subroutine smoothz(f)
     
    use parameter_mod

    real :: f(nz),f0(nz+2)
     
    do i = 1,nz
      f0(i+1) = f(i)
    enddo
     
!   zero-gradient in z
     
    f0(1)     = f0(2)
    f0(nz+2)  = f0(nz+1)

!   sweep in z (forward)
        
    do i = 1,nz+1
      f0 (i) = 0.5 * ( f0(i) + f0(i+1) )
    enddo

!   sweep in z (backward)
     
    do i = nz+2,2,-1
      f0 (i) = 0.5 * ( f0(i) + f0(i-1) )
    enddo
     
!       now get f
     
    do i = 1,nz
      f(i)  = f0(i+1)
    enddo
     
    return
    end subroutine smoothz
    
    


pure real(8) function kcc1(nx,ny,x,y,f,xp,yp)

    implicit none

! Input variables

    integer,intent(in)     :: nx
    integer,intent(in)     :: ny
    real(8),intent(in)     :: x(nx)
    real(8),intent(in)     :: y(ny)
    real(8),intent(in)     :: f(nx,ny)

    real(8),intent(in)     :: xp
    real(8),intent(in)     :: yp

! Internal variables

    real(8)     :: x1,x2,x3,x4
    real(8)     :: y1,y2,y3,y4
    integer     :: i1,i2,i3,i4
    integer     :: j1,j2,j3,j4
    integer     :: i,j

! Apply the Key's Cubic Convolution interpolation algorithm

! Calculate the interpolation kernel for the x direction

    i1 = int(xp) + 1
    i2 = i1 - 1
    i3 = i1 + 1
    i4 = i1 + 2
    x1 = abs(xp - x(i1))
    x2 = x1 + 1.0d0
    x3 = 1.0d0 - x1
    x4 = 2.0d0 - x1
    x1 =  1.5d0*x1*x1*x1 - 2.5d0*x1*x1 + 1.0d0
    x2 = -0.5d0*x2*x2*x2 + 2.5d0*x2*x2 - 4.0d0*x2 + 2.0d0
    x3 =  1.5d0*x3*x3*x3 - 2.5d0*x3*x3 + 1.0d0
    x4 = -0.5d0*x4*x4*x4 + 2.5d0*x4*x4 - 4.0d0*x4 + 2.0d0

    ! Calculate the interpolation kernel for the y direction

    j1 = int(yp) + 1
    j2 = j1 - 1
    j3 = j1 + 1
    j4 = j1 + 2
    y1 = abs(yp - y(j1))
    y2 = y1 + 1.0d0
    y3 = 1.0d0 - y1
    y4 = 2.0d0 - y1
    y1 =  1.5d0*y1*y1*y1 - 2.5d0*y1*y1 + 1.0d0
    y2 = -0.5d0*y2*y2*y2 + 2.5d0*y2*y2 - 4.0d0*y2 + 2.0d0
    y3 =  1.5d0*y3*y3*y3 - 2.5d0*y3*y3 + 1.0d0
    y4 = -0.5d0*y4*y4*y4 + 2.5d0*y4*y4 - 4.0d0*y4 + 2.0d0

        ! Sum the contributions from the eight adjacent data points

    kcc1 = y1*(f(i1,j1)*x1 + f(i2,j1)*x2 + f(i3,j1)*x3 + f(i4,j1)*x4) + &
           y2*(f(i1,j2)*x1 + f(i2,j2)*x2 + f(i3,j2)*x3 + f(i4,j2)*x4) + &
           y3*(f(i1,j3)*x1 + f(i2,j3)*x2 + f(i3,j3)*x3 + f(i4,j3)*x4) + &
           y4*(f(i1,j4)*x1 + f(i2,j4)*x2 + f(i3,j4)*x3 + f(i4,j4)*x4)

    return

end function kcc1
