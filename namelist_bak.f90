!     namelist data

module namelist

    use param3

    logical :: hall,restart
    logical :: lmadala,lcr,lvs,lweimer,lhwm93,lhwm14
    real :: snn(nneut)
    integer :: psmooth

    integer :: maxstep,mmass,nz1,nz2,nz3,nz4

    real :: hrmax,dthr,hrpr,dt0, &
            grad_in,glat_in,glon_in, &
            rmin,rmax, &
            altmin, &
            fbar,f10p7,ap, &
            year,day, &
            nion1,nion2,hrinit,tvn0,tvexb0,ver,veh,vw, &
            gams1,gams1m,gamp1, &
            gams2,gams2m,gamp2, &
            gams3,gams3m,gamp3, &
            gams4,gams4m,gamp4, &
            r_min1,r_max1, &
            r_max2, &
            blat_max3,blat_max4, &
            stn,denmin,alt_crit,cqe,plat,plon, &
            dellon, &
            storm_ti,storm_tf,vexb_max, &
            decay_time,pcrit, &
            vsi0,delta_vsio,anu_drag0

end module namelist

