!     namelist data

module namelist_mod

  use parameter_mod

    logical :: hall,restart
    logical :: lmadala,lcr,lvs,lweimer,lhwm93,lhwm14
    logical :: lfism

    integer :: psmooth,nion1,nion2
    integer :: maxstep,mmass,kp(8)
    integer :: nday(8),iyear

    real :: snn(nneut),fbar(8),f10p7(8),ap(8)
    real :: hrmax, dthr, hrpr, dt0, &
            rmin, altmin, &                            
            hrinit, tvn0, tvexb0, ver, veh, vw,&
            gams, gamp, alt_crit, cqe, alt_crit_avg
    real :: storm_ti, storm_tf, vexb_max, &
            decay_time, pcrit, anu_drag0, &
            blat_max, stn,denmin,blat_min,xhardy,tphi,&
            tmax,euv_fac

end module namelist_mod
