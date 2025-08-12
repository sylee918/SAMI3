!******************************************

module RT_mod

    use parameter_mod

    real :: a_rt(nf,nl),b_rt(nf,nl),c_rt(nf,nl),cg_rt(nf,nl)
    real :: fs1(nf,nl),fs2(nf,nl),fs3(nf,nl)
    real :: Ln_inv(nz,nf,nl)

!    real :: a_rtt(nf,nlt),b_rtt(nf,nlt),c_rtt(nf,nlt)
!    real :: fs1t(nf,nlt),fs2t(nf,nlt),fs3t(nf,nlt)

end module RT_mod
