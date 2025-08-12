
module photo_production_mod

  use parameter_mod

!     zenith data 
 
    real :: cosbdec(nz,nf,nl),sinbdec(nz,nf,nl),cx(nz,nf,nl), & 
            coschicrit(nz,nf,nl) 
 
!     photodeposition rates 
!     used 3 (absorption) and 7 (nion) explicitly 
!     used 4 (number of angles in nighttime deposition) 
 
    real :: sigabsdt_fism(linesuv_fism,3)
    real :: sigidt_fism(linesuv_fism,nneut) 
    real :: flux_fism(linesuv_fism)

    real :: sigabsdt_euvac(linesuv_euvac,3)
    real :: sigidt_euvac(linesuv_euvac,nneut) 
    real :: flux_euvac(linesuv_euvac)
    real :: fluxdat_euvac(linesuv_euvac,2)

    real :: sigint(linesnt,nneut),fluxnt(nz,nf,nl,91,linesnt) 
    real :: thetant(linesnt,4),zaltnt(linesnt,2) 

 
end module photo_production_mod
