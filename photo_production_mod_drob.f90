
module photo_production_mod

  use parameter_mod

!     zenith data 
! 
    real :: cosbdec(nz,nf,nl),sinbdec(nz,nf,nl),cx(nz,nf,nl), & 
            coschicrit(nz,nf,nl) 
!     photodeposition rates 
!     used 3 (absorption) and 7 (nion) explicitly 
!     used 4 (number of angles in nighttime deposition) 

!!$    real :: sigabsdt(linesuv,3),flux(linesuv),sigidt(linesuv,nneut) 
    real :: sigint(linesnt,nneut)
  
    real :: thetant(linesnt,4),zaltnt(linesnt,2),fluxnt(nz,nf,nl,91,linesnt) 
 
    real :: sigabsdt_fism(linesuv_fism,3)
    real :: sigidt_fism(linesuv_fism,nneut) 
    real :: flux_fism(linesuv_fism)

    real :: sigabsdt_euvac(linesuv_euvac,3)
    real :: sigidt_euvac(linesuv_euvac,nneut) 
    real :: flux_euvac(linesuv_euvac)
    real :: fluxdat_euvac(linesuv_euvac,2)

!    real,allocatable :: cosbdec(:,:,:),sinbdec(:,:,:),cx(:,:,:), &
!            coschicrit(:,:,:)
!
!     photodeposition rates
!     used 3 (absorption) and 7 (nion) explicitly
!     used 4 (number of angles in nighttime deposition)

!!$    real,allocatable :: sigabsdt(:,:),flux(:),sigidt(:,:)
!!$    real,allocatable :: sigint(:,:),fluxnt(:,:,:,:,:)
!!$    real,allocatable :: thetant(:,:),zaltnt(:,:)

!    real :: emask(nz) = 1.0
    
    integer,parameter :: ncfx = 2000
    integer,parameter :: ncfchi = 1800
    real(8)           :: cfx(ncfx),cfchi(ncfchi),chapfunc(ncfx,ncfchi)
        

contains
   
    real(8) function chapman(x,chi)
    
        implicit none
        
        real(8),intent(in) :: x
        real(8),intent(in) :: chi
        real(8),external   :: kcc1  
        real(8),external   :: atm8_chapman
        real(8)            :: arg
        real(8)             :: xb,chib
          
        xb = max(x,1.0d0)
!        xb = min(xb,1998.0d0)
        xb = min(xb,1997.0d0)
        chib = max(chi,1.0d0)
!        chib = min(chib,178.0d0)    
        chib = min(chib,177.0d0)    
        arg = chib/0.1d0
        chapman = kcc1(ncfx,ncfchi,cfx,cfchi,chapfunc,xb,arg)
        chapman = exp(chapman)
        chapman = min(chapman,1.0e22)
    
        return
        
    end function chapman



    subroutine allocphoto(nz,nf,nl,linesuv,nneut,linesnt)

        integer,intent(in) :: nz,nf,nl,linesuv,nneut,linesnt

!        allocate(cosbdec(nz,nf,nl),sinbdec(nz,nf,nl),cx(nz,nf,nl), &
!                coschicrit(nz,nf,nl))
!        allocate(sigabsdt(linesuv,3),flux(linesuv),sigidt(linesuv,nneut))
!        allocate(sigint(linesnt,nneut),fluxnt(linesnt,91,nz,nf,nl))
!        allocate(thetant(linesnt,4),zaltnt(linesnt,2))
        return

    end subroutine allocphoto

end module photo_production_mod
