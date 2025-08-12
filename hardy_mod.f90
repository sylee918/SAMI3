module hardy_mod

  use parameter_mod

  integer,parameter:: no_lamzor = 250
  integer,parameter:: no_atm_ps = 250

  real :: lambda(no_lamzor)
  real :: zor(no_lamzor)
  real :: height(no_atm_ps)

  real :: z(no_atm_ps),nNn(no_atm_ps),nM(no_atm_ps),rho(no_atm_ps)

  real :: preciprn(nz,nf,nl),preciprs(nz,nf,nl)
  real :: tpn(nz,nf,nl),tps(nz,nf,nl)
  integer:: iin(nf,nl),iis(nf,nl)

end module hardy_mod
