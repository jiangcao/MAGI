module parameters_mod
    
    implicit none 

    !constants
    integer, parameter :: dp=8
    REAL(kind=dp), PARAMETER :: pi=3.14159265359d0
    REAL(kind=dp), PARAMETER :: twopi = 3.14159265359d0*2.0d0
    REAL(kind=dp), PARAMETER :: e_charge=1.6d-19            ! charge of an electron (C)
    REAL(kind=dp), PARAMETER :: epsilon0=8.85e-12    ! Permittivity of free space (m^-3 kg^-1 s^4 A^2)    
    REAL(kind=dp), PARAMETER :: light_speed=2.998d8           ! m/s
    REAL(kind=dp), PARAMETER :: m0_charge=5.6856D-16        ! eV s2 / cm2
    REAL(kind=dp), PARAMETER :: hbar=1.0546d-34     ! value of hbar=h/2pi (J s)
    REAL(kind=dp), PARAMETER :: hbar_eV=hbar/e_charge ! eV s    
    COMPLEX(kind=dp), PARAMETER :: cone = dcmplx(1.0d0,0.0d0)
    COMPLEX(kind=dp), PARAMETER :: czero  = dcmplx(0.0d0,0.0d0)
    COMPLEX(kind=dp), PARAMETER :: c1i  = dcmplx(0.0d0,1.0d0)     
    REAL(kind=dp), PARAMETER  :: BOLTZ = 8.61734d-05 !eV K-1
   
end module parameters_mod