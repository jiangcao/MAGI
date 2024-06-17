
module parameters_mod   
    implicit none 
    !constants
    integer, parameter :: dp=8
    REAL(kind=dp), PARAMETER :: pi=3.14159265359_dp
    REAL(kind=dp), PARAMETER :: twopi = 6.2831853072_dp
    REAL(kind=dp), PARAMETER :: e_charge=1.60217663e-19_dp      ! charge of an electron (C)
    REAL(kind=dp), PARAMETER :: epsilon0=8.8541878188e-12_dp    ! Permittivity of free space (m^-3 kg^-1 s^4 A^2)    
    REAL(kind=dp), PARAMETER :: light_speed=2.99792458e8_dp     ! m/s
    REAL(kind=dp), PARAMETER :: m0_charge=5.6856e-16_dp         ! electron mass eV s2 / cm2
    REAL(kind=dp), PARAMETER :: m0_kg=9.109e-31_dp         ! electron mass kg 
    REAL(kind=dp), PARAMETER :: hbar=1.05457182e-34_dp     ! value of hbar=h/2pi (J s)
    REAL(kind=dp), PARAMETER :: hbar_eV=6.582295486e-16_dp ! eV s    
    COMPLEX(kind=dp), PARAMETER :: cone = dcmplx(1.0_dp,0.0_dp)
    COMPLEX(kind=dp), PARAMETER :: czero  = dcmplx(0.0_dp,0.0_dp)
    COMPLEX(kind=dp), PARAMETER :: c1i  = dcmplx(0.0_dp,1.0_dp)    
    real(dp), parameter :: zero= 0.0_dp
    real(dp), parameter :: one= 1.0_dp
    real(dp), parameter :: tol6= 0.000001_dp
    real(dp), parameter :: tol12=0.000000000001_dp 
    REAL(kind=dp), PARAMETER  :: BOLTZ = 8.617333262e-05_dp !eV K-1 
end module parameters_mod


