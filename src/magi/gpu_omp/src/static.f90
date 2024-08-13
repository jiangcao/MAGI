MODULE static
    !! constants

    IMPLICIT NONE

    integer, parameter :: dp = 8
    complex(dp), parameter :: cone = dcmplx(1.0d0, 0.0d0)
    complex(dp), parameter :: czero = dcmplx(0.0d0, 0.0d0)
    complex(dp), parameter :: c1i = dcmplx(0.0d0, 1.0d0)
    real(dp), parameter :: hbar = 1.0546d-34 ! m^2 kg / s
    real(dp), parameter :: m0 = 9.109d-31 ! kg
    real(dp), parameter :: eps0 = 8.854d-12 ! C/V/m
    real(dp), parameter :: c0 = 2.998d8 ! m/s  v light
    real(dp), parameter :: e0 = 1.6022d-19 ! C electron charge
    REAL(dp), PARAMETER :: pi = 3.14159265359d0
    REAL(dp), PARAMETER :: twopi = 3.14159265359d0*2.0d0

    REAL(dp), PARAMETER  :: m0_ev = 5.6856D-16 !eV s^2 / cm^2   rest mass of electron
    real(dp), PARAMETER  :: kg2eV = m0_ev/m0    ! eV s^2 / cm^2 / gram
    REAL(dp), PARAMETER  :: hbar_ev = 6.58211899E-16 !eV s
    REAL(dp), PARAMETER  :: BOLTZ = 8.61734d-05 !eV K-1
    real(dp), parameter :: hb2m = 7.6305d-16              ! eV*cm^2

END MODULE static
