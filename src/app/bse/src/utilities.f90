module utilities

  contains

    FUNCTION STRING(inn)
        INTEGER, PARAMETER :: POS= 4
        INTEGER, INTENT(IN) :: inn
        CHARACTER(LEN=POS) :: STRING
        !............................................................
        INTEGER :: cifra, np, mm, num  
        IF (inn > (10**POS)-1) stop "ERROR: (inn > (10**3)-1)  in STRING"
        num= inn
        DO np= 1, POS
            mm= pos-np
            cifra= num/(10**mm)            
            STRING(np:np)= ACHAR(48+cifra)
            num= num - cifra*(10**mm)
        END DO
    END FUNCTION STRING


end module utilities
