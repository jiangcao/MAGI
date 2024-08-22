
module epw_interface

implicit none 

contains

! copied from QE 
SUBROUTINE davcio( vect, nword, unit, nrec, io )
  !----------------------------------------------------------------------------
  !! Direct-access vector input/output.  
  !! read/write \(\text{nword}\) words starting from the address specified by
  !! \(\text{vect}\).
  !
  USE parameters_mod ,     ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nword
  !! the dimension of vect
  INTEGER, INTENT(IN) :: unit
  !! the unit where to read/write
  INTEGER, INTENT(IN) :: nrec
  !! the record where to read/write
  INTEGER, INTENT(IN) :: io
  !! flag if < 0 reading if > 0 writing
  REAL(DP), INTENT(INOUT) :: vect(nword)
  !! input/output: the vector to read/write
  !
  INTEGER :: ios
    ! integer variable for I/O control
  LOGICAL :: opnd
  CHARACTER*256 :: name
  !  
  !
  IF ( unit  <= 0 ) print *, 'davcio', 'wrong unit', 1 
  IF ( nrec  <= 0 ) print *, 'davcio', 'wrong record number', 2 
  IF ( nword <= 0 ) print *, 'davcio', 'wrong record length', 3 
  IF ( io    == 0 ) print *, 'davcio', 'nothing to do?' 
  !
  INQUIRE( UNIT = unit, OPENED = opnd, NAME = name )
  !
  IF ( .NOT. opnd ) &
     print *,  'davcio', 'unit is not opened', unit 
  !
  ios = 0
  !
  IF ( io < 0 ) THEN
     !
     READ( UNIT = unit, REC = nrec, IOSTAT = ios ) vect
     IF ( ios /= 0 ) print *, 'davcio', &
         & 'error reading file "' // TRIM(name) // '"', unit 
     !
  ELSE IF ( io > 0 ) THEN
     !
     WRITE( UNIT = unit, REC = nrec, IOSTAT = ios ) vect
     IF ( ios /= 0 ) print *, 'davcio', &
         & 'error writing file "' // TRIM(name) // '"', unit 
     !
  END IF  
  !
  RETURN
  !
END SUBROUTINE davcio


end module epw_interface