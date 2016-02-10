MODULE mod_random
  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC :: Random_init, Random_close
  
CONTAINS
  !-----------------------------------------------------------
  SUBROUTINE random_init (nom, rang)
    
    IMPLICIT NONE
    CHARACTER (LEN=*), INTENT(IN)      :: nom
    
    INTEGER                            :: k, i, rang
    LOGICAL                            :: est_la
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed
    
    ! On lit la taille de la "graine"
    CALL RANDOM_SEED (SIZE=k)
    ALLOCATE (seed(k))
    INQUIRE (FILE=nom, EXIST=est_la)
    
    ! Si le fichier "seed" existe, on le lit
    IF (est_la) THEN
       OPEN (UNIT=10, FILE=nom, STATUS="old")
       READ (10,*) seed
    ! Sinon, on le cree avec une valeur par defaut
    ELSE
       seed = (/ (200*i+1, i=1,k) /)
       OPEN (UNIT=10, FILE=nom, STATUS="new")
       WRITE (10,*) seed
    END IF
    
    ! On impose le point de depart du generateur
    CALL RANDOM_SEED (put=seed + rang)
    CLOSE (10)
    DEALLOCATE (seed)
    
  END SUBROUTINE random_init
  
  !----------------------------------------------
  
  SUBROUTINE Random_close (nom)
    
    IMPLICIT NONE
    CHARACTER (LEN=*), INTENT(IN)      :: nom

    INTEGER                            :: k
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed

    ! On recupere la valeur courante de la "graine"
    CALL RANDOM_SEED (SIZE=k)
    ALLOCATE (seed(k))
    CALL RANDOM_SEED (get=seed)

    ! On la sauvegarde a la place de la precedente
    ! Inutile de tester l existence de "nom", il a ete cree par Random_init si necessaire
    OPEN (UNIT=10, FILE=nom, STATUS="old")
    WRITE (10,*) seed
    CLOSE (10)
    DEALLOCATE (seed)

  END SUBROUTINE Random_close
END MODULE mod_random


