module mod_variables
  use mod_random_n
  use mod_constants
  private
  public :: init, newton, leapfrog, kinetic
contains
  !------------------------------------------------------
  subroutine init(masse, position)
    implicit none

    character(len = 40)                                       :: jnk = 'jnk'
    integer                                                   :: i
    real(kind=xp), intent(inout), dimension(nb_particules)    :: masse
    real(kind=xp), intent(inout), dimension(3, nb_particules) :: position
    real(kind=xp)                                             :: u1, u2, u3
    real(kind=xp)                                             :: mu
 
    call random_init(jnk)
    call random_number(masse)
    mu    = sum(masse)
    masse = masse / mu
    
    
    i = 1
    open(13, file = "position.dat", position = 'append')
    do while (i <= nb_particules)
       call random_number(u1)
       call random_number(u2)
       call random_number(u3)
       u1 = 2._xp * u1 - 1
       u2 = 2._xp * u2 - 1
       u3 = 2._xp * u3 - 1
       if (sqrt(u1**2 + u2**2 + u3**2) < 1) then 
          position(1,i) = u1
          position(2,i) = u2
          position(3,i) = u3
          write(13,*)position(1,i),position(2,i),position(3,i) 
          i = i + 1
       endif
    enddo
    close(13)
    call random_close(jnk)
  end subroutine init

  subroutine newton(position, masse, acceleration, potentiel)
    implicit none
    integer                                                 :: i, j
    real(kind=xp), intent(in), dimension(nb_particules)     :: masse
    real(kind=xp), intent(in), dimension(3, nb_particules)  :: position
    real(kind=xp), intent(out), dimension(3, nb_particules) :: acceleration
    real(kind=xp), intent(out)                              :: potentiel
    real(kind=xp)                                           :: distance

    potentiel = 0._xp
    acceleration = 0._xp

    !$omp parallel private(distance, j)reduction(+:potentiel)
    !$omp do schedule(dynamic,50)
    do i=1, nb_particules
       do j=1, nb_particules
          if (i .ne. j) then
             distance  = sqrt( (position(1,j) - position(1,i))**2 + &
                  (position(2,j) - position(2,i))**2 + (position(3,j) -  &
                  position(3,i))**2 + epsilon**2)

             acceleration(:,i) = acceleration(:,i) +  masse (j) / &
                  abs(distance**3) * (position(:,j) - position(:,i))
           
             potentiel  = potentiel + ( - (masse(i) * masse(j)) &
                  / abs(distance) / 2._xp )
          endif
       enddo
    enddo
    !$omp end do
    !$omp end parallel
  end subroutine newton

  subroutine leapfrog (acceleration, dt, vitesse, position)
    implicit none
    real(kind=xp), intent(in), dimension(3, nb_particules)    :: acceleration
    real(kind=xp), intent(inout), dimension(3, nb_particules) :: position
    real(kind=xp), intent(inout), dimension(3, nb_particules) :: vitesse
    real(kind=xp), intent(in)                                 :: dt

    vitesse  = vitesse + acceleration * dt
    position = position + vitesse * dt
  end subroutine leapfrog

  subroutine kinetic(vitesse, masse, cinetique)
    implicit none
    integer                                                   :: i
    real(kind=xp), intent(in), dimension(nb_particules)       :: masse
    real(kind=xp), intent(inout), dimension(3, nb_particules) :: vitesse
    real(kind=xp),intent(out)                                 :: cinetique

    cinetique = 0._xp
    !$omp parallel reduction(+:cinetique)
    !$omp do schedule(dynamic,10)
    do i=1, nb_particules
       cinetique = cinetique + 0.5_xp * masse(i) * (vitesse(1,i)**2 + &
            vitesse(2,i)**2 + vitesse(3,i)**2)
    enddo
    !$omp end do
    !$omp end parallel
  end subroutine kinetic
end module mod_variables
