module mod_variables
  use mod_constants
  use mod_random
  use mpi
  private
  public :: init, newton, leapfrog, kinetic
contains
  !--------------------------------------------------------------------
  subroutine init(masse, position)
    implicit none
    character(len = 40)                                       :: jnk = 'jnk'
    integer                                                   :: i
    real(kind=xp), intent(inout), dimension(nb_particules)    :: masse
    real(kind=xp), intent(inout), dimension(3, nb_particules) :: position
    real(kind=xp)                                             :: u1, u2, u3
    real(kind=xp)                                             :: mu
 
   ! call random_init(jnk)
    call random_number(masse)
    mu    = sum(masse)
    masse = masse / mu
    
    i = 1
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
          i = i + 1
       endif
    enddo
   ! call random_close(jnk)
  end subroutine init

  subroutine newton(position, masse, acceleration, potentiel, nb_proces, rang)      
    implicit none
    integer                                                 :: i, j
    real(kind=xp), intent(in), dimension(nb_particules)     :: masse
    real(kind=xp), intent(in), dimension(3, nb_particules)  :: position
    real(kind=xp), intent(out), dimension(3, nb_particules) :: acceleration
    real(kind=xp), dimension(3,nb_particules)               :: acc_tmp
    real(kind=xp), intent(inout)                            :: potentiel
    real(kind=xp)                                           :: V_tmp
    real(kind=xp)                                           :: distance

    integer                                                 :: err, rang, nb_proces
    integer                                                 :: count
  
    potentiel = 0._xp
    V_tmp = 0._xp
    acceleration = 0._xp
    acc_tmp = 0._xp

    count = nb_particules/nb_proces
    if (mod(nb_particules,nb_proces) .ne. 0) then
       stop 'modulo(N,n_proces) diff 0'
    endif
    
    do i=1+(rang*count), (rang+1)*count
       do j=1 ,nb_particules
          if (i .ne. j) then
             distance  = sqrt( (position(1,j) - position(1,i))**2 + &
                  (position(2,j) - position(2,i))**2 + (position(3,j) -  &
                  position(3,i))**2 + epsilon**2)

             acceleration(:,i) = acceleration(:,i) +  masse (j) / &
                  abs(distance**3) * (position(:,j) - position(:,i))
           
             V_tmp  = V_tmp + ( - (masse(i) * masse(j)) &
                  / abs(distance) / 2._xp )
          endif
       enddo
       call MPI_REDUCE(V_tmp, potentiel, 1, mpi_double_precision, &
            MPI_SUM, 0, mpi_comm_world, err)
       acc_tmp(:,i) = acceleration(:,i) 
    enddo
    call mpi_allgather(acc_tmp(1,1+(rang*count)), 3*count, & 
         mpi_double_precision, acceleration(1,1), 3*count,& 
         mpi_double_precision, mpi_comm_world, err)
  end subroutine newton
  !--------------------------------------------------------------------------
  subroutine leapfrog (acceleration, dt, vitesse, position)
    implicit none
    real(kind=xp), intent(in), dimension(3, nb_particules)    :: acceleration
    real(kind=xp), intent(inout), dimension(3, nb_particules) :: position
    real(kind=xp), intent(inout), dimension(3, nb_particules) :: vitesse
    real(kind=xp), intent(in)                                 :: dt

    vitesse  = vitesse + acceleration * dt
    position = position + vitesse * dt
  end subroutine leapfrog

  subroutine kinetic(vitesse, masse, cinetique, nb_proces, rang)
    implicit none
    integer                                                   :: i
    real(kind=xp), intent(in), dimension(nb_particules)       :: masse
    real(kind=xp), intent(inout), dimension(3, nb_particules) :: vitesse
    real(kind=xp), intent(inout)                              :: cinetique
    real(kind=xp)                                             :: cin_tmp
    integer                                                   :: nb_proces, rang, err
    integer                                                   :: count

    cin_tmp = 0._xp
    count = nb_particules / nb_proces
    cinetique = 0._xp
    do i=1+(rang*count), (rang+1)*count
       cin_tmp = cin_tmp + 0.5_xp * masse(i) * (vitesse(1,i)**2 + &
            vitesse(2,i)**2 + vitesse(3,i)**2)
    enddo
    call MPI_REDUCE(cin_tmp, cinetique, 1, mpi_double_precision, &
            MPI_SUM, 0, mpi_comm_world, err)
  end subroutine kinetic
 
end module mod_variables
