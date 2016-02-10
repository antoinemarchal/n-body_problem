program projet_para
  use mod_constants
  use mod_variables
  implicit none

  integer                                    :: i, j, ios = 0 
  real(kind=xp), dimension(nb_particules)    :: masse
  real(kind=xp), dimension(3, nb_particules) :: position
  real(kind=xp), dimension(3, nb_particules) :: vitesse
  real(kind=xp), dimension(3, nb_particules) :: acceleration
  real(kind=xp)                              :: potentiel 
  real(kind=xp)                              :: cinetique
  real(kind=xp)                              :: energie_tot, E_tot_save
  real(kind=xp)                              :: distance

  !------------Ouverture du fichier output.dat--------------
  open(11, file="output.dat", status="replace", iostat=ios)
  if (ios /= 0) then
     stop "Error while opening output file."
  end if
  write(11, fmt='(4(a16))')'i', 'V', 'U', 'Etot'

  !--------Initialisation masse - position - vitesse--------
  call init(masse, position)
  vitesse(1,:) =  - w * position(2,:)
  vitesse(2,:) =  w * position(1,:)
  vitesse(3,:) = 0._xp

  position(1,:) = position(1,:) - dt / 2._xp * vitesse(1,:)
  position(2,:) = position(2,:) - dt / 2._xp * vitesse(2,:)
  position(3,:) = position(3,:) - dt / 2._xp * vitesse(3,:)
!  vitesse = vitesse * acceleration * dt / 2._xp

  !--------------Initialisation du potentiel----------------
  !$omp parallel private(distance, j)reduction(+:potentiel)
  !$omp do schedule(dynamic,25)
  do i=1, nb_particules
     do j=1, nb_particules
        if (i .ne. j) then
          distance  = (position(1,j) - position(1,i))**2 + &
               (position(2,j) - position(2,i))**2 + (position(3,j) -  &
               position(3,i))**2 + epsilon**2
          potentiel = potentiel + ( - (masse(i) * masse(j)) &
               / abs(sqrt(distance)) / 2._xp )
       endif
     enddo
  enddo
  !$omp end do
  !$omp end parallel


  !------------Initialisation energie cinetique-------------
  call kinetic(vitesse, masse, cinetique)

  !------------------Debut des iterations-------------------
  do i=1, n_iteration-1
     E_tot_save  = potentiel + cinetique
     
     call newton (position, masse, acceleration, potentiel)
     call leapfrog(acceleration, dt, vitesse, position)
     call kinetic(vitesse, masse, cinetique)
     energie_tot = potentiel + cinetique
     
    ! write(*,*)i, potentiel, cinetique, energie_tot
     write(11,*)i, real(potentiel), real(cinetique), real(energie_tot)
     if (i >= 2 .and. abs(energie_tot - E_tot_save)/E_tot_save >= limit) then
        stop 'pb'
     endif
  enddo
  close(11)
  
  
end program projet_para
