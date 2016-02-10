program projet_para
  use mod_constants
  use mod_variables
  use mpi
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
  
  integer                                    :: nb_proces, err, rang
 
  call mpi_init(err)
  
  call mpi_comm_rank(mpi_comm_world, rang, err)
  call mpi_comm_size(mpi_comm_world, nb_proces, err)
  if (rang == 0) then
     open(11, file="output.dat", status="replace", iostat=ios)
     if (ios /= 0) then
        stop "Error while opening output file."
     end if
     write(11, fmt='(4(a16))')'i', 'V', 'U', 'Etot'
  endif
  !--------Initialisation masse - position - vitesse--------
  call init(masse, position)
  vitesse(1,:) =  - w * position(2,:)
  vitesse(2,:) =  w * position(1,:)
  vitesse(3,:) = 0._xp

  position(1,:) = position(1,:) - dt / 2._xp * vitesse(1,:)
  position(2,:) = position(2,:) - dt / 2._xp * vitesse(2,:)
  position(3,:) = position(3,:) - dt / 2._xp * vitesse(3,:)
  !--------------Initialisation du potentiel----------------
 ! do i=1, nb_particules
 !    do j=1, nb_particules
 !       if (i .ne. j) then
 !         distance  = (position(1,j) - position(1,i))**2 + &
 !              (position(2,j) - position(2,i))**2 + (position(3,j) -  &
 !              position(3,i))**2 + epsilon**2
 !         potentiel = potentiel + ( - (masse(i) * masse(j)) &
 !              / abs(sqrt(distance)) / 2._xp )
 !      endif
 !    enddo
 ! enddo

  !------------Initialisation energie cinetique-------------
  call kinetic(vitesse, masse, cinetique, nb_proces, rang)

  !------------------Debut des iterations-------------------
  acceleration = 0._xp
  if (rang == 0) then
  !  open (12, file = 'acceleration.d', status = 'replace')
    open (12, file = 'position_open.d', status = 'replace')
    open (14, file = 'position_close.d', status = 'replace')
    write(12,*)position(1,:), position(2,:), position(3,:)
  endif
  do i=1, n_iteration-1
     E_tot_save  = potentiel + cinetique
     
     call newton (position, masse, acceleration, potentiel, nb_proces, rang)
     call leapfrog(acceleration, dt, vitesse, position)
     call kinetic(vitesse, masse, cinetique, nb_proces, rang)
     energie_tot = potentiel + cinetique
     
    ! write(*,*)i, potentiel, cinetique, energie_tot
     if (rang == 0) then 
        write(11,"(3F13.7)")potentiel, cinetique, energie_tot
     endif
     if (i >= 2 .and. abs((energie_tot - E_tot_save)/E_tot_save) >= limit) then
        stop 'pb'
     endif
     if (rang == 0) then
      !   write(12,"(1F13.7)") acceleration(1,1)
     endif
  enddo
  if (rang == 0) then
    write(14,*)position(1,:), position(2,:), position(3,:)
  endif
  if (rang == 0) then
    close(12)
  endif
  call mpi_finalize(err)
  if (rang == 0) then
     close(11)
  endif
end program projet_para
