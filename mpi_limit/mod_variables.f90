module mod_variables
  use mod_constants
  use mod_random
  use mpi
  private
  public :: init, newton, leapfrog, kinetic
contains
  !--------------------------------------------------------------------
  subroutine init(masse, position, nb_proces, rang)
    implicit none
    character(len = 40)                                       :: jnk = 'jnk'
    integer                                                   :: i
    integer                                                   :: count

    real(kind=xp)                                             :: u1, u2, u3
    real(kind=xp)                                             :: mu
    real(kind=xp), intent(inout), dimension(:,:)              :: position
    real(kind=xp), intent(inout), dimension(:)                :: masse
    integer                                                   :: nb_proces, rang


    call random_init(jnk, rang)

    count = nb_particules/nb_proces
   ! if (mod(nb_particules,nb_proces) .ne. 0) then
   !    stop 'modulo(N,n_proces) diff 0'
   ! endif

    call random_number(masse)
    mu    = sum(masse)
    masse = masse / mu / (nb_proces+1)
        
    i = 1
    do while (i <= count)
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
    call random_close(jnk)
  end subroutine init


  subroutine newton(position, masse, acceleration, potentiel, &
       nb_proces, rang)      
    implicit none
    integer                                                 :: i, j,k
    real(kind=xp), intent(in), dimension(:,:)               :: position
    real(kind=xp), intent(in), dimension(:)                 :: masse
    real(kind=xp), dimension(:),allocatable                 :: masse_tmp
    real(kind=xp), dimension(:,:),allocatable               :: position_tmp
    real(kind=xp), intent(inout), dimension(:,:)            :: acceleration
    real(kind=xp), intent(inout)                            :: potentiel
    real(kind=xp)                                           :: V_tmp
    real(kind=xp)                                           :: distance
    integer                                                 :: err, rang, nb_proces
    integer                                                 :: count
  
    V_tmp = 0._xp
  !  position_tmp = 0._p
  !  masse_tmp = 0._xp
    count = nb_particules/nb_proces

    allocate(position_tmp(3,count))
    allocate(masse_tmp(count))

   
    
   do k=0, nb_proces - 1
     position_tmp = position
     masse_tmp = masse
     ! V_tmp = 0._xp
      call mpi_bcast(position_tmp(1,1),3*count,mpi_double_precision,&
            k ,mpi_comm_world,err)
      call mpi_bcast(masse_tmp(1),count,mpi_double_precision,&
           k ,mpi_comm_world,err)
      do i=1, count
         do j=1, count
            if (i .ne. j) then 
               distance  = sqrt( (position_tmp(1,j) - position(1,i))**2 &
                    + (position_tmp(2,j) - position(2,i))**2 + &
                    (position_tmp(3,j) - position(3,i))**2 + epsilon**2)
               
               acceleration(:,i) = acceleration(:,i) +  masse_tmp(j) / &
                    abs(distance**3) * (position_tmp(:,j) - position(:,i))
               
               V_tmp  = V_tmp + ( - (masse(i) * masse_tmp(j)) &
                    / abs(distance) / 2._xp )
            endif
         enddo
         call MPI_REDUCE(V_tmp, potentiel, 1, mpi_double_precision, &
        MPI_SUM, 0, mpi_comm_world, err)
      enddo
   enddo
   
   deallocate(masse_tmp)
   deallocate(position_tmp)   
  end subroutine newton
  !--------------------------------------------------------------------------
  subroutine leapfrog (acceleration, dt, vitesse, position, nb_proces, rang)
    implicit none
    real(kind=xp), intent(in), dimension(:,:)                 :: acceleration
    real(kind=xp), intent(inout), dimension(:,:)              :: position
    real(kind=xp), intent(inout), dimension(:,:)              :: vitesse
    real(kind=xp), intent(in)                                 :: dt
    integer                                                   :: nb_proces, rang
    integer                                                   :: count

    count = nb_particules/nb_proces

    vitesse  = vitesse + acceleration * dt
    position = position + vitesse * dt
  end subroutine leapfrog
  !--------------------------------------------------------------------------
  subroutine kinetic(vitesse, masse, cinetique, nb_proces, rang)
    implicit none
    integer                                                   :: i,count, nb_proces
    integer                                                   :: rang
    real(kind=xp), intent(in), dimension(:)                   :: masse
    real(kind=xp), intent(inout), dimension(:,:)              :: vitesse
    real(kind=xp),intent(inout)                               :: cinetique
    real(kind=xp)                                             :: cin_tmp
    integer                                                   :: err

    cin_tmp = 0._xp
    count = nb_particules/nb_proces
    do i=1, count
       cin_tmp = cin_tmp + 0.5_xp * masse(i) * (vitesse(1,i)**2 + &
            vitesse(2,i)**2 + vitesse(3,i)**2)
    enddo
  !  cin_tmp = cinetique
  !  if (rang .ne. 0) then    
       call MPI_REDUCE(cin_tmp, cinetique, 1, mpi_double_precision, &
            MPI_SUM, 0, mpi_comm_world, err)
  !  endif
  end subroutine kinetic 
end module mod_variables
