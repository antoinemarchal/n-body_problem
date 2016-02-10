module mod_constants
  
  implicit none
  integer      , parameter, public :: xp = selected_real_kind(8)
  real(kind=xp), parameter, public :: pi = 4.*datan(1.d0)
  real(kind=xp), parameter, public :: c = 2.99792458e8_xp

  integer, parameter               :: nb_particules = 7680
  integer, parameter               :: n_iteration = 1000
  real(kind=xp), parameter         :: dt = 1.e-2_xp
  real(kind=xp), parameter         :: w = 1._xp, epsilon = 1.e-1_xp
  real(kind=xp), parameter         :: limit = 2.e-2_xp
end module mod_constants
