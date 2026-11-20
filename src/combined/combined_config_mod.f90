! Configuration of kicked_top_combined, read from a namelist file (default
! config_combined.nml). One parameter set drives the classical kicked top and
! the qubit chain for every rotation angle h in the list h(:).
module combined_config_mod
  use kinds_mod, only: wp
  use io_utils_mod, only: check_nml_read
  implicit none
  private
  public :: rconfig_t, read_combined_config, str_len

  integer, parameter :: str_len = 256   ! length of string parameters
  integer, parameter :: max_h = 20      ! maximum number of h values

  type :: rconfig_t
    integer :: n_q                 ! number of qubits, j = n_q/2
    real(wp) :: eps                ! coupling: coefficient of H_zz (eps = k/(2j))
    real(wp), allocatable :: h(:)  ! kick: rotation angles about x, one run per value
    integer :: s0                  ! initial state: -1 all up, 0 all down (along pol)
    character(len=1) :: pol        ! axis of the initial polarisation: 'z', 'x', 'y'
    logical :: has_angles          ! theta0 and phi0 given: they replace s0 and pol
    real(wp) :: theta0, phi0       ! initial product coherent state |theta0, phi0>
    integer :: n_kicks             ! classical: kicks per orbit
    integer :: n_theta             ! classical: grid points in Z = cos(theta)
    integer :: n_phi               ! classical: grid points in phi
    integer :: t_r                 ! chain: forward steps before the reversal
    logical :: butterfly_effect    ! chain: flip qubit i_flip at t_r, before the backward leg
    integer :: i_flip              ! chain: index (1 .. n_q) of the flipped qubit
    character(len=1) :: flip_axis  ! chain: rotation axis a = 'x', 'y', 'z' of qubit i_flip
    real(wp) :: flip_angle         ! chain: rotation angle about a (pi = full flip sigma^a)
    integer :: hus_n_theta         ! Husimi grid points in Z
    integer :: hus_n_phi           ! Husimi grid points in phi
    character(len=str_len) :: out_dir   ! output directory
  end type rconfig_t

contains

  subroutine read_combined_config(path, cfg)
    character(len=*), intent(in) :: path
    type(rconfig_t), intent(out) :: cfg

    integer :: n_q, s0, n_kicks, n_theta, n_phi, t_r, i_flip
    logical :: butterfly_effect
    real(wp) :: eps, h(max_h), flip_angle, theta0, phi0
    character(len=1) :: pol, flip_axis
    character(len=str_len) :: out_dir
    integer :: u, ios

    namelist /system/ n_q, eps, h
    namelist /initial/ s0, pol, theta0, phi0
    namelist /classical/ n_kicks, n_theta, n_phi
    namelist /chain/ t_r, butterfly_effect, i_flip, flip_axis, flip_angle
    namelist /husimi/ n_theta, n_phi
    namelist /output/ out_dir

    n_q = 10
    eps = 2.0_wp
    h = huge(1.0_wp)
    h(1:3) = [0.02_wp, 0.1_wp, 1.0_wp]
    s0 = -1
    pol = 'z'
    theta0 = huge(1.0_wp)
    phi0 = huge(1.0_wp)
    t_r = 30
    butterfly_effect = .false.
    i_flip = 1
    flip_axis = 'x'
    flip_angle = acos(-1.0_wp)
    out_dir = 'output_combined'

    open(newunit=u, file=path, status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'combined_config_mod: cannot open namelist file '//trim(path)

    read(u, nml=system, iostat=ios)
    call check_nml_read(ios, 'system', 'combined_config_mod')
    rewind(u)
    read(u, nml=initial, iostat=ios)
    call check_nml_read(ios, 'initial', 'combined_config_mod')
    rewind(u)
    n_kicks = 500
    n_theta = 20
    n_phi = 20
    read(u, nml=classical, iostat=ios)
    call check_nml_read(ios, 'classical', 'combined_config_mod')
    cfg%n_kicks = n_kicks
    cfg%n_theta = n_theta
    cfg%n_phi = n_phi
    rewind(u)
    read(u, nml=chain, iostat=ios)
    call check_nml_read(ios, 'chain', 'combined_config_mod')
    rewind(u)
    n_theta = 100
    n_phi = 200
    read(u, nml=husimi, iostat=ios)
    call check_nml_read(ios, 'husimi', 'combined_config_mod')
    cfg%hus_n_theta = n_theta
    cfg%hus_n_phi = n_phi
    rewind(u)
    read(u, nml=output, iostat=ios)
    call check_nml_read(ios, 'output', 'combined_config_mod')
    close(u)

    cfg%n_q = n_q
    cfg%eps = eps
    cfg%h = pack(h, h < huge(1.0_wp))
    cfg%s0 = s0
    cfg%pol = pol
    cfg%has_angles = theta0 < huge(1.0_wp) .and. phi0 < huge(1.0_wp)
    cfg%theta0 = theta0
    cfg%phi0 = phi0
    cfg%t_r = t_r
    cfg%butterfly_effect = butterfly_effect
    cfg%i_flip = i_flip
    cfg%flip_axis = flip_axis
    cfg%flip_angle = flip_angle
    cfg%out_dir = out_dir

    call validate(cfg)
  end subroutine read_combined_config

  subroutine validate(cfg)
    type(rconfig_t), intent(in) :: cfg

    if (cfg%n_q < 1 .or. cfg%n_q > 26) error stop 'combined_config_mod: n_q must be in 1 .. 26'
    if (size(cfg%h) < 1) error stop 'combined_config_mod: give at least one value of h'
    if (cfg%n_kicks < 1) error stop 'combined_config_mod: n_kicks must be >= 1'
    if (cfg%n_theta < 1 .or. cfg%n_phi < 1) error stop 'combined_config_mod: n_theta and n_phi must be >= 1'
    if (cfg%t_r < 1) error stop 'combined_config_mod: t_r must be >= 1'
    if (cfg%i_flip < 1 .or. cfg%i_flip > cfg%n_q) error stop 'combined_config_mod: i_flip must be in 1 .. n_q'
    if (index('xyz', cfg%flip_axis) == 0) error stop "combined_config_mod: flip_axis must be 'x', 'y' or 'z'"
    if (cfg%flip_angle < 0.0_wp) error stop 'combined_config_mod: flip_angle must be >= 0'
    if (cfg%hus_n_theta < 1 .or. cfg%hus_n_phi < 1) error stop 'combined_config_mod: Husimi grid must be >= 1'
  end subroutine validate

end module combined_config_mod
