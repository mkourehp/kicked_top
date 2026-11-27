! Configuration of kicked_top_scan (make scan), read from a namelist file (default
! config_scan.nml): the classical map with (h, k), the spin-j top with j = n_q/2 and
! eps = k/(2j) = k/n_q, and the n_q-qubit chain that runs the echo.
module scan_config_mod
  use kinds_mod, only: wp
  use io_utils_mod, only: check_nml_read
  implicit none
  private
  public :: sconfig_t, read_scan_config, str_len

  integer, parameter :: str_len = 256

  type :: sconfig_t
    ! &system
    real(wp) :: h                  ! kick: rotation angle about x
    real(wp) :: k                  ! torsion strength; eps = k/n_q
    integer :: n_q                 ! number of qubits, j = n_q/2
    ! &scan
    integer :: n_theta, n_phi      ! equal-area grid of initial coherent states (theta0, phi0)
    integer :: n_lyap              ! kicks of the finite-time Lyapunov exponent
    integer :: n_renorm            ! rescale the tangent vector every n_renorm kicks
    integer :: t_scan              ! kicks of the spin-j evolution of every candidate
    real(wp) :: lam_regular        ! island candidates: Lambda < lam_regular
    real(wp) :: lam_chaotic        ! sea candidates: Lambda > lam_chaotic
    real(wp) :: sat_tol            ! sea candidates: |PR_late/PR_early - 1| < sat_tol
    ! &classical
    integer :: n_kicks             ! kicks per orbit of the phase portrait
    integer :: sec_n_theta, sec_n_phi   ! grid of orbits of the phase portrait
    ! &chain
    integer :: t_r                 ! forward kicks before the reversal
    integer :: i_flip              ! index (1 .. n_q) of the rotated qubit
    character(len=1) :: flip_axis  ! rotation axis a = 'x', 'y', 'z'
    real(wp) :: flip_angle         ! rotation angle theta_f of exp(-i theta_f sigma^a/2)
    integer :: n_angle             ! flip-angle scan: angles pi*l/(n_angle-1) for all axes (0: off)
    integer :: tr_step, tr_max     ! t_r scan: t_r = tr_step, 2 tr_step, ..., tr_max (0: off)
    ! &husimi
    integer :: hus_n_theta, hus_n_phi
    ! &output
    character(len=str_len) :: out_dir
  end type sconfig_t

contains

  subroutine read_scan_config(path, cfg)
    character(len=*), intent(in) :: path
    type(sconfig_t), intent(out) :: cfg
    real(wp) :: h, k, lam_regular, lam_chaotic, sat_tol, flip_angle
    integer :: n_q, n_theta, n_phi, n_lyap, n_renorm, t_scan, n_kicks, t_r, i_flip, n_angle, tr_step, tr_max
    character(len=1) :: flip_axis
    character(len=str_len) :: out_dir
    integer :: u, ios

    namelist /system/ h, k, n_q
    namelist /scan/ n_theta, n_phi, n_lyap, n_renorm, t_scan, lam_regular, lam_chaotic, sat_tol
    namelist /classical/ n_kicks, n_theta, n_phi
    namelist /chain/ t_r, i_flip, flip_axis, flip_angle, n_angle, tr_step, tr_max
    namelist /husimi/ n_theta, n_phi
    namelist /output/ out_dir

    h = 0.5_wp*acos(-1.0_wp)
    k = 2.4_wp
    n_q = 12
    open(newunit=u, file=path, status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'scan_config_mod: cannot open namelist file '//trim(path)
    read(u, nml=system, iostat=ios)
    call check_nml_read(ios, 'system', 'scan_config_mod')
    cfg%h = h
    cfg%k = k
    cfg%n_q = n_q

    n_theta = 40
    n_phi = 80
    n_lyap = 2000
    n_renorm = 10
    t_scan = 200
    lam_regular = 0.01_wp
    lam_chaotic = 0.1_wp
    sat_tol = 0.05_wp
    rewind(u)
    read(u, nml=scan, iostat=ios)
    call check_nml_read(ios, 'scan', 'scan_config_mod')
    cfg%n_theta = n_theta
    cfg%n_phi = n_phi
    cfg%n_lyap = n_lyap
    cfg%n_renorm = n_renorm
    cfg%t_scan = t_scan
    cfg%lam_regular = lam_regular
    cfg%lam_chaotic = lam_chaotic
    cfg%sat_tol = sat_tol

    n_kicks = 1000
    n_theta = 12
    n_phi = 12
    rewind(u)
    read(u, nml=classical, iostat=ios)
    call check_nml_read(ios, 'classical', 'scan_config_mod')
    cfg%n_kicks = n_kicks
    cfg%sec_n_theta = n_theta
    cfg%sec_n_phi = n_phi

    t_r = 30
    i_flip = 1
    flip_axis = 'y'
    flip_angle = acos(-1.0_wp)
    n_angle = 0
    tr_step = 0
    tr_max = 0
    rewind(u)
    read(u, nml=chain, iostat=ios)
    call check_nml_read(ios, 'chain', 'scan_config_mod')
    cfg%t_r = t_r
    cfg%i_flip = i_flip
    cfg%flip_axis = flip_axis
    cfg%flip_angle = flip_angle
    cfg%n_angle = n_angle
    cfg%tr_step = tr_step
    cfg%tr_max = tr_max

    n_theta = 100
    n_phi = 200
    rewind(u)
    read(u, nml=husimi, iostat=ios)
    call check_nml_read(ios, 'husimi', 'scan_config_mod')
    cfg%hus_n_theta = n_theta
    cfg%hus_n_phi = n_phi

    out_dir = 'output_scan'
    rewind(u)
    read(u, nml=output, iostat=ios)
    call check_nml_read(ios, 'output', 'scan_config_mod')
    cfg%out_dir = out_dir
    close(u)

    if (cfg%n_q < 2) error stop 'scan_config_mod: n_q must be >= 2'
    if (cfg%t_scan < 4) error stop 'scan_config_mod: t_scan must be >= 4'
    if (cfg%t_r < 2) error stop 'scan_config_mod: t_r must be >= 2'
    if (cfg%i_flip < 1 .or. cfg%i_flip > cfg%n_q) error stop 'scan_config_mod: i_flip must be in 1 .. n_q'
    if (index('xyz', cfg%flip_axis) == 0) error stop "scan_config_mod: flip_axis must be 'x', 'y' or 'z'"
    if (cfg%n_angle == 1) error stop 'scan_config_mod: n_angle must be 0 or >= 2'
  end subroutine read_scan_config

end module scan_config_mod
