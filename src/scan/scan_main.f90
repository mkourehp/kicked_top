! kicked_top_scan: selection of two initial coherent states and their butterfly echo.
!  1. Equal-area grid of (theta0, phi0): classical finite-time Lyapunov exponent Lambda
!     (lyapunov_mod) and spin-j evolution for t_scan kicks (scan_mod).
!  2. Island state: Lambda < lam_regular with the smallest late-window PR.
!     Sea state: Lambda > lam_chaotic, saturated spread |PR_late/PR_early - 1| < sat_tol,
!     with the largest late-window PR.
!  3. Classical phase portrait and the orbits of both states.
!  4. Butterfly echo of both states on the n_q-qubit chain (chain_run_mod), with and without
!     the flip, and the optional flip-angle/axis and t_r scans.
! Usage: kicked_top_scan [config_scan.nml]
program scan_main
  use kinds_mod, only: wp
  use io_utils_mod, only: make_dir
  use sphere_mod, only: equal_area_grid, spherical_to_cartesian
  use config_mod, only: config_t
  use evolve_mod, only: evolve
  use lyapunov_mod, only: lyapunov
  use io_mod, only: write_results
  use floquet_mod, only: diagonalise_jx, floquet_operator
  use scan_config_mod, only: sconfig_t, read_scan_config
  use scan_mod, only: spin_series, spin_scan
  use combined_config_mod, only: rconfig_t
  use chain_state_mod, only: product_coherent_state, dicke_amplitudes
  use chain_floquet_mod, only: uniform_coupling, zz_energies, rotate_qubit, forward_step, backward_step
  use chain_observables_mod, only: fidelity, site_spins
  use chain_run_mod, only: run_chain
  implicit none

  character(len=*), parameter :: prog = 'kicked_top_scan'
  character(len=6), parameter :: names(2) = ['island', 'sea   ']
  type(sconfig_t) :: cfg
  type(config_t) :: c
  type(rconfig_t) :: r
  character(len=256) :: path
  real(wp), allocatable :: theta(:), phi(:), s0(:,:), lam(:), stats(:,:), orbit(:,:,:), e_zz(:)
  real(wp), allocatable :: hus_theta(:), hus_phi(:), f(:), pr(:), mu(:)
  complex(wp), allocatable :: v(:,:), u(:,:)
  real(wp) :: j, eps
  integer :: sel(2), ic, n, i

  path = 'config_scan.nml'
  if (command_argument_count() >= 1) call get_command_argument(1, path)
  call read_scan_config(trim(path), cfg)
  j = 0.5_wp*cfg%n_q
  eps = cfg%k/cfg%n_q
  call make_dir(cfg%out_dir)

  ! 1. scan grid
  n = cfg%n_theta*cfg%n_phi
  allocate(theta(n), phi(n), s0(3, n))
  call equal_area_grid(cfg%n_theta, cfg%n_phi, theta, phi)
  do ic = 1, n
    s0(:, ic) = spherical_to_cartesian(theta(ic), phi(ic))
  end do

  c%h = cfg%h
  c%k = cfg%k
  c%n_kicks = cfg%n_lyap
  c%n_transient = 0
  c%save_every = 1
  c%mode = 'grid'
  c%n_renorm = cfg%n_renorm
  c%lyapunov = .true.
  call lyapunov(c, s0, lam)

  allocate(v(cfg%n_q + 1, cfg%n_q + 1), mu(cfg%n_q + 1))
  call diagonalise_jx(j, cfg%n_q + 1, v, mu)
  u = floquet_operator(j, cfg%h, eps, v, mu)
  allocate(stats(5, n))
  call spin_scan(u, j, theta, phi, cfg%t_scan, stats)

  ! 2. selection
  sel = 0
  do ic = 1, n
    if (lam(ic) < cfg%lam_regular) then
      if (sel(1) == 0) then
        sel(1) = ic
      else if (stats(5, ic) < stats(5, sel(1))) then
        sel(1) = ic
      end if
    end if
    if (lam(ic) > cfg%lam_chaotic .and. abs(stats(5, ic)/stats(4, ic) - 1.0_wp) < cfg%sat_tol) then
      if (sel(2) == 0) then
        sel(2) = ic
      else if (stats(5, ic) > stats(5, sel(2))) then
        sel(2) = ic
      end if
    end if
  end do
  if (sel(1) == 0) error stop prog//': no island candidate (Lambda < lam_regular)'
  if (sel(2) == 0) error stop prog//': no saturated sea candidate (Lambda > lam_chaotic)'
  call write_scan()

  ! 3. classical phase portrait and orbits of the selected states
  c%n_kicks = cfg%n_kicks
  c%n_theta = cfg%sec_n_theta
  c%n_phi = cfg%sec_n_phi
  c%out_dir = trim(cfg%out_dir)//'/section'
  block
    real(wp) :: th(cfg%sec_n_theta*cfg%sec_n_phi), ph(cfg%sec_n_theta*cfg%sec_n_phi)
    real(wp) :: sg(3, cfg%sec_n_theta*cfg%sec_n_phi)
    call equal_area_grid(cfg%sec_n_theta, cfg%sec_n_phi, th, ph)
    do ic = 1, size(th)
      sg(:, ic) = spherical_to_cartesian(th(ic), ph(ic))
    end do
    call evolve(c, sg, orbit)
    call write_results(c, orbit)
  end block
  c%mode = 'single'
  do i = 1, 2
    c%out_dir = trim(cfg%out_dir)//'/'//trim(names(i))//'/classical'
    call evolve(c, s0(:, sel(i):sel(i)), orbit)
    call write_results(c, orbit)
  end do

  ! 4. chain echo
  allocate(hus_theta(cfg%hus_n_theta*cfg%hus_n_phi), hus_phi(cfg%hus_n_theta*cfg%hus_n_phi))
  call equal_area_grid(cfg%hus_n_theta, cfg%hus_n_phi, hus_theta, hus_phi)
  e_zz = zz_energies(cfg%n_q, uniform_coupling(cfg%n_q))
  r%n_q = cfg%n_q
  r%eps = eps
  r%h = [cfg%h]
  r%s0 = -1
  r%pol = 'z'
  r%has_angles = .true.
  r%n_kicks = cfg%n_kicks
  r%n_theta = cfg%sec_n_theta
  r%n_phi = cfg%sec_n_phi
  r%t_r = cfg%t_r
  r%i_flip = cfg%i_flip
  r%flip_axis = cfg%flip_axis
  r%flip_angle = cfg%flip_angle
  r%hus_n_theta = cfg%hus_n_theta
  r%hus_n_phi = cfg%hus_n_phi
  r%out_dir = cfg%out_dir
  do i = 1, 2
    r%theta0 = theta(sel(i))
    r%phi0 = phi(sel(i))
    r%butterfly_effect = .true.
    call run_chain(r, cfg%h, r%theta0, r%phi0, e_zz, hus_theta, hus_phi, &
                   trim(cfg%out_dir)//'/'//trim(names(i))//'/chain', prog)
    r%butterfly_effect = .false.
    call run_chain(r, cfg%h, r%theta0, r%phi0, e_zz, hus_theta, hus_phi, &
                   trim(cfg%out_dir)//'/'//trim(names(i))//'/chain_no_flip', prog)
    if (cfg%n_angle >= 2) call flip_scan(i)
    if (cfg%tr_step > 0) call tr_scan(i)
  end do

contains

  ! scan.dat (all grid points), selected.dat, and <name>/scan_series.dat.
  subroutine write_scan()
    integer :: uu, ii, t
    real(wp) :: t_e

    open(newunit=uu, file=trim(cfg%out_dir)//'/scan.dat', status='replace', action='write')
    write(uu, '(a)') '# '//title()
    write(uu, '(a)') '# ic theta0 phi0 Z0 Lambda F(t_scan) PR(t_scan) F_late PR_early PR_late'
    do ii = 1, n
      write(uu, '(i0,9(1x,es23.15e3))') ii, theta(ii), phi(ii), cos(theta(ii)), lam(ii), stats(:, ii)
    end do
    close(uu)

    open(newunit=uu, file=trim(cfg%out_dir)//'/selected.dat', status='replace', action='write')
    write(uu, '(a)') '# '//title()
    write(uu, '(a)') '# t_E = ln(j)/(2 Lambda) (Ehrenfest estimate), t_scan/t_E; F_random = 1/(2j+1)'
    write(uu, '(a)') '# name ic theta0 phi0 Z0 Lambda F(t_scan) PR(t_scan) F_late PR_early PR_late t_E t_scan/t_E F_random'
    do ii = 1, 2
      t_e = huge(1.0_wp)
      if (lam(sel(ii)) > 0.0_wp) t_e = log(j)/(2.0_wp*lam(sel(ii)))
      write(uu, '(a,1x,i0,13(1x,es23.15e3))') trim(names(ii)), sel(ii), theta(sel(ii)), phi(sel(ii)), &
        cos(theta(sel(ii))), lam(sel(ii)), stats(:, sel(ii)), t_e, cfg%t_scan/t_e, 1.0_wp/(2.0_wp*j + 1.0_wp)
      write(*, '(a,a,a,f0.4,a,f0.4,a,es10.3,a,f0.4,a,f0.3)') prog//': ', trim(names(ii)), ' theta0 = ', &
        theta(sel(ii)), ', phi0 = ', phi(sel(ii)), ', Lambda = ', lam(sel(ii)), ', F_late = ', stats(3, sel(ii)), &
        ', PR_late = ', stats(5, sel(ii))
    end do
    close(uu)

    allocate(f(0:cfg%t_scan), pr(0:cfg%t_scan))
    do ii = 1, 2
      call make_dir(trim(cfg%out_dir)//'/'//trim(names(ii)))
      call spin_series(u, j, theta(sel(ii)), phi(sel(ii)), cfg%t_scan, f, pr)
      open(newunit=uu, file=trim(cfg%out_dir)//'/'//trim(names(ii))//'/scan_series.dat', status='replace', action='write')
      write(uu, '(a)') '# '//title()
      write(uu, '(a)') '# t F PR'
      do t = 0, cfg%t_scan
        write(uu, '(i0,2(1x,es23.15e3))') t, f(t), pr(t)
      end do
      close(uu)
    end do
  end subroutine write_scan

  ! Echo F(2t_r) for flip angles pi*l/(n_angle-1) about x, y, z (one forward leg, then one
  ! backward leg per angle), with the prediction cos^2 + sin^2 <sigma^a_i>^2 and the symmetric weight.
  subroutine flip_scan(i)
    integer, intent(in) :: i
    complex(wp), allocatable :: psi0(:), psi_tr(:), psi(:)
    real(wp) :: sig(3, cfg%n_q), ang
    integer :: uu, t, a, l

    psi0 = product_coherent_state(cfg%n_q, theta(sel(i)), phi(sel(i)))
    psi_tr = psi0
    do t = 1, cfg%t_r
      call forward_step(psi_tr, cfg%n_q, cfg%h, eps, e_zz)
    end do
    sig = site_spins(psi_tr, cfg%n_q)
    open(newunit=uu, file=trim(cfg%out_dir)//'/'//trim(names(i))//'/flip_scan.dat', status='replace', action='write')
    write(uu, '(a,i0,a,i0)') '# '//title()//', t_r = ', cfg%t_r, ', i_flip = ', cfg%i_flip
    write(uu, '(a)') '# axis flip_angle F(2t_r) F_predicted symmetric_weight(2t_r)'
    do a = 1, 3
      do l = 0, cfg%n_angle - 1
        ang = acos(-1.0_wp)*real(l, wp)/real(cfg%n_angle - 1, wp)
        psi = psi_tr
        call rotate_qubit(psi, cfg%i_flip, 'xyz'(a:a), ang)
        do t = 1, cfg%t_r
          call backward_step(psi, cfg%n_q, cfg%h, eps, e_zz)
        end do
        write(uu, '(a,4(1x,es23.15e3))') 'xyz'(a:a), ang, fidelity(psi0, psi), &
          cos(0.5_wp*ang)**2 + sin(0.5_wp*ang)**2*sig(a, cfg%i_flip)**2, sum(abs(dicke_amplitudes(psi, cfg%n_q))**2)
      end do
      write(uu, '(a)') ''
    end do
    close(uu)
    write(*, '(a)') prog//': flip scan -> '//trim(cfg%out_dir)//'/'//trim(names(i))//'/flip_scan.dat'
  end subroutine flip_scan

  ! Echo F(2t_r) with the configured flip for t_r = tr_step, 2 tr_step, ..., tr_max.
  subroutine tr_scan(i)
    integer, intent(in) :: i
    complex(wp), allocatable :: psi0(:), psi_fw(:), psi(:)
    real(wp) :: sig(3, cfg%n_q), ang
    integer :: uu, t, tr, a

    a = index('xyz', cfg%flip_axis)
    ang = cfg%flip_angle
    psi0 = product_coherent_state(cfg%n_q, theta(sel(i)), phi(sel(i)))
    psi_fw = psi0
    open(newunit=uu, file=trim(cfg%out_dir)//'/'//trim(names(i))//'/tr_scan.dat', status='replace', action='write')
    write(uu, '(a,i0,a,a,a,f0.6)') '# '//title()//', i_flip = ', cfg%i_flip, ', flip_axis = ', cfg%flip_axis, &
      ', flip_angle = ', ang
    write(uu, '(a)') '# t_r F(t_r) F(2t_r) F_predicted symmetric_weight(2t_r)'
    t = 0
    do tr = cfg%tr_step, cfg%tr_max, cfg%tr_step
      do while (t < tr)
        call forward_step(psi_fw, cfg%n_q, cfg%h, eps, e_zz)
        t = t + 1
      end do
      sig = site_spins(psi_fw, cfg%n_q)
      psi = psi_fw
      call rotate_qubit(psi, cfg%i_flip, cfg%flip_axis, ang)
      do t = 1, tr
        call backward_step(psi, cfg%n_q, cfg%h, eps, e_zz)
      end do
      t = tr
      write(uu, '(i0,4(1x,es23.15e3))') tr, fidelity(psi0, psi_fw), fidelity(psi0, psi), &
        cos(0.5_wp*ang)**2 + sin(0.5_wp*ang)**2*sig(a, cfg%i_flip)**2, sum(abs(dicke_amplitudes(psi, cfg%n_q))**2)
    end do
    close(uu)
    write(*, '(a)') prog//': t_r scan -> '//trim(cfg%out_dir)//'/'//trim(names(i))//'/tr_scan.dat'
  end subroutine tr_scan

  function title() result(s)
    character(len=:), allocatable :: s
    character(len=160) :: buf

    write(buf, '(a,f0.6,a,f0.6,a,i0,a,f0.6,a,i0)') 'kicked_top_scan: h = ', cfg%h, ', k = ', cfg%k, &
      ', n_q = ', cfg%n_q, ', eps = k/n_q = ', eps, ', t_scan = ', cfg%t_scan
    s = trim(buf)
  end function title

end program scan_main
