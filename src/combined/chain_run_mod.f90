! Time evolution of the qubit chain for one rotation angle h, shared by
! kicked_top_combined (make combined) and kicked_top_spin_chain (make spin_chain):
! t_r forward steps U, then (with butterfly_effect) a rotation of qubit i_flip
! about a = flip_axis by flip_angle (pi = flip sigma^a), then t_r backward
! steps (coupling 2 pi - eps, kick -h).
! Writes the files of combined_io_mod in dir.
module chain_run_mod
  use kinds_mod, only: wp
  use io_utils_mod, only: make_dir
  use husimi_mod, only: husimi, husimi_norm, participation_ratio
  use combined_config_mod, only: rconfig_t
  use combined_io_mod, only: write_chain, write_husimi, write_summary
  use chain_state_mod, only: product_coherent_state, dicke_amplitudes
  use chain_floquet_mod, only: rotate_qubit, forward_step, backward_step
  use chain_observables_mod, only: fidelity, site_spins, total_spin
  implicit none
  private
  public :: run_chain

contains

  ! Forward evolution for t_r steps, then (after the optional flip) backward evolution for t_r steps.
  ! obs(:, t_r) and the Husimi snapshot at t_r are taken before the flip.
  subroutine run_chain(cfg, h, theta0, phi0, e_zz, hus_theta, hus_phi, dir, prog)
    type(rconfig_t), intent(in) :: cfg
    real(wp), intent(in) :: h
    real(wp), intent(in) :: theta0, phi0      ! direction of the initial polarisation
    real(wp), intent(in) :: e_zz(0:)          ! diagonal of H_zz
    real(wp), intent(in) :: hus_theta(:), hus_phi(:)
    character(len=*), intent(in) :: dir
    character(len=*), intent(in) :: prog      ! program name for the progress line
    complex(wp), allocatable :: psi0(:), psi(:)
    real(wp), allocatable :: obs(:,:), q(:,:)
    real(wp) :: sig(3, cfg%n_q), j, sym_dev, sym_end, c2, s2, sig_a
    integer :: steps(5), t, l, a

    j = 0.5_wp*cfg%n_q
    steps = [0, cfg%t_r/2, cfg%t_r, 3*cfg%t_r/2, 2*cfg%t_r]
    allocate(obs(9 + 3*cfg%n_q, 0:2*cfg%t_r), q(size(hus_theta), size(steps)))
    psi0 = product_coherent_state(cfg%n_q, theta0, phi0)
    psi = psi0
    sym_dev = 0.0_wp

    do t = 0, 2*cfg%t_r
      if (t >= 1 .and. t <= cfg%t_r) call forward_step(psi, cfg%n_q, h, cfg%eps, e_zz)
      if (t == cfg%t_r + 1 .and. cfg%butterfly_effect) call rotate_qubit(psi, cfg%i_flip, cfg%flip_axis, cfg%flip_angle)
      if (t > cfg%t_r) call backward_step(psi, cfg%n_q, h, cfg%eps, e_zz)

      sig = site_spins(psi, cfg%n_q)
      obs(1, t) = fidelity(psi0, psi)
      obs(2, t) = sum(abs(psi)**2)
      obs(3:9, t) = total_spin(psi, cfg%n_q)
      obs(10:, t) = [sig(1, :), sig(2, :), sig(3, :)]

      do l = 1, size(steps)
        if (t == steps(l)) then
          sym_dev = max(sym_dev, abs(1.0_wp - sum(abs(dicke_amplitudes(psi, cfg%n_q))**2)))
          call husimi(dicke_amplitudes(psi, cfg%n_q), j, hus_theta, hus_phi, q(:, l))
        end if
      end do
    end do

    ! With R = exp(-i flip_angle sigma^a_i/2) = cos(flip_angle/2) - i sin(flip_angle/2) sigma^a_i and exact
    ! time reversal, F(2t_r) = |<psi(t_r)|R|psi(t_r)>|**2 = c2 + s2 <sigma^a_i>**2 at t_r (any state), and
    ! for a symmetric psi(t_r) the symmetric weight after the flip is c2 + s2 <J_a**2>/j**2.
    sym_end = sum(abs(dicke_amplitudes(psi, cfg%n_q))**2)
    a = index('xyz', cfg%flip_axis)
    c2 = cos(0.5_wp*cfg%flip_angle)**2
    s2 = sin(0.5_wp*cfg%flip_angle)**2
    if (.not. cfg%butterfly_effect) then
      c2 = 1.0_wp
      s2 = 0.0_wp
    end if
    sig_a = obs(9 + (a - 1)*cfg%n_q + cfg%i_flip, cfg%t_r)

    call make_dir(dir)
    call write_chain(cfg, h, dir, obs)
    call write_husimi(cfg, h, dir, steps, hus_theta, hus_phi, q)
    call write_summary(cfg, h, dir, &
      [character(len=24) :: 'F(t_r)', 'F(2t_r)', 'max|norm-1|', 'max|J^2-j(j+1)|', &
       'max(1-symmetric_weight)', 'husimi_norm(0)', 'husimi_norm(t_r/2)', 'husimi_norm(t_r)', &
       'husimi_norm(3t_r/2)', 'husimi_norm(2t_r)', 'symmetric_weight(2t_r)', &
       'F(2t_r)_predicted', 'sym_weight_predicted', 'participation(0)', 'participation(t_r/2)', &
       'participation(t_r)', 'participation(3t_r/2)', 'participation(2t_r)'], &
      [obs(1, cfg%t_r), obs(1, 2*cfg%t_r), maxval(abs(obs(2, :) - 1.0_wp)), &
       maxval(abs(obs(9, :) - j*(j + 1.0_wp))), sym_dev, &
       (husimi_norm(q(:, l), j), l = 1, size(steps)), sym_end, c2 + s2*sig_a**2, &
       c2 + s2*obs(5 + a, cfg%t_r)/j**2, (participation_ratio(q(:, l), j), l = 1, size(steps))])

    write(*, '(a,f0.4,a,f0.4,a,es10.3,a)') prog//': h = ', h, ', eps = ', cfg%eps, &
      ', F(2t_r) - 1 = ', obs(1, 2*cfg%t_r) - 1.0_wp, ' -> '//dir
  end subroutine run_chain

end module chain_run_mod
