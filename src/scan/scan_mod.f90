! Spin-j part of the initial-condition scan: each coherent state |theta0, phi0> is
! evolved by the Floquet matrix U = exp(-i eps J_z^2) exp(-i h J_x) for t_scan kicks,
! recording the fidelity F(t) = |<psi(0)|psi(t)>|^2 and the Husimi participation
! ratio PR(t) = (2j+1) (int Q)^2/(4 pi int Q^2) (husimi_mod).
module scan_mod
  use kinds_mod, only: wp
  use spin_mod, only: coherent_state
  use husimi_mod, only: participation_ratio_state
  implicit none
  private
  public :: spin_series, spin_scan

contains

  ! F(t) and PR(t), t = 0 .. t_scan, of one initial coherent state.
  subroutine spin_series(u, j, theta0, phi0, t_scan, f, pr)
    complex(wp), intent(in) :: u(:,:)
    real(wp), intent(in) :: j, theta0, phi0
    integer, intent(in) :: t_scan
    real(wp), intent(out) :: f(0:t_scan), pr(0:t_scan)
    complex(wp) :: psi0(size(u, 1)), psi(size(u, 1))
    integer :: t

    psi0 = coherent_state(j, size(u, 1), theta0, phi0)
    psi = psi0
    do t = 0, t_scan
      if (t > 0) psi = matmul(u, psi)
      f(t) = abs(dot_product(psi0, psi))**2
      pr(t) = participation_ratio_state(psi, j)
    end do
  end subroutine spin_series

  ! For every initial state: F(t_scan), PR(t_scan), and the means of F and PR over the
  ! early window t_scan/2 < t <= 3 t_scan/4 and the late window 3 t_scan/4 < t <= t_scan.
  ! stats(:, ic) = [F(t_scan), PR(t_scan), F_late, PR_early, PR_late].
  subroutine spin_scan(u, j, theta0, phi0, t_scan, stats)
    complex(wp), intent(in) :: u(:,:)
    real(wp), intent(in) :: j, theta0(:), phi0(:)
    integer, intent(in) :: t_scan
    real(wp), intent(out) :: stats(5, size(theta0))
    real(wp) :: f(0:t_scan), pr(0:t_scan)
    integer :: ic, t1, t2

    t1 = t_scan/2
    t2 = (3*t_scan)/4
    !$omp parallel do private(f, pr) schedule(dynamic)
    do ic = 1, size(theta0)
      call spin_series(u, j, theta0(ic), phi0(ic), t_scan, f, pr)
      stats(:, ic) = [f(t_scan), pr(t_scan), sum(f(t2 + 1:t_scan))/real(t_scan - t2, wp), &
                      sum(pr(t1 + 1:t2))/real(t2 - t1, wp), sum(pr(t2 + 1:t_scan))/real(t_scan - t2, wp)]
    end do
    !$omp end parallel do
  end subroutine spin_scan

end module scan_mod
