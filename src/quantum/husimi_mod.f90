! Husimi distribution Q(theta, phi) = |<theta,phi|psi>|**2 of a spin-j state psi,
! with |theta,phi> the spin coherent state of spin_mod. Q >= 0, and
! (2j+1)/(4 pi) * integral of Q over the sphere = 1 for a normalised psi.
module husimi_mod
  use kinds_mod, only: wp
  use spin_mod, only: coherent_state
  implicit none
  private
  public :: husimi, husimi_norm, participation_ratio, participation_ratio_state

contains

  ! Q at the points (theta(i), phi(i)), in parallel over the points (OpenMP).
  subroutine husimi(psi, j, theta, phi, q)
    complex(wp), intent(in) :: psi(:)
    real(wp), intent(in) :: j
    real(wp), intent(in) :: theta(:)
    real(wp), intent(in) :: phi(:)
    real(wp), intent(out) :: q(:)
    integer :: i

    !$omp parallel do schedule(static)
    do i = 1, size(theta)
      q(i) = abs(dot_product(coherent_state(j, size(psi), theta(i), phi(i)), psi))**2
    end do
    !$omp end parallel do
  end subroutine husimi

  ! (2j+1)/(4 pi) * integral of Q, evaluated on an equal-area grid of n points
  ! (each of area 4 pi/n): (2j+1) * mean(Q). Should be 1.
  pure function husimi_norm(q, j) result(norm)
    real(wp), intent(in) :: q(:)
    real(wp), intent(in) :: j
    real(wp) :: norm

    norm = (2.0_wp*j + 1.0_wp)*sum(q)/real(size(q), wp)
  end function husimi_norm

  ! Participation ratio (2j+1) (integral of Q)**2 / (4 pi * integral of Q**2) on an
  ! equal-area grid (each point carries the area 4 pi/n): (2j+1) mean(Q)**2/mean(Q**2).
  ! (4j+1)/(2j+1) ~ 2 for a coherent state, j+1 on average for a random state, 2j+1 for Q uniform.
  pure function participation_ratio(q, j) result(pr)
    real(wp), intent(in) :: q(:)
    real(wp), intent(in) :: j
    real(wp) :: pr

    pr = (2.0_wp*j + 1.0_wp)*(sum(q)/real(size(q), wp))**2/(sum(q**2)/real(size(q), wp))
  end function participation_ratio

  ! The same participation ratio from the amplitudes psi_m (ordered as in spin_mod), without a grid.
  ! |theta,phi> (x) |theta,phi> is the spin-2j coherent state, so Q**2 = |<theta,phi|_{2j} chi|**2 with
  ! chi_M = sum_{m1+m2=M} c(m1,m2) psi_m1 psi_m2, c = sqrt(binom(2j,j-m1) binom(2j,j-m2)/binom(4j,2j-M)),
  ! the projection of psi (x) psi on spin 2j. The resolution of the identity for spin 2j gives
  ! integral of Q**2 = 4 pi/(4j+1) |chi|**2, and integral of Q = 4 pi/(2j+1) |psi|**2, so
  ! pr = (4j+1) |psi|**4/((2j+1) |chi|**2).
  pure function participation_ratio_state(psi, j) result(pr)
    complex(wp), intent(in) :: psi(:)
    real(wp), intent(in) :: j
    real(wp) :: pr
    real(wp) :: lb(0:size(psi) - 1), lb2(0:2*(size(psi) - 1))
    complex(wp) :: chi(0:2*(size(psi) - 1))
    integer :: n, r1, r2

    n = size(psi) - 1                      ! 2j; index i = r + 1 holds m = j - r
    lb = [(log_binom(n, r1), r1 = 0, n)]
    lb2 = [(log_binom(2*n, r1), r1 = 0, 2*n)]
    chi = (0.0_wp, 0.0_wp)
    do r1 = 0, n
      do r2 = 0, n
        chi(r1 + r2) = chi(r1 + r2) + exp(0.5_wp*(lb(r1) + lb(r2) - lb2(r1 + r2)))*psi(r1 + 1)*psi(r2 + 1)
      end do
    end do
    pr = (4.0_wp*j + 1.0_wp)*sum(abs(psi)**2)**2/((2.0_wp*j + 1.0_wp)*sum(abs(chi)**2))
  end function participation_ratio_state

  ! log of the binomial coefficient binom(n, r).
  pure function log_binom(n, r) result(lb)
    integer, intent(in) :: n, r
    real(wp) :: lb

    lb = log_gamma(real(n + 1, wp)) - log_gamma(real(r + 1, wp)) - log_gamma(real(n - r + 1, wp))
  end function log_binom

end module husimi_mod
