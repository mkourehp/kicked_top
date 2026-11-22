! One-period (Floquet) operator of the quantum kicked top,
!   U = exp(-i eps J_z**2) exp(-i h J_x)      (first the kick, then the torsion),
! the quantum counterpart of the classical map, Sec. 9 [sec:quantum].
module floquet_mod
  use kinds_mod, only: wp
  use spin_mod, only: m_values, jx_matrix
  implicit none
  private
  public :: diagonalise_jx, exp_jx, floquet_operator

contains

  ! Eigen-decomposition J_x = V diag(mu) V^dagger with LAPACK zheev.
  ! The eigenvalues mu are returned in ascending order (-j, -j+1, ..., j).
  subroutine diagonalise_jx(j, n, v, mu)
    real(wp), intent(in) :: j
    integer, intent(in) :: n
    complex(wp), intent(out) :: v(n, n)
    real(wp), intent(out) :: mu(n)
    complex(wp), allocatable :: work(:)
    complex(wp) :: work_query(1)
    real(wp) :: rwork(max(1, 3*n - 2))
    integer :: info
    external :: zheev

    v = jx_matrix(j, n)
    call zheev('V', 'U', n, v, n, mu, work_query, -1, rwork, info)
    allocate(work(max(1, int(real(work_query(1), wp)))))
    call zheev('V', 'U', n, v, n, mu, work, size(work), rwork, info)
    if (info /= 0) error stop 'floquet_mod: zheev failed'
  end subroutine diagonalise_jx

  ! exp(-i alpha J_x) = V diag(exp(-i alpha mu)) V^dagger, a rotation about x by alpha.
  pure function exp_jx(v, mu, alpha) result(r)
    complex(wp), intent(in) :: v(:,:)
    real(wp), intent(in) :: mu(:)
    real(wp), intent(in) :: alpha
    complex(wp) :: r(size(mu), size(mu))
    complex(wp) :: w(size(mu), size(mu))
    integer :: l

    do l = 1, size(mu)
      w(:, l) = v(:, l)*cmplx(cos(alpha*mu(l)), -sin(alpha*mu(l)), wp)
    end do
    r = matmul(w, conjg(transpose(v)))
  end function exp_jx

  ! U = T R with the kick R = exp(-i h J_x) and the diagonal torsion
  ! T = diag(exp(-i eps m**2)). v, mu: eigen-decomposition of J_x.
  pure function floquet_operator(j, h, eps, v, mu) result(u)
    real(wp), intent(in) :: j
    real(wp), intent(in) :: h
    real(wp), intent(in) :: eps
    complex(wp), intent(in) :: v(:,:)
    real(wp), intent(in) :: mu(:)
    complex(wp) :: u(size(mu), size(mu))
    real(wp) :: m(size(mu)), phase
    integer :: i

    m = m_values(j, size(mu))
    u = exp_jx(v, mu, h)
    do i = 1, size(mu)
      phase = eps*m(i)**2
      u(i, :) = u(i, :)*cmplx(cos(phase), -sin(phase), wp)
    end do
  end function floquet_operator

end module floquet_mod
