! Floquet step of the qubit chain, applied to the state vector only (no
! Hamiltonian or 2**n_q x 2**n_q matrix is built). Basis as in chain_state_mod.
!
! The chain realises the kicked top with J = (1/2) sum_i sigma_i, j = n_q/2:
!   H(t) = eps H_zz + h J_x sum_n delta(t-n),   H_zz = sum_{i,l} V_il sigma^z_i sigma^z_l,
! with i and l running over all qubits (i = l included) and h J_x = (h/2) sum_i sigma^x_i.
! The couplings act continuously; at integer times a pulse rotates every qubit about x by h.
! For V_il = 1/4, H_zz = (1/2 sum_i sigma^z_i)**2 = J_z**2. Sampled just before
! each kick, one period is
!   U = exp(-i eps H_zz) exp(-i h J_x),
! the kicked-top Floquet operator of floquet_mod, lecture notes eq. (39) [eq:floquet].
!
! Time reversal (as in spin_chain): the sign of the couplings cannot be flipped,
! so the coefficient eps -> 2 pi - eps, h -> -h and the two factors are
! swapped. The eigenvalues m**2 of H_zz (m all integer or all half-integer) give
! exp(-i 2 pi H_zz) = global phase, so backward_step = U^dagger up to that phase.
!
! Static impurity (as imp_z, imp_x in spin_chain): drawn once per run, the same in
! the forward and the backward leg (not reversed), so with it the echo is not exact.
!   coupling: coefficient of sigma^z_i sigma^z_l is eps V_il + xi_il, xi_il = xi_li
!             uniform in [-deltaz/2, deltaz/2] for i < l, xi_ii = 0;
!   kick:     qubit k is rotated about x by h + eta_k, eta_k uniform in [-deltax/2, deltax/2].
! With e_imp the diagonal of sum_{i,l} xi_il sigma^z_i sigma^z_l (zz_energies of xi):
!   forward:  exp(-i (eps H_zz + H_xi)) prod_k exp(-i (h + eta_k)/2 sigma^x_k),
!   backward: prod_k exp(-i (-h + eta_k)/2 sigma^x_k) exp(-i ((2 pi - eps) H_zz + H_xi)).
module chain_floquet_mod
  use kinds_mod, only: wp, two_pi
  implicit none
  private
  public :: uniform_coupling, static_impurity, zz_energies, zz_coupling, rotate_x, rotate_qubit, &
            forward_step, backward_step

contains

  ! V_il = 1/4 for all i, l: H_zz = J_z**2 (kicked-top limit).
  pure function uniform_coupling(n_q) result(v)
    integer, intent(in) :: n_q
    real(wp) :: v(n_q, n_q)

    v = 0.25_wp
  end function uniform_coupling

  ! Static impurity, drawn once from the seed: xi(i, l) = xi(l, i) uniform in
  ! [-deltaz/2, deltaz/2] for i < l, xi(i, i) = 0; eta(k) uniform in [-deltax/2, deltax/2].
  ! deltaz = deltax = 0 gives xi = eta = 0 (clean chain).
  subroutine static_impurity(n_q, deltaz, deltax, seed, xi, eta)
    integer, intent(in) :: n_q
    real(wp), intent(in) :: deltaz
    real(wp), intent(in) :: deltax
    integer, intent(in) :: seed
    real(wp), intent(out) :: xi(n_q, n_q)
    real(wp), intent(out) :: eta(n_q)
    integer, allocatable :: seed_arr(:)
    integer :: ss, k, i, l
    real(wp) :: r

    call random_seed(size=ss)
    seed_arr = [(seed + 37*k, k = 1, ss)]
    call random_seed(put=seed_arr)
    xi = 0.0_wp
    do i = 1, n_q
      do l = i + 1, n_q
        call random_number(r)
        xi(i, l) = deltaz*(r - 0.5_wp)
        xi(l, i) = xi(i, l)
      end do
    end do
    do k = 1, n_q
      call random_number(r)
      eta(k) = deltax*(r - 0.5_wp)
    end do
  end subroutine static_impurity

  ! Diagonal of H_zz: e(b) = sum_{i,l} V_il z_i z_l, z_i = +1 if bit i-1 of b
  ! is 1 (qubit i up) and -1 otherwise. Computed once per run.
  function zz_energies(n_q, v) result(e)
    integer, intent(in) :: n_q
    real(wp), intent(in) :: v(n_q, n_q)
    real(wp), allocatable :: e(:)
    real(wp) :: z(n_q)
    integer :: b, i

    allocate(e(0:2**n_q - 1))
    !$omp parallel do private(z, i) schedule(static)
    do b = 0, 2**n_q - 1
      do i = 1, n_q
        z(i) = merge(1.0_wp, -1.0_wp, btest(b, i - 1))
      end do
      e(b) = dot_product(z, matmul(v, z))
    end do
    !$omp end parallel do
  end function zz_energies

  ! psi -> exp(-i (t H_zz + H_xi)) psi, with e and e_imp the diagonals of H_zz and
  ! of the static impurity H_xi from zz_energies.
  subroutine zz_coupling(psi, e, t, e_imp)
    complex(wp), intent(inout) :: psi(0:)
    real(wp), intent(in) :: e(0:)
    real(wp), intent(in) :: t
    real(wp), intent(in) :: e_imp(0:)
    real(wp) :: w
    integer :: b

    !$omp parallel do private(w) schedule(static)
    do b = 0, size(psi) - 1
      w = t*e(b) + e_imp(b)
      psi(b) = psi(b)*cmplx(cos(w), -sin(w), wp)
    end do
    !$omp end parallel do
  end subroutine zz_coupling

  ! psi -> prod_i exp(-i (alpha + eta(i))/2 sigma^x_i) psi, the x case of rotate_qubit
  ! applied to every qubit. For eta = 0 this is exp(-i alpha J_x), the rotation of the
  ! spin about x by alpha, lecture notes eq. (17) [eq:rotation-x].
  subroutine rotate_x(psi, n_q, alpha, eta)
    complex(wp), intent(inout) :: psi(0:)
    integer, intent(in) :: n_q
    real(wp), intent(in) :: alpha
    real(wp), intent(in) :: eta(:)
    integer :: i

    do i = 1, n_q
      call rotate_qubit(psi, i, 'x', alpha + eta(i))
    end do
  end subroutine rotate_x

  ! Rotation of qubit i (1 .. n_q) about a = axis ('x', 'y' or 'z') by theta:
  ! psi -> exp(-i theta/2 sigma^a_i) psi. theta = pi is the spin flip sigma^a_i
  ! (up to the global phase -i). With c = cos(theta/2), s = sin(theta/2), in the
  ! basis (up, down): x: [[c, -i s], [-i s, c]], y: [[c, -s], [s, c]],
  ! z: diag(exp(-i theta/2), exp(+i theta/2)).
  subroutine rotate_qubit(psi, i, axis, theta)
    complex(wp), intent(inout) :: psi(0:)
    integer, intent(in) :: i
    character(len=1), intent(in) :: axis
    real(wp), intent(in) :: theta
    complex(wp), parameter :: im = (0.0_wp, 1.0_wp)
    real(wp) :: c, s
    complex(wp) :: a_up, a_dn
    integer :: b, mask

    c = cos(theta/2.0_wp)
    s = sin(theta/2.0_wp)
    mask = ishft(1, i - 1)
    !$omp parallel do private(a_up, a_dn) schedule(static)
    do b = 0, size(psi) - 1
      if (iand(b, mask) /= 0) cycle
      a_dn = psi(b)
      a_up = psi(ior(b, mask))
      select case (axis)
      case ('x')
        psi(ior(b, mask)) = c*a_up - im*s*a_dn
        psi(b) = -im*s*a_up + c*a_dn
      case ('y')
        psi(ior(b, mask)) = c*a_up - s*a_dn
        psi(b) = s*a_up + c*a_dn
      case default
        psi(ior(b, mask)) = cmplx(c, -s, wp)*a_up
        psi(b) = cmplx(c, s, wp)*a_dn
      end select
    end do
    !$omp end parallel do
  end subroutine rotate_qubit

  ! One period: U = exp(-i eps H_zz) exp(-i h J_x)
  ! (first the kick, then the couplings for one period), with the static impurity e_imp, eta.
  subroutine forward_step(psi, n_q, h, eps, e, e_imp, eta)
    complex(wp), intent(inout) :: psi(0:)
    integer, intent(in) :: n_q
    real(wp), intent(in) :: h
    real(wp), intent(in) :: eps
    real(wp), intent(in) :: e(0:)
    real(wp), intent(in) :: e_imp(0:)
    real(wp), intent(in) :: eta(:)

    call rotate_x(psi, n_q, h, eta)
    call zz_coupling(psi, e, eps, e_imp)
  end subroutine forward_step

  ! Backward period: exp(+i h J_x) exp(-i (2 pi - eps) H_zz) = U^dagger (up to a global phase)
  ! without impurity; the static impurity e_imp, eta enters unchanged (not reversed).
  subroutine backward_step(psi, n_q, h, eps, e, e_imp, eta)
    complex(wp), intent(inout) :: psi(0:)
    integer, intent(in) :: n_q
    real(wp), intent(in) :: h
    real(wp), intent(in) :: eps
    real(wp), intent(in) :: e(0:)
    real(wp), intent(in) :: e_imp(0:)
    real(wp), intent(in) :: eta(:)

    call zz_coupling(psi, e, two_pi - eps, e_imp)
    call rotate_x(psi, n_q, -h, eta)
  end subroutine backward_step

end module chain_floquet_mod
