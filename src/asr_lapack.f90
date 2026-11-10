!> Explicit interfaces for the LAPACK routines the code calls (standard LAPACK 3 signatures, double precision).
!> With them, every call is checked at compile time (types, number of arguments), and compilers that reject implicit interfaces
!> (gfortran under fpm: -Werror=implicit-interface) build the code unchanged.  The routines come from any LAPACK library (MKL,
!> OpenBLAS, reference LAPACK); nothing here depends on the vendor.  Array dummies are assumed-size, so callers pass arrays of
!> any rank and contiguous sections by sequence association, as with the implicit interface.
module asr_lapack
  implicit none
  private
  public :: dgesv, zgesv, zgetrf, zgetri, zgetrs, dsterf, dsyev, dsygv, zheev, zhegv, zgeev

  integer, parameter :: dp = kind(1.0d0)

  interface
    subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
      import :: dp
      integer, intent(in) :: n, nrhs, lda, ldb
      real(dp), intent(inout) :: a(lda, *), b(ldb, *)
      integer, intent(out) :: ipiv(*), info
    end subroutine dgesv

    subroutine zgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
      import :: dp
      integer, intent(in) :: n, nrhs, lda, ldb
      complex(dp), intent(inout) :: a(lda, *), b(ldb, *)
      integer, intent(out) :: ipiv(*), info
    end subroutine zgesv

    subroutine zgetrf(m, n, a, lda, ipiv, info)
      import :: dp
      integer, intent(in) :: m, n, lda
      complex(dp), intent(inout) :: a(lda, *)
      integer, intent(out) :: ipiv(*), info
    end subroutine zgetrf

    subroutine zgetri(n, a, lda, ipiv, work, lwork, info)
      import :: dp
      integer, intent(in) :: n, lda, lwork
      complex(dp), intent(inout) :: a(lda, *)
      integer, intent(in) :: ipiv(*)
      complex(dp), intent(out) :: work(*)
      integer, intent(out) :: info
    end subroutine zgetri

    subroutine zgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
      import :: dp
      character(len=1), intent(in) :: trans
      integer, intent(in) :: n, nrhs, lda, ldb
      complex(dp), intent(in) :: a(lda, *)
      integer, intent(in) :: ipiv(*)
      complex(dp), intent(inout) :: b(ldb, *)
      integer, intent(out) :: info
    end subroutine zgetrs

    subroutine dsterf(n, d, e, info)
      import :: dp
      integer, intent(in) :: n
      real(dp), intent(inout) :: d(*), e(*)
      integer, intent(out) :: info
    end subroutine dsterf

    subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
      import :: dp
      character(len=1), intent(in) :: jobz, uplo
      integer, intent(in) :: n, lda, lwork
      real(dp), intent(inout) :: a(lda, *)
      real(dp), intent(out) :: w(*), work(*)
      integer, intent(out) :: info
    end subroutine dsyev

    subroutine dsygv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, info)
      import :: dp
      integer, intent(in) :: itype, n, lda, ldb, lwork
      character(len=1), intent(in) :: jobz, uplo
      real(dp), intent(inout) :: a(lda, *), b(ldb, *)
      real(dp), intent(out) :: w(*), work(*)
      integer, intent(out) :: info
    end subroutine dsygv

    subroutine zheev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info)
      import :: dp
      character(len=1), intent(in) :: jobz, uplo
      integer, intent(in) :: n, lda, lwork
      complex(dp), intent(inout) :: a(lda, *)
      real(dp), intent(out) :: w(*), rwork(*)
      complex(dp), intent(out) :: work(*)
      integer, intent(out) :: info
    end subroutine zheev

    subroutine zhegv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, info)
      import :: dp
      integer, intent(in) :: itype, n, lda, ldb, lwork
      character(len=1), intent(in) :: jobz, uplo
      complex(dp), intent(inout) :: a(lda, *), b(ldb, *)
      real(dp), intent(out) :: w(*), rwork(*)
      complex(dp), intent(out) :: work(*)
      integer, intent(out) :: info
    end subroutine zhegv

    subroutine zgeev(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info)
      import :: dp
      character(len=1), intent(in) :: jobvl, jobvr
      integer, intent(in) :: n, lda, ldvl, ldvr, lwork
      complex(dp), intent(inout) :: a(lda, *)
      complex(dp), intent(out) :: w(*), vl(ldvl, *), vr(ldvr, *), work(*)
      real(dp), intent(out) :: rwork(*)
      integer, intent(out) :: info
    end subroutine zgeev
  end interface
end module asr_lapack
