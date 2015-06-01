      module maindef

      use aAdjKeep
      use common

      implicit none

      integer so, sf, count, icntu, icntp, ni,nj,nk, kk,
     &     icntpp, irst, nrest, solf, exist,
     &     restf, StepCount, i, j, k

      integer, dimension(:), allocatable ::
     &     rowp,  colp,
     &     rowpp, colpp,
     &     rowu,  colu

      character*30 fname

      integer lenseg, is, itask
      integer ierr, knotf, lcp

      logical:: exts, debug=.true.

                                ! MEAN COMPUTATION VARIABLES
      integer
     &     NPTSu, NPTSv, NPTSw,NPTSu1, NPTSv1, NPTSw1,
     &     PFLAGu, PFLAGv, PFLAGw,
     &     avg, iavg, iavg0, iavgfreq!, varcount

      real*8, allocatable, dimension(:)   :: PTSu, PTSv, PTSw
      real*4, allocatable, dimension(:,:,:,:) :: SOL_AR_AVG
      real*4, allocatable, dimension(:,:,:,:) :: STRAIN_AVG	  

      real*8 divfac, T_eddy

      end module
