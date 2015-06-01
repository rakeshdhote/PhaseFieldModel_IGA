c------------------------------------------------------------------------c
c                                                                        c
c        Module for storing arrays and allocation routines               c
c                                                                        c
c                                                                        c
c------------------------------------------------------------------------c
      
      module aAdjKeep

      implicit none

c...  Knot vectors in u and v directions
      real*8, allocatable :: U_KNOTu(:),V_KNOTu(:),W_KNOTu(:)
      real*8, allocatable :: U_KNOTp(:),V_KNOTp(:),W_KNOTp(:)

c...  Control Net
      real*8, allocatable :: B_NETu(:,:,:,:)
      real*8, allocatable :: B_NETini(:,:,:,:)
      real*8, allocatable :: B_NETp(:,:,:,:)

c...  The Load applied to each edge and element face
      real*8, allocatable :: LD_FACE(:,:), LD_EL(:,:)

c...  The right hand side load vector G and left hand stiffness matrix K
      real*8, allocatable :: RHSGu(:,:), LHSK(:,:)

c...  IEN matches element number a local node number with global node number
      integer, allocatable :: IENu(:,:), IENp(:,:)

c...  INN relate global node number to the (i,j) "NURBS coordinates"
      integer, allocatable :: INNu(:,:), INNp(:,:)

c...  IPER contains all master/slave relationships.
      integer, allocatable, dimension(:) ::
     &     IPERu, IPERp, EL_CON, FACE_CON

c...  Matrix relating face number and local node to global node number
      integer, allocatable :: FACE_IENu(:,:), FACE_IENp(:,:)

c...  Face orientation matrix
      integer, allocatable :: FACE_OR(:)

c...  Boundary condition indicator for global nodes and edges respectively
      integer, allocatable :: IBC(:,:), IBC_FACE(:,:)

c...  Prescribed displacements of Dirichlet nodes.
      real*8, allocatable :: DIR_BC(:,:)

c...  Array containing Young's modulus in each element
      real*8, allocatable :: DENS(:)

c...  Array containing Poisson's ration in each element
      real*8, allocatable :: VISC(:)

c...  Array containing Energy at each time
      real*8, allocatable :: VENER(:)

c...  Array containing time step at each time
      real*8, allocatable :: VTIME(:)

c...  Array containing results
      real*8, allocatable :: VDELT(:)
	  
c...  Array containing results
      real*8, allocatable :: VTEMP(:)	  
      real*8, allocatable :: VSxx(:)
      real*8, allocatable :: VSxy(:)
      real*8, allocatable :: VSyy(:)

      real*8, allocatable :: VExx(:)
      real*8, allocatable :: VExy(:)
      real*8, allocatable :: VEyy(:)

c...  Array containing velocity
      real*8, allocatable :: xor(:)

c...  Solution vector (primitive variables)
      real*8, allocatable :: acg(:,:), acgold(:,:) ! velocities
      real*8, allocatable :: ygBE(:,:)     ! backward euler sol vector
      real*8, allocatable :: yg(:,:), ygold(:,:) ! solution vector

c...  Restriction operator
      real*8, dimension(:,:), allocatable ::
     &     TuMu, TuSu,
     &     TvMu, TvSu,
     &     TwMu, TwSu,
     &     TuMp, TuSp,
     &     TvMp, TvSp,
     &     TwMp, TwSp

      end module

