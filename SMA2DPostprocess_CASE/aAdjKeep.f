
     
c------------------------------------------------------------------------c
c                                                                        c
c        Module for storing arrays and allocation routines               c
c                                                                        c
c                                                                        c
c------------------------------------------------------------------------c
      
      module aAdjKeep

c...  Knot vectors in u and v directions
      real*8, allocatable :: U_KNOT(:), V_KNOT(:),W_KNOT(:)

c...  Control Net
      real*8, allocatable :: B_NET(:,:,:,:)

c...  Projective Coordinate control net
!      real*8, allocatable :: BW(:,:,:,:)

c...  The Load applied to each edge and element face
      real*8, allocatable :: LD_FACE(:,:), LD_EL(:,:)

c...  The right hand side load vector G and left hand stiffness matrix K
      real*8, allocatable :: RHSG(:,:), LHSK(:,:), LHSKV(:)

c...  IEN matches element number a local node number with global node number
      integer, allocatable :: IEN(:,:),  IEN_GLOB_LIN(:,:)

c...  INN relate global node number to the (i,j) "NURBS coordinates"
      integer, allocatable :: INN(:,:)

c...  IPER contains all master/slave relationships.
      integer, allocatable :: IPER(:)

c...  Matrix relating face number and local node to global node number
      integer, allocatable :: FACE_IEN(:,:)

c...  Face orientation matrix
      integer, allocatable :: FACE_OR(:)

c...  Boundary condition indicator for global nodes and edges respectively
      integer, allocatable :: IBC(:,:), IBC_FACE(:,:)

c...  Prescribed displacements of Dirichlet nodes.
      real*8, allocatable :: DIR_BC(:,:)

c...  Array containing Young's modulus in each element
!      real*8, allocatable :: EY(:)

c...  Array containing Poisson's ration in each element
!      real*8, allocatable :: NU(:)

c...  Solution vector
      real*8, allocatable :: yg(:,:), x_glob_lin(:,:),
     &  y_glob_lin(:,:), glob_mat_dat(:)

      real*8, allocatable :: U_KNOTu(:),V_KNOTu(:),W_KNOTu(:)	 

      end module

c-------------------------------------------------------------------------c







