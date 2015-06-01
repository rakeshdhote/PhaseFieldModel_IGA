      
      subroutine input_3D_u
      
c---------------------------------------------------------------------c
c     c
c     Subroutine to open file and read vertex, edge                   c
c     and element numbers. Initialze parameters of interest.          c
c     c
c     3D Version                                                      c
c     c
c---------------------------------------------------------------------c

      use aAdjKeep
!      use mpi
      use common

      implicit none
      
      integer :: i, j, k,l, meshf, ink, ios,IBC1(NNODZu,NDOF)
      character*30 fname
      character*10 cname
      real*8 jnk,DIRBC1(NNODZu,NDOF)

      meshf = 11
      
      fname = trim('meshu' // cname(myid+1)) // '.dat'
!      fname = 'meshu.1.dat'
c...  Read in preliminary information
c     
      open (meshf, file=fname, status='old', IOSTAT=ios)
      IF (ios.NE.0) PRINT*, ios

      ! number of spatial dimensions
      read (meshf,*)  NSD

      ! degree of curves in u, v, and w direction
      read (meshf,*)  Pu, Qu, Ru

      ! number of control points in each direction
      read (meshf,*)  MCPu, NCPu, OCPu
c
c...  Allocate arrays for knot vectors and for control net
c
      allocate(U_KNOTu(MCPu+Pu+1))
      allocate(V_KNOTu(NCPu+Qu+1))
      allocate(W_KNOTu(OCPu+Ru+1))
      allocate(B_NETu(MCPu,NCPu,OCPu,NSD+1))

!      print *,'NSD                  =', NSD

      do i = 1,MCPu+Pu+1
         if (i < MCPu+Pu+1) then
            read (meshf,"(F13.9)", ADVANCE = "NO") U_KNOTu(i) ! read U_KNOT
         else
            read (meshf,"(F13.9)", ADVANCE = "YES") U_KNOTu(i)
         endif
      enddo

      do i = 1,NCPu+Qu+1
         if (i < NCPu+Qu+1) then
            read (meshf,"(F13.9)", ADVANCE = "NO") V_KNOTu(i) ! read V_KNOT
         else
            read (meshf,"(F13.9)", ADVANCE = "YES") V_KNOTu(i)
         endif
      enddo

      do i = 1,OCPu+Ru+1
         if (i < OCPu+Ru+1) then
            read (meshf,"(F13.9)", ADVANCE = "NO") W_KNOTu(i) ! read W_KNOT
         else
            read (meshf,"(F13.9)", ADVANCE = "YES") W_KNOTu(i)
         endif
      enddo


      do k = 1,(OCPu)
         do j = 1,(NCPu)
            do i = 1,(MCPu)
               do l = 1,NSD+1
                  if (l < NSD+1) then
                     read (meshf,"(F13.9)", ADVANCE = "NO")
     &                    B_NETu(i,j,k,l)
                  else
                     read (meshf,"(F13.9)", ADVANCE = "YES")
     &                    B_NETu(i,j,k,l)
                  endif
               enddo
            enddo
         enddo
      enddo

!      print *,'B_NETini                  =', B_NETini

c
c...  Set number of global nodes

      NNODZu = MCPu*NCPu*OCPu

c...  Set number of elements

      NELu = (MCPu-Pu)*(NCPu-Qu)*(OCPu-Ru)

c...  Read in Master/Slave relationships. These arise from periodic boundary
c     conditions, corners in the geometry, or any other case where two
c     control points must always remain equal.
      allocate(IPERu(NNODZu))
      do i = 1, NNODZu
         read (meshf,*)  IPERu(i)
      enddo


      deallocate(IPERu)
c...  Read in Boundary/Face information
      read (meshf,*) CLOSED_U_flag, CLOSED_V_flag, CLOSED_W_flag
                                ! flag determines if surface is closed
                                !  in some direction

c...  Set number of faces
      NFACEu =
     &     2*(MCPu-Pu)*(NCPu-Qu)*(1-CLOSED_W_flag) +
     &     2*(NCPu-Qu)*(OCPu-Ru)*(1-CLOSED_U_flag) +
     &     2*(MCPu-Pu)*(OCPu-Ru)*(1-CLOSED_V_flag)

      !write(*,*)"NNODZu, NELu, NFACEu = ", NNODZu, NELu, NFACEu

c...  Read in boundary condition indicator at each node, in each direction
c     0 - nothing, 1 - Dir., 2 - Neu., 3 - Periodic! do not trust

      allocate(IBC(NNODZu,NDOF))! This gets overrided in the
      IBC = 0 					! next statement
      do i = 1,NNODZu
         read (meshf,*) IBC(i,1:3)
      enddo

c...  Read in Nodal displacements for inhomogeneous Dirichlet nodes.
      allocate(DIR_BC(NNODZu,NDOF))	! This gets overrided in the
      DIR_BC = 0d0					! next statement	
      do i = 1, NNODZu
         read (meshf,*) DIR_BC(i,1:3)
      enddo

      l = 0
      IBC=0
      DIR_BC=0d0

!      if((periodicx.eq.0).or.(periodicy.eq.0).or.(periodicz.eq.0)) then !non-periodic
      do k = 1, OCPu
         do j = 1, NCPu
            do i = 1, MCPu
               l = l + 1

!      if((periodicx==0)) then ! periodicx =0 => non-periodic
!     X-direction constraint
               if((i==MCPu)) then !xplus
                  IBC(l,1)    = pullxplusu1 ! IBC = 2 pull; =1 constraint
                  IBC(l,2)    = pullxplusu2
                  IBC(l,3)    = pullxplusu3
                  IBC(l,4)    = pullxplusu4
                  IBC(l,5)    = pullxplusu5
                  DIR_BC(l,1) = xplusu1
                  DIR_BC(l,2) = xplusu2
                  DIR_BC(l,3) = xplusu3
                  DIR_BC(l,4) = xplusu4
                  DIR_BC(l,5) = xplusu5
                  print *,'NNoDzu_Right            =', l			   
               endif

               if((i==1)) then !xminus
                  IBC(l,1)    = pullxminusu1 ! IBC = 2 pull; =1 constraint
                  IBC(l,2)    = pullxminusu2
                  IBC(l,3)    = pullxminusu3
                  IBC(l,4)    = pullxminusu4
                  IBC(l,5)    = pullxminusu5
                  DIR_BC(l,1) = xminusu1
                  DIR_BC(l,2) = xminusu2
                  DIR_BC(l,3) = xminusu3
                  DIR_BC(l,4) = xminusu4
                  DIR_BC(l,5) = xminusu5
                  print *,'NNoDzu_left             =', l
               endif
!      endif  ! periodic-X

!      if((periodicy==0)) then ! periodicy =0 => non-periodic
!     Y- direction constraint
               if((j==NCPu)) then
                  IBC(l,1)    = pullyplusu1 ! IBC = 2 pull; =1 constraint
                  IBC(l,2)    = pullyplusu2
                  IBC(l,3)    = pullyplusu3		 !yplus
                  IBC(l,4)    = pullyplusu4
                  IBC(l,5)    = pullyplusu5
                  DIR_BC(l,1) = yplusu1
                  DIR_BC(l,2) = yplusu2
                  DIR_BC(l,3) = yplusu3
                  DIR_BC(l,4) = yplusu4
                  DIR_BC(l,5) = yplusu5
               endif

               if((j==1)) then !yminus
                  IBC(l,1)    = pullyminusu1 ! IBC = 2 pull; =1 constraint
                  IBC(l,2)    = pullyminusu2
                  IBC(l,3)    = pullyminusu3
                  IBC(l,4)    = pullyminusu4
                  IBC(l,5)    = pullyminusu5
                  DIR_BC(l,1) = yminusu1
                  DIR_BC(l,2) = yminusu2
                  DIR_BC(l,3) = yminusu3
                  DIR_BC(l,4) = yminusu4
                  DIR_BC(l,5) = yminusu5
               endif
!      endif  ! periodic-Y

            enddo
         enddo
      enddo
!      endif ! non-periodic X-Y condition
	  
      IBC1 = IBC
      DIRBC1 = DIR_BC

c...  Read in boundary condition indicator on each face, in each direction
c     0 - nothing, 1 - constraint, 2 - load
      allocate(IBC_FACE(NFACEu,NDOF))
      IBC_FACE = 0
      do i = 1,NFACEu
         read (meshf,*) IBC_FACE(i,1:3)
      enddo

c...  Read in actual load data on each face
c     This assumes (for the time being) that the load is contant
c     across each face. This will be modified in the future!
      allocate (LD_FACE(NFACEu,NDOF))
      LD_FACE = 0d0
      do i = 1,NFACEu
         read (meshf,*) LD_FACE(i,1:3)
      enddo

c
c...  Read in body force loads on elements. This assumes constant
c     force across each element. We can modify this in the future.
      allocate(LD_EL(NELu,NDOF))
      LD_EL = 0d0
      do i = 1,NELu
         read (meshf,*) LD_EL(i,1:3)
      enddo

c
c...  Read in material data
c

c...  Density in each element
      allocate(DENS(NELu))
      do i = 1,NELu
         read (meshf,*) DENS(i)
      enddo


c...  Viscosity in each element
      allocate(VISC(NELu))
      do i = 1,NELu
         read (meshf,*) VISC(i)
      enddo

!         print *,'DENS             =', DENS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      ! Rescaling the B_NETu matrix
      close(meshf)

!      open (58, file='ActualGeom.csv', status='unknown')
!      do k = 1,(OCPu)
!         do j = 1,(NCPu)
!            do i = 1,(MCPu)
!                     write(58,*) B_NETu(i,j,k,1), B_NETu(i,j,k,2),
!     &                B_NETu(i,j,k,3), B_NETu(i,j,k,4)
!            enddo
!         enddo
!      enddo

!      close(58)

      B_NETu(:,:,:,1) = B_NETu(:,:,:,1)*scalex
      B_NETu(:,:,:,2) = B_NETu(:,:,:,2)*scaley
      B_NETu(:,:,:,3) = B_NETu(:,:,:,3)*scalez

      allocate(B_NETini(MCPu,NCPu,OCPu,NSD+1))
      B_NETini = B_NETu

!      print *,'B_NETini                  =', B_NETini(:,:,:,2)

!      open (59, file='RescaledGeom.csv', status='unknown')
!      do k = 1,(OCPu)
!         do j = 1,(NCPu)
!            do i = 1,(MCPu)
!                     write(59,*) B_NETu(i,j,k,1), B_NETu(i,j,k,2),
!     &                B_NETu(i,j,k,3), B_NETu(i,j,k,4)
!            enddo
!         enddo
!      enddo
!
!      close(59)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      return

      end
