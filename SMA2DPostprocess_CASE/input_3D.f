
      subroutine input_3D(iproc,scalex,scaley,scalez)

      use aAdjKeep
      
      include "common.h"
      
      integer :: i, j, k,l, meshf, iproc, istep
      character*30 fname
      character*10 cname

      meshf = 11
!      scalex = 16.0
!      scaley = 16.0
!      scalez = 1.0

c
c...  Read in preliminary information
c
      
      fname = trim('meshu' // cname(iproc)) // '.dat'
      
      open (meshf, file=fname, status='old')
      read (meshf,*)  NSD           ! number of spatial dimensions
      read (meshf,*)  P, Q, R       ! degree of curves in u, v, and w direction
      read (meshf,*)  MCP, NCP, OCP ! number of control points in each direction

c
c...  Allocate arrays for knot vectors and for control net
c
      allocate(U_KNOT(MCP+P+1))    
      allocate(V_KNOT(NCP+Q+1))
      allocate(W_KNOT(OCP+R+1))
      allocate(B_NET(MCP,NCP,OCP,NSD+1))
!      allocate(BW(MCP,NCP,OCP,NSD+1))

c
c...  Read knot vectors and control net
c
      do i = 1,MCP+P+1
         if (i < MCP+P+1) then
            read (meshf,"(F13.9)", ADVANCE = "NO") U_KNOT(i) ! read U_KNOT
         else
            read (meshf,"(F13.9)", ADVANCE = "YES") U_KNOT(i) 
         endif
      enddo

      do i = 1,NCP+Q+1
         if (i < NCP+Q+1) then
            read (meshf,"(F13.9)", ADVANCE = "NO") V_KNOT(i) ! read V_KNOT
         else
            read (meshf,"(F13.9)", ADVANCE = "YES") V_KNOT(i) 
         endif
      enddo

      do i = 1,OCP+R+1
         if (i < OCP+R+1) then
            read (meshf,"(F13.9)", ADVANCE = "NO") W_KNOT(i) ! read W_KNOT
         else
            read (meshf,"(F13.9)", ADVANCE = "YES") W_KNOT(i) 
         endif
      enddo

      do k = 1,(OCP)
         do j = 1,(NCP)
            do i = 1,(MCP)
               do l = 1,NSD+1
                  if (l < NSD+1) then
                     read (meshf,"(F13.9)", ADVANCE = "NO")
     &                    B_NET(i,j,k,l)
                  else
                     read (meshf,"(F13.9)", ADVANCE = "YES") 
     &                    B_NET(i,j,k,l)
                  endif
               enddo
            enddo
         enddo
      enddo

      
c     
c...  Set number of global nodes
      NNODZ = MCP*NCP*OCP
c...  Set number of elements
      NEL = (MCP-P)*(NCP-Q)*(OCP-R)

   
c
c...  Read in Master/Slave relationships. These arise from periodic boundary
c       conditions, corners in the geometry, or any other case where two
c       control points must always remain equal.
c
      allocate(IPER(NNODZ))
      do i = 1, NNODZ
         read (meshf,*)  IPER(i) 
      enddo

c
c...  Read in Boundary/Face information
c
      read (meshf,*) CLOSED_U_flag, CLOSED_V_flag, CLOSED_W_flag
                                ! flag determines if surface is closed
                                !  in some direction

c...  Set number of faces
      NFACE = 2*(MCP-P)*(NCP-Q)*(1-CLOSED_W_flag) +
     &     2*(NCP-Q)*(OCP-R)*(1-CLOSED_U_flag) +
     &     2*(MCP-P)*(OCP-R)*(1-CLOSED_V_flag) 

!      write(*,*)"NNODZ, NEL, NFACE = ", NNODZ, NEL, NFACE

c...  Read in boundary condition indicator at each node, in each direction
c     0 - nothing, 1 - Dir., 2 - Neu., 3 - Periodic
      allocate(IBC(NNODZ,NSD))
      IBC = 0
      do i = 1,NNODZ
         read (meshf,'(I2)', ADVANCE = "NO") IBC(i,1)
         read (meshf,'(I2)', ADVANCE = "NO") IBC(i,2)
         read (meshf,'(I2)', ADVANCE = "YES") IBC(i,3)
      enddo


c...  Read in Nodal displacements for inhomogeneous Dirichlet nodes.
      allocate(DIR_BC(NNODZ,NSD))
      DIR_BC = 0d+0 
      do i = 1, NNODZ
         read (meshf,"(F13.9)", ADVANCE = "NO") DIR_BC(i,1)
         read (meshf,"(F13.9)", ADVANCE = "NO") DIR_BC(i,2)
         read (meshf,"(F13.9)", ADVANCE = "YES") DIR_BC(i,3)
      enddo


c...  Read in boundary condition indicator on each face, in each direction
c     0 - nothing, 1 - constraint, 2 - load
      allocate(IBC_FACE(NFACE,NSD))
      IBC_FACE = 0
      do i = 1,NFACE
         read (meshf,'(I2)', ADVANCE = "NO") IBC_FACE(i,1)
         read (meshf,'(I2)', ADVANCE = "NO") IBC_FACE(i,2)
         read (meshf,'(I2)', ADVANCE = "YES") IBC_FACE(i,3)
      enddo


c...  Read in actual load data on each face
c      This assumes (for the time being) that the load is contant
c      across each face. This will be modified in the future!
      allocate (LD_FACE(NFACE,NSD))
      LD_FACE = 0d+0          
      do i = 1,NFACE
         read (meshf,"(F13.9)", ADVANCE = "NO") LD_FACE(i,1)
         read (meshf,"(F13.9)", ADVANCE = "NO") LD_FACE(i,2)
         read (meshf,"(F13.9)", ADVANCE = "YES") LD_FACE(i,3)
      enddo


     
c
c...  Read in body force loads on elements. This assumes constant
c       force across each element. We can modify this in the future.
      allocate(LD_EL(NEL,NSD))
      LD_EL = 0d+0
      do i = 1,NEL
         read (meshf,"(F13.9)", ADVANCE = "NO") LD_EL(i,1)
         read (meshf,"(F13.9)", ADVANCE = "NO") LD_EL(i,2)
         read (meshf,"(F13.9)", ADVANCE = "YES") LD_EL(i,3)
      enddo

      
c
c...  Read in material data
c

c...  Young's Modulus in each element
!      allocate(EY(NEL))
!      do i = 1,NEL
!         read (meshf,*) EY(i)
!      enddo
!
!c... Poisson's ratio in each element
!      allocate(NU(NEL))
!      do i = 1,NEL
!         read (meshf,*) NU(i)
!      enddo
         
      close(meshf)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      ! Rescaling the B_NETu matrix
!

!      write(*,*) " -- Writing Actual geometry file -- "
!      open (58, file='ActualGeom.csv', status='new')
!      do k = 1,(OCP)
!         do j = 1,(NCP)
!            do i = 1,(MCP)
!                     write(58,*) B_NET(i,j,k,1), B_NET(i,j,k,2),
!     &                B_NET(i,j,k,3), B_NET(i,j,k,4)
!            enddo
!         enddo
!      enddo
!      close(58)

      B_NET(:,:,:,1) = B_NET(:,:,:,1)*scalex
      B_NET(:,:,:,2) = B_NET(:,:,:,2)*scaley
      B_NET(:,:,:,3) = B_NET(:,:,:,3)*scalez

!      allocate(B_NETini(MCPu,NCPu,OCPu,NSD+1))
!      B_NETini = B_NETu

!      print *,'B_NETini                  =', B_NETini(:,:,:,2)

!      write(*,*) " -- Writing Rescaled geometry file -- "
!      open (59, file='RescaledGeom.csv', status='new')
!      do k = 1,(OCP)
!         do j = 1,(NCP)
!            do i = 1,(MCP)
!                     write(59,*) B_NET(i,j,k,1), B_NET(i,j,k,2),
!     &                B_NET(i,j,k,3), B_NET(i,j,k,4)
!            enddo
!         enddo
!      enddo
!
!      close(59)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      return
      end
