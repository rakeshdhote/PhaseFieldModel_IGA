c       This subroutine generates the IEN matrix, which relates element numbers
c       and local node numbers to the appropriate global node numbers. The
c       routine also generates the INN matrix, which relates global node
c       number to the "NURBS coordinates" of the node. This routine is for
c       the 3D code
c
c       Nov 15, 2003
c
c       J. Austin Cottrell
c       CAM Graduate Student
c       Institute for Computational Engineering Science
c       The University of Texas at Austin
	
      subroutine genIEN_INN_3D
      use aAdjKeep        ! this f90 module contains all allocatables
c       
      include "common.h"	! common file defines all variables
c       
c...    Local variables
	  integer :: i, j, k, loop1,loop2,loop3,g, e, gtemp, ln
	
      allocate(INN(NNODZ,3))
      allocate(IEN(NEL,(P+1)*(Q+1)*(R+1)))
c       
c...    Initialize matrices and variables
c       
      IEN = 0
      INN = 0
      g = 0
      e = 0

c       
c...    Loop through control points assigning global node
c       numbers and filling out IEN and INN as we go
      do k = 1,OCP            ! loop through control points in W direction
       do j = 1,NCP		! loop through control points in V direction
	        do i = 1,MCP	! loop through control points in U direction
             g = g+1
             INN(g,1) = i
             INN(g,2) = j
             INN(g,3) = k
             if (((i .ge. (P+1)).and.(j .ge. (Q+1))).and.
     &                (k.ge.(R+1))) then
                e = e +1
                do loop1 = 0,r
                   do loop2 = 0,q
                     do loop3 = 0,p
                 gtemp = g - loop1*MCP*NCP -MCP*loop2-
     &                        loop3
                 ln = (P+1)*(Q+1)*loop1+(P+1)*loop2 +
     &                        loop3 + 1
                 IEN(e,ln) = gtemp
			         enddo
		          enddo
		        enddo
		     endif
	        enddo
         enddo
      enddo
	   
      return
      end
	
