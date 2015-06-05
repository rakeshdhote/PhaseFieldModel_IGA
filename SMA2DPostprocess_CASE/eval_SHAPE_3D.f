c
c  Subroutine eval_SHAPE.f consumes an element number and the coordinates in the
c     parent element of an integration point and returns the vector of all local
c     basis functions evaluated at the point and the matrix of gradients for all
c     nonzero bais functions with respect to parameters u and v and with respect
c     to x and y. 
c
c     Future modificatons should allow for the 3D case...here NSD is assumed 2.
c
c
c
c  January 20, 2004
c
c
c  J. Austin Cottrell
c  Jurijs Bazilevs
c
c  CAM Graduate Students
c  Institute for Computational Engineering Science
c  The University of Texas at Austin



      subroutine eval_SHAPE_3D(e,u_hat,v_hat,w_hat,shl,shgradl,shgradg)

      use aAdjKeep

      include "common.h"
      
c     --------------VARIABLE DECLARATIONS--------------------------------
c...  Element number
      integer e
c...  u and v coordinates of integration point in parent element
      real*8 u_hat, v_hat, w_hat
      
c...  Vector of Local basis function values at (u_hat, v_hat), local and 
c          global gradients.
      real*8 shl(NSHL),shgradl(NSHL,NSD),shgradg(NSHL,NSD),
     &  dxdxi(NSD,NSD), dxidx(NSD,NSD)

c...  Local Variables
c    
c     1D nonrational basis functions and derivs in u and v
      real*8 N(2,P+1), M(2,Q+1), O(2,R+1)
c     u and v coordinates of integration point, denominator and derivative sums
      real*8 u, v, w, denom_sum, derv_sum_U, derv_sum_V,
     &  derv_sum_W
c     NURBS coordinates, counters for loops
      integer ni, nj, nk, i, j, k, icount
c     temporary variables
      real*8 tmp
c ------------------------------------------------------------------


c     Get NURBS coordinates for local node 1
      ni = INN(IEN(e,1),1)
      nj = INN(IEN(e,1),2)
      nk = INN(IEN(e,1),3)

c     Get u and v coordinates of integration point
      u = ((U_KNOT(ni+1)-U_KNOT(ni))*u_hat + 
     &     U_KNOT(ni+1) + U_KNOT(ni))/2d+0
      v = ((V_KNOT(nj+1)-V_KNOT(nj))*v_hat + 
     &     V_KNOT(nj+1) + V_KNOT(nj))/2d+0
      w = ((W_KNOT(nk+1)-W_KNOT(nk))*w_hat + 
     &     W_KNOT(nk+1) + W_KNOT(nk))/2d+0
      
      
c     Evaluate 1D shape functions and derivatives each direction
      call dersbasisfuns(ni,P,MCP,u,U_KNOT,N)    ! calculate in u direction
      call dersbasisfuns(nj,Q,NCP,v,V_KNOT,M)    ! calculate in v direction
      call dersbasisfuns(nk,R,OCP,w,W_KNOT,O)    ! calculate in v direction

c     Form basis functions and derivatives dR/du and dR/dv
      icount = 0
      denom_sum = 0d+0
      derv_sum_U = 0d+0
      derv_sum_V = 0d+0
      derv_sum_W = 0d+0
      shgradl = 0d+0


      do k = 0, R
        do j = 0,Q
          do i = 0,P
            icount = icount+1
            
c...        basis functions
            shl(icount) = N(1,P+1-i)*M(1,Q+1-j)*O(1,R+1-k)*
     &        B_NET(ni-i,nj-j,nk-k,NSD+1)
            denom_sum = denom_sum + shl(icount)
            
c...        derivatives
            shgradl(icount,1) = 
     &        N(2,P+1-i)*M(1,Q+1-j)*O(1,R+1-k)*
     &        B_NET(ni-i,nj-j,nk-k,NSD+1) ! u
            derv_sum_U = derv_sum_U + shgradl(icount,1)
            shgradl(icount,2) = 
     &        N(1,P+1-i)*M(2,Q+1-j)*O(1,R+1-k)*
     &        B_NET(ni-i,nj-j,nk-k,NSD+1) ! v
            derv_sum_V = derv_sum_V + shgradl(icount,2)
            shgradl(icount,3) = 
     &        N(1,P+1-i)*M(1,Q+1-j)*O(2,R+1-k)*
     &        B_NET(ni-i,nj-j,nk-k,NSD+1) ! w
            derv_sum_W = derv_sum_W + shgradl(icount,3)
            
          enddo
        enddo
      enddo
      
c     Divide through by denominator
      
      do i = 1,NSHL
        shgradl(i,1) = shgradl(i,1)/denom_sum - 
     &    (shl(i)*derv_sum_U)/(denom_sum**2)
        shgradl(i,2) = shgradl(i,2)/denom_sum - 
     &    (shl(i)*derv_sum_V)/(denom_sum**2)
        shgradl(i,3) = shgradl(i,3)/denom_sum - 
     &    (shl(i)*derv_sum_W)/(denom_sum**2)
        shl(i) = shl(i)/denom_sum
      enddo
      
      
c     
c...  Now calculate gradients.
c

c     calculate dx/dxi
      


      dxdxi = 0d+0
      icount = 0
      
      do k = 0, R
        do j = 0, Q
          do i = 0, P
            icount = icount + 1
            
            dxdxi(1,1) = dxdxi(1,1) + B_NET(ni-i,nj-j,nk-k,1) *
     &        shgradl(icount,1)
            dxdxi(1,2) = dxdxi(1,2) + B_NET(ni-i,nj-j,nk-k,1) *
     &        shgradl(icount,2)
            dxdxi(1,3) = dxdxi(1,3) + B_NET(ni-i,nj-j,nk-k,1) *
     &        shgradl(icount,3)
            dxdxi(2,1) = dxdxi(2,1) + B_NET(ni-i,nj-j,nk-k,2) *
     &        shgradl(icount,1)
            dxdxi(2,2) = dxdxi(2,2) + B_NET(ni-i,nj-j,nk-k,2) *
     &        shgradl(icount,2)
            dxdxi(2,3) = dxdxi(2,3) + B_NET(ni-i,nj-j,nk-k,2) *
     &      shgradl(icount,3)
            dxdxi(3,1) = dxdxi(3,1) + B_NET(ni-i,nj-j,nk-k,3) *
     &        shgradl(icount,1)
            dxdxi(3,2) = dxdxi(3,2) + B_NET(ni-i,nj-j,nk-k,3) *
     &        shgradl(icount,2)
            dxdxi(3,3) = dxdxi(3,3) + B_NET(ni-i,nj-j,nk-k,3) *
     &        shgradl(icount,3)
            
          enddo
        enddo
      enddo

      
c
c.... compute the inverse of deformation gradient
c
      
      dxidx(1,1) =   dxdxi(2,2) * dxdxi(3,3) 
     &     - dxdxi(3,2) * dxdxi(2,3)
      dxidx(1,2) =   dxdxi(3,2) * dxdxi(1,3) 
     &     - dxdxi(1,2) * dxdxi(3,3)
      dxidx(1,3) =  dxdxi(1,2) * dxdxi(2,3) 
     &     - dxdxi(1,3) * dxdxi(2,2)
      tmp          = 1d+0 / ( dxidx(1,1) * dxdxi(1,1) 
     &     + dxidx(1,2) * dxdxi(2,1)  
     &     + dxidx(1,3) * dxdxi(3,1) )
      dxidx(1,1) = dxidx(1,1) * tmp
      dxidx(1,2) = dxidx(1,2) * tmp
      dxidx(1,3) = dxidx(1,3) * tmp
      dxidx(2,1) = (dxdxi(2,3) * dxdxi(3,1) 
     &     - dxdxi(2,1) * dxdxi(3,3)) * tmp
      dxidx(2,2) = (dxdxi(1,1) * dxdxi(3,3) 
     &     - dxdxi(3,1) * dxdxi(1,3)) * tmp
      dxidx(2,3) = (dxdxi(2,1) * dxdxi(1,3) 
     &     - dxdxi(1,1) * dxdxi(2,3)) * tmp
      dxidx(3,1) = (dxdxi(2,1) * dxdxi(3,2) 
     &  - dxdxi(2,2) * dxdxi(3,1)) * tmp
      dxidx(3,2) = (dxdxi(3,1) * dxdxi(1,2) 
     &     - dxdxi(1,1) * dxdxi(3,2)) * tmp
      dxidx(3,3) = (dxdxi(1,1) * dxdxi(2,2) 
     &     - dxdxi(1,2) * dxdxi(2,1)) * tmp
      
      DetJ = 1d+0/tmp           ! Note that DetJ resides in common


      

      do i = 1, NSHL

        shgradg(i,1) = shgradl(i,1) * dxidx(1,1) + 
     &    shgradl(i,2) * dxidx(2,1) +
     &    shgradl(i,3) * dxidx(3,1)
        shgradg(i,2) = shgradl(i,1) * dxidx(1,2) + 
     &    shgradl(i,2) * dxidx(2,2) +
     &    shgradl(i,3) * dxidx(3,2) 
        shgradg(i,3) = shgradl(i,1) * dxidx(1,3) + 
     &    shgradl(i,2) * dxidx(2,3) +
     &    shgradl(i,3) * dxidx(3,3) 
        
      enddo
      
      
      
      if (DetJ < 0d+0) DetJ = -1d+0*DetJ
      
      return
      end
      
      


