c
c  Subroutine eval_FACE_3D.f is modified from eval_SHAPE.f in the 2D version.
c    Here, in the 3D version, we need to incorporate the boundary conditions
c    over the faces. 
c
c
c  January 21, 2004
c
c
c  J. Austin Cottrell
c  CAM Graduate Student
c  Institute for Computational Engineering Science
c  The University of Texas at Austin



      subroutine eval_FACE_3D(f, temp1, temp2, ni, nj, nk, nip,
     &  njp, nkp, shb, shbg, shbp, dxidx, nor)

      use aAdjKeep

      use common
      implicit none

c     --------------VARIABLE DECLARATIONS--------------------------------
c...  Face number
      integer f
c...  u, v, w coordinates of integration point in parent element and temps
      real*8 u_hat, v_hat, w_hat, temp1, temp2, tmp,
     &  nor_0(NSD), nor(NSD), cofF(NSD,NSD), du, dv, dw
      
c...  Vector of Local basis function values at (u_hat, v_hat, w_hat).
      real*8 shb(NSHLu), shbg(NSHLu,NSD), shbp(NSHLp), dxdxi(NSD,NSD),
     &  ee,ff,gg, dxidx(NSD,NSD), xloc(2)
      
c...  Local Variables
c     
c     1D nonrational basis functions and derivs in u and v
      real*8 N(2,Pu+1), M(2,Qu+1), O(2,Ru+1), shbl(NSHLu,NSD)
      real*8 Np(2,Pp+1), Mp(2,Qp+1), Op(2,Rp+1)
c     u,v,w coordinates of integration point, denominator
      real*8 u, v, w, denom_sum, derv_sum_U, derv_sum_V,derv_sum_W

c     NURBS coordinates, counters for loops
      integer ni, nj, nk, icount, i, j, k, nip, njp, nkp
c ------------------------------------------------------------------

      

c... The details are dependent of face orientation.

      nor_0 = 0d+0
      
      if (FACE_OR(f).eq.1) then
        
        u_hat = temp1
        v_hat = temp2
        w_hat = -1d+0           ! -1.0 ensure w will be at beginning
                                !  knot vector
        nor_0(3) = -1d+0
        
      else if (FACE_OR(f).eq.2) then
        
        u_hat = temp1
        v_hat = -1d+0
        w_hat = temp2
        
        nor_0(2) = -1d+0
        
      else if (FACE_OR(f).eq.3) then
        
        u_hat = 1d+0            ! 1.0 ensures u will be at end of knot
        v_hat = temp1           !  vector
        w_hat = temp2
        
        nor_0(1) = 1d+0
        
      else if (FACE_OR(f).eq.4) then
        
        u_hat = temp1  
        v_hat = 1d+0
        w_hat = temp2
        
        nor_0(2) = 1d+0
        
      else if (FACE_OR(f).eq.5) then
        
        u_hat = -1d+0      
        v_hat = temp1
        w_hat = temp2   
        
        nor_0(1) = -1d+0
        
      else 
        
        u_hat = temp1
        v_hat = temp2
        w_hat = 1d+0 
        
        nor_0(3) = 1d+0
        
      endif

      du = U_KNOTu(ni+1)-U_KNOTu(ni)
      dv = V_KNOTu(nj+1)-V_KNOTu(nj)
      dw = W_KNOTu(nk+1)-W_KNOTu(nk)


c...  Find integration point in parameter space

      u = ((U_KNOTu(ni+1)-U_KNOTu(ni))*u_hat + 
     &  U_KNOTu(ni+1) + U_KNOTu(ni))/2d+0
      
      v = ((V_KNOTu(nj+1)-V_KNOTu(nj))*v_hat + 
     &  V_KNOTu(nj+1) + V_KNOTu(nj))/2d+0
      
      w = ((W_KNOTu(nk+1)-W_KNOTu(nk))*w_hat + 
     &  W_KNOTu(nk+1) + W_KNOTu(nk))/2d+0
      

      
c     Evaluate 1D shape functions for velocity
         
      call dersbasisfuns(ni,Pu,MCPu,u,1,U_KNOTu,N) ! calculate in u direction
      call dersbasisfuns(nj,Qu,NCPu,v,1,V_KNOTu,M) ! calculate in v direction
      call dersbasisfuns(nk,Ru,OCPu,w,1,W_KNOTu,O) ! calculate in w direction

c     Evaluate 1D shape functions for pressure
         
      call dersbasisfuns(nip,Pp,MCPp,u,1,U_KNOTp,Np) ! calculate in u direction
      call dersbasisfuns(njp,Qp,NCPp,v,1,V_KNOTp,Mp) ! calculate in v direction
      call dersbasisfuns(nkp,Rp,OCPp,w,1,W_KNOTp,Op) ! calculate in w direction
      
c     Form basis functions for pressure      
      
      icount = 0
      denom_sum = 0d+0
      
      shbp = 0d+0
      
      do k = 0, Rp
        do j = 0, Qp
          do i = 0, Pp
            icount = icount+1
            
c...        basis functions
            shbp(icount) = Np(1,Pp+1-i)*Mp(1,Qp+1-j)*Op(1,Rp+1-k)*
     &        B_NETu(nip-i,njp-j,nkp-k,NSD+1)
            denom_sum = denom_sum + shbp(icount)
            
          enddo
        enddo
      enddo
      
      
c     Divide through by denominator
      do i = 1,NSHLp
        shbp(i) = shbp(i)/denom_sum
      enddo      
      
      
      
c     Form basis functions for velocity
      icount = 0
      denom_sum = 0d+0

      shb = 0d+0
      derv_sum_U = 0d+0
      derv_sum_V = 0d+0
      derv_sum_W = 0d+0
      shbl = 0d+0
      shbg = 0d+0
      
      do k = 0,Ru
        do j = 0,Qu
          do i = 0,Pu
            icount = icount+1
            
c...        basis functions
            shb(icount) = N(1,Pu+1-i)*M(1,Qu+1-j)*O(1,Ru+1-k)*
     &        B_NETu(ni-i,nj-j,nk-k,NSD+1)
            denom_sum = denom_sum + shb(icount)
            
            shbl(icount,1) = 
     &        N(2,Pu+1-i)*M(1,Qu+1-j)*O(1,Ru+1-k)*
     &        B_NETu(ni-i,nj-j,nk-k,NSD+1) ! u
            derv_sum_U = derv_sum_U + shbl(icount,1)
            shbl(icount,2) = 
     &        N(1,Pu+1-i)*M(2,Qu+1-j)*O(1,Ru+1-k)*
     &        B_NETu(ni-i,nj-j,nk-k,NSD+1) ! v
            derv_sum_V = derv_sum_V + shbl(icount,2)
            shbl(icount,3) = 
     &        N(1,Pu+1-i)*M(1,Qu+1-j)*O(2,Ru+1-k)*
     &        B_NETu(ni-i,nj-j,nk-k,NSD+1) ! w
            derv_sum_W = derv_sum_W + shbl(icount,3)
            
          enddo
        enddo
      enddo
      
      
c     Divide through by denominator
      do i = 1,NSHLu
        shbl(i,1) = shbl(i,1)/denom_sum - 
     &    (shb(i)*derv_sum_U)/(denom_sum**2)
        shbl(i,2) = shbl(i,2)/denom_sum - 
     &    (shb(i)*derv_sum_V)/(denom_sum**2)
        shbl(i,3) = shbl(i,3)/denom_sum - 
     &    (shb(i)*derv_sum_W)/(denom_sum**2)
        shb(i) = shb(i)/denom_sum
      enddo      

c...  Now we calculate the face Jacobian

      dxdxi(:,:) = 0d+0
      icount = 0
      
      do k = 0, Ru
        do j = 0, Qu
          do i = 0, Pu
            icount = icount + 1
            
            dxdxi(1,1) = dxdxi(1,1) + B_NETu(ni-i,nj-j,nk-k,1) *
     &        shbl(icount,1)
            dxdxi(1,2) = dxdxi(1,2) + B_NETu(ni-i,nj-j,nk-k,1) *
     &        shbl(icount,2)
            dxdxi(1,3) = dxdxi(1,3) + B_NETu(ni-i,nj-j,nk-k,1) *
     &        shbl(icount,3)
            dxdxi(2,1) = dxdxi(2,1) + B_NETu(ni-i,nj-j,nk-k,2) *
     &        shbl(icount,1)
            dxdxi(2,2) = dxdxi(2,2) + B_NETu(ni-i,nj-j,nk-k,2) *
     &        shbl(icount,2)
            dxdxi(2,3) = dxdxi(2,3) + B_NETu(ni-i,nj-j,nk-k,2) *
     &        shbl(icount,3)
            dxdxi(3,1) = dxdxi(3,1) + B_NETu(ni-i,nj-j,nk-k,3) *
     &        shbl(icount,1)
            dxdxi(3,2) = dxdxi(3,2) + B_NETu(ni-i,nj-j,nk-k,3) *
     &        shbl(icount,2)
            dxdxi(3,3) = dxdxi(3,3) + B_NETu(ni-i,nj-j,nk-k,3) *
     &        shbl(icount,3)
               
          enddo
        enddo
      enddo

c... Calculate Cofactor

      cofF(1,1) = dxdxi(2,2)*dxdxi(3,3) - dxdxi(2,3)*dxdxi(3,2)
      cofF(1,2) = -(dxdxi(2,1)*dxdxi(3,3) - dxdxi(2,3)*dxdxi(3,1))
      cofF(1,3) = dxdxi(2,1)*dxdxi(3,2) - dxdxi(2,2)*dxdxi(3,1)
      
      cofF(2,1) = -(dxdxi(1,2)*dxdxi(3,3) - dxdxi(1,3)*dxdxi(3,2))
      cofF(2,2) = dxdxi(1,1)*dxdxi(3,3) - dxdxi(1,3)*dxdxi(3,1)
      cofF(2,3) = -(dxdxi(1,1)*dxdxi(3,2) - dxdxi(1,2)*dxdxi(3,1))
      
      cofF(3,1) = dxdxi(1,2)*dxdxi(2,3) - dxdxi(1,3)*dxdxi(2,2)
      cofF(3,2) = -(dxdxi(1,1)*dxdxi(2,3) - dxdxi(1,3)*dxdxi(2,1))
      cofF(3,3) = dxdxi(1,1)*dxdxi(2,2) - dxdxi(1,2)*dxdxi(2,1)
      
c...  Calculate Normal Based on Nansen's Formula

c*      nor = 0d+0
c*      
c*      do i = 1, NSD
c*        do j = 1, NSD
c*          nor(i) = nor(i) + cofF(i,j)*nor_0(j)
c*        enddo
c*      enddo
c*      
c*      tmp = sqrt(nor(1)**2+nor(2)**2+nor(3)**2)
c*      
c*      nor(:) = nor(:)/tmp

cccccc      nor(1) = xloc(1)
cccccc      nor(2) = xloc(2)
      
c...  Jacobian of Face Mapping
      
      if ((FACE_OR(f).eq.1).or.(FACE_OR(f).eq.6)) then
        ee = dxdxi(1,1)**2+dxdxi(2,1)**2+dxdxi(3,1)**2
        ff = dxdxi(1,1)*dxdxi(1,2)+dxdxi(2,1)*dxdxi(2,2)+ 
     &    dxdxi(3,1)*dxdxi(3,2)
        gg = dxdxi(1,2)**2+dxdxi(2,2)**2+dxdxi(3,2)**2
        
      else if ((FACE_OR(f).eq.2).or.(FACE_OR(f).eq.4)) then
        ee = dxdxi(1,1)**2+dxdxi(2,1)**2+dxdxi(3,1)**2
        ff = dxdxi(1,1)*dxdxi(1,3)+dxdxi(2,1)*dxdxi(2,3)+ 
     &    dxdxi(3,1)*dxdxi(3,3)
        gg = dxdxi(1,3)**2+dxdxi(2,3)**2+dxdxi(3,3)**2
        
      else if ((FACE_OR(f).eq.3).or.(FACE_OR(f).eq.5)) then
        ee = dxdxi(1,2)**2+dxdxi(2,2)**2+dxdxi(3,2)**2
        ff = dxdxi(1,2)*dxdxi(1,3)+dxdxi(2,2)*dxdxi(2,3)+ 
     &    dxdxi(3,2)*dxdxi(3,3)
        gg = dxdxi(1,3)**2+dxdxi(2,3)**2+dxdxi(3,3)**2
      endif
      
      DetJb = sqrt(ee*gg-ff**2) ! Jacobian of face mapping
      
c------------------------------------------------------------------------

c...  Computation of normal
      
      if (FACE_OR(f).eq.1) then
        nor(1) = dxdxi(2,2)*dxdxi(3,1) - dxdxi(3,2)*dxdxi(2,1)
        nor(2) = dxdxi(3,2)*dxdxi(1,1) - dxdxi(1,2)*dxdxi(3,1)
        nor(3) = dxdxi(1,2)*dxdxi(2,1) - dxdxi(2,2)*dxdxi(1,1)
      else if(FACE_OR(f).eq.2) then
        nor(1) = dxdxi(2,1)*dxdxi(3,3) - dxdxi(3,1)*dxdxi(2,3)
        nor(2) = dxdxi(3,1)*dxdxi(1,3) - dxdxi(1,1)*dxdxi(3,3)
        nor(3) = dxdxi(1,1)*dxdxi(2,3) - dxdxi(2,1)*dxdxi(1,3)
      else if(FACE_OR(f).eq.3) then
        nor(1) = dxdxi(2,2)*dxdxi(3,3) - dxdxi(3,2)*dxdxi(2,3)
        nor(2) = dxdxi(3,2)*dxdxi(1,3) - dxdxi(1,2)*dxdxi(3,3)
        nor(3) = dxdxi(1,2)*dxdxi(2,3) - dxdxi(2,2)*dxdxi(1,3)
      else if(FACE_OR(f).eq.4) then  
        nor(1) = -dxdxi(2,1)*dxdxi(3,3) + dxdxi(3,1)*dxdxi(2,3)
        nor(2) = -dxdxi(3,1)*dxdxi(1,3) + dxdxi(1,1)*dxdxi(3,3)
        nor(3) = -dxdxi(1,1)*dxdxi(2,3) + dxdxi(2,1)*dxdxi(1,3)
      else if(FACE_OR(f).eq.5) then  
        nor(1) = -dxdxi(2,2)*dxdxi(3,3) + dxdxi(3,2)*dxdxi(2,3)
        nor(2) = -dxdxi(3,2)*dxdxi(1,3) + dxdxi(1,2)*dxdxi(3,3)
        nor(3) = -dxdxi(1,2)*dxdxi(2,3) + dxdxi(2,2)*dxdxi(1,3)
      else if(FACE_OR(f).eq.6) then  
        nor(1) = -dxdxi(2,2)*dxdxi(3,1) + dxdxi(3,2)*dxdxi(2,1)
        nor(2) = -dxdxi(3,2)*dxdxi(1,1) + dxdxi(1,2)*dxdxi(3,1)
        nor(3) = -dxdxi(1,2)*dxdxi(2,1) + dxdxi(2,2)*dxdxi(1,1)
      endif
      
      tmp = sqrt( nor(1)**2 + nor(2)**2 + nor(3)**2)
      nor(:) = nor(:)/tmp
      
c*      write(*,*) xloc(1), xloc(2)
c*      write(*,*) nor(1), nor(2), nor(3)
c*
c*      write(*,*) dxdxi(1,1), dxdxi(1,2), dxdxi(1,3)
c*      write(*,*) dxdxi(2,1), dxdxi(2,2), dxdxi(2,3)
c*      write(*,*) dxdxi(3,1), dxdxi(3,2), dxdxi(3,3)
      
c-----------------------------------------------------------------------c      
      
      
c
c.... compute the inverse of deformation gradient
c
      
                                ! 1/ Det
      
      tmp = 1d+0/(dxdxi(1,1)*dxdxi(2,2)*dxdxi(3,3) -
     &  dxdxi(1,1)*dxdxi(2,3)*dxdxi(3,2) -
     &  dxdxi(2,1)*dxdxi(1,2)*dxdxi(3,3) +
     &  dxdxi(2,1)*dxdxi(1,3)*dxdxi(3,2) +
     &  dxdxi(3,1)*dxdxi(1,2)*dxdxi(2,3) -
     &  dxdxi(3,1)*dxdxi(1,3)*dxdxi(2,2))
      
      
                                ! F^(-1) = 1/Det * (cofF)^T
      do i = 1, NSD
        do j = 1, NSD
          dxidx(i,j) = cofF(j,i)
        enddo
      enddo
      
      dxidx(:,:) = dxidx(:,:)*tmp
      
      
c*      dxidx(1,1) =   dxdxi(2,2) * dxdxi(3,3) 
c*     &     - dxdxi(3,2) * dxdxi(2,3)
c*      dxidx(1,2) =   dxdxi(3,2) * dxdxi(1,3) 
c*     &     - dxdxi(1,2) * dxdxi(3,3)
c*      dxidx(1,3) =  dxdxi(1,2) * dxdxi(2,3) 
c*     &     - dxdxi(1,3) * dxdxi(2,2)
c*      tmp          = 1d+0 / ( dxidx(1,1) * dxdxi(1,1) 
c*     &     + dxidx(1,2) * dxdxi(2,1)  
c*     &     + dxidx(1,3) * dxdxi(3,1) )
c*      dxidx(1,1) = dxidx(1,1) * tmp
c*      dxidx(1,2) = dxidx(1,2) * tmp
c*      dxidx(1,3) = dxidx(1,3) * tmp
c*      dxidx(2,1) = (dxdxi(2,3) * dxdxi(3,1) 
c*     &     - dxdxi(2,1) * dxdxi(3,3)) * tmp
c*      dxidx(2,2) = (dxdxi(1,1) * dxdxi(3,3) 
c*     &     - dxdxi(3,1) * dxdxi(1,3)) * tmp
c*      dxidx(2,3) = (dxdxi(2,1) * dxdxi(1,3) 
c*     &     - dxdxi(1,1) * dxdxi(2,3)) * tmp
c*      dxidx(3,1) = (dxdxi(2,1) * dxdxi(3,2) 
c*     &     - dxdxi(2,2) * dxdxi(3,1)) * tmp
c*      dxidx(3,2) = (dxdxi(3,1) * dxdxi(1,2) 
c*     &     - dxdxi(1,1) * dxdxi(3,2)) * tmp
c*      dxidx(3,3) = (dxdxi(1,1) * dxdxi(2,2) 
c*     &     - dxdxi(1,2) * dxdxi(2,1)) * tmp

      do i = 1, NSHLu
        
        shbg(i,1) = shbl(i,1) * dxidx(1,1) + 
     &    shbl(i,2) * dxidx(2,1) +
     &    shbl(i,3) * dxidx(3,1)
        shbg(i,2) = shbl(i,1) * dxidx(1,2) + 
     &    shbl(i,2) * dxidx(2,2) +
     &    shbl(i,3) * dxidx(3,2) 
        shbg(i,3) = shbl(i,1) * dxidx(1,3) + 
     &    shbl(i,2) * dxidx(2,3) +
     &    shbl(i,3) * dxidx(3,3) 
        
      enddo

      dxidx(1,:) = dxidx(1,:)*2d+0/du
      dxidx(2,:) = dxidx(2,:)*2d+0/dv
      dxidx(3,:) = dxidx(3,:)*2d+0/dw
      
      return
      end
