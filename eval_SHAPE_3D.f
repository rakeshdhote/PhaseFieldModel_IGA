      subroutine eval_SHAPE_3D(
     &     e,
     &     u_hat,v_hat,w_hat,
     &     shlu,shgradgu,shhessg,
     &     dxidx)

      use aAdjKeep

      use common
      implicit none

c     --------------VARIABLE DECLARATIONS--------------------------------
c...  Element number
      integer e, ep, nders
c...  u and v coordinates of integration point in parent element
      real*8 u_hat, v_hat, w_hat, du, dv, dw, ds, ds2, ds3
      
c...  Vector of Local basis function values at (u_hat, v_hat), local and 
c          global gradients.
      real*8
     &     shlu(NSHLu),shgradlu(NSHLu,NSD),shgradgu(NSHLu,NSD),
     &     dxdxi(NSD,NSD), dxidx(NSD,NSD), shhessg(NSHLu,NSD,NSD),
     &     lapgu(NSHLu),
     &     tempshl(NSHLu),tempshgradl(NSHLu,NSD),shhessl(NSHLu,6),
     &     tempshhessl(NSHLu,6), dxdxixj(NSD,6), locLHS(6,6)

c...  Local Variables
c    
c     1D nonrational basis functions and derivs in u and v

      real*8 
     &     Nnu(3,Pu+1), ! Nu = Nnu, since nu is viscosity
     &     Mmu(3,Qu+1), ! Mu = Mmu, since mu is viscosity
     &     Ou (3,Ru+1)
      
c     u and v coordinates of integration point, denominator and derivative sums
      real*8 u, v, w, denom_sum, derv_sum_U, derv_sum_V,
     &     derv_sum_W, derv_sum_UU,
     &     derv_sum_UV, derv_sum_UW,derv_sum_VV, derv_sum_VW,derv_sum_WW
c     NURBS coordinates, counters for loops
      integer niu, nju, nku, i, j, k, icount, aa, bb, cc
c     temporary variables
      real*8 tmp, tmprow(6)
c ------------------------------------------------------------------


c     Get NURBS coordinates for local node 1
      
      niu = INNu(IENu(e,1),1)
      nju = INNu(IENu(e,1),2)
      nku = INNu(IENu(e,1),3)

c     Get u and v coordinates of integration point
      
      u = ( (U_KNOTu(niu+1) - U_KNOTu(niu))*u_hat
     &     + U_KNOTu(niu+1) + U_KNOTu(niu)  )*0.5d0
      v = ( (V_KNOTu(nju+1) - V_KNOTu(nju))*v_hat
     &     + V_KNOTu(nju+1) + V_KNOTu(nju)  )*0.5d0
      w = ( (W_KNOTu(nku+1) - W_KNOTu(nku))*w_hat
     &     + W_KNOTu(nku+1) + W_KNOTu(nku)  )*0.5d0

c     Get knot span sizes
      
      du = U_KNOTu(niu+1)-U_KNOTu(niu)
      dv = V_KNOTu(nju+1)-V_KNOTu(nju)
      dw = W_KNOTu(nku+1)-W_KNOTu(nku)
      
c     Evaluate 1D shape functions and derivatives each direction
      
      nders = 2
      
      call dersbasisfuns(niu,Pu,MCPu,u,nders,U_KNOTu,Nnu) !  u direction
      call dersbasisfuns(nju,Qu,NCPu,v,nders,V_KNOTu,Mmu) !  v direction
      call dersbasisfuns(nku,Ru,OCPu,w,nders,W_KNOTu,Ou)  !   w direction


c...  DISPLACEMENT
      
      icount = 0
      denom_sum = 0d0
      derv_sum_U = 0d0
      derv_sum_V = 0d0
      derv_sum_W = 0d0
      shgradlu = 0d0
      shhessl = 0d0
      shhessg = 0d0
      tempshl = 0d0
      tempshgradl = 0d0
      tempshhessl = 0d0
      
      
      do k = 0, Ru
        do j = 0,Qu
          do i = 0,Pu
            icount = icount+1
            
c...        basis functions
            shlu(icount) = Nnu(1,Pu+1-i)*Mmu(1,Qu+1-j)*Ou(1,Ru+1-k)*
     &        B_NETu(niu-i,nju-j,nku-k,NSD+1)
            denom_sum = denom_sum + shlu(icount)
            
c...        derivatives
            shgradlu(icount,1) = 
     &        Nnu(2,Pu+1-i)*Mmu(1,Qu+1-j)*Ou(1,Ru+1-k)*
     &        B_NETu(niu-i,nju-j,nku-k,NSD+1) ! u
            derv_sum_U = derv_sum_U + shgradlu(icount,1)
            shgradlu(icount,2) = 
     &        Nnu(1,Pu+1-i)*Mmu(2,Qu+1-j)*Ou(1,Ru+1-k)*
     &        B_NETu(niu-i,nju-j,nku-k,NSD+1) ! v
            derv_sum_V = derv_sum_V + shgradlu(icount,2)
            shgradlu(icount,3) = 
     &        Nnu(1,Pu+1-i)*Mmu(1,Qu+1-j)*Ou(2,Ru+1-k)*
     &        B_NETu(niu-i,nju-j,nku-k,NSD+1) ! w
            derv_sum_W = derv_sum_W + shgradlu(icount,3)
            
          enddo
        enddo
      enddo
      
c     Divide through by denominator
      
      tempshl = shlu
      tempshgradl = shgradlu
      
      ds  = 1d0/denom_sum
      ds2 = ds*ds

      shgradlu(:,1) = shgradlu(:,1)*ds - (shlu*derv_sum_U)*ds2
      shgradlu(:,2) = shgradlu(:,2)*ds - (shlu*derv_sum_V)*ds2
      shgradlu(:,3) = shgradlu(:,3)*ds - (shlu*derv_sum_W)*ds2
      shlu = shlu*ds


c...  Now calculate gradients.

c     calculate dx/dxi
      
      dxdxi = 0d0
      icount = 0
      
      do k = 0, Ru
      do j = 0, Qu
      do i = 0, Pu
         icount = icount + 1
            
         dxdxi(1,1) = dxdxi(1,1) + B_NETu(niu-i,nju-j,nku-k,1) *
     &        shgradlu(icount,1)
         dxdxi(1,2) = dxdxi(1,2) + B_NETu(niu-i,nju-j,nku-k,1) *
     &        shgradlu(icount,2)
         dxdxi(1,3) = dxdxi(1,3) + B_NETu(niu-i,nju-j,nku-k,1) *
     &        shgradlu(icount,3)
         dxdxi(2,1) = dxdxi(2,1) + B_NETu(niu-i,nju-j,nku-k,2) *
     &        shgradlu(icount,1)
         dxdxi(2,2) = dxdxi(2,2) + B_NETu(niu-i,nju-j,nku-k,2) *
     &        shgradlu(icount,2)
         dxdxi(2,3) = dxdxi(2,3) + B_NETu(niu-i,nju-j,nku-k,2) *
     &        shgradlu(icount,3)
         dxdxi(3,1) = dxdxi(3,1) + B_NETu(niu-i,nju-j,nku-k,3) *
     &        shgradlu(icount,1)
         dxdxi(3,2) = dxdxi(3,2) + B_NETu(niu-i,nju-j,nku-k,3) *
     &        shgradlu(icount,2)
         dxdxi(3,3) = dxdxi(3,3) + B_NETu(niu-i,nju-j,nku-k,3) *
     &        shgradlu(icount,3)

      enddo
      enddo
      enddo

      
c.... compute the inverse of deformation gradient
      
      dxidx(1,1) =   dxdxi(2,2) * dxdxi(3,3) 
     &     - dxdxi(3,2) * dxdxi(2,3)
      dxidx(1,2) =   dxdxi(3,2) * dxdxi(1,3) 
     &     - dxdxi(1,2) * dxdxi(3,3)
      dxidx(1,3) =  dxdxi(1,2) * dxdxi(2,3) 
     &     - dxdxi(1,3) * dxdxi(2,2)
      tmp = 1d0 /
     &     ( dxidx(1,1) * dxdxi(1,1) 
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
      
      DetJ = 1d0/tmp           ! Note that DetJ resides in common

      shgradgu(:,1) =
     &     + shgradlu(:,1) * dxidx(1,1)
     &     + shgradlu(:,2) * dxidx(2,1)
     &     + shgradlu(:,3) * dxidx(3,1)
      shgradgu(:,2) =
     &     + shgradlu(:,1) * dxidx(1,2)
     &     + shgradlu(:,2) * dxidx(2,2)
     &     + shgradlu(:,3) * dxidx(3,2) 
      shgradgu(:,3) =
     &     + shgradlu(:,1) * dxidx(1,3)
     &     + shgradlu(:,2) * dxidx(2,3)
     &     + shgradlu(:,3) * dxidx(3,3) 
        
        
      icount = 0
      derv_sum_UU = 0d0
      derv_sum_UV = 0d0
      derv_sum_UW = 0d0
      derv_sum_VV = 0d0
      derv_sum_VW = 0d0
      derv_sum_WW = 0d0
      do k = 0, Ru
      do j = 0,Qu
      do i = 0,Pu
         icount = icount+1
              
c...  2nd derivatives
         tempshhessl(icount,1) = 
     &        Nnu(3,Pu+1-i)*Mmu(1,Qu+1-j)*Ou(1,Ru+1-k)*
     &        B_NETu(niu-i,nju-j,nku-k,NSD+1) ! ,uu
         derv_sum_UU = derv_sum_UU + tempshhessl(icount,1)
         tempshhessl(icount,2) = 
     &        Nnu(2,Pu+1-i)*Mmu(2,Qu+1-j)*Ou(1,Ru+1-k)*
     &        B_NETu(niu-i,nju-j,nku-k,NSD+1) ! ,uv
         derv_sum_UV = derv_sum_UV + tempshhessl(icount,2)
         tempshhessl(icount,3) = 
     &        Nnu(2,Pu+1-i)*Mmu(1,Qu+1-j)*Ou(2,Ru+1-k)*
     &        B_NETu(niu-i,nju-j,nku-k,NSD+1) ! ,uw
         derv_sum_UW = derv_sum_UW + tempshhessl(icount,3)
         tempshhessl(icount,4) = 
     &        Nnu(1,Pu+1-i)*Mmu(3,Qu+1-j)*Ou(1,Ru+1-k)*
     &        B_NETu(niu-i,nju-j,nku-k,NSD+1) ! ,vv
         derv_sum_VV = derv_sum_VV + tempshhessl(icount,4)
         tempshhessl(icount,5) = 
     &        Nnu(1,Pu+1-i)*Mmu(2,Qu+1-j)*Ou(2,Ru+1-k)*
     &        B_NETu(niu-i,nju-j,nku-k,NSD+1) ! ,vw
         derv_sum_VW = derv_sum_VW + tempshhessl(icount,5)
         tempshhessl(icount,6) = 
     &        Nnu(1,Pu+1-i)*Mmu(1,Qu+1-j)*Ou(3,Ru+1-k)*
     &        B_NETu(niu-i,nju-j,nku-k,NSD+1) ! ,ww
         derv_sum_WW = derv_sum_WW + tempshhessl(icount,6)
      enddo
      enddo
      enddo
        

c...  local Hessian
      ds  = 1d0/denom_sum
      ds2 = ds*ds
      ds3 = ds*ds2

      shhessl(:,1) =
     &     + tempshhessl(:,1)*ds
     &     - tempshgradl(:,1)*derv_sum_U*ds2
     &     - ((tempshgradl(:,1)*derv_sum_U+tempshl*derv_sum_UU)*ds2
     &     -   2d0*tempshl*derv_sum_U*derv_sum_U*ds3)
      shhessl(:,2) =
     &     + tempshhessl(:,2)*ds
     &     - tempshgradl(:,1)*derv_sum_V*ds2
     &     - ((tempshgradl(:,2)*derv_sum_U+tempshl*derv_sum_UV)*ds2
     &     -   2d0*tempshl*derv_sum_U*derv_sum_V*ds3)
      shhessl(:,3) =
     &     + tempshhessl(:,3)*ds
     &     - tempshgradl(:,1)*derv_sum_W*ds2
     &     - ((tempshgradl(:,3)*derv_sum_U+tempshl*derv_sum_UW)*ds2
     &     -   2d0*tempshl*derv_sum_U*derv_sum_W*ds3)
      shhessl(:,4) =
     &     + tempshhessl(:,4)*ds
     &     - tempshgradl(:,2)*derv_sum_V*ds2
     &     - ((tempshgradl(:,2)*derv_sum_V+tempshl*derv_sum_VV)*ds2
     &     -   2d0*tempshl*derv_sum_V*derv_sum_V*ds3)
      shhessl(:,5) =
     &     + tempshhessl(:,5)*ds
     &     - tempshgradl(:,2)*derv_sum_W*ds2
     &     - ((tempshgradl(:,3)*derv_sum_V+tempshl*derv_sum_VW)*ds2
     &     -   2d0*tempshl*derv_sum_V*derv_sum_W*ds3)
      shhessl(:,6) =
     &     + tempshhessl(:,6)*ds
     &     - tempshgradl(:,3)*derv_sum_W*ds2
     &     - ((tempshgradl(:,3)*derv_sum_W+tempshl*derv_sum_WW)*ds2
     &     -   2d0*tempshl*derv_sum_W*derv_sum_W*ds3)

c...    global Hessian

c...    Second derivatives of the geometrical map

      dxdxixj = 0d0
      icount = 0
        
      do k = 0, Ru
      do j = 0, Qu
      do i = 0, Pu
         icount = icount + 1
              
         dxdxixj(1,:) = dxdxixj(1,:)+B_NETu(niu-i,nju-j,nku-k,1)*
     &        shhessl(icount,:)
              
         dxdxixj(2,:) = dxdxixj(2,:)+B_NETu(niu-i,nju-j,nku-k,2)*
     &        shhessl(icount,:)
              
         dxdxixj(3,:) = dxdxixj(3,:)+B_NETu(niu-i,nju-j,nku-k,3)*
     &        shhessl(icount,:)

      enddo
      enddo
      enddo
        
c...  RHS of the matrix equation for the second derivatives of bases.
c...  Reuse local hess. array        

        
        
      shhessl(:,1) = shhessl(:,1) - shgradgu(:,1)*dxdxixj(1,1) -
     &     shgradgu(:,2)*dxdxixj(2,1) - shgradgu(:,3)*dxdxixj(3,1)
        
      shhessl(:,2) = shhessl(:,2) - shgradgu(:,1)*dxdxixj(1,2) -
     &     shgradgu(:,2)*dxdxixj(2,2) - shgradgu(:,3)*dxdxixj(3,2)
        
      shhessl(:,3) = shhessl(:,3) - shgradgu(:,1)*dxdxixj(1,3) -
     &     shgradgu(:,2)*dxdxixj(2,3) - shgradgu(:,3)*dxdxixj(3,3)
        
      shhessl(:,4) = shhessl(:,4) - shgradgu(:,1)*dxdxixj(1,4) -
     &     shgradgu(:,2)*dxdxixj(2,4) - shgradgu(:,3)*dxdxixj(3,4)
        
      shhessl(:,5) = shhessl(:,5) - shgradgu(:,1)*dxdxixj(1,5) -
     &     shgradgu(:,2)*dxdxixj(2,5) - shgradgu(:,3)*dxdxixj(3,5)
        
      shhessl(:,6) = shhessl(:,6) - shgradgu(:,1)*dxdxixj(1,6) -
     &     shgradgu(:,2)*dxdxixj(2,6) - shgradgu(:,3)*dxdxixj(3,6)
        
        
c...  LHS (6x6, same for every basis function)

      locLHS(1,1) = dxdxi(1,1)*dxdxi(1,1)
      locLHS(1,2) = 2d0*dxdxi(1,1)*dxdxi(2,1)
      locLHS(1,3) = 2d0*dxdxi(1,1)*dxdxi(3,1)
      locLHS(1,4) = dxdxi(2,1)*dxdxi(2,1)
      locLHS(1,5) = 2d0*dxdxi(2,1)*dxdxi(3,1)
      locLHS(1,6) = dxdxi(3,1)*dxdxi(3,1)

      locLHS(2,1) = dxdxi(1,1)*dxdxi(1,2)
      locLHS(2,2) = dxdxi(1,1)*dxdxi(2,2) + dxdxi(1,2)*dxdxi(2,1)
      locLHS(2,3) = dxdxi(1,1)*dxdxi(3,2) + dxdxi(1,2)*dxdxi(3,1)
      locLHS(2,4) = dxdxi(2,1)*dxdxi(2,2)
      locLHS(2,5) = dxdxi(2,1)*dxdxi(3,2) + dxdxi(2,2)*dxdxi(3,1)
      locLHS(2,6) = dxdxi(3,1)*dxdxi(3,2)

      locLHS(3,1) = dxdxi(1,1)*dxdxi(1,3)
      locLHS(3,2) = dxdxi(1,1)*dxdxi(2,3) + dxdxi(1,3)*dxdxi(2,1)
      locLHS(3,3) = dxdxi(1,1)*dxdxi(3,3) + dxdxi(1,3)*dxdxi(3,1)
      locLHS(3,4) = dxdxi(2,1)*dxdxi(2,3)
      locLHS(3,5) = dxdxi(2,1)*dxdxi(3,3) + dxdxi(2,3)*dxdxi(3,1)
      locLHS(3,6) = dxdxi(3,1)*dxdxi(3,3)

      locLHS(4,1) = dxdxi(1,2)*dxdxi(1,2)
      locLHS(4,2) = 2d0*dxdxi(1,2)*dxdxi(2,2)
      locLHS(4,3) = 2d0*dxdxi(1,2)*dxdxi(3,2)
      locLHS(4,4) = dxdxi(2,2)*dxdxi(2,2)
      locLHS(4,5) = 2d0*dxdxi(2,2)*dxdxi(3,2)
      locLHS(4,6) = dxdxi(3,2)*dxdxi(3,2)

      locLHS(5,1) = dxdxi(1,2)*dxdxi(1,3)
      locLHS(5,2) = dxdxi(1,2)*dxdxi(2,3) + dxdxi(1,3)*dxdxi(2,2)
      locLHS(5,3) = dxdxi(1,2)*dxdxi(3,3) + dxdxi(1,3)*dxdxi(3,2)
      locLHS(5,4) = dxdxi(2,2)*dxdxi(2,3)
      locLHS(5,5) = dxdxi(2,2)*dxdxi(3,3) + dxdxi(2,3)*dxdxi(3,2)
      locLHS(5,6) = dxdxi(3,2)*dxdxi(3,3)

      locLHS(6,1) = dxdxi(1,3)*dxdxi(1,3)
      locLHS(6,2) = 2d0*dxdxi(1,3)*dxdxi(2,3)
      locLHS(6,3) = 2d0*dxdxi(1,3)*dxdxi(3,3)
      locLHS(6,4) = dxdxi(2,3)*dxdxi(2,3)
      locLHS(6,5) = 2d0*dxdxi(2,3)*dxdxi(3,3)
      locLHS(6,6) = dxdxi(3,3)*dxdxi(3,3)

c...    (6x6) - Gaussian elimination

      tmprow = 0d0
        
      do k = 1,6

                                ! BEGIN PIVOT
         tmp = locLHS(k,k)
         cc = k
         do aa = k+1,6
            if (abs(tmp).le.abs(locLHS(aa,k))) then
               cc = aa          ! Store row number with bigger pivot
               tmp = locLHS(aa,k)
            endif
         enddo

         tmprow(k:6) = locLHS(k,k:6) ! Swap rows of matrix
         locLHS(k,k:6) = locLHS(cc,k:6)
         locLHS(cc,k:6) = tmprow(k:6) 

         do bb = 1, NSHLu
                                ! Swap RHS entries
            tmp = shhessl(bb,k)
            shhessl(bb,k) = shhessl(bb,cc)
            shhessl(bb,cc) = tmp

         enddo

                                ! END PIVOT
          
         do i = k+1,6
            tmp = locLHS(i,k)/locLHS(k,k)
            do j = k+1,6
               locLHS(i,j) = locLHS(i,j) - tmp*locLHS(k,j)
            enddo
            shhessl(:,i) = shhessl(:,i) - tmp*shhessl(:,k)
         enddo
      enddo
        
      do i = 6,1,-1
         do j = i+1,6
            shhessl(:,i) = shhessl(:,i) - locLHS(i,j)*shhessl(:,j)
         enddo
         shhessl(:,i) = shhessl(:,i)/locLHS(i,i)
      enddo

c...  Assign to global hessian of basis functions
        
      shhessg(:,1,1) = shhessl(:,1)
      shhessg(:,1,2) = shhessl(:,2)
      shhessg(:,1,3) = shhessl(:,3)
        
      shhessg(:,2,1) = shhessl(:,2)
      shhessg(:,2,2) = shhessl(:,4)
      shhessg(:,2,3) = shhessl(:,5)
        
      shhessg(:,3,1) = shhessl(:,3)
      shhessg(:,3,2) = shhessl(:,5)
      shhessg(:,3,3) = shhessl(:,6)

c$$$      lapgu = shhessg(:,1,1)+shhessg(:,2,2)+shhessg(:,3,3)
     
      dxidx(1,:) = dxidx(1,:)*2d0/du
      dxidx(2,:) = dxidx(2,:)*2d0/dv
      dxidx(3,:) = dxidx(3,:)*2d0/dw
      
      if (DetJ < 0d0) DetJ = -1d0*DetJ
      
      return
      end



!###########################################################################

      subroutine eval_SHAPE_3D_mean(e,u,v,w,shl,shgradl,shgradg)

      use aAdjKeep
      use common

      implicit none
      
c     --------------VARIABLE DECLARATIONS--------------------------------
c...  Element number
      integer e
c...  u and v coordinates of integration point in parent element
      real*8 u_hat, v_hat, w_hat
      
c...  Vector of Local basis function values at (u_hat, v_hat), local and 
c          global gradients.
      real*8 shl(NSHLu),shgradl(NSHLu,NSD),shgradg(NSHLu,NSD),
     &  dxdxi(NSD,NSD), dxidx(NSD,NSD)

c...  Local Variables
c    
c     1D nonrational basis functions and derivs in u and v
      real*8 N(2,Pu+1), M(2,Qu+1), O(2,Ru+1)
c     u and v coordinates of integration point, denominator and derivative sums
      real*8 u, v, w, denom_sum, derv_sum_U, derv_sum_V,
     &  derv_sum_W
c     NURBS coordinates, counters for loops
      integer ni, nj, nk, i, j, k, icount
c     temporary variables
      real*8 tmp, ds, ds2
c ------------------------------------------------------------------


c     Get NURBS coordinates for local node 1
      ni = INNu(IENu(e,1),1)
      nj = INNu(IENu(e,1),2)
      nk = INNu(IENu(e,1),3)
      
c     Evaluate 1D shape functions and derivatives each direction
      call dersbasisfuns(ni,Pu,MCPu,u,1,U_KNOTu,N) ! calculate in u direction
      call dersbasisfuns(nj,Qu,NCPu,v,1,V_KNOTu,M) ! calculate in v direction
      call dersbasisfuns(nk,Ru,OCPu,w,1,W_KNOTu,O) ! calculate in v direction

c     Form basis functions and derivatives dR/du and dR/dv
      icount     = 0
      denom_sum  = 0d0
      derv_sum_U = 0d0
      derv_sum_V = 0d0
      derv_sum_W = 0d0
      shgradl    = 0d0

      
      do k = 0, Ru
        do j = 0,Qu
          do i = 0,Pu
            icount = icount+1
            
c...        basis functions
            shl(icount) = N(1,Pu+1-i)*M(1,Qu+1-j)*O(1,Ru+1-k)*
     &           B_NETu(ni-i,nj-j,nk-k,NSD+1)
            denom_sum = denom_sum + shl(icount)
            
c...        derivatives
            shgradl(icount,1) = 
     &        N(2,Pu+1-i)*M(1,Qu+1-j)*O(1,Ru+1-k)*
     &        B_NETu(ni-i,nj-j,nk-k,NSD+1) ! u
            derv_sum_U = derv_sum_U + shgradl(icount,1)
            shgradl(icount,2) = 
     &        N(1,Pu+1-i)*M(2,Qu+1-j)*O(1,Ru+1-k)*
     &        B_NETu(ni-i,nj-j,nk-k,NSD+1) ! v
            derv_sum_V = derv_sum_V + shgradl(icount,2)
            shgradl(icount,3) = 
     &        N(1,Pu+1-i)*M(1,Qu+1-j)*O(2,Ru+1-k)*
     &        B_NETu(ni-i,nj-j,nk-k,NSD+1) ! w
            derv_sum_W = derv_sum_W + shgradl(icount,3)
            
          enddo
        enddo
      enddo
      
c     Divide through by denominator
      
      ds  = 1d0/denom_sum
      ds2 = ds*ds

      shgradl(:,1) = shgradl(:,1)*ds - (shl*derv_sum_U)*ds2
      shgradl(:,2) = shgradl(:,2)*ds - (shl*derv_sum_V)*ds2
      shgradl(:,3) = shgradl(:,3)*ds - (shl*derv_sum_W)*ds2
      shl = shl*ds
      
      
c...  Now calculate gradients.
      
c     calculate dx/dxi

      dxdxi  = 0d0
      icount = 0
      
      do k = 0, Ru
      do j = 0, Qu
      do i = 0, Pu
         icount = icount + 1
            
         dxdxi(1,1) = dxdxi(1,1) + B_NETu(ni-i,nj-j,nk-k,1) *
     &        shgradl(icount,1)
         dxdxi(1,2) = dxdxi(1,2) + B_NETu(ni-i,nj-j,nk-k,1) *
     &        shgradl(icount,2)
         dxdxi(1,3) = dxdxi(1,3) + B_NETu(ni-i,nj-j,nk-k,1) *
     &        shgradl(icount,3)
         dxdxi(2,1) = dxdxi(2,1) + B_NETu(ni-i,nj-j,nk-k,2) *
     &        shgradl(icount,1)
         dxdxi(2,2) = dxdxi(2,2) + B_NETu(ni-i,nj-j,nk-k,2) *
     &        shgradl(icount,2)
         dxdxi(2,3) = dxdxi(2,3) + B_NETu(ni-i,nj-j,nk-k,2) *
     &        shgradl(icount,3)
         dxdxi(3,1) = dxdxi(3,1) + B_NETu(ni-i,nj-j,nk-k,3) *
     &        shgradl(icount,1)
         dxdxi(3,2) = dxdxi(3,2) + B_NETu(ni-i,nj-j,nk-k,3) *
     &        shgradl(icount,2)
         dxdxi(3,3) = dxdxi(3,3) + B_NETu(ni-i,nj-j,nk-k,3) *
     &        shgradl(icount,3)
            
      enddo
      enddo
      enddo
      
      
c.... compute the inverse of deformation gradient
      
      dxidx(1,1) = 
     &     +  dxdxi(2,2) * dxdxi(3,3) 
     &     -  dxdxi(3,2) * dxdxi(2,3)
      dxidx(1,2) = 
     &     +  dxdxi(3,2) * dxdxi(1,3) 
     &     -  dxdxi(1,2) * dxdxi(3,3)
      dxidx(1,3) = 
     &     +  dxdxi(1,2) * dxdxi(2,3) 
     &     -  dxdxi(1,3) * dxdxi(2,2)
      tmp =
     &     + dxidx(1,1) * dxdxi(1,1)
     &     + dxidx(1,2) * dxdxi(2,1)
     &     + dxidx(1,3) * dxdxi(3,1)
      tmp = 1d0/tmp
      dxidx(1,1) = dxidx(1,1) * tmp
      dxidx(1,2) = dxidx(1,2) * tmp
      dxidx(1,3) = dxidx(1,3) * tmp
      dxidx(2,1) =
     &  + (dxdxi(2,3) * dxdxi(3,1) 
     &  -  dxdxi(2,1) * dxdxi(3,3)) * tmp
      dxidx(2,2) =
     &  + (dxdxi(1,1) * dxdxi(3,3) 
     &  -  dxdxi(3,1) * dxdxi(1,3)) * tmp
      dxidx(2,3) =
     &  + (dxdxi(2,1) * dxdxi(1,3) 
     &  -  dxdxi(1,1) * dxdxi(2,3)) * tmp
      dxidx(3,1) =
     &  + (dxdxi(2,1) * dxdxi(3,2) 
     &  -  dxdxi(2,2) * dxdxi(3,1)) * tmp
      dxidx(3,2) =
     &  + (dxdxi(3,1) * dxdxi(1,2) 
     &  -  dxdxi(1,1) * dxdxi(3,2)) * tmp
      dxidx(3,3) =
     &  + (dxdxi(1,1) * dxdxi(2,2) 
     &  -  dxdxi(1,2) * dxdxi(2,1)) * tmp
      
      DetJ = 1d0/tmp           ! Note that DetJ resides in common
      
      shgradg(:,1) =
     &     + shgradl(:,1) * dxidx(1,1) 
     &     + shgradl(:,2) * dxidx(2,1)
     &     + shgradl(:,3) * dxidx(3,1)
      shgradg(:,2) =
     &     + shgradl(:,1) * dxidx(1,2) 
     &     + shgradl(:,2) * dxidx(2,2)
     &     + shgradl(:,3) * dxidx(3,2) 
      shgradg(:,3) =
     &     + shgradl(:,1) * dxidx(1,3)
     &     + shgradl(:,2) * dxidx(2,3)
     &     + shgradl(:,3) * dxidx(3,3)      
      
      if (DetJ < 0d0) DetJ = -1d0*DetJ
      
      return
      end



!###########################################################################

      subroutine eval_SHAPE_3D_int(e,u_hat,v_hat,w_hat,shl,shgradl)

      use aAdjKeep
      use common

      implicit none
      
c     --------------VARIABLE DECLARATIONS--------------------------------
c...  Element number
      integer e
c...  u and v coordinates of integration point in parent element
      real*8 u_hat, v_hat, w_hat
      
c...  Vector of Local basis function values at (u_hat, v_hat), local and 
c          global gradients.
      real*8 shl(NSHLu),shgradl(NSHLu,NSD),
     &  dxdxi(NSD,NSD), dxidx(NSD,NSD)

c...  Local Variables
c    
c     1D nonrational basis functions and derivs in u and v
      real*8 N(2,Pu+1), M(2,Qu+1), O(2,Ru+1)
c     u and v coordinates of integration point, denominator and derivative sums
      real*8 u, v, w, 
     &     denom_sum,
     &     derv_sum_U, derv_sum_V, derv_sum_W
c     NURBS coordinates, counters for loops
      integer ni, nj, nk, i, j, k, icount
c     temporary variables
      real*8 tmp, ds, ds2
c ------------------------------------------------------------------


c     Get NURBS coordinates for local node 1
      ni = INNu(IENu(e,1),1)
      nj = INNu(IENu(e,1),2)
      nk = INNu(IENu(e,1),3)
      
c     Get u and v coordinates of integration point
      
      u = u_hat
      v = v_hat
      w = w_hat

c$$$      u = ( (U_KNOTu(ni+1) - U_KNOTu(ni))*u_hat
c$$$     &     + U_KNOTu(ni+1) + U_KNOTu(ni)  )*0.5d0
c$$$      v = ( (V_KNOTu(nj+1) - V_KNOTu(nj))*v_hat
c$$$     &     + V_KNOTu(nj+1) + V_KNOTu(nj)  )*0.5d0
c$$$      w = ( (W_KNOTu(nk+1) - W_KNOTu(nk))*w_hat
c$$$     &     + W_KNOTu(nk+1) + W_KNOTu(nk)  )*0.5d0

c     Evaluate 1D shape functions and derivatives each direction
      call dersbasisfuns(ni,Pu,MCPu,u,1,U_KNOTu,N) ! calculate in u direction
      call dersbasisfuns(nj,Qu,NCPu,v,1,V_KNOTu,M) ! calculate in v direction
      call dersbasisfuns(nk,Ru,OCPu,w,1,W_KNOTu,O) ! calculate in w direction

c     Form basis functions and derivatives dR/du and dR/dv
      icount     = 0
      denom_sum  = 0d0
      derv_sum_U = 0d0
      derv_sum_V = 0d0
      derv_sum_W = 0d0
      shgradl    = 0d0
      
      do k = 0, Ru
        do j = 0,Qu
          do i = 0,Pu
            icount = icount+1
            
c...        basis functions
            shl(icount) = N(1,Pu+1-i)*M(1,Qu+1-j)*O(1,Ru+1-k)*
     &           B_NETu(ni-i,nj-j,nk-k,NSD+1)
            denom_sum = denom_sum + shl(icount)
            
c...        derivatives
            shgradl(icount,1) = 
     &        N(2,Pu+1-i)*M(1,Qu+1-j)*O(1,Ru+1-k)*
     &        B_NETu(ni-i,nj-j,nk-k,NSD+1) ! u
            derv_sum_U = derv_sum_U + shgradl(icount,1)
            shgradl(icount,2) = 
     &        N(1,Pu+1-i)*M(2,Qu+1-j)*O(1,Ru+1-k)*
     &        B_NETu(ni-i,nj-j,nk-k,NSD+1) ! v
            derv_sum_V = derv_sum_V + shgradl(icount,2)
            shgradl(icount,3) = 
     &        N(1,Pu+1-i)*M(1,Qu+1-j)*O(2,Ru+1-k)*
     &        B_NETu(ni-i,nj-j,nk-k,NSD+1) ! w
            derv_sum_W = derv_sum_W + shgradl(icount,3)
            
          enddo
        enddo
      enddo
      
c     Divide through by denominator
      
      ds  = 1d0/denom_sum
      ds2 = ds*ds

      shgradl(:,1) = shgradl(:,1)*ds - (shl*derv_sum_U)*ds2
      shgradl(:,2) = shgradl(:,2)*ds - (shl*derv_sum_V)*ds2
      shgradl(:,3) = shgradl(:,3)*ds - (shl*derv_sum_W)*ds2
      shl = shl*ds
      
      
c...  Now calculate gradients.
      
c     calculate dx/dxi

      dxdxi  = 0d0
      icount = 0
      
      do k = 0, Ru
        do j = 0, Qu
          do i = 0, Pu
            icount = icount + 1
            
            dxdxi(1,1) = dxdxi(1,1) + B_NETu(ni-i,nj-j,nk-k,1) *
     &        shgradl(icount,1)
            dxdxi(1,2) = dxdxi(1,2) + B_NETu(ni-i,nj-j,nk-k,1) *
     &        shgradl(icount,2)
            dxdxi(1,3) = dxdxi(1,3) + B_NETu(ni-i,nj-j,nk-k,1) *
     &        shgradl(icount,3)
            dxdxi(2,1) = dxdxi(2,1) + B_NETu(ni-i,nj-j,nk-k,2) *
     &        shgradl(icount,1)
            dxdxi(2,2) = dxdxi(2,2) + B_NETu(ni-i,nj-j,nk-k,2) *
     &        shgradl(icount,2)
            dxdxi(2,3) = dxdxi(2,3) + B_NETu(ni-i,nj-j,nk-k,2) *
     &        shgradl(icount,3)
            dxdxi(3,1) = dxdxi(3,1) + B_NETu(ni-i,nj-j,nk-k,3) *
     &        shgradl(icount,1)
            dxdxi(3,2) = dxdxi(3,2) + B_NETu(ni-i,nj-j,nk-k,3) *
     &        shgradl(icount,2)
            dxdxi(3,3) = dxdxi(3,3) + B_NETu(ni-i,nj-j,nk-k,3) *
     &        shgradl(icount,3)
            
          enddo
        enddo
      enddo
      
      
c.... compute the inverse of deformation gradient
      
      dxidx(1,1) = 
     &     +  dxdxi(2,2) * dxdxi(3,3) 
     &     -  dxdxi(3,2) * dxdxi(2,3)
      dxidx(1,2) = 
     &     +  dxdxi(3,2) * dxdxi(1,3) 
     &     -  dxdxi(1,2) * dxdxi(3,3)
      dxidx(1,3) = 
     &     +  dxdxi(1,2) * dxdxi(2,3) 
     &     -  dxdxi(1,3) * dxdxi(2,2)
      tmp =
     &     + dxidx(1,1) * dxdxi(1,1)
     &     + dxidx(1,2) * dxdxi(2,1)
     &     + dxidx(1,3) * dxdxi(3,1)
      
      DetJ = tmp           ! Note that DetJ resides in common
      
      if (DetJ < 0d0) DetJ = -1d0*DetJ
      
      return
      end




!###########################################################################

      subroutine eval_SHAPE_3D_pres(e,u_hat,v_hat,w_hat,shl)

      use aAdjKeep
      use common

      implicit none
      
c     --------------VARIABLE DECLARATIONS--------------------------------
c...  Element number
      integer e, ep
c...  u and v coordinates of integration point in parent element
      real*8 u_hat, v_hat, w_hat
      
c...  Vector of Local basis function values at (u_hat, v_hat), local and 
c          global gradients.
      real*8 shl(NSHLp),shgradl(NSHLp,NSD),
     &  dxdxi(NSD,NSD), dxidx(NSD,NSD)

c...  Local Variables
c    
c     1D nonrational basis functions and derivs in u and v
      real*8 N(2,Pp+1), M(2,Qp+1), O(2,Rp+1)
c     u and v coordinates of integration point, denominator and derivative sums
      real*8 u, v, w, 
     &     denom_sum,
     &     derv_sum_U, derv_sum_V, derv_sum_W
c     NURBS coordinates, counters for loops
      integer ni, nj, nk, i, j, k, icount
c     temporary variables
      real*8 tmp, ds, ds2
c ------------------------------------------------------------------


c     Get NURBS coordinates for local node 1
      
      ni = INNp(IENp(e,1),1)
      nj = INNp(IENp(e,1),2)
      nk = INNp(IENp(e,1),3)

c     Get u and v coordinates of integration point
      
      u = ( (U_KNOTp(ni+1) - U_KNOTp(ni))*u_hat
     &     + U_KNOTp(ni+1) + U_KNOTp(ni)  )*0.5d0
      v = ( (V_KNOTp(nj+1) - V_KNOTp(nj))*v_hat
     &     + V_KNOTp(nj+1) + V_KNOTp(nj)  )*0.5d0
      w = ( (W_KNOTp(nk+1) - W_KNOTp(nk))*w_hat
     &     + W_KNOTp(nk+1) + W_KNOTp(nk)  )*0.5d0

      
c     Get u and v coordinates of integration point
      
c     Evaluate 1D shape functions and derivatives each direction
      call dersbasisfuns(ni,Pp,MCPp,u,1,U_KNOTp,N) ! calculate in u direction
      call dersbasisfuns(nj,Qp,NCPp,v,1,V_KNOTp,M) ! calculate in v direction
      call dersbasisfuns(nk,Rp,OCPp,w,1,W_KNOTp,O) ! calculate in v direction

c     Form basis functions and derivatives dR/du and dR/dv
      icount     = 0
      denom_sum  = 0d0
      derv_sum_U = 0d0
      derv_sum_V = 0d0
      derv_sum_W = 0d0
      shgradl    = 0d0

      
      do k = 0, Rp
        do j = 0, Qp
          do i = 0, Pp
            icount = icount+1
            
c...        basis functions
            shl(icount) = N(1,Pp+1-i)*M(1,Qp+1-j)*O(1,Rp+1-k)*
     &           B_NETp(ni-i,nj-j,nk-k,NSD+1)
            denom_sum = denom_sum + shl(icount)
            
c...        derivatives
            shgradl(icount,1) = 
     &        N(2,Pp+1-i)*M(1,Qp+1-j)*O(1,Rp+1-k)*
     &        B_NETp(ni-i,nj-j,nk-k,NSD+1) ! u
            derv_sum_U = derv_sum_U + shgradl(icount,1)
            shgradl(icount,2) = 
     &        N(1,Pp+1-i)*M(2,Qp+1-j)*O(1,Rp+1-k)*
     &        B_NETp(ni-i,nj-j,nk-k,NSD+1) ! v
            derv_sum_V = derv_sum_V + shgradl(icount,2)
            shgradl(icount,3) = 
     &        N(1,Pp+1-i)*M(1,Qp+1-j)*O(2,Rp+1-k)*
     &        B_NETp(ni-i,nj-j,nk-k,NSD+1) ! w
            derv_sum_W = derv_sum_W + shgradl(icount,3)
            
          enddo
        enddo
      enddo
      
c     Divide through by denominator
      
      ds  = 1d0/denom_sum
      ds2 = ds*ds

      shgradl(:,1) = shgradl(:,1)*ds - (shl*derv_sum_U)*ds2
      shgradl(:,2) = shgradl(:,2)*ds - (shl*derv_sum_V)*ds2
      shgradl(:,3) = shgradl(:,3)*ds - (shl*derv_sum_W)*ds2
      shl = shl*ds
      
      
c...  Now calculate gradients.
      
c     calculate dx/dxi

      dxdxi  = 0d0
      icount = 0
      
      do k = 0, Rp
        do j = 0, Qp
          do i = 0, Pp
            icount = icount + 1
            
            dxdxi(1,1) = dxdxi(1,1) + B_NETp(ni-i,nj-j,nk-k,1) *
     &        shgradl(icount,1)
            dxdxi(1,2) = dxdxi(1,2) + B_NETp(ni-i,nj-j,nk-k,1) *
     &        shgradl(icount,2)
            dxdxi(1,3) = dxdxi(1,3) + B_NETp(ni-i,nj-j,nk-k,1) *
     &        shgradl(icount,3)
            dxdxi(2,1) = dxdxi(2,1) + B_NETp(ni-i,nj-j,nk-k,2) *
     &        shgradl(icount,1)
            dxdxi(2,2) = dxdxi(2,2) + B_NETp(ni-i,nj-j,nk-k,2) *
     &        shgradl(icount,2)
            dxdxi(2,3) = dxdxi(2,3) + B_NETp(ni-i,nj-j,nk-k,2) *
     &        shgradl(icount,3)
            dxdxi(3,1) = dxdxi(3,1) + B_NETp(ni-i,nj-j,nk-k,3) *
     &        shgradl(icount,1)
            dxdxi(3,2) = dxdxi(3,2) + B_NETp(ni-i,nj-j,nk-k,3) *
     &        shgradl(icount,2)
            dxdxi(3,3) = dxdxi(3,3) + B_NETp(ni-i,nj-j,nk-k,3) *
     &        shgradl(icount,3)
            
          enddo
        enddo
      enddo
      
      
c.... compute the inverse of deformation gradient
      
      dxidx(1,1) = 
     &     +  dxdxi(2,2) * dxdxi(3,3) 
     &     -  dxdxi(3,2) * dxdxi(2,3)
      dxidx(1,2) = 
     &     +  dxdxi(3,2) * dxdxi(1,3) 
     &     -  dxdxi(1,2) * dxdxi(3,3)
      dxidx(1,3) = 
     &     +  dxdxi(1,2) * dxdxi(2,3) 
     &     -  dxdxi(1,3) * dxdxi(2,2)
      tmp =
     &     + dxidx(1,1) * dxdxi(1,1)
     &     + dxidx(1,2) * dxdxi(2,1)
     &     + dxidx(1,3) * dxdxi(3,1)
      
      DetJ = tmp           ! Note that DetJ resides in common
      
      if (DetJ < 0d0) DetJ = -1d0*DetJ
      
      return
      end



