
      subroutine e3STAB_3D(dxidx, duidxi, ui, tauM)
      
      use aAdjKeep
      
      use common
      
      implicit none

      real*8 dxidx(NSD,NSD), tauM, ui, uiold
      real*8 sigma, forc, forcold, temp(NSD,NSD)
      real*8 duidxi(NSD), duidxiold(NSD)
      integer i,j,k,l

c...  build remaining diagonal entries of tauM  
    
      temp = 0d0
      sigma = 0d0
      forc  = 0d0   ! quadratic form (grad ^ T c) g (grad c)
      
      do i = 1, NSD             ! gij
         do j = 1, NSD
            temp(i,j) = 
     &           + dxidx(1,i)*dxidx(1,j)
     &           + dxidx(2,i)*dxidx(2,j)
     &           + dxidx(3,i)*dxidx(3,j)
     &           + dxidx(4,i)*dxidx(4,j)
         enddo
      enddo

      do i = 1, NSD
         do j = 1, NSD
            forc  = forc + duidxi(i)*temp(i,j)*duidxi(j)
            forcold  = forcold + duidxiold(i)*temp(i,j)*duidxiold(j)
            sigma = sigma + temp(i,j)*temp(i,j)
         enddo
      enddo
      
      sigma = sigma/4d1

c... Build tauM

      tauM = 0d0

      return
      end

      
!########################################################################
      
      subroutine e3bSTAB_3D(kappa_l, nor, dxidx, kapoh)
      
      use aAdjKeep
      
      use common
      
      implicit none

      real*8 dxidx(NSD,NSD), kappa_l, kapoh
      real*8 m_k, temp1(NSD), kappa_n(NSD), nor(NSD)
      integer i,j,k,l, aa, bb
      
      
      kappa_n = 0d0
      
      do aa = 1, NSD
        kappa_n(aa) = kappa_l*nor(aa)
      enddo
      
      temp1 = 0d0
      
      do i = 1, NSD
        temp1(1) = temp1(1)+dxidx(1,i)*kappa_n(i)
        temp1(2) = temp1(2)+dxidx(2,i)*kappa_n(i)
        temp1(3) = temp1(3)+dxidx(3,i)*kappa_n(i)
      enddo
      
      m_k = max(Pu**2,Qu**2,Ru**2)*2d0 ! Polynomial order dependent
      
      kapoh = m_k*sqrt(temp1(1)*temp1(1) +
     &  temp1(2)*temp1(2) +
     &  temp1(3)*temp1(3))


      return
      end
      

      
      
