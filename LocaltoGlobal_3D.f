c     This routine takes the local load vector and places it into the
c     global vector.

      subroutine LocaltoGlobal_3D_u (iel, Rhs)
      
      use aAdjKeep
      
      use common
      implicit none

      integer iel
      integer aa, bb, idof
      real*8 Rhs(NDOF,NSHLu)
 
      do aa = 1, NDOF     
       do bb = 1, NSHLu
            
         RHSGu(IENu(iel,bb),aa) = RHSGu(IENu(iel,bb),aa) + Rhs(aa,bb)

       enddo           
      enddo
      
      
      return
      end
