
      subroutine FillSparseMat_3D_u(iel, col, row, xKebe) 
      
      use aAdjKeep
      
      use common
      implicit none
      
      integer iel, row(nnodzu*8*(Pu+1)*(Qu+1)*(Ru+1)), col(nnodzu+1)
      
      integer a, b, c, d, ee, n, k, locn, i, idof
      
      real*8 xKebe(NDOF*NDOF,NSHLu,NSHLu)
      
      do a = 1, NSHLu
         i = IENu(iel,a)
         c = col(i)
         n = col(i+1) - c
         do b = 1, NSHLu
            
            call SparseMatLoc_3D (row(c), n, IENu(iel,b), locn) 
            
            k = locn + c - 1
            
            do idof = 1, NDOF*NDOF
             LHSK(idof,k) = LHSK(idof,k) + xKebe(idof,a,b)
            enddo
            
         enddo
         
      enddo
      
c      print*, 'aqui ok', lhsk(1)
c      print*, 'xKebe', xKebe(1,1)

      return
      end
