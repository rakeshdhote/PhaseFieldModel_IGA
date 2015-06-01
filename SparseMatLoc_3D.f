c
c modified from original code of Jurijs Bazilevs, who stole it from someone else
c  (this code is the village bicycle)
c
c July 1, 2003 
c J. Austin Cottrell, III
c
c renamed SparseMatLoc_3D to be in keeping with 3D suffix convention, but no
c    changes made.


      subroutine SparseMatLoc_3D( list, n, target, locat )
      implicit none

      integer	n
      integer	list(n),	target
      integer	rowvl,	rowvh,	rowv, locat
c     
c.... Initialize
c     
      rowvl = 1
      rowvh = n + 1
c     
c.... do a binary search
c     
 100  if ( rowvh-rowvl .gt. 1 ) then
         rowv = ( rowvh + rowvl ) / 2
         if ( list(rowv) .gt. target ) then
            rowvh = rowv
         else
            rowvl = rowv
         endif
         goto 100
      endif
c     
c.... return
c     
      locat = rowvl
c     
      return
      end
