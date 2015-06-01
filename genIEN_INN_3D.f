c       This subroutine generates the IEN matrix, which relates element numbers
c       and local node numbers to the appropriate global node numbers. The
c       routine also generates the INN matrix, which relates global node
c       number to the "NURBS coordinates" of the node. This routine is for
c       the 3D code


      subroutine genIEN_INN_3D_u


      use aAdjKeep		! this f90 module contains all allocatables
c     
      use common		! common file defines all variables

!      use mpi

      implicit none

c     
c...  Local variables
      integer :: i, j, k, loop1,loop2,loop3, g, e, gtemp, ln
      
      allocate(INNu(NNODZu,3))
      allocate(IENu(NELu,(Pu+1)*(Qu+1)*(Ru+1)))

c...  Initialize matrices and variables
      IENu = 0
      INNu = 0
      g = 0
      e = 0


c...  Loop through control points assigning global node
c     numbers and filling out IEN and INN as we go
      do k = 1,OCPu		! loop through control points in W direction
         do j = 1,NCPu          ! loop through control points in V direction
            do i = 1,MCPu	! loop through control points in U direction
               g = g+1
               INNu(g,1) = i
               INNu(g,2) = j
               INNu(g,3) = k

c$$$                write(1000+myid,'(a,20i7)')
c$$$     &              'INNu',INNu(g,:),g

               if (((i .ge. (Pu+1)).and.
     &              (j .ge. (Qu+1))).and.
     &              (k .ge. (Ru+1))) then
                  e = e + 1
                  do loop1 = 0,Ru
                     do loop2 = 0,Qu
                        do loop3 = 0,Pu
                           gtemp = g 
     &                          - loop1*MCPu*NCPu
     &                          - loop2*MCPu
     &                          - loop3
                           ln = + loop1*(Pu+1)*(Qu+1)
     &                          + loop2*(Pu+1)
     &                          + loop3
     &                          + 1
                           IENu(e,ln) = gtemp
                        enddo
                     enddo
                  enddo

c$$$                  write(1000+myid,'(a,20i7)')'e, IENu = ',e,IENu(e,:)

               endif
            enddo
         enddo
      enddo	      
      
      return
      end


!#########################################################################

      subroutine genIEN_INN_3D_p


      use aAdjKeep		! this f90 module contains all allocatables
c     
      use common                ! common file defines all variables
c     
c...  Local variables
      integer :: i, j, k, loop1,loop2,loop3,g, e, gtemp, ln
      
      allocate(INNp(NNODZp,3))
      allocate(IENp(NELp,(Pp+1)*(Qp+1)*(Rp+1)))
c     
c...  Initialize matrices and variables
c     
      IENp = 0
      INNp = 0
      g = 0
      e = 0

c     
c...  Loop through control points assigning global node
c     numbers and filling out IEN and INN as we go
      do k = 1,OCPp             ! loop through control points in W direction
         do j = 1,NCPp		! loop through control points in V direction
            do i = 1,MCPp	! loop through control points in U direction
               g = g+1
               INNp(g,1) = i
               INNp(g,2) = j
               INNp(g,3) = k
               if (((i .ge. (Pp+1)).and.(j .ge. (Qp+1))).and.
     &              (k.ge.(Rp+1))) then
                  e = e +1
                  do loop1 = 0,Rp
                     do loop2 = 0,Qp
                        do loop3 = 0,Pp
                           gtemp = g - loop1*MCPp*NCPp -MCPp*loop2-
     &                          loop3
                           ln = (Pp+1)*(Qp+1)*loop1+(Pp+1)*loop2 + 
     &                          loop3 + 1
                           IENp(e,ln) = gtemp
                        enddo
                     enddo
                  enddo
               endif
            enddo
         enddo
      enddo	      
      
      return
      end
