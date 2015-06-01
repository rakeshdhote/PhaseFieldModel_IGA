      subroutine BCLhs_3D(iel, xKebe, Rhs)
      
      use aAdjKeep
      use common

      implicit none
      
      integer iel, aa, bb, cc
      real*8 xKebe(NDOF*NDOF,NSHLu,NSHLu), Rhs(NDOF,NSHLu)
      

c...  loop through local nodes
         
      do aa = 1, NSHLu
         
         cc = IENu(iel,aa)      ! global node number of local node aa
         ! iel - total # of elements
         
         if ((IBC(cc,1).eq.1).or.(IBC(cc,1).eq.2)) then !----------------------------
                                                              !
            do bb = 1, NSHLu    ! project away row (u1)    !
                                                              !
               xKebe(1,aa,bb) = 0d+0                          !
               xKebe(2,aa,bb) = 0d+0                          !
               xKebe(3,aa,bb) = 0d+0                          !
               xKebe(4,aa,bb) = 0d+0                          !
               xKebe(5,aa,bb) = 0d+0                          !
                                                              !
            enddo                                             !
                                                              !
            do bb = 1, NSHLu    !project away column (u1)  !
                                                              !
               xKebe(1,bb,aa) = 0d+0                          !
               xKebe(6,bb,aa) = 0d+0                          !
               xKebe(11,bb,aa) = 0d+0                          !
               xKebe(16,bb,aa) = 0d+0                         !
               xKebe(21,bb,aa) = 0d+0                          !
                                                              !                                                                             !
            enddo                                             !
                                                              !
            xKebe(1,aa,aa) = 1d+0                             !
            Rhs(1,aa) = 0d+0    ! rhs                         !
                                                              !
         endif  ! cc = 1                    !-----------------------
         
         if ((IBC(cc,2).eq.1).or.(IBC(cc,2).eq.2)) then
            
            do bb = 1, NSHLu    ! project away row (u2)
               
               xKebe(6,aa,bb) = 0d+0
               xKebe(7,aa,bb) = 0d+0
               xKebe(8,aa,bb) = 0d+0
               xKebe(9,aa,bb) = 0d+0
               xKebe(10,aa,bb) = 0d+0
               
            enddo
            
            do bb = 1, NSHLu    !project away column (u2)
               
               xKebe(2,bb,aa) = 0d+0
               xKebe(7,bb,aa) = 0d+0
               xKebe(12,bb,aa) = 0d+0
               xKebe(17,bb,aa) = 0d+0
               xKebe(22,bb,aa) = 0d+0
               
            enddo
            
            xKebe(7,aa,aa) = 1d+0
            Rhs(2,aa) = 0d+0    ! rhs
            
         endif  ! cc = 2
         
         
         if ((IBC(cc,3).eq.1).or.(IBC(cc,3).eq.2)) then
            
            do bb = 1, NSHLu    ! project away row (v1)
               
               xKebe(11,aa,bb) = 0d+0
               xKebe(12,aa,bb) = 0d+0
               xKebe(13,aa,bb) = 0d+0
               xKebe(14,aa,bb) = 0d+0
               xKebe(15,aa,bb) = 0d+0
               
            enddo
            
            do bb = 1, NSHLu    !project away column (v1)
               
               xKebe(3,bb,aa) = 0d+0
               xKebe(8,bb,aa) = 0d+0
               xKebe(13,bb,aa) = 0d+0
               xKebe(18,bb,aa) = 0d+0
               xKebe(23,bb,aa) = 0d+0
               
            enddo
            
            xKebe(13,aa,aa) = 1d+0
            Rhs(3,aa) = 0d+0    ! rhs
            
         endif  ! cc = 3
         
         
         if ((IBC(cc,4).eq.1).or.(IBC(cc,4).eq.2)) then
            
            do bb = 1, NSHLu    ! project away row (v2)
               
               xKebe(16,aa,bb) = 0d+0
               xKebe(17,aa,bb) = 0d+0
               xKebe(18,aa,bb) = 0d+0
               xKebe(19,aa,bb) = 0d+0
               xKebe(20,aa,bb) = 0d+0
               
            enddo
               
            do bb = 1, NSHLu    !project away column (v2)
               
               xKebe(4,bb,aa) = 0d+0
               xKebe(9,bb,aa) = 0d+0
               xKebe(14,bb,aa) = 0d+0
               xKebe(19,bb,aa) = 0d+0
               xKebe(24,bb,aa) = 0d+0
               
            enddo
            
            xKebe(19,aa,aa) = 1d+0
            Rhs(4,aa) = 0d+0    ! rhs
            
         endif   ! cc = 4

         if ((IBC(cc,5).eq.1).or.(IBC(cc,5).eq.2)) then

            do bb = 1, NSHLu    ! project away row (v2)

               xKebe(21,aa,bb) = 0d+0
               xKebe(22,aa,bb) = 0d+0
               xKebe(23,aa,bb) = 0d+0
               xKebe(24,aa,bb) = 0d+0
               xKebe(25,aa,bb) = 0d+0

            enddo

            do bb = 1, NSHLu    !project away column (v2)

               xKebe(5,bb,aa) = 0d+0
               xKebe(10,bb,aa) = 0d+0
               xKebe(15,bb,aa) = 0d+0
               xKebe(20,bb,aa) = 0d+0
               xKebe(25,bb,aa) = 0d+0

            enddo

            xKebe(25,aa,aa) = 1d+0
            Rhs(5,aa) = 0d+0    ! rhs

         endif   ! cc = 5
            
      enddo ! Over all local nodes
      
      return
      end
      
