      subroutine SparseProdK_scalar (lhsM, colvec, rowvec,
     &     rhstmp, prodtmp, NNODZ, Pu, Qu, Ru, icntu)
      implicit none

      integer aa, bb, cc, NNODZ, Pu, Qu, Ru, icntu, NSD
      integer colvec(NNODZ+1), rowvec(NNODZ*8*(Pu+1)*(Qu+1)*(Ru+1))
      real*8 rhstmp(NNODZ), prodtmp(NNODZ)
      real*8 tmp1, tmp2, tmp3, pisave
      real*8 lhsM(icntu)
      
c.... clear the vector
      prodtmp = 0d0
      
      do aa = 1, NNODZ
         
         tmp1 = 0d+0
         
         do bb = colvec(aa), colvec(aa+1)-1
            cc = rowvec(bb)
            
            tmp1 = tmp1 + LHSM(bb)*rhstmp(cc)
            
         enddo
         
         prodtmp(aa) = prodtmp(aa) + tmp1
         
      enddo
      
      return
      end
      
!#########################################################################

      subroutine SparseProdFull_3D (lhsK, 
     &     colu, rowu,
     &     rhstmpu, prodtmpu,
     &     NNODZu, Pu, Qu, Ru,
     &     NDOF, icntu)
      implicit none

      integer aa, bb, cc, NNODZu, Pu, Qu, Ru, icntu, NDOF
      integer colu(NNODZu+1), rowu(NNODZu*8*(Pu+1)*(Qu+1)*(Ru+1))
      real*8 rhstmpu(NNODZu,NDOF), prodtmpu(NNODZu,NDOF)
      real*8 tmp1, tmp2, tmp3, tmp4, tmp5, pisave
      real*8 lhsK(NDOF*NDOF,icntu) 
      
c.... clear the vector
      
      prodtmpu = 0d+0
      
      do aa = 1, NNODZu         ! K*u
         
         tmp1 = 0d+0
         tmp2 = 0d+0
         tmp3 = 0d+0
         tmp4 = 0d+0
         tmp5 = 0d+0
         
         do bb = colu(aa), colu(aa+1)-1
            cc = rowu(bb)

            tmp1 = tmp1 + LHSK(1,bb)*rhstmpu(cc,1) + 
     &           LHSK(2,bb)*rhstmpu(cc,2) + LHSK(3,bb)*rhstmpu(cc,3) +
     &           LHSK(4,bb)*rhstmpu(cc,4) + LHSK(5,bb)*rhstmpu(cc,5)
            tmp2 = tmp2 + LHSK(6,bb)*rhstmpu(cc,1) +
     &           LHSK(7,bb)*rhstmpu(cc,2) + LHSK(8,bb)*rhstmpu(cc,3) +
     &           LHSK(9,bb)*rhstmpu(cc,4) + LHSK(10,bb)*rhstmpu(cc,5)
            tmp3 = tmp3 + LHSK(11,bb)*rhstmpu(cc,1) +
     &           LHSK(12,bb)*rhstmpu(cc,2) + LHSK(13,bb)*rhstmpu(cc,3) +
     &           LHSK(14,bb)*rhstmpu(cc,4) + LHSK(15,bb)*rhstmpu(cc,5)
            tmp4 = tmp4 + LHSK(16,bb)*rhstmpu(cc,1) +
     &           LHSK(17,bb)*rhstmpu(cc,2) + LHSK(18,bb)*rhstmpu(cc,3) +
     &           LHSK(19,bb)*rhstmpu(cc,4) + LHSK(20,bb)*rhstmpu(cc,5)
            tmp5 = tmp5 + LHSK(21,bb)*rhstmpu(cc,1) +
     &           LHSK(22,bb)*rhstmpu(cc,2) + LHSK(23,bb)*rhstmpu(cc,3) +
     &           LHSK(24,bb)*rhstmpu(cc,4) + LHSK(25,bb)*rhstmpu(cc,5)
         enddo
         
         prodtmpu(aa,1) = prodtmpu(aa,1) + tmp1
         prodtmpu(aa,2) = prodtmpu(aa,2) + tmp2
         prodtmpu(aa,3) = prodtmpu(aa,3) + tmp3
         prodtmpu(aa,4) = prodtmpu(aa,4) + tmp4
         prodtmpu(aa,5) = prodtmpu(aa,5) + tmp5

      enddo

      do aa = 1, NNODZu         ! G*p
         
         tmp1 = 0d+0
         tmp2 = 0d+0
         tmp3 = 0d+0

c$$$         do bb = colp(aa), colp(aa+1)-1
c$$$            cc = rowp(bb)
c$$$            
c$$$            tmp1 = tmp1 + LHSP(1,bb)*rhstmpp(cc)
c$$$            tmp2 = tmp2 + LHSP(2,bb)*rhstmpp(cc)
c$$$            tmp3 = tmp3 + LHSP(3,bb)*rhstmpp(cc)
c$$$            
c$$$         enddo
         
         prodtmpu(aa,1) = prodtmpu(aa,1) + tmp1
         prodtmpu(aa,2) = prodtmpu(aa,2) + tmp2
         prodtmpu(aa,3) = prodtmpu(aa,3) + tmp3

      enddo

      do aa = 1, NNODZu         ! G^t*u
         
         tmp1 = rhstmpu(aa,1)
         tmp2 = rhstmpu(aa,2)
         tmp3 = rhstmpu(aa,3)

c$$$         do bb = colp(aa), colp(aa+1)-1
c$$$            cc = rowp(bb)
c$$$            
c$$$            prodtmpp(cc) = prodtmpp(cc) + LHSPt(1,bb)*tmp1
c$$$     &           + LHSPt(2,bb)*tmp2 + LHSPt(3,bb)*tmp3
c$$$            
c$$$         enddo
         
      enddo

c$$$      do aa = 1, NNODZp         ! M*p
c$$$         
c$$$         tmp1 = 0d+0
c$$$         
c$$$         do bb = colpp(aa), colpp(aa+1)-1
c$$$            cc = rowpp(bb)
c$$$            
c$$$            tmp1 = tmp1 + LHSM(bb)*rhstmpp(cc)
c$$$            
c$$$         enddo
c$$$         
c$$$         prodtmpp(aa) = prodtmpp(aa) + tmp1
c$$$         
c$$$      enddo
      
      return
      end

      
