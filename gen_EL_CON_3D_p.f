c     This routine links each element in the p-refined mesh with its
c        parent element in the unrefined mesh.
c
c     EL_CON(i) = j implies that element j of the OLD mesh is the parent
c           of element i in the new mesh.
c
c      
c     May 02, 2004
c       J. Austin Cottrell
c       CAM Graduate Student
c       Institute for Computational Engineering Science
c       The University of Texas at Austin

      
      subroutine gen_EL_CON_3D_p(new_NEL,new_MCP,new_NCP,
     &     new_OCP,new_P,new_Q,new_R, MCP, NCP, OCP, P, Q, R)
      
      use aAdjKeep
      use common
      implicit none

      integer e,new_NEL,i,j,k,l,ll,lll,icount,jcount,new_MCP,new_NCP
      integer new_OCP,kcount,tu,tv,tw
      integer new_P,new_Q,new_R, P,Q,R,MCP,NCP,OCP
      
      tu = new_P-P
      tv = new_Q-Q
      tw = new_R-R
      
      EL_CON = 0
      e = 0
      kcount = 0
      do k = 1,OCP-R
        jcount = 0
        do j = 1,NCP-Q
          icount = 0
          e = (MCP-P)*(NCP-Q)*(k-1)+(MCP-P)*(j-1)
          do i = 1,MCP-P
            e = e+1
            EL_CON(i+icount+(j-1+jcount)*(new_MCP-new_P)+
     &        (k-1+kcount)*(new_MCP-new_P)*(new_NCP-new_Q)) = e
            if (((U_KNOTp(i+P).ne.U_KNOTp(i+P+1)).and.(tu.gt.0)).and.
     &        (i.lt.(MCP-P))) then
              do l = 1,tu
                icount = icount + 1
                EL_CON(i+icount+(j-1+jcount)*(new_MCP-new_P)+
     &            (k-1+kcount)*(new_MCP-new_P)*(new_NCP-new_Q)) = e
              enddo
            endif
          enddo
          if (((V_KNOTp(j+Q).ne.V_KNOTp(j+Q+1)).and.(tv.gt.0)).and.
     &      (j.lt.(NCP-Q))) then
            do l = 1,tv
              icount = 0
              jcount = jcount+1
              e = (MCP-P)*(NCP-Q)*(k-1)+
     &          (j-1)*(MCP-P)
              do i = 1,MCP-P
                e = e+1
                EL_CON(i+icount+(j-1+jcount)*(new_MCP-new_P)+
     &            (k-1+kcount)*(new_MCP-new_P)*(new_NCP-new_Q)) = e
               if (((U_KNOTp(i+P).ne.U_KNOTp(i+P+1)).and.(tu.gt.0)).and.
     &            (i.lt.(MCP-P))) then
                  do ll = 1,tu
                    icount = icount + 1
                    EL_CON(i+icount+(j-1+jcount)*(new_MCP-new_P)+
     &                (k-1+kcount)*(new_MCP-new_P)*(new_NCP-new_Q)) = e
                  enddo
                endif
              enddo
            enddo
          endif
        enddo
        if (((W_KNOTp(k+R).ne.W_KNOTp(k+R+1)).and.(tw.gt.0)).and.
     &    (k.lt.(OCP-R))) then
          do l = 1,tw
            kcount = kcount+1
            jcount = 0
            do j = 1,NCP-Q
              icount = 0
              e = (MCP-P)*(NCP-Q)*(k-1)+(MCP-P)*(j-1)
              do i = 1,MCP-P
                e = e+1
                EL_CON(i+icount+(j-1+jcount)*(new_MCP-new_P)+
     &            (k-1+kcount)*(new_MCP-new_P)*(new_NCP-new_Q)) = e
               if (((U_KNOTp(i+P).ne.U_KNOTp(i+P+1)).and.(tu.gt.0)).and.
     &            (i.lt.(MCP-P))) then
                  do ll = 1,tu
                    icount = icount + 1
                    EL_CON(i+icount+(j-1+jcount)*(new_MCP-new_P)+
     &                (k-1+kcount)*(new_MCP-new_P)*(new_NCP-new_Q)) = e
                  enddo
                endif
              enddo
              if (((V_KNOTp(j+Q).ne.V_KNOTp(j+Q+1)).and.(tv.gt.0)).and.
     &          (j.lt.(NCP-Q))) then
                do ll = 1,tv
                  icount = 0
                  jcount = jcount+1
                  e = (MCP-P)*(NCP-Q)*(k-1)+
     &              (j-1)*(MCP-P)
                  do i = 1,MCP-P
                    e = e+1
                    EL_CON(i+icount+(j-1+jcount)*(new_MCP-new_P)+
     &                (k-1+kcount)*(new_MCP-new_P)*(new_NCP-new_Q)) = e
                    if (((U_KNOTp(i+P).ne.U_KNOTp(i+P+1)).and.(tu.gt.0))
     &                .and.(i.lt.(MCP-P))) then
                      do lll = 1,tu
                        icount = icount + 1
                        EL_CON(i+icount+(j-1+jcount)*(new_MCP-new_P)+
     &                    (k-1+kcount)*(new_MCP-new_P)*(new_NCP-new_Q))
     &                    = e
                      enddo
                    endif
                  enddo
                enddo
              endif
            enddo
          enddo
        endif
      enddo
      
      
      return
      end
