c     This routine links each face in a p-refined mesh with its parent face
c        in the unrefined mesh.
c
c     FACE_CON(i) = j indicates that element j in the OLD mesh is the parent
c        of element i in the new mesh.
c
c      
c       May 01, 2004
c       J. Austin Cottrell
c       CAM Graduate Student
c       Institute for Computational Engineering Science
c       The University of Texas at Austin        

      subroutine gen_FACE_CON_3D_p(new_NFACE,new_MCP,new_NCP,
     &     new_OCP,new_P,new_Q,new_R, MCP, NCP, OCP, P, Q, R)

      use aAdjKeep
      use common
      implicit none

      integer e, icount,i,j,k,l,jcount,new_MCP,new_NCP,new_NFACE,
     &  new_OCP, kcount,
     &  tu,tv,tw,dummy,new_P,new_Q,new_R,P,Q,R,MCP,NCP,OCP


      tu = new_P-P
      tv = new_Q-Q
      tw = new_R-R
      
      e = 0    
      
      if (CLOSED_W_flag.ne.1) then
c...    Loop through orientation 1
        
        jcount = 0
        do j= 1,NCP-Q
          icount = 0
          do i = 1,MCP-P
            e = e+1
            FACE_CON(i+icount + (new_MCP-new_P)*(j-1+jcount)) = e
          if (((U_KNOTp(i+P).ne.U_KNOTp(i+P+1)).and.(tu.gt.0)).and.
     &        (i.lt.(MCP-P))) then
              do k = 1,tu
                icount = icount+1
                FACE_CON(i+icount + (new_MCP-new_P)*(j-1+jcount)) = e
              enddo
          endif
          enddo         
        if (((V_KNOTp(j+Q).ne.V_KNOTp(j+1+Q)).and.(tv.gt.0)).and.
     &      (j.lt.(NCP-Q))) then
            do l = 1,tv
              jcount = jcount+1
              icount = 0
              e = (j-1)*(MCP-P)
              do i = 1,MCP-P
                e = e+1
                FACE_CON(i+icount + (new_MCP-new_P)*(j-1+jcount)) = e
              if (((U_KNOTp(i+P).ne.U_KNOTp(i+P+1)).and.(tu.gt.0)).and.
     &            (i.lt.(MCP-P))) then
                  do k = 1,tu
                    icount = icount+1
                    FACE_CON(i+icount + (new_MCP-new_P)*(j-1+jcount))
     &                = e
                  enddo
              endif
              enddo
            enddo
        endif
        enddo
      endif
      
      
      
      if (CLOSED_V_flag.ne.1) then
c...    Loop Through orientation 2
        
        kcount = 0
        do k= 1,OCP-R
          icount = 0
          do i = 1,MCP-P
            e = e+1
            FACE_CON(i+icount +
     &        (1-CLOSED_W_flag)*(new_MCP-new_P)*(new_NCP-new_Q)+
     &        (k-1+kcount)*(new_MCP-new_P)) = e
          if (((U_KNOTp(i+P).ne.U_KNOTp(i+P+1)).and.(tu.gt.0)).and.
     &        (i.lt.(MCP-P))) then
              do j = 1,tu
                icount = icount+1
                FACE_CON(i+icount +
     &            (1-CLOSED_W_flag)*(new_MCP-new_P)*(new_NCP-new_Q)+
     &            (k-1+kcount)*(new_MCP-new_P)) = e
              enddo
          endif
          enddo
        if (((W_KNOTp(k+R).ne.W_KNOTp(k+1+R)).and.(tw.gt.0)).and.
     &        (k.lt.(OCP-R))) then
            do j = 1,tw
              kcount = kcount+1
              icount = 0
              e = (k-1)*(MCP-P)+(1-CLOSED_W_flag)*(MCP-P)*(NCP-Q)
              do i = 1,MCP-P
                e = e+1
                FACE_CON(i+icount +
     &            (1-CLOSED_W_flag)*(new_MCP-new_P)*(new_NCP-new_Q)+
     &            (k-1+kcount)*(new_MCP-new_P)) = e
              if (((U_KNOTp(i+P).ne.U_KNOTp(i+P+1)).and.(tu.gt.0)).and.
     &            (i.lt.(MCP-P))) then
                  do l = 1,tu
                    icount = icount+1
                    FACE_CON(i+icount +
     &                (1-CLOSED_W_flag)*(new_MCP-new_P)*(new_NCP-new_Q)+
     &                (k-1+kcount)*(new_MCP-new_P)) = e
                  enddo
              endif
              enddo
            enddo
        endif
        enddo
      endif
      
      
      
      
      
      if (CLOSED_U_flag.ne.1) then
c...    Loop Through orientation 3
        
        kcount = 0
        do k= 1,OCP-R
          jcount = 0            
          do j = 1,NCP-Q
            e = e+1
            FACE_CON(j+jcount +
     &        (1-CLOSED_W_flag)*(new_MCP-new_P)*(new_NCP-new_Q)+
     &        (1-CLOSED_V_flag)*(new_MCP-new_P)*(new_OCP-new_R)+
     &        (k-1+kcount)*(new_NCP-new_Q)) = e
          if (((V_KNOTp(j+Q).ne.V_KNOTp(j+1+Q)).and.(tv.gt.0)).and.
     &        (j.lt.(NCP-Q))) then
              do i = 1,tv
                jcount = jcount+1
                FACE_CON(j+jcount +
     &            (1-CLOSED_W_flag)*(new_MCP-new_P)*(new_NCP-new_Q)+
     &            (1-CLOSED_V_flag)*(new_MCP-new_P)*(new_OCP-new_R)+
     &            (k-1+kcount)*(new_NCP-new_Q)) = e
              enddo
          endif
          enddo
          if (((W_KNOTp(k+R).ne.W_KNOTp(k+1+R)).and.(tw.gt.0)).and.
     &      (k.lt.(OCP-R))) then
            do i = 1,tw
              kcount = kcount+1
              jcount = 0
              e = (k-1)*(NCP-Q)+(1-CLOSED_W_flag)*(MCP-P)*(NCP-Q)+
     &          (1-CLOSED_V_flag)*(MCP-P)*(OCP-R)
              do j = 1,NCP-Q
                e = e+1
                FACE_CON(j+jcount +
     &            (1-CLOSED_W_flag)*(new_MCP-new_P)*(new_NCP-new_Q)+
     &            (1-CLOSED_V_flag)*(new_MCP-new_P)*(new_OCP-new_R)+
     &            (k-1+kcount)*(new_NCP-new_Q)) = e
              if (((V_KNOTp(j+Q).ne.V_KNOTp(j+1+Q)).and.(tv.gt.0)).and.
     &            (j.lt.(NCP-Q))) then
                  do l = 1,tv
                    jcount = jcount+1
                    FACE_CON(j+jcount +
     &                (1-CLOSED_W_flag)*(new_MCP-new_P)*(new_NCP-new_Q)+
     &                (1-CLOSED_V_flag)*(new_MCP-new_P)*(new_OCP-new_R)+
     &                (k-1+kcount)*(new_NCP-new_Q)) = e
                  enddo
              endif
              enddo
            enddo
          endif
        enddo
      endif
      
      
      
      
      if (CLOSED_V_flag.ne.1) then
c...    Loop Through orientation 4
        
        kcount = 0
        do k= 1,OCP-R
          icount = 0
          do i = 1,MCP-P
            e = e+1
            FACE_CON(i+icount +
     &        (1-CLOSED_W_flag)*(new_MCP-new_P)*(new_NCP-new_Q)+
     &        (1-CLOSED_V_flag)*(new_MCP-new_P)*(new_OCP-new_R)+
     &        (1-CLOSED_U_flag)*(new_NCP-new_Q)*(new_OCP-new_R)+
     &        (k-1+kcount)*(new_MCP-new_P)) = e
          if (((U_KNOTp(i+P).ne.U_KNOTp(i+P+1)).and.(tu.gt.0)).and.
     &        (i.lt.(MCP-P))) then
              do j = 1,tu
                icount = icount+1
                FACE_CON(i+icount +
     &            (1-CLOSED_W_flag)*(new_MCP-new_P)*(new_NCP-new_Q)+
     &            (1-CLOSED_V_flag)*(new_MCP-new_P)*(new_OCP-new_R)+
     &            (1-CLOSED_U_flag)*(new_NCP-new_Q)*(new_OCP-new_R)+
     &            (k-1+kcount)*(new_MCP-new_P)) = e
              enddo
          endif
          enddo
          if (((W_KNOTp(k+R).ne.W_KNOTp(k+1+R)).and.(tw.gt.0)).and.
     &      (k.lt.(OCP-R))) then
            do j = 1,tw
              kcount = kcount+1
              icount = 0
              e = (k-1)*(MCP-P)+(1-CLOSED_W_flag)*(MCP-P)*(NCP-Q)+
     &          (1-CLOSED_V_flag)*(MCP-P)*(OCP-R)+
     &          (1-CLOSED_U_flag)*(NCP-Q)*(OCP-R)
              do i = 1,MCP-P
                e = e+1
                FACE_CON(i+icount +
     &            (1-CLOSED_W_flag)*(new_MCP-new_P)*(new_NCP-new_Q)+
     &            (1-CLOSED_V_flag)*(new_MCP-new_P)*(new_OCP-new_R)+
     &            (1-CLOSED_U_flag)*(new_NCP-new_Q)*(new_OCP-new_R)+
     &            (k-1+kcount)*(new_MCP-new_P)) = e
              if (((U_KNOTp(i+P).ne.U_KNOTp(i+P+1)).and.(tu.gt.0)).and.
     &            (i.lt.(MCP-P))) then
                  do l = 1,tu
                    icount = icount+1
                    FACE_CON(i+icount +
     &                (1-CLOSED_W_flag)*(new_MCP-new_P)*(new_NCP-new_Q)+
     &                (1-CLOSED_V_flag)*(new_MCP-new_P)*(new_OCP-new_R)+
     &                (1-CLOSED_U_flag)*(new_NCP-new_Q)*(new_OCP-new_R)+
     &                (k-1+kcount)*(new_MCP-new_P)) = e
                  enddo
              endif
              enddo
            enddo
          endif
        enddo
      endif
      
      

      
      if (CLOSED_U_flag.ne.1) then
c...    Loop Through orientation 5
        
        kcount = 0
        do k= 1,OCP-R
          jcount = 0            
          do j = 1,NCP-Q
            e = e+1
            FACE_CON(j+jcount +
     &        (1-CLOSED_W_flag)*(new_MCP-new_P)*(new_NCP-new_Q)+
     &        2*(1-CLOSED_V_flag)*(new_MCP-new_P)*(new_OCP-new_R)+
     &        (1-CLOSED_U_flag)*(new_NCP-new_Q)*(new_OCP-new_R)+
     &        (k-1+kcount)*(new_NCP-new_Q)) = e
            if (((V_KNOTp(j+Q).ne.V_KNOTp(j+1+Q)).and.(tv.gt.0)).and.
     &        (j.lt.(NCP-Q))) then
              do i = 1,tv
                jcount = jcount+1
                FACE_CON(j+jcount +
     &            (1-CLOSED_W_flag)*(new_MCP-new_P)*(new_NCP-new_Q)+
     &            2*(1-CLOSED_V_flag)*(new_MCP-new_P)*(new_OCP-new_R)+
     &            (1-CLOSED_U_flag)*(new_NCP-new_Q)*(new_OCP-new_R)+
     &            (k-1+kcount)*(new_NCP-new_Q)) = e
              enddo
            endif
          enddo
          if (((W_KNOTp(k+R).ne.W_KNOTp(k+1+R)).and.(tw.gt.0)).and.
     &      (k.lt.(OCP-R))) then
            do i = 1,tw
              kcount = kcount+1
              jcount = 0
              e = (k-1)*(NCP-Q)+(1-CLOSED_W_flag)*(MCP-P)*(NCP-Q)+
     &          2*(1-CLOSED_V_flag)*(MCP-P)*(OCP-R)+
     &          (1-CLOSED_U_flag)*(NCP-Q)*(OCP-R)
              do j = 1,NCP-P
                e = e+1
                FACE_CON(j+jcount +
     &            (1-CLOSED_W_flag)*(new_MCP-new_P)*(new_NCP-new_Q)+
     &            2*(1-CLOSED_V_flag)*(new_MCP-new_P)*(new_OCP-new_R)+
     &            (1-CLOSED_U_flag)*(new_NCP-new_Q)*(new_OCP-new_R)+
     &            (k-1+kcount)*(new_NCP-new_Q)) = e
              if (((V_KNOTp(j+Q).ne.V_KNOTp(j+1+Q)).and.(tv.gt.0)).and.
     &            (j.lt.(NCP-Q))) then
                  do l = 1,tv
                    jcount = jcount+1
                    FACE_CON(j+jcount +
     &                (1-CLOSED_W_flag)*(new_MCP-new_P)*(new_NCP-new_Q)+
     &                2*(1-CLOSED_V_flag)*(new_MCP-new_P)*
     &                (new_OCP-new_R)+
     &                (1-CLOSED_U_flag)*(new_NCP-new_Q)*(new_OCP-new_R)+
     &                (k-1+kcount)*(new_NCP-new_Q)) = e
                  enddo
              endif
              enddo
            enddo
          endif
        enddo
      endif


      
      if (CLOSED_W_flag.ne.1) then
c...    Loop through orientation 6
        
        jcount = 0
        do j= 1,NCP-Q
          icount = 0
          do i = 1,MCP-P
            e = e+1
            FACE_CON(i+icount +
     &        (1-CLOSED_W_flag)*(new_MCP-new_P)*(new_NCP-new_Q)+
     &        2*(1-CLOSED_U_flag)*(new_NCP-new_Q)*(new_OCP-new_R) +
     &        2*(1-CLOSED_V_flag)*(new_MCP-new_P)*(new_OCP-new_R) +
     &        (new_MCP-new_P)*(j-1+jcount)) = e
            if (((U_KNOTp(i+P).ne.U_KNOTp(i+P+1)).and.(tu.gt.0)).and.
     &        (i.lt.(MCP-P))) then
              do k = 1,tu
                icount = icount+1
                FACE_CON(i+icount +
     &            (1-CLOSED_W_flag)*(new_MCP-new_P)*(new_NCP-new_Q)+
     &            2*(1-CLOSED_U_flag)*(new_NCP-new_Q)*(new_OCP-new_R) +
     &            2*(1-CLOSED_V_flag)*(new_MCP-new_P)*(new_OCP-new_R) +
     &            (new_MCP-new_P)*(j-1+jcount)) = e
              enddo
            endif
          enddo
          if (((V_KNOTp(j+Q).ne.V_KNOTp(j+1+Q)).and.(tv.gt.0)).and.
     &      (j.lt.(NCP-Q))) then
            do k = 1,tv
              jcount = jcount+1
              icount = 0
              e = (1-CLOSED_W_flag)*(MCP-P)*(NCP-Q)+
     &          2*(1-CLOSED_U_flag)*(NCP-Q)*(OCP-R) +
     &          2*(1-CLOSED_V_flag)*(MCP-P)*(OCP-R) +
     &          (MCP-P)*(j-1)
              do i = 1,MCP-P
                e = e+1
                FACE_CON(i+icount +
     &            (1-CLOSED_W_flag)*(new_MCP-new_P)*(new_NCP-new_Q)+
     &            2*(1-CLOSED_U_flag)*(new_NCP-new_Q)*(new_OCP-new_R) +
     &            2*(1-CLOSED_V_flag)*(new_MCP-new_P)*(new_OCP-new_R) +
     &            (new_MCP-new_P)*(j-1+jcount)) = e
              if (((U_KNOTp(i+P).ne.U_KNOTp(i+P+1)).and.(tu.gt.0)).and.
     &            (i.lt.(MCP-P))) then
                  do l = 1,tu
                    icount = icount+1
                    FACE_CON(i+icount +
     &                (1-CLOSED_W_flag)*(new_MCP-new_P)*(new_NCP-new_Q)+
     &                2*(1-CLOSED_U_flag)*(new_NCP-new_Q)*
     &                (new_OCP-new_R) +
     &                2*(1-CLOSED_V_flag)*(new_MCP-new_P)*
     &                (new_OCP-new_R) +
     &                (new_MCP-new_P)*(j-1+jcount)) = e
                  enddo
              endif
              enddo
            enddo
          endif
        enddo
      endif
      
      
      
      return
      end
