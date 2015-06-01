
      subroutine genSparStruc_3D_u (colm, rowp, icnt)
      
      
      use aAdjKeep

      use common
      implicit none

      integer tmpr(NNODZu) , colm(NNODZu+1),
     &  rowp(NNODZu*8*(Pu+1)*(Qu+1)*(Ru+1))
      integer adjcnt(NNODZu),  mloc(1)
      integer i, j, k, imin, icnt, ibig, ncol 
      
      integer, allocatable :: row_fill_list(:,:)
      
      allocate(row_fill_list(NNODZu,6*8*(Pu+1)*(Qu+1)*(Ru+1)))
      
      row_fill_list = 0 
      
      adjcnt = 0
      
c.... compute sparse matrix data structures
c     
      call Asadj_u (row_fill_list,  adjcnt)

c     build the colm array
c     
      colm(1)=1
      
      do i=1,NNODZu
        colm(i+1)=colm(i)+adjcnt(i)
      enddo
      
c     sort the rowp into increasing order
c     
      ibig=10*NNODZu
      icnt=0
      do i=1,NNODZu
        ncol=adjcnt(i)
        tmpr(1:ncol)=row_fill_list(i,1:ncol)
        do j=1,ncol
          icnt=icnt+1
          imin=minval(tmpr(1:ncol))
          mloc=minloc(tmpr(1:ncol))
          rowp(icnt)=imin
          tmpr(mloc(1))=ibig
        enddo
      enddo
      
      deallocate(row_fill_list)
      
      return
      end



      subroutine Asadj_u (row_fill_list, adjcnt)
      
      use aAdjKeep
      
      use common
      implicit none

      integer row_fill_list(NNODZu,6*8*(Pu+1)*(Qu+1)*(Ru+1)),
     &     k, i, j, ni, nj, nk,
     &     l, ibroke, knd,
     &     adjcnt(NNODZu), ndlist((Pu+1)*(Qu+1)*(Ru+1)),
     &     jnd, jlngth
      
      
      do i=1,NELu
        
c...    Check to see if current element has nonzero area
        ni = INNu(IENu(i,1),1)  ! get NURB coordinates
        nj = INNu(IENu(i,1),2)
        nk = INNu(IENu(i,1),3)
        
        if ((U_KNOTu(ni).ne.U_KNOTu(ni+1)).and.
     &    (V_KNOTu(nj).ne.V_KNOTu(nj+1)).and.
     &    (W_KNOTu(nk).ne.W_KNOTu(nk+1))) then
          
          nshlu = (Pu+1)*(Qu+1)*(Ru+1)
          
          do j=1,nshlu
            ndlist(j)=IENu(i,j) ! gen list of global "nodes" for this element
          enddo
          do j=1,nshlu
            jnd=ndlist(j)       ! jnd is the global "node" we are working on
            jlngth=adjcnt(jnd)  ! current length of j's list
            do k=1,nshlu 
              knd=ndlist(k)
              ibroke=0
              
              do l= 1,jlngth    ! row_fill_list is, for each node, the
                                ! list of nodes that I have already
                                ! detected interaction with
                if(row_fill_list(jnd,l).eq. knd) then
                  ibroke=1
                  exit
                endif
              enddo

c             to get here k was not in  j's list so add it
              if(ibroke.eq.0) then
                jlngth=jlngth+1 ! lengthen list
                row_fill_list(jnd,jlngth)=knd ! add unique entry to list
              endif
            enddo               ! finished checking all the k's for this j
            adjcnt(jnd)=jlngth  ! update the counter
          enddo                 ! done with j's
          
        endif
        
      enddo                     ! done with elements in this block
        
      
      
      
      return
      end
      
      
      subroutine genSparStruc_3D_p (colm, rowp, icnt)

 
      use aAdjKeep
      
      use common
      implicit none

      integer tmpr(NNODZu) , colm(NNODZu+1),
     &  rowp(NNODZu*8*(Pu+1)*(Qu+1)*(Ru+1))
      integer adjcnt(NNODZu),  mloc(1)
      integer i, j, k, imin, icnt, ibig, ncol 
      
      integer, allocatable :: row_fill_list(:,:)
      
      allocate(row_fill_list(NNODZu,6*8*(Pu+1)*(Qu+1)*(Ru+1)))
      
      row_fill_list = 0 
      
      adjcnt = 0
      
c.... compute sparse matrix data structures
c     
      call Asadj_p (row_fill_list,  adjcnt)
      
c     build the colm array
c     
      colm(1)=1
      
      do i=1,NNODZu
        colm(i+1)=colm(i)+adjcnt(i)
      enddo
      
c     sort the rowp into increasing order
c     
      ibig=10*NNODZu
      icnt=0
      do i=1,NNODZu
        ncol=adjcnt(i)
        tmpr(1:ncol)=row_fill_list(i,1:ncol)
        do j=1,ncol
          icnt=icnt+1
          imin=minval(tmpr(1:ncol))
          mloc=minloc(tmpr(1:ncol))
          rowp(icnt)=imin
          tmpr(mloc(1))=ibig
        enddo
      enddo
      
      deallocate(row_fill_list)
      
      return
      end
      
      

      subroutine Asadj_p(row_fill_list, adjcnt)
      
      use aAdjKeep
      
      use common
      implicit none

      integer row_fill_list(NNODZu,6*8*(Pu+1)*(Qu+1)*(Ru+1)),
     &     k, i, j, ni, nj, nk,
     &     l, ibroke, knd,
     &     adjcnt(NNODZu), ndlistu((Pu+1)*(Qu+1)*(Ru+1)),
     &     jnd, jlngth, ndlistp((Pp+1)*(Qp+1)*(Rp+1))
      
      
      do i=1,NELu
        
c...    Check to see if current element has nonzero area
        ni = INNu(IENu(i,1),1)  ! get NURB coordinates
        nj = INNu(IENu(i,1),2)
        nk = INNu(IENu(i,1),3)
        
        if ((U_KNOTu(ni).ne.U_KNOTu(ni+1)).and.
     &    (V_KNOTu(nj).ne.V_KNOTu(nj+1)).and.
     &    (W_KNOTu(nk).ne.W_KNOTu(nk+1))) then
          
          nshlu = (Pu+1)*(Qu+1)*(Ru+1)
          nshlp = (Pp+1)*(Qp+1)*(Rp+1)
          
          do j=1,nshlu
            ndlistu(j)=IENu(i,j) ! gen list of global "nodes" for this element
          enddo
          
          do j=1,nshlp
            ndlistp(j)=IENp(EL_CON(i),j) ! gen list of global "nodes" for this element
          enddo
          
          
          do j=1,nshlu
            jnd=ndlistu(j)      ! jnd is the global "node" we are working on
            jlngth=adjcnt(jnd)  ! current length of j's list
            do k=1,nshlp 
              knd=ndlistp(k)
              ibroke=0
              
              do l= 1,jlngth    ! row_fill_list is, for each node, the
                                ! list of nodes that I have already
                                ! detected interaction with
                if(row_fill_list(jnd,l).eq. knd) then
                  ibroke=1
                  exit
                endif
              enddo
              
c             
c             to get here k was not in  j's list so add it
c             
              if(ibroke.eq.0) then
                jlngth=jlngth+1 ! lengthen list

                if(jlngth.gt.48*nshlu) then
                  write(*,*) 'increase overflow factor in genadj'
                  stop
                endif

                row_fill_list(jnd,jlngth)=knd ! add unique entry to list
              endif
            enddo               ! finished checking all the k's for this j
            adjcnt(jnd)=jlngth  ! update the counter
          enddo                 ! done with j's
          
        endif
        
      enddo                     ! done with elements in this block
      
      
      return
      end
      
      subroutine genSparStruc_3D_pp (colm, rowp, icnt)
      
      
      use aAdjKeep

      use common
      implicit none

      integer tmpr(NNODZp) , colm(NNODZp+1),
     &  rowp(NNODZp*8*(Pp+1)*(Qp+1)*(Rp+1))
      integer adjcnt(NNODZp),  mloc(1)
      integer i, j, k, imin, icnt, ibig, ncol 
      
      integer, allocatable :: row_fill_list(:,:)
      
      allocate(row_fill_list(NNODZp,6*8*(Pp+1)*(Qp+1)*(Rp+1)))
      
      row_fill_list = 0 
      
      adjcnt = 0
      
c.... compute sparse matrix data structures
c     
      call Asadj_pp (row_fill_list,  adjcnt)

c     build the colm array
c     
      colm(1)=1
      
      do i=1,NNODZp
        colm(i+1)=colm(i)+adjcnt(i)
      enddo
      
c     sort the rowp into increasing order
c     
      ibig=10*NNODZp
      icnt=0
      do i=1,NNODZp
        ncol=adjcnt(i)
        tmpr(1:ncol)=row_fill_list(i,1:ncol)
        do j=1,ncol
          icnt=icnt+1
          imin=minval(tmpr(1:ncol))
          mloc=minloc(tmpr(1:ncol))
          rowp(icnt)=imin
          tmpr(mloc(1))=ibig
        enddo
      enddo
      
      deallocate(row_fill_list)
      
      return
      end



      subroutine Asadj_pp (row_fill_list, adjcnt)
      
      use aAdjKeep
      
      use common
      implicit none

      integer row_fill_list(NNODZp,6*8*(Pp+1)*(Qp+1)*(Rp+1)),
     &     k, i, j, ni, nj, nk,
     &     l, ibroke, knd,
     &     adjcnt(NNODZp), ndlist((Pp+1)*(Qp+1)*(Rp+1)),
     &     jnd, jlngth
      
      
      do i=1,NELp
        
c...    Check to see if current element has nonzero area
        ni = INNp(IENp(i,1),1)  ! get NURB coordinates
        nj = INNp(IENp(i,1),2)
        nk = INNp(IENp(i,1),3)
        
        if ((U_KNOTp(ni).ne.U_KNOTp(ni+1)).and.
     &    (V_KNOTp(nj).ne.V_KNOTp(nj+1)).and.
     &    (W_KNOTp(nk).ne.W_KNOTp(nk+1))) then
          
          nshlp = (Pp+1)*(Qp+1)*(Rp+1)
          
          do j=1,nshlp
            ndlist(j)=IENp(i,j) ! gen list of global "nodes" for this element
          enddo
          do j=1,nshlp
            jnd=ndlist(j)       ! jnd is the global "node" we are working on
            jlngth=adjcnt(jnd)  ! current length of j's list
            do k=1,nshlp 
              knd=ndlist(k)
              ibroke=0
              
              do l= 1,jlngth    ! row_fill_list is, for each node, the
                                ! list of nodes that I have already
                                ! detected interaction with
                if(row_fill_list(jnd,l).eq. knd) then
                  ibroke=1
                  exit
                endif
              enddo
              
c             
c             to get here k was not in  j's list so add it
c             
              if(ibroke.eq.0) then
                jlngth=jlngth+1 ! lengthen list
                row_fill_list(jnd,jlngth)=knd ! add unique entry to list
              endif
            enddo               ! finished checking all the k's for this j
            adjcnt(jnd)=jlngth  ! update the counter
          enddo                 ! done with j's
          
        endif
        
      enddo                     ! done with elements in this block
        
      
      
      
      return
      end
