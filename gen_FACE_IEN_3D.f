      subroutine gen_FACE_IEN_3D_u

      use aAdjKeep		! this f90 module contains all allocatables
c     
      use common                ! common file defines all variables

      implicit none

c     
c...  Local variables
      integer :: i,j,face,e
      
      allocate(FACE_IENu(NFACEu,NSHLu))
      allocate(FACE_OR(NFACEu))
c     
c...  Initialize matrices and variables
c     
      FACE_IENu = 0
      FACE_OR = 0
      face = 0    
      
      
c     This paragraph doesn't really make sense. FIXME/    
c...  Loop through faces assigning numbers and filling out matrices.
c     Face numbers will be assigned sequentially starting with the
c     u-v surface at w = 1, denoted (u,v,1). From there we count
c     the (u,1,w), then the (MCP,v,w), the (u,NCP,w),the (1,v,w),
c     and lastly (u,v,OCP) surface. Once the face has been labled,
c     we find the number of the element the face belongs to and
c     assign the local node numbers.
c     
c     
c     Face Orientation scheme: 
c     1 - surface (u,v,1)
c     2 - surface (u,1,w)
c     3 - surface (MCP,v,w)
c     4 - surface (u,NCP,w)
c     5 - surface (1,v,w)
c     6 - surface (u,v,OCP)



c...  Orientation 1
      e = 0
      if (CLOSED_W_flag.ne.1) then
         do j = 1,NCPu
	    do i = 1,MCPu
               if ((i .ge. (Pu+1)).and.(j .ge. (Qu+1))) then
                  face = face + 1 ! face number
                  FACE_OR(face) = 1 ! face orientation
                  e = e + 1	! element number
                  FACE_IENu(face,:) = IENu(e,:)
               endif
	    enddo
         enddo
      endif
      
      
      
c...  Orientation 2
      if (CLOSED_V_flag.ne.1) then
         do j = Ru+1,OCPu
	    e = (MCPu-Pu)*(NCPu-Qu)*(j-(Ru+1))
	    do i = Pu+1,MCPu
               e = e + 1        ! element number
               face = face+1
               FACE_OR(face) = 2
               FACE_IENu(face,:) = IENu(e,:)
	    enddo
         enddo
      endif
      
      
      
c...  Orientation 3
      if (CLOSED_U_flag.ne.1) then
         do j = Ru+1,OCPu
	    e = (MCPu-Pu)*(NCPu-Qu)*(j-(Ru+1))
	    do i = Qu+1,NCPu
               face = face+1
               FACE_OR(face) = 3
               e = e + (MCPu-Pu)
               FACE_IENu(face,:) = IENu(e,:)
	    enddo
         enddo
      endif
      
      
      
c...  Orientation 4
      if (CLOSED_V_flag.ne.1) then
         do j = Ru+1,OCPu
	    e = (MCPu-Pu)*(NCPu-Qu)*(j-Ru)+1
	    do i = Pu+1,MCPu
               face = face+1
               FACE_OR(face) = 4
               e = e - 1 
               FACE_IENu(face,:) = IENu(e,:)
	    enddo
         enddo
      endif
      
      
      
c...  Orientation 5
      if (CLOSED_U_flag.ne.1) then
         do j = Ru+1,OCPu
	    e = (MCPu-Pu)*(NCPu-Qu)*(j-Ru)+1
	    do i = Qu+1,NCPu
               face = face+1
               FACE_OR(face) = 5
               e = e - (MCPu-Pu)
               FACE_IENu(face,:) = IENu(e,:)
	    enddo
         enddo
      endif
      
      
c...  Orientation 6
      if (CLOSED_W_flag.ne.1) then
         e = (MCPu-Pu)*(NCPu-Qu)*(OCPu-Ru-1)
         do j = 1,NCPu
	    do i = 1,MCPu
               if ((i .ge. (Pu+1)).and.(j .ge. (Qu+1))) then
                  face = face + 1
                  FACE_OR(face) = 6
                  e = e + 1
                  FACE_IENu(face,:) = IENu(e,:)
               endif
	    enddo
         enddo
      endif
      
      return
      end
	
	

      subroutine gen_FACE_IEN_3D_p

      use aAdjKeep		! this f90 module contains all allocatables
c       
      use common	! common file defines all variables
      implicit none
c       
c...    Local variables
      integer :: i,j,face,e
	
      allocate(FACE_IENp(NFACEp,NSHLp))
c       
c...    Initialize matrices and variables
c       
      FACE_IENp = 0
        face = 0    
	
	
c       This paragraph doesn't really make sense. FIXME/    
c...    Loop through faces assigning numbers and filling out matrices.
c          Face numbers will be assigned sequentially starting with the
c          u-v surface at w = 1, denoted (u,v,1). From there we count
c          the (u,1,w), then the (MCP,v,w), the (u,NCP,w),the (1,v,w),
c          and lastly (u,v,OCP) surface. Once the face has been labled,
c          we find the number of the element the face belongs to and
c          assign the local node numbers.
c
c
c       Face Orientation scheme: 
c           1 - surface (u,v,1)
c           2 - surface (u,1,w)
c           3 - surface (MCP,v,w)
c           4 - surface (u,NCP,w)
c           5 - surface (1,v,w)
c           6 - surface (u,v,OCP)



c...    Orientation 1
      e = 0
      if (CLOSED_W_flag.ne.1) then
	  do j = 1,NCPp
	    do i = 1,MCPp
	      if ((i .ge. (Pp+1)).and.(j .ge. (Qp+1))) then
		face = face + 1	! face number
		e = e + 1	! element number
		FACE_IENp(face,:) = IENp(e,:)
	      endif
	    enddo
	  enddo
      endif
	
	
	
c...    Orientation 2
      if (CLOSED_V_flag.ne.1) then
	  do j = Rp+1,OCPp
	    e = (MCPp-Pp)*(NCPp-Qp)*(j-(Rp+1))
	    do i = Pp+1,MCPp
	      e = e + 1		! element number
	      face = face+1
	      FACE_IENp(face,:) = IENp(e,:)
	    enddo
	  enddo
      endif
	
	
	
c...    Orientation 3
      if (CLOSED_U_flag.ne.1) then
	  do j = Rp+1,OCPp
	    e = (MCPp-Pp)*(NCPp-Qp)*(j-(Rp+1))
	    do i = Qp+1,NCPp
	      face = face+1
	      e = e + (MCPp-Pp)
	      FACE_IENp(face,:) = IENp(e,:)
	    enddo
	  enddo
      endif
	
	
	
c...    Orientation 4
      if (CLOSED_V_flag.ne.1) then
	  do j = Rp+1,OCPp
	    e = (MCPp-Pp)*(NCPp-Qp)*(j-Rp)+1
	    do i = Pp+1,MCPp
	      face = face+1
	      e = e - 1 
	      FACE_IENp(face,:) = IENp(e,:)
	    enddo
	  enddo
      endif
	
	
	
c...    Orientation 5
      if (CLOSED_U_flag.ne.1) then
	  do j = Rp+1,OCPp
	    e = (MCPp-Pp)*(NCPp-Qp)*(j-Rp)+1
	    do i = Qp+1,NCPp
	      face = face+1
	      e = e - (MCPp-Pp)
	      FACE_IENp(face,:) = IENp(e,:)
	    enddo
	  enddo
      endif
	
	
c...    Orientation 6
      if (CLOSED_W_flag.ne.1) then
      e = (MCPp-Pp)*(NCPp-Qp)*(OCPp-Rp-1)
      do j = 1,NCPp
	    do i = 1,MCPp
	      if ((i .ge. (Pp+1)).and.(j .ge. (Qp+1))) then
		face = face + 1
		e = e + 1
		FACE_IENp(face,:) = IENp(e,:)
	      endif
	    enddo
      enddo
      endif
	
      return
      end


