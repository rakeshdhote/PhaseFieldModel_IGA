c     This subroutine generates the prolongation/restriction operators
      
      subroutine gen_constraint(
     *     W_KNOT_M_, W_KNOT_S_,  !Master, Slave knot positions
     *     Tw,R)

      implicit none
 
      integer knotf,i,R,pflag,dir, ierr
      real*8, dimension(2*R,R) :: Tw
      real*8, allocatable, dimension(:) :: 
     &     w_knot_f,w_knot_c,w_knot_m,w_knot_s
      real*8  w_knot_m_(R+1), w_knot_s_(2*R)
      logical :: stopper = .false.

      allocate(w_knot_m (R+1));  w_knot_m = w_knot_m_
      allocate(w_knot_s (2*R));  w_knot_s = w_knot_s_
      allocate(w_knot_f(3*R+1)); w_knot_f = 0.d0
      allocate(w_knot_c(2*R+1)); w_knot_c = 0.d0

      do i = 1,R
         if ( W_KNOT_S(R+1) /= W_KNOT_S(R+i) ) then
            print*, "Knot vector does not have the required order"
            stopper = .true.
         endif
      enddo
                                ! Assemble knot arrays to send to T_mat
      if ( w_knot_s(R+1) == W_KNOT_M(1)) then ! Verify if patches match
         continue
      elseif (W_KNOT_M(1)==0.0) then ! Otherwise connecting last to first
         W_KNOT_M = W_KNOT_M + W_KNOT_S(R+1)
      else
         print*, "Something is wrong with the patch alignment"
         stopper = .true.
      endif

      if (stopper)
     &     stop

      w_knot_c(    1:R  ) = w_knot_s(1:R)
      w_knot_c(  R+1:   ) = w_knot_m( : )
      w_knot_f(    1:2*R) = w_knot_s( : )
      w_knot_f(2*R+1:   ) = w_knot_m( : )

      call T_mat(R,w_knot_c,w_knot_f,Tw)

      return
      end


!#######################################################################      
c     This subroutine generates the coefficient matrix corresponding to
c     knot refinement. It is used to enforce interface conditions between
c     two patches with different levels of refinement.

c     This is for the multi-patch code.
c     
c     Oct 22, 2004, J. Austin Cottrell
c     
c     Has been modified to impose pl-continuity across patch boundaries

      subroutine T_mat(pl,xc,xf,Tu)

      implicit none
                                !  Local variables
      integer :: pl,mlc,mlf,i,j,k
      real*8, dimension(2*pl,2*pl+1) :: T, temp
      real*8, dimension(2*pl,pl) :: Tu
      real*8, dimension(3*pl+1) :: xf
      real*8, dimension(2*pl+1) :: xc
      real*8 :: tiny=1d-12

      T=0.0                     ! initialize T matrix
      do i=1,pl+1
         T(i,i) = 1.d0
      enddo
      T(pl+1:,pl+1) = 1.d0

      do k = 1,pl               ! Compute T matrix
         temp = T;
         T=0.d0;
         do i = 1,2*pl
            do j = 1,pl
               if ( abs(xc(j+k)-xc(j))>tiny ) 
     &              T(i,j) = T(i,j)+(xf(i+k)-xc(j))/
     &              (xc(j+k)-xc(j))*temp(i,j);
               if ( abs(xc(j+1+k)-xc(j+1))>tiny ) 
     &              T(i,j) = T(i,j)+(xc(j+k+1)-xf(i+k))/
     &              (xc(j+1+k)-xc(j+1))*temp(i,j+1);
            enddo
         enddo
      enddo

      Tu = T(:,1:pl)

      return
      end

!#######################################################################      
      subroutine constraint_input(
     *     TwM,TwS,             !Master, Slave restriction operators
     *     R,iknot, lprt,var,dir)

 
      implicit none

      integer knotf,i,j, R, ri,rj, iknot, lprt, ierr
      real*8, dimension(2*R,R) :: TwM,TwS
      character*30 fname
      character*10 cname
      logical exts, stopper
      character*1 var,dir

      stopper = .false.
      knotf   = 14

      fname = trim(dir//'knotf'//var//cname(iknot))//'.dat'
      inquire(file=fname,exist=exts)
      if (exts) then
         open (knotf,file=fname, status='old')

         read (knotf,"(3i5)") ri,rj,lprt
         if (ri/=2*R .or. rj/=R) then
            print*,'Incompatible restriction operators'
            stopper = .true.
         else
            if ( lprt==1 .or. lprt==2 ) then
               do i = 1,2*R
                  do j = 1,R
                     if (j < R) then
                        read (knotf,"(F13.9)", ADVANCE = "NO") TwM(i,j)
                     else
                        read (knotf,"(F13.9)", ADVANCE = "YES")TwM(i,j)
                     endif
                  enddo
               enddo
            endif
            if ( lprt==1 .or. lprt==0 ) then
               do i = 1,2*R
                  do j = 1,R
                     if (j < R) then
                        read (knotf,"(F13.9)", ADVANCE = "NO") TwS(i,j)
                     else
                        read (knotf,"(F13.9)", ADVANCE = "YES")TwS(i,j)
                     endif
                  enddo
               enddo
            endif

            close (knotf)

         endif

      else
         stopper = .true.
         print*," Missing ",fname
      endif

      return
      end



!###########################################################################

      subroutine shuffle_u(
     &     global,
     &     TuM,TuS,             ! Restriction operators
     &     lprt,
     &     P, NNODZ, MCP, NCP, OCP, NSD,
     &     cflag,var)


      implicit none

      integer ni, nj, nk, counti, countj
      integer P, NNODZ, NSD, MCP, NCP, OCP,lprt
      real*8 global(NNODZ,NSD)
      real*8 tmpu(OCP*NCP,NSD)
      real*8 TuM(2*P,P), TuS(2*P,P)
      character*3 cflag
      character var


      if (cflag.eq.'in ') then

!     U-Direction: Control points on plane are spaced MCP apart
         do ni = P,1,-1         ! First P planes (master)
            tmpu = 0d+0
            do nj = 1, ni
               tmpu = tmpu
     &              + TuM(nj+P,ni)
     &              * global(nj:nj+MCP*(NCP*OCP-1):MCP,:) 
            enddo
            global(ni:ni+(NCP*OCP-1)*MCP:MCP,:) = tmpu
         enddo

         do ni = 1, P           ! Last P planes (slave)
            tmpu = 0d+0
            do nj = ni, P
               tmpu = tmpu
     &              + TuS(nj,ni)
     &              * global(MCP-P+nj:MCP*NCP*OCP-P+nj:MCP,:)
            enddo
            global(MCP-P+ni:MCP*NCP*OCP-P+ni:MCP,:) = tmpu
         enddo

         do ni = 1, P           ! Master = Slave+Master
            global       (ni:ni+(NCP*OCP-1)*MCP:MCP,:) =
     &           + global(ni:ni+(NCP*OCP-1)*MCP:MCP,:)
     &           + global(MCP-P+ni:MCP*NCP*OCP-P+ni:MCP,:)
            global(MCP-P+ni:MCP*NCP*OCP-P+ni:MCP,:) = 0d0 ! zero out Slave
         enddo

      elseif (cflag.eq.'out') then

         do ni = 1, P           ! Slave = Master
            global(MCP-P+ni:MCP*NCP*OCP-P+ni:MCP,:) =
     &           global(ni:ni+(NCP*OCP-1)*MCP:MCP,:)
         enddo

!     U-Direction: Control points on plane are spaced MCP apart
         do ni = 1,P            ! First P planes (master)
            tmpu = 0d+0
            do nj = ni, P
               tmpu = tmpu
     &              + TuM(P+ni,nj)
     &              * global(nj:nj+MCP*(NCP*OCP-1):MCP,:) 
            enddo
            global(ni:ni+(NCP*OCP-1)*MCP:MCP,:) = tmpu
         enddo

         do ni = P,1,-1         ! Last P planes (slave)
            tmpu = 0d+0
            do nj = 1,ni
               tmpu = tmpu
     &              + TuS(ni,nj)
     &              * global(MCP-P+nj:MCP*NCP*OCP-P+nj:MCP,:)
            enddo
            global(MCP-P+ni:MCP*NCP*OCP-P+ni:MCP,:) = tmpu
         enddo

      endif


      return
      end



!###########################################################################

      subroutine shuffle_v(     ! V - direction
     &     global,
     &     TvM,TvS,             ! Restriction operators
     &     lprt,
     &     Q, NNODZ, MCP, NCP, OCP, NSD,
     &     cflag,var)


      implicit none

      integer ni, nj, nk, counti, countj
      integer Q, NNODZ, NSD, MCP, NCP, OCP,lprt
      real*8 global(NNODZ,NSD)
      real*8 tmpv(MCP,NSD)
      real*8 TvM(2*Q,Q), TvS(2*Q,Q)
      character*3 cflag
      character var


      if (cflag.eq.'in ') then

         do ni = Q,1,-1         ! First Q planes (master)
            do nk = 1, OCP
               tmpv = 0d+0
               do nj = 1, ni
                  tmpv = tmpv
     &                 + TvM(nj+Q,ni)
     &                 * global
     &                 (      MCP*(NCP * (nk-1) + nj-1 ) + 1
     &                 :      MCP*(NCP * (nk-1) + nj),:) 
               enddo
               global( MCP*(NCP * (nk-1) + ni-1 ) + 1
     &              :  MCP*(NCP * (nk-1) + ni),:) = tmpv
            enddo
         enddo

         do ni = 1,Q            ! Last Q planes (slave)
            do nk = 1, OCP
               tmpv = 0d+0
               do nj = ni, Q
                  tmpv = tmpv
     &                 + TvS(nj,ni)
     &                 * global
     &                 (    MCP*(nk * NCP - Q + nj -1) +1
     &                 :    MCP*(nk * NCP - Q + nj),:)
               enddo
               global
     &              (    MCP*(nk * NCP - Q + ni -1) +1
     &              :    MCP*(nk * NCP - Q + ni),:) = tmpv
            enddo
         enddo

         do ni = 1, Q           ! Master = Slave + Master
            do nk = 1, OCP
               global
     &              (  MCP*(NCP * (nk-1) + ni-1 ) + 1
     &              :  MCP*(NCP * (nk-1) + ni),:) =
     &              + global
     &              (  MCP*(NCP * (nk-1) + ni-1 ) + 1
     &              :  MCP*(NCP * (nk-1) + ni),:) 
     &	            + global
     &              (    MCP*(nk * NCP - Q + ni -1) +1
     &              :    MCP*(nk * NCP - Q + ni),:)
               global
     &              (    MCP*(nk * NCP - Q + ni -1) +1
     &              :    MCP*(nk * NCP - Q + ni),:) = 0d0
            enddo
         enddo

      elseif (cflag.eq.'out') then

         do ni = 1, Q           ! Master = Slave + Master
            do nk = 1, OCP
               global
     &              (   MCP*(nk * NCP - Q + ni -1) +1
     &              :   MCP*(nk * NCP - Q + ni),:) = 
     &              global
     &              (   MCP*(NCP * (nk-1) + ni-1 ) + 1
     &              :   MCP*(NCP * (nk-1) + ni),:)
            enddo
         enddo

         do ni = 1,Q            ! First Q planes (slave)
            do nk = 1, OCP
               tmpv = 0d+0
               do nj = ni,Q
                  tmpv = tmpv
     &                 + TvM(ni+Q,nj)
     &                 * global
     &                 (      MCP*(NCP * (nk-1) + nj-1 ) + 1
     &                 :      MCP*(NCP * (nk-1) + nj),:) 
               enddo
               global(MCP*(NCP * (nk-1) + ni-1) + 1
     &              : MCP*(NCP * (nk-1) + ni),:)  = tmpv
            enddo
         enddo

         do ni = Q,1,-1         ! Last Q planes (master)
            do nk = 1, OCP
               tmpv = 0d+0
               do nj = 1,ni
                  tmpv = tmpv
     &                 + TvS(ni,nj)
     &                 * global
     &                 (    MCP*(nk * NCP - Q + nj -1) +1
     &                 :    MCP*(nk * NCP - Q + nj),:)
               enddo
               global
     &              (MCP*(nk * NCP - Q + ni -1) +1
     &              :MCP*(nk * NCP - Q + ni),:) = tmpv
            enddo

         enddo

      endif

      return
      end



!###########################################################################

      subroutine shuffle_w(
     &     global,
     &     TwM,TwS,             ! Restriction operators
     &     lprt,                ! 0: S, 1: M+S, 2: M
     &     R, NNODZ, MCP, NCP, OCP, NSD,
     &     cflag,var)

      implicit none

      integer ni, nj, counti, countj
      integer lprt, R, NNODZ, NSD, MCP, NCP, OCP
      real*8 global(NNODZ,NSD)
      real*8 tmp2(MCP*NCP,NSD)
      real*8, dimension(2*R,R) :: TwM, TwS
      character*3 cflag
      character var

!     Control points on plane are contiguous
      if (cflag.eq.'in ') then
         do ni = R,1,-1         ! First R planes (master)
            tmp2 = 0d+0
            do nj = 1, ni
               tmp2 = tmp2
     &              + TwM(nj+R,ni)
     &              * global(MCP*NCP*(nj-1)+1:MCP*NCP*nj,:) 
            enddo
            global(MCP*NCP*(ni-1)+1:MCP*NCP*ni,:) = tmp2
         enddo

         do ni = 1,R            ! Last R planes (slave)
            tmp2 = 0d+0
            do nj = ni,R
               tmp2 = tmp2
     &              + TwS(nj,ni)
     &              * global(MCP*NCP*(OCP-R+nj-1)+1
     &              :        MCP*NCP*(OCP-R+nj),:) 
            enddo
            global(MCP*NCP*(OCP-R+ni-1)+1
     &           : MCP*NCP*(OCP-R+ni),:) = tmp2
         enddo

         do ni = 1, R
            global(MCP*NCP*(ni-1)+1:MCP*NCP*ni,:) = 
     &           + global(MCP*NCP*(ni-1)+1:MCP*NCP*ni,:) 
     &           + global(MCP*NCP*(OCP-R+ni-1)+1
     &           :        MCP*NCP*(OCP-R+ni),:)
            global(MCP*NCP*(OCP-R+ni-1)+1
     &           : MCP*NCP*(OCP-R+ni),:) = 0d0
         enddo

      elseif (cflag.eq.'out') then

         do ni = 1, R
            global(MCP*NCP*(OCP-R+ni-1)+1
     &           : MCP*NCP*(OCP-R+ni),:) = 
     &           global(MCP*NCP*(ni-1)+1:MCP*NCP*ni,:) 
         enddo

         do ni = 1,R            ! First R planes (master)
            tmp2 = 0d+0
            do nj = ni,R
               tmp2 = tmp2
     &              + TwM(ni+R,nj)
     &              * global(MCP*NCP*(nj-1)+1:MCP*NCP*nj,:) 
            enddo
            global(MCP*NCP*(ni-1)+1:MCP*NCP*ni,:) = tmp2
         enddo

         do ni = R,1,-1         ! Last R planes (slave)
            tmp2 = 0d+0
            do nj = 1, ni
               tmp2 = tmp2
     &              + TwS(ni,nj)
     &              * global(MCP*NCP*(OCP-R+nj-1)+1
     &              :        MCP*NCP*(OCP-R+nj),:) 
            enddo
            global(MCP*NCP*(OCP-R+ni-1)+1
     &           : MCP*NCP*(OCP-R+ni),:) = tmp2
         enddo

      endif

      return
      end


