      subroutine solflowGMRES(
     &     rowu,  colu, icntu)
      
      use aAdjKeep
      use common

      implicit none
      
      real*8
     &     ygAlpha (NNODZu,NDOF),
     &     acgAlpha(NNODZu,NDOF),
     &     fact1, fact2,
     &     timeAlpha, timeFinal,RHSGu1(NNODZu,NDOF)

      real*8
     &     momresl, momres, resm0, resd0,
     &     Tempresl, Tempres, rest0,
     &     res, resl, dumresl, dumres

      real*8      ! rotation
     &     VEL_NET(NNODZu,NSD)

      integer
     &     rowu(NNODZu*8*(Pu+1)*(Qu+1)*(Ru+1)),
     &     colu(NNODZu+1),
     &     inewt, icntu, icntp, icntpp,
     &     i, j, itask, is, lenseg, n, solf,
     &     inod, idof, timealpload, timealpunload

      character*30 fname
	  
      timeAlpha = time + alfi*Delt
      timeFinal = time + Delt

      fact1 = gami*Delt
      fact2 = alfi*gami*Delt

      VEL_NET = 0d0

      newton:do inewt = 1, Nnewt   ! Newton iteration
         
         RHSGu = 0d0               ! initialize
         LHSK  = 0d0

!!!!!!!!!!!!!!!!!!!!!!!

      if((timeAlpha.le.tloadT)) then
       timealpload = 1
       timealpunload = 0
      endif
      if((timeAlpha>tloadT)) then
       timealpload = 0
       timealpunload = 1
      endif

      if((pullxplusu1==2)) then ! if pulled in X-direction
        VXnp1=1d0*(velT*timeFinal*timealpload+velT
     &         *tloadbyunloadT*(ttotT-timeFinal)*timealpunload)
       endif

            ! ADJUST for Strong BCs on n+1 levels
            where (IBC==2) ! IBC = 2 pull; =1 constraint
              yg  = VXnp1*DIR_BC                  
              acg = (yg - ygold)/(gami*Delt) + (gami-1d0)/gami*acgold
            end where

            where (IBC==1) ! IBC = 2 pull; =1 constraint
              yg  = DIR_BC
              acg = 0d0
            end where

        ygAlpha  = ygold  + alfi*( yg -  ygold)
        acgAlpha = acgold + almi*(acg - acgold)

!!!!!!!!!!!!!!!!!!!!!!!
      call IntElmAss_3D(
     &        ygAlpha, acgAlpha,
     &        VEL_NET,
     &       colu,  rowu)
	 
      RHSGu1 = RHSGu
	  
      if (numnodes.gt.1) then
            resl = Ener
      endif


      ! Communicate Res - GLLG
!----------------------------------------------------------------------------
      ! Master for MPI collects RHS
         call shuffle_u(
     &        RHSGu,
     &        TuMu,TuSu,lpu, Pu,
     &        NNODZu, MCPu, NCPu, OCPu, NDOF,
     &        'in ','u')


	  
         call shuffle_v(
     &        RHSGu,
     &        TvMu,TvSu,lpv, Qu,
     &        NNODZu, MCPu, NCPu, OCPu, NDOF,
     &        'in ','u')


	  
         call shuffle_w(
     &        RHSGu,
     &        TwMu,TwSu,lpw, Ru,
     &        NNODZu, MCPu, NCPu, OCPu, NDOF,
     &        'in ','u')

 
	  
!----------------------------------------------------------------------------
      ! Compute  Dummy variable Residual Norm
         dumresl = 0d0
         dumresl = sum(RHSGu(:,1)*RHSGu(:,1))
     &           + sum(RHSGu(:,2)*RHSGu(:,2))

      ! Compute  Momentum Residual Norm
         momresl = 0d0
         momresl = sum(RHSGu(:,3)*RHSGu(:,3))
     &           + sum(RHSGu(:,4)*RHSGu(:,4))

      ! Temperature Norm
         Tempresl = 0d0
         Tempresl = sum(RHSGu(:,5)*RHSGu(:,5))

      ! Field variables DUMMY u1, u2 residual
         if (numnodes.gt.1) then
            resl = dumresl
            dumres = sqrt(res) !I think this should be resl
         else
            dumres = sqrt(dumresl)
         endif

      ! Velocity v1, v2 residual
         if (numnodes.gt.1) then
            resl = momresl
            momres = sqrt(res)
         else
            momres = sqrt(momresl)
         endif

      ! Temperature residual
         if (numnodes.gt.1) then
            resl = Tempresl
            Tempres = sqrt(res)
         else
            Tempres = sqrt(Tempresl)
         endif

         if (inewt == 1) resm0 = 1d2/momres
         if (inewt == 1) resd0 = 1d2/dumres
         if (inewt == 1) rest0 = 1d2/Tempres

            write(*,'(a,i10)')  "Newton Iteration Number  = ", inewt
            write(*,'(a,e11.4)')"Dummy Residual Norm      = ", dumres
            write(*,'(a,e11.4)')"Mom. Residual Norm       = ", momres
            write(*,'(a,e11.4)')"Temp. Residual Norm      = ", Tempres

            if (inewt > 1) then
               write(*,'(a,f10.4)') "   DUMMY Var.: It. Red. (%) = ",
     &              dumres*resd0
               write(*,'(a,f10.4)') "     MOMENTUM: It. Red. (%) = ",
     &              momres*resm0
               write(*,'(a,f10.4)') "     TEMP CON: It. Red. (%) = ",
     &              Tempres*rest0
            endif

         if ( inewt > 1) then
            if ((momres*resm0<NLtol).and.(Tempres*rest0<NLtol).and.
     &           (dumres*resd0<NLtol)) then
               print*,'- Newton steps converged -'
               return
            endif
         endif
         
!---------------------------------------------------------------------------
         ygAlpha  = 0d0         ! Use these for increments
         acgAlpha = 0d0

		! output fo GMRES is acg i.e. velocities
         call SparseGMRES_up(
     &        lhsK,
     &        Utol, colu, rowu,
     &        rhsGu,
     &        acgAlpha,
     &        Kspaceu, Kspaceu_mn, Ksp_usd,
     &        icntu,
     &        NNODZu,
     &        Pu, Qu, Ru, NDOF,
     &        TuMu, TuSu, TvMu, TvSu, TwMu, TwSu,
     &        MCPu, NCPu, OCPu,
     &        lpu,lpv,lpw
     &        )

                                ! Update Current Values         
         acg = acg + acgAlpha
         yg  = yg  + fact1*acgAlpha

      enddo newton


      return
      end

!###################################################################

      subroutine rot2D(Tim, rpm, xor,
     &                 B_NETini, B_NETu, VEL_NET,
     &                 MCPu, NCPu, OCPu, NSD)

      implicit none 

      integer i, j, k, l, m, count

      integer MCPu, NCPu, OCPu, NSD

      real*8 Tim, rpm, xor(2), xl(NSD), beta, pi
      
      real*8 B_NETini(MCPu,NCPu,OCPu,NSD+1),
     &       VEL_NET (MCPu*NCPu*OCPu,NSD  ),
     &       B_NETu  (MCPu,NCPu,OCPu,NSD+1)

      real*8 B_TMP1(NSD), 
     &       B_TMP2(NSD)

      real*8 RotMat (NSD,NSD), 
     &       DRotMat(NSD,NSD)

      pi = acos(-1d0)
      
      beta = 2d0*pi*rpm*Tim   ! \beta = 2*pi*f*t_alpha
                              !  rpm  = revolutions per second
      
      RotMat(:,:) =  0d0      ! Rotation matrix R(t_{n+\alpha_f})
      RotMat(1,1) =  cos(beta)
      RotMat(1,2) =  sin(beta)
      RotMat(2,1) = -sin(beta)
      RotMat(2,2) =  cos(beta)
      RotMat(3,3) =  1d0
      
      
      DRotMat(:,:) = 0d0       ! \dot{R}
      DRotMat(1,1) = -2d0*pi*rpm*sin(beta)
      DRotMat(1,2) =  2d0*pi*rpm*cos(beta)
      DRotMat(2,1) = -2d0*pi*rpm*cos(beta)
      DRotMat(2,2) = -2d0*pi*rpm*sin(beta)
      
      VEL_NET = 0d0
      
      count = 0
      
      do k = 1, OCPu
        do j = 1, NCPu
          do i = 1, MCPu
            
            count = count+1
            
            B_TMP1(:) = 0d0    ! Temp
            B_TMP2(:) = 0d0    ! Temp
            
            xl(1) = B_NETini(i,j,k,1)-xor(1)
            xl(2) = B_NETini(i,j,k,2)-xor(2)
            xl(3) = B_NETini(i,j,k,3)
            
            do l = 1,NSD
              do m = 1,NSD
                B_TMP1(l) = B_TMP1(l) +  RotMat(l,m)*xl(m)
                B_TMP2(l) = B_TMP2(l) + DRotMat(l,m)*xl(m)
              enddo
            enddo
            
            B_NETu (i,j,k,1) = B_TMP1(1) ! Mesh position at alfi
            B_NETu (i,j,k,2) = B_TMP1(2) ! Mesh position at alfi
            B_NETu (i,j,k,3) = B_TMP1(3) ! Mesh position at alfi

            VEL_NET(count,:) = B_TMP2(:) ! Mesh velocity at alfi
            
          enddo
        enddo
      enddo

      return
      end


!###################################################################

      subroutine rotBC2D(Tim, rpm, xor,
     &                   B_NETini,
     &                   MCPu, NCPu, OCPu,
     &                   NSD, NNODZu, NDOF,
     &                   IBC, vel, ace)

      implicit none

      integer i, j, k, l, m, count, idof
      integer MCPu, NCPu, OCPu, 
     &        NNODZu, NDOF, NSD
      real*8 Tim, rpm, beta, pi, xor(2), xl(NSD) 
      real*8
     &     B_NETini (MCPu,NCPu,OCPu,NSD+1)
      real*8
     &     B_TMP1(NSD),
     &     B_TMP2(NSD)
      real*8 
     &     RotMat(NSD,NSD),
     &     DRotMat(NSD,NSD),
     &     DDRotMat(NSD,NSD)

      real*8 vel (NNODZu,NDOF), 
     &       ace(NNODZu,NDOF)
      integer*4 IBC(NNODZu,NDOF)

      
      pi = acos(-1d0)
      
      beta = 2d0*pi*rpm*Tim     ! \beta = 2*pi*f*t_alpha
      
      RotMat(:,:) =  0d0        ! Rotation matrix R(t_{n+1})
      RotMat(1,1) =  cos(beta)
      RotMat(1,2) =  sin(beta)
      RotMat(2,1) = -sin(beta)
      RotMat(2,2) =  cos(beta)
      RotMat(3,3) =  1d0      
      
      DRotMat(:,:) =  0d0       ! \dot{R}
      DRotMat(1,1) = -2d0*pi*rpm*sin(beta)
      DRotMat(1,2) =  2d0*pi*rpm*cos(beta)
      DRotMat(2,1) = -2d0*pi*rpm*cos(beta)
      DRotMat(2,2) = -2d0*pi*rpm*sin(beta)
      
      DDRotMat(:,:) =  0d0      ! \ddot{R}
      DDRotMat(1,1) = -(2d0*pi*rpm)*(2d0*pi*rpm)*cos(beta)
      DDRotMat(1,2) = -(2d0*pi*rpm)*(2d0*pi*rpm)*sin(beta)
      DDRotMat(2,1) =  (2d0*pi*rpm)*(2d0*pi*rpm)*sin(beta)
      DDRotMat(2,2) = -(2d0*pi*rpm)*(2d0*pi*rpm)*cos(beta) 

      count = 0
      
      do k = 1, OCPu
         do j = 1, NCPu
            do i = 1, MCPu
               
               count = count+1
               
               B_TMP1(:) = 0d+0 ! Temp
               B_TMP2(:) = 0d+0 ! Temp
               
               xl(1) = B_NETini(i,j,k,1)-xor(1)
               xl(2) = B_NETini(i,j,k,2)-xor(2)
               xl(3) = B_NETini(i,j,k,3)
               
               do l = 1,NSD
                  do m = 1,NSD
                     B_TMP1(l) = B_TMP1(l) +  DRotMat(l,m)*xl(m) ! velocidad
                     B_TMP2(l) = B_TMP2(l) + DDRotMat(l,m)*xl(m) ! aceleracion
                  enddo
               enddo

               
               if(j==NCPu) then
                  do idof = 2, NDOF
                     if (IBC(count,idof).eq.1) then
                        vel(count,idof) = 0d0 
                        ace(count,idof) = 0d0 
                     endif
                  enddo
               endif
               
               if(j==1) then
                  do idof = 2, NDOF
                     if (IBC(count,idof).eq.1) then
                        vel(count,idof) = B_TMP1(idof-1)
                        ace(count,idof) = B_TMP2(idof-1)
                     endif
                  enddo
               endif
               
            enddo
         enddo
      enddo  


      return
      end
      
      
