      subroutine SparseGMRES_up(
     &     lhsK,
     &     Utol, colu, rowu,
     &     rhsGu,
     &     solu,
     &     Kspaceu, Kspaceu_mn, Ksp_usd,
     &     icntu,
     &     NNODZu,
     &     Pu, Qu, Ru, NDOF,
     &     TuMu, TuSu, TvMu, TvSu, TwMu, TwSu,
     &     MCPu, NCPu, OCPu,
     &     lpu,lpv,lpw
     &     )


      implicit none

      real*8 HBrg(Kspaceu+1,Kspaceu)
      real*8 Rcos(Kspaceu), Rsin(Kspaceu)
      real*8 lhsKBdiag(NNODZu,NDOF)
      real*8 uBrg1(NNODZu,NDOF,Kspaceu+1)

      real unorm_ref

      integer
     &     n, i, j, k, iK, iKs, jK, lK, kK, Kspaceu, count, idof,
     &     icntu, Kspaceu_mn,Ksp_usd, k_f, k_f2,
     &     NNODZu, Pu, Qu, Ru, NSD, NDOF,
     &     is, lenseg,
     &     MCPu, NCPu, OCPu,
     &     lpu,lpv,lpw

      integer, dimension(NNODZu+1):: colu
      integer, dimension(NNODZu*8*(Pu+1)*(Qu+1)*(Ru+1)):: rowu

      real*8
     &     rhstmp1(NNODZu,NDOF),
     &     eBrg(Kspaceu+1),
     &     temp1u(NNODZu,NDOF),
     &     yBrg(Kspaceu),
     &     rhsGu(NNODZu,NDOF), solu(NNODZu,NDOF),
     &     rr, unorm, 
     &     rru, rrul,
     &     epsnrm, beta, rrglob,
     &     ercheck, tmp, tmp1, Utol,
     &     lhsK(NDOF*NDOF,icntu),
     &     Binv(NDOF,NDOF)

      real*8
     &     TuMu(2*Pu,Pu), TuSu(2*Pu,Pu),
     &     TvMu(2*Qu,Qu), TvSu(2*Qu,Qu),
     &     TwMu(2*Ru,Ru), TwSu(2*Ru,Ru)

      rhstmp1 = RHSGu     ! Algorithm from Shakib/Johan's theses

!------------- zero out common values --------------
      lhsKBdiag = 0d0
      uBrg1     = 0d0
      Rcos  = 0d0
      Rsin  = 0d0
      HBrg  = 0d0

!---------------   Block Diagonal Preconditioner   ---------------

      ! Form Block Diagonal components
      lhsKBdiag = 0d0
      do i = 1, NNODZu
         do j = colu(i), colu(i+1)-1
            n = rowu(j)
            if ( n==i ) then

               lhsKBdiag(i,:) = LHSK(1:NDOF*NDOF:NDOF+1,j)
               
            endif
         enddo

      enddo

      ! communicate the block-diagonal ! GLLG
      call shuffle_u(
     &     lhsKBdiag,
     &     TuMu,TuSu, lpu, Pu,
     &     NNODZu, MCPu, NCPu, OCPu, NDOF,
     &     'in ','u')

      call shuffle_v(
     &     lhsKBdiag,
     &     TvMu,TvSu,lpv,Qu,
     &     NNODZu, MCPu, NCPu, OCPu, NDOF,
     &     'in ','u')

      call shuffle_w(
     &     lhsKBdiag,
     &     TwMu,TwSu,lpw,Ru,
     &     NNODZu, MCPu, NCPu, OCPu, NDOF,
     &     'in ','u')

      call shuffle_w(
     &     lhsKBdiag,
     &     TwMu,TwSu,lpw,Ru,
     &     NNODZu, MCPu, NCPu, OCPu, NDOF,
     &     'out','u')

      call shuffle_v(
     &     lhsKBdiag,
     &     TvMu,TvSu,lpv, Qu,
     &     NNODZu, MCPu, NCPu, OCPu, NDOF,
     &     'out','u')

      call shuffle_u(
     &     lhsKBdiag,
     &     TuMu,TuSu,lpu,Pu,
     &     NNODZu, MCPu, NCPu, OCPu, NDOF,
     &     'out','u')


      do i= 1, NNODZu
        lhsKBdiag(i,:) = 1d0/lhsKBdiag(i,:)
c        print*, 'lhs', lhsKBdiag(i,:)
        if (dabs(sum(lhsKBdiag(i,:)*lhsKBdiag(i,:)))<1d-8) then
         print*,'diag is almost zero',lhsKBdiag(i,:)
        endif
      enddo

      !  Precondition momentum residual
      uBrg1(:,1,1) = lhsKBdiag(:,1)*rhstmp1(:,1)
      uBrg1(:,2,1) = lhsKBdiag(:,2)*rhstmp1(:,2)
      uBrg1(:,3,1) = lhsKBdiag(:,3)*rhstmp1(:,3)
      uBrg1(:,4,1) = lhsKBdiag(:,4)*rhstmp1(:,4)
      uBrg1(:,5,1) = lhsKBdiag(:,5)*rhstmp1(:,5)

!--------------------------------------------------------------------------
      rhstmp1 = uBrg1(:,:,1)

      rr =   sum(rhstmp1(:,1)*rhstmp1(:,1)) ! calculate norm of rhs
     &     + sum(rhstmp1(:,2)*rhstmp1(:,2))
     &     + sum(rhstmp1(:,3)*rhstmp1(:,3))
     &     + sum(rhstmp1(:,4)*rhstmp1(:,4))
     &     + sum(rhstmp1(:,5)*rhstmp1(:,5))

c      if (numnodes .gt. 1) then
         rrglob = rr
c      endif

      unorm = sqrt(rr)

      epsnrm    = Utol * unorm !  set up tolerance
                                ! set up RHS of the Hessenberg's problem
      eBrg    = 0d0
      eBrg(1) = unorm

      unorm = 1d0/unorm         ! normalize the first Krylov vector
      uBrg1(:,:,1) = uBrg1(:,:,1) * unorm

      unorm_ref = unorm*1d2

      ! loop through GMRES iterations
      Kspu:do iK = 1, Kspaceu

         iKs = iK

      !  matrix-vector product
         temp1u = uBrg1(:,:,iKs)

      ! Periodicity (Slave = Master)  - LG
         call shuffle_w(
     &        temp1u,
     &        TwMu,TwSu,lpw,Ru,
     &        NNODZu, MCPu, NCPu, OCPu, NDOF,
     &        'out','u')

         call shuffle_v(
     &        temp1u,
     &        TvMu,TvSu,lpv, Qu,
     &        NNODZu, MCPu, NCPu, OCPu, NDOF,
     &        'out','u')

         call shuffle_u(
     &        temp1u,
     &        TuMu,TuSu,lpu,Pu,
     &        NNODZu, MCPu, NCPu, OCPu, NDOF,
     &        'out','u')

         call SparseProdFull_3D( !  Product (forms new Krylov vector)
     &        lhsK, colu, rowu,
     &        temp1u, uBrg1(:,:,iKs+1),
     &        NNODZu, Pu, Qu, Ru, NDOF, icntu)

!     Communicate to Masters, Zero out Slaves - GLLG
!----------------------------------------------------------------------

      !  VELOCITY
         call shuffle_u(
     &        uBrg1(:,:,iKs+1),
     &        TuMu,TuSu,lpu,Pu,
     &        NNODZu, MCPu, NCPu, OCPu, NDOF,
     &        'in ','u')

         call shuffle_v(
     &        uBrg1(:,:,iKs+1),
     &        TvMu,TvSu,lpv,Qu,
     &        NNODZu, MCPu, NCPu, OCPu, NDOF,
     &        'in ','u')

         call shuffle_w(
     &        uBrg1(:,:,iKs+1),
     &        TwMu,TwSu,lpw,Ru,
     &        NNODZu, MCPu, NCPu, OCPu, NDOF,
     &        'in ','u')


c-------------------------------------------------------------------------
      !  Precondition momentum residual
         rhstmp1 = uBrg1(:,:,iKs+1)

         uBrg1(:,1,iKs+1) = lhsKBdiag(:,1)*rhstmp1(:,1)
         uBrg1(:,2,iKs+1) = lhsKBdiag(:,2)*rhstmp1(:,2)
         uBrg1(:,3,iKs+1) = lhsKBdiag(:,3)*rhstmp1(:,3)
         uBrg1(:,4,iKs+1) = lhsKBdiag(:,4)*rhstmp1(:,4)
         uBrg1(:,5,iKs+1) = lhsKBdiag(:,5)*rhstmp1(:,5)

!--------------------------------------------------------------------------
         do jK = 1, iKs+1       ! orthogonalize and get the norm

            if (jK .eq. 1) then

               rr = 0d0
               do i = 1, NDOF    !{u_{i+1}*u_1} vector
                  rr = rr + sum(uBrg1(:,i,iKs+1)*uBrg1(:,i,1))
               enddo

c               if (numnodes .gt. 1) then
                  rrglob = rr
c               endif

               beta = rr

            else
                                !     project off jK-1 vector
               uBrg1(:,:,iKs+1) =
     &              uBrg1(:,:,iKs+1) - beta * uBrg1(:,:,jK-1)

               rr = 0d0
               do i = 1, NDOF    !{u_{i+1}*u_1} vector
                  rr = rr + sum(uBrg1(:,i,iKs+1)*uBrg1(:,i,jK))
               enddo

c               if (numnodes .gt. 1) then
                  rrglob = rr
c               endif

               beta = rr

            endif

            HBrg(jK,iKs) = beta ! put this in the Hessenberg Matrix

         enddo

         unorm           = sqrt(beta)
         HBrg(iKs+1,iKs) = unorm ! this fills the 1 sub diagonal band

!     normalize the Krylov vector
         unorm = 1d0/unorm      ! normalize the next Krylov vector
         uBrg1(:,:,iKs+1) = uBrg1(:,:,iKs+1) * unorm

c.... construct and reduce the Hessenberg Matrix
c     since there is only one subdiagonal we can use a Givens rotation
c     to rotate off each subdiagonal AS IT IS FORMED. We do this because it
c     allows us to check progress of solution and quit when satisfied.  Note
c     that all future K vects will put a subdiagonal in the next column so
c     there is no penalty to work ahead as  the rotation for the next vector
c     will be unaffected by this rotation.

c     H Y = E ========>   R_i H Y = R_i E

         do jK = 1, iKs-1
            tmp =  Rcos(jK) * HBrg(jK,iKs) + Rsin(jK) * HBrg(jK+1,iKs)
            HBrg(jK+1,iKs) = 
     &           - Rsin(jK) * HBrg(jK,  iKs)
     &           + Rcos(jK) * HBrg(jK+1,iKs)
            HBrg(jK,  iKs) =  tmp
         enddo

         tmp  = sqrt(HBrg(iKs,iKs)**2 + HBrg(iKs+1,iKs)**2)
         tmp1 = 1d0/tmp
         Rcos(iKs) = HBrg(iKs,  iKs) * tmp1
         Rsin(iKs) = HBrg(iKs+1,iKs) * tmp1
         HBrg(iKs,  iKs) = tmp
         HBrg(iKs+1,iKs) = 0d0

      ! rotate eBrg    R_i E
         tmp         = + Rcos(iKs) * eBrg(iKs) + Rsin(iKs) * eBrg(iKs+1)
         eBrg(iKs+1) = - Rsin(iKs) * eBrg(iKs) + Rcos(iKs) * eBrg(iKs+1)
         eBrg(iKs)   = tmp

      ! check for convergence
         ercheck=abs(eBrg(iKs+1))

         if (mod(iKs,100).eq.0) then
c            if (myid.eq.mpi_master) 
            write(*,8000) iKs, ercheck,epsnrm, ercheck*unorm_ref
         endif 

         if (ercheck .le. epsnrm .and. iK >= Kspaceu_mn) exit


      enddo Kspu ! end of GMRES iteration loop

      Ksp_usd = iK 


!     ------------------------->   Solution   <------------------------
!     if converged or end of Krylov space

      ! solve for yBrg
      do jK = iKs, 1, -1
         yBrg(jK)     = eBrg(jK) / HBrg(jK,jK)
         eBrg(1:jK-1) = eBrg(1:jK-1) - yBrg(jK) * HBrg(1:jK-1,jK)
      enddo

      ! update Dy
      do jK = 1, iKs
         solu = solu + yBrg(jK) * uBrg1(:,:,jK)
      enddo
     
      !  communicate solution - LG
      !  VELOCITY
         call shuffle_w(
     &        solu,
     &        TwMu,TwSu,lpw,Ru,
     &        NNODZu, MCPu, NCPu, OCPu, NDOF,
     &        'out','u')

         call shuffle_v(
     &        solu,
     &        TvMu,TvSu,lpv,Qu,
     &        NNODZu, MCPu, NCPu, OCPu, NDOF,
     &        'out','u')

         call shuffle_u(
     &        solu,
     &        TuMu,TuSu,lpu,Pu,
     &        NNODZu, MCPu, NCPu, OCPu, NDOF,
     &        'out','u')
    

      write(*,9000) iKs, ercheck, ercheck*unorm_ref

      return

 2000 format(i10, f17.12)
 8000 format(
     &     10x,' number of GMRES iterations:', 2x,i10,/,
     &     10x,'              residual goal:', 2x,e11.4,/,
     &     10x,'             residual value:', 2x,e11.4,/,
     &     10x,'       reduction percentage:', 2x,f10.6,/ )
 9000 format(
     &     10x,' number of GMRES iterations:', 2x,i10,/,
     &     10x,'       residual final value:', 2x,e11.4,/,
     &     10x,'       reduction percentage:', 2x,f10.6,/ )
      
      end


