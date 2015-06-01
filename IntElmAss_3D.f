      subroutine IntElmAss_3D(
     &     ygAlpha, acgAlpha, 
     &     VEL_NET,
     &     colu, rowu)
      
      use aAdjKeep
      use common
!      use mpi

      implicit none

c...  Local variables
      integer iel, igauss, jgauss, kgauss, i, j, k, icount

      integer
     &     rowu (NNODZu*8*(Pu+1)*(Qu+1)*(Ru+1)),
     &     colu (NNODZu+1),
     &     ni, nj, nk

      real*8
     &     shgu(NSHLu),
     &     shgradgu(NSHLu,NSD),
     &     shhessgu(NSHLu,NSD,NSD)

      real*8 gp(NGAUSS), gw(NGAUSS), gwt, da

      real*8
     &     xl(NSHLu,NSD),
     &     ul(NSHLu,NDOF),
     &     acl(NSHLu,NDOF),
     &     chil(NSHLu,NSD),
     &     dxidx(NSD,NSD)

      real*8 
     &     Rhsu(NDOF,NSHLu),
     &     xKebe(NDOF*NDOF,NSHLu,NSHLu)

      real*8 
     &     ygAlpha (NNODZu,NDOF),
     &     acgAlpha(NNODZu,NDOF)

      real*8 
     &     xi(NSD),
     &     ui(NDOF),
     &     aci(NDOF),
     &     chii(NSD),
     &     duidxi(NDOF,NSD),
     &     duidxixj(NDOF,NSD,NSD),
     &     rLi(NDOF),
     &     tauM(NDOF,NDOF)

      real*8 VEL_NET(NNODZu,NSD)

      real*8 grade2energy, lapu1, lapu2
      real*8 u1, u2, v1, v2, th, tau
      integer solf
      character*30 fname

c...  get Gaussian points and weights
      gp    = 0d0
      gw    = 0d0
      call genGPandGW(gp, gw, NGAUSS)

!---------------------------------------------------------------
!     Compute element matrices and residuals
!---------------------------------------------------------------

      Ener = 0d0
      AvgTemp = 0d0	  

      AvgSxx = 0d0
      AvgSxy = 0d0
      AvgSyy = 0d0

      AvgExx = 0d0
      AvgExy = 0d0
      AvgEyy = 0d0
      havg = 0d0
      hmax =-1d99
      hmin = 1d99

c...  loop over elements
      do iel = 1, NELu
         
c...  Check to see if current element has nonzero area
         ni = INNu(IENu(iel,1),1) ! get NURB coordinates
         nj = INNu(IENu(iel,1),2)
         nk = INNu(IENu(iel,1),3)
         
         if ( (U_KNOTu(ni).ne.U_KNOTu(ni+1)).and.
     &        (V_KNOTu(nj).ne.V_KNOTu(nj+1)).and.
     &        (W_KNOTu(nk).ne.W_KNOTu(nk+1))) then
            
            da = (U_KNOTu(ni+1) - U_KNOTu(ni))*
     &           (V_KNOTu(nj+1) - V_KNOTu(nj))*
     &           (W_KNOTu(nk+1) - W_KNOTu(nk))/8d0 
            ! used in calculating
            ! quadrature points. The factor of 8d0
            ! comes from mapping from the [-1,1]
            ! line onto a real segment...
            
            do i = 1, NSHLu
               ul  (i,:) = ygAlpha (IENu(iel,i),:)
               acl (i,:) = acgAlpha(IENu(iel,i),:)
               chil(i,:) = VEL_NET (IENu(iel,i),:)
            enddo
            
            icount = 0
            
            do k = 0, Ru
               do j = 0, Qu
                  do i = 0, Pu
                     icount = icount + 1
                     xl(icount,:) = B_NETu(ni-i,nj-j,nk-k,1:NSD)
                  enddo
               enddo
            enddo
            
            xKebe  = 0d0   ! initialize local stiffness matrix
            Rhsu = 0d0     ! initialize local load vector
            
c...  Loop over integration points (NGAUSS in each direction)
            do igauss = 1, NGAUSS
               do jgauss = 1, NGAUSS
                  do kgauss = 1, 1 !NGAUSS
                     
c...  Get Element Shape functions and their gradients
                     
                     shgu     = 0d0 ! initialize
                     shgradgu = 0d0
                     shhessgu = 0d0
                     
                     call eval_SHAPE_3D(
     &                    iel,
     &                    gp(igauss),gp(jgauss),gp(kgauss),
     &                    shgu, shgradgu, shhessgu, dxidx)
                     
      ! Compute values at integration points
                     call e3int(
     &                    ul, acl, xl, chil,
     &                    shgu, shgradgu, shhessgu,
     &                    ui, aci, xi, chii, grade2energy,
     &                    duidxi, duidxixj, rLi)
                     
      ! Compute intrinsic-time scales
c                    call e3STAB_3D(dxidx, duidxi, ui, tauM)
                     tauM = 0d0
                     
c...  Calculate given element stiffness matrix and force vector
                     gwt = gw(igauss)*gw(jgauss)*gw(kgauss)*da

                     call e3Rhs_3D (
     &                      xi, ui, aci, rLi, chii,
     &                      duidxi, duidxixj, gwt,
     &                      shgu, shgradgu, shhessgu,
     &                      grade2energy,lapu1, lapu2,
     &                      Rhsu)
                     
                     call Energy (ui, duidxi, grade2energy, gwt)
c                     call Elliptic(ui, laprho, gwt)
                     
                     call e3LHS_3D (
     &                     ui, rLi, chii,
     &                     duidxi, duidxixj, gwt,
     &                     shgu, shgradgu, shhessgu,
     &                     grade2energy,lapu1, lapu2,
     &                     u1, u2, v1, v2, th,
     &                     xKebe)
                     
                     
                  enddo
               enddo
            enddo

  
c...  Modify local matrices and vectors for Dirichlet boundary conditions
!      if((periodicx.eq.0).or.(periodicy.eq.0).or.(periodicz.eq.0)) then !non-periodic
            call BCLhs_3D(iel, xKebe, Rhsu)
!      endif

  
!      print*, 'xKebe'
!            write(*,*) xKebe
!      print*, 'Rhsu'
!            write(*,*) Rhsu

c...  assemble into the Sparse Global Stiffness Matrix and Rhs Vector
                                ! Assemble stiffness
            call FillSparseMat_3D_u (iel, colu,  rowu,  xKebe)
                                ! Assemble load vector
            call LocaltoGlobal_3D_u (iel, Rhsu)
            
         endif
         
      enddo
      
      
      return
      end

