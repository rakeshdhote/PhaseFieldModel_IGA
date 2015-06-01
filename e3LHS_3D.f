      subroutine e3LHS_3D (ui, rLi, chii, 
     &                     duidxi, duidxixj, gwt,
     &                     shgu, shgradgu, shhessgu, 
     &                     grade2energy,lapu1, lapu2,
     &                     u1, u2, v1, v2, th,
     &                     xKebe)  
      
      use aAdjKeep
      use common
      implicit none

      integer aa, bb, idof, jdof, ipos

      ! Tangent Matrix xKebe
      real*8 xKebe (NDOF*NDOF,NSHLu,NSHLu)

      ! scalar parameters
      real*8 gwt, DetJgwt,
     &       fact1, fact2

      ! Shape functions and derivatives
      real*8 shgu(NSHLu), shgradgu(NSHLu,NSD), shhessgu(NSHLu,NSD,NSD)
      real*8 grade2energy

      ! Variables and derivatives at Gauss points
      real*8 ui(NDOF), rLi(NDOF), chii(NSD), duidxi(NDOF,NSD), 
     &       duidxixj(NDOF,NSD,NSD)

      real*8 A1(NDOF,NDOF), A2(NDOF,NDOF)  !,
c     &       A2(NDOF,NDOF), A3(NDOF,NDOF)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Variables
      real*8 u1, u2, v1, v2, th, tau
!        real*8 u1ad, u2ad, v1ad, v2ad, thad

      real*8 du1dx, du2dx, dv1dx, dv2dx, dthdx
      real*8 du1dy, du2dy, dv1dy, dv2dy, dthdy
      real*8 du1dxixj(NSD,NSD), du2dxixj(NSD,NSD), dv1dxixj(NSD,NSD),
     &       dv2dxixj(NSD,NSD), dthdxixj(NSD,NSD)
      real*8 e1T, e2T, e3T, se1T, se2T, se3T
      real*8 e1Tdot, e2Tdot, e3Tdot, se2Tdot
      real*8 de2dx, de2dy
      real*8 lapu1, lapu2
      real*8 S11, S12, S21, S22
      real*8 f2f3(NDOF), f4x(NDOF),f4y(NDOF), f5xx(NDOF),f5xy(NDOF),
     &         f5yy(NDOF)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      real*8 delta(NSD,NSD)

!      real*8 u, v, w, p, dpdrho, rho, divu, dote
!      real*8 uad, vad, wad
!
!      real*8 drdx, drdy, drdz
!      real*8 dudx, dudy, dudz
!      real*8 dvdx, dvdy, dvdz
!      real*8 dwdx, dwdy, dwdz

      DetJgwt = DetJ*gwt


c      zero = 0d0   ! zero
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Variables
      u1 = ui(1)
      u2 = ui(2)
      v1 = ui(3)
      v2 = ui(4)
      th = ui(5)

!      u1ad = u1 - chii(1)
!      u2ad = u2 - chii(2)
!      v1ad = v1 - chii(3)
!      v2ad = v2 - chii(4)
!      thad = th - chii(5)

!     Derivatives of ui in x direction
      du1dx = duidxi(1,1)
      du2dx = duidxi(2,1)
      dv1dx = duidxi(3,1)
      dv2dx = duidxi(4,1)
      dthdx = duidxi(5,1)

!     Derivatives of ui in y direction
      du1dy = duidxi(1,2)
      du2dy = duidxi(2,2)
      dv1dy = duidxi(3,2)
      dv2dy = duidxi(4,2)
      dthdy = duidxi(5,2)

!     Derivatives of ui in z direction
!      du1dz = duidxi(1,3)
!      du2dz = duidxi(2,3)
!      dv1dz = duidxi(3,3)
!      dv2dz = duidxi(4,3)
!      dthdz = duidxi(5,3)

!     Second derivatives of ui
      du1dxixj(:,:) = duidxixj (1,:,:)
      du2dxixj(:,:) = duidxixj (2,:,:)
      dv1dxixj(:,:) = duidxixj (3,:,:)
      dv2dxixj(:,:) = duidxixj (4,:,:)
      dthdxixj(:,:) = duidxixj (5,:,:)

!     Laplacians in u1 and u2 dofs
      lapu1 = du1dxixj(1,1)+du1dxixj(2,2)
      lapu2 = du2dxixj(1,1)+du2dxixj(2,2)

!     Strain definations
      e1T = (du1dx+du2dy)/sqrt(2d0)   ! Strains
      e2T = (du1dx-du2dy)/sqrt(2d0)
      e3T = (du1dy+du2dx)/2d0
      se1T = (du1dx+du2dy)            ! rescaled strains
      se2T = (du1dx-du2dy)
      se3T = (du1dy+du2dx)

      e1Tdot = (dv1dx+dv2dy)/sqrt(2d0)
      e2Tdot = (dv1dx-dv2dy)/sqrt(2d0)
      e3Tdot = (dv1dy+dv2dx)/2d0
      se2Tdot = (dv1dx-dv2dy)

!     Temperature Variable
      tau = (th-Tm)/Tm

!     Stress definations
!      S11 =  aa1*se1T+aa2*Tau*se2T-aa4*se2T**3+aa6*se2T**5
!      S12 =  aa3*se3T
!      S21 =  aa3*se3T
!      S22 =  aa1*se1T-aa2*Tau*se2T+aa4*se2T**3-aa6*se2T**5

!      de2dx = (duidxixj(1,1,1)-duidxixj(2,1,2))/sqrt(2d0)
!      de2dy = (duidxixj(1,1,2)-duidxixj(2,2,2))/sqrt(2d0)
!
!      grade2energy = de2dx**2+de2dy**2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Time dependent matrix (initialization)
      A1   = 0d0
!     Populating the A1 matrix cells
      A1(1,1) = 1d0
      A1(2,2) = 1d0
      A1(3,3) = rhoT
      A1(4,4) = rhoT
      A1(5,5) = cvT

!     New variables and forcing functions (initialization)
      A2   = 0d0
!     Populating the A1 matrix cells
      A2(1,3) = -1d0
      A2(2,4) = -1d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Loop over local shape functions in each direction

      fact1  = almi
      fact2  = alfi*gami*Delt

!     Tangent matrix xKebe - time dependent part
!                            + New variables and forcing functions

      do bb = 1, NSHLu
         do aa = 1, NSHLu
            do idof = 1, NDOF
               do jdof = 1, NDOF
                  
                  ipos = NDOF*(idof-1) + jdof
            
        xKebe(ipos,aa,bb) = xKebe(ipos,aa,bb) + 
     &       (
     &          fact1*shgu(aa)*shgu(bb)*A1(idof,jdof) +
     &          fact2*shgu(aa)*shgu(bb)*A2(idof,jdof)
     &       )*DetJgwt

               enddo ! jdof
            enddo ! idof
         enddo ! aa
      enddo ! bb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Tangent matrix xKebe - Non-linear coupling term

      do bb = 1, NSHLu
         do aa = 1, NSHLu

            xKebe(21,aa,bb) = xKebe(21,aa,bb) +
     &           fact2*( shgu(aa)*Ca2th*th*se2Tdot*shgradgu(bb,1)
     &                 )*DetJgwt

            xKebe(22,aa,bb) = xKebe(22,aa,bb) +
     &           fact2*( -shgu(aa)*Ca2th*th*se2Tdot*shgradgu(bb,2)
     &                 )*DetJgwt

            xKebe(23,aa,bb) = xKebe(23,aa,bb) +
     &           fact2*( shgu(aa)*Ca2th*th*se2T*shgradgu(bb,1)
     &                 )*DetJgwt

            xKebe(24,aa,bb) = xKebe(24,aa,bb) +
     &           fact2*( -shgu(aa)*Ca2th*th*se2T*shgradgu(bb,2)
     &                 )*DetJgwt

            xKebe(25,aa,bb) = xKebe(25,aa,bb) +
     &           fact2*( shgu(aa)*Ca2th*se2Tdot*se2T*shgu(bb)
     &                 )*DetJgwt

         enddo ! aa
      enddo ! bb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Tangent matrix xKebe - Stress flux part - X direction A4_x matrix

      do bb = 1, NSHLu
         do aa = 1, NSHLu

            xKebe(11,aa,bb) = xKebe(11,aa,bb) + 
     &           fact2*( shgradgu(aa,1)*shgradgu(bb,1)*
     &                   (aa1+aa2*tau-3d0*aa4*se2T**2+5d0*aa6*se2T**4)
     &                 )*DetJgwt

            xKebe(12,aa,bb) = xKebe(12,aa,bb) + 
     &           fact2*( shgradgu(aa,1)*shgradgu(bb,2)*
     &                   (aa1-aa2*tau+3d0*aa4*se2T**2-5d0*aa6*se2T**4)
     &                 )*DetJgwt

            xKebe(13,aa,bb) = xKebe(13,aa,bb) +
     &           fact2*( shgradgu(aa,1)*shgradgu(bb,1)*etaT
     &                 )*DetJgwt

            xKebe(15,aa,bb) = xKebe(15,aa,bb) + 
     &           fact2*( shgradgu(aa,1)*se2T*shgu(bb)*aa2/Tm
     &                 )*DetJgwt

            xKebe(16,aa,bb) = xKebe(16,aa,bb) +
     &           fact2*( shgradgu(aa,1)*shgradgu(bb,2)*aa3
     &                 )*DetJgwt

            xKebe(17,aa,bb) = xKebe(17,aa,bb) +
     &           fact2*( shgradgu(aa,1)*shgradgu(bb,1)*aa3
     &                 )*DetJgwt

            xKebe(19,aa,bb) = xKebe(19,aa,bb) +
     &           fact2*( shgradgu(aa,1)*shgradgu(bb,1)*etaT
     &                 )*DetJgwt

            xKebe(25,aa,bb) = xKebe(25,aa,bb) +
     &           fact2*( shgradgu(aa,1)*shgradgu(bb,1)*kappaT
     &                 )*DetJgwt

         enddo ! aa
      enddo ! bb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Tangent matrix xKebe - Stress flux part - Y direction A4_y matrix

      do bb = 1, NSHLu
         do aa = 1, NSHLu

            xKebe(11,aa,bb) = xKebe(11,aa,bb) +
     &           fact2*( shgradgu(aa,2)*shgradgu(bb,2)*aa3
     &                 )*DetJgwt

            xKebe(12,aa,bb) = xKebe(12,aa,bb) +
     &           fact2*( shgradgu(aa,2)*shgradgu(bb,1)*aa3
     &                 )*DetJgwt

            xKebe(13,aa,bb) = xKebe(13,aa,bb) +
     &           fact2*( shgradgu(aa,2)*shgradgu(bb,2)*etaT
     &                 )*DetJgwt

            xKebe(16,aa,bb) = xKebe(16,aa,bb) + 
     &           fact2*( shgradgu(aa,2)*shgradgu(bb,1)*
     &                   (aa1-aa2*tau+3d0*aa4*se2T**2-5d0*aa6*se2T**4)
     &                 )*DetJgwt
 
            xKebe(17,aa,bb) = xKebe(17,aa,bb) +
     &           fact2*( shgradgu(aa,2)*shgradgu(bb,2)*
     &                   (aa1+aa2*tau-3d0*aa4*se2T**2+5d0*aa6*se2T**4)
     &                 )*DetJgwt

            xKebe(19,aa,bb) = xKebe(19,aa,bb) +
     &           fact2*( shgradgu(aa,2)*shgradgu(bb,2)*etaT
     &                 )*DetJgwt

            xKebe(20,aa,bb) = xKebe(20,aa,bb) +
     &           fact2*( -shgradgu(aa,2)*se2T*shgu(bb)*aa2/Tm
     &                 )*DetJgwt

            xKebe(25,aa,bb) = xKebe(25,aa,bb) +
     &           fact2*( shgradgu(aa,2)*shgradgu(bb,2)*kappaT
     &                 )*DetJgwt

         enddo ! aa
      enddo ! bb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Tangent matrix xKebe - Phase field part

      do bb = 1, NSHLu
         do aa = 1, NSHLu

            xKebe(11,aa,bb) = xKebe(11,aa,bb) +
     &           fact2*( kgTS*shhessgu(aa,1,1)
     &                   *( shhessgu(bb,1,1)+shhessgu(bb,2,2) )
     &                 )*DetJgwt

            xKebe(12,aa,bb) = xKebe(12,aa,bb) +
     &           fact2*( -kgTS*shhessgu(aa,1,2)
     &                   *( shhessgu(bb,1,1)+shhessgu(bb,2,2) )
     &                 )*DetJgwt

            xKebe(16,aa,bb) = xKebe(16,aa,bb) +
     &           fact2*( -kgTS*shhessgu(aa,2,1)
     &                   *( shhessgu(bb,1,1)+shhessgu(bb,2,2) )
     &                 )*DetJgwt

            xKebe(17,aa,bb) = xKebe(17,aa,bb) +
     &           fact2*( kgTS*shhessgu(aa,2,2)
     &                   *( shhessgu(bb,1,1)+shhessgu(bb,2,2) )
     &                 )*DetJgwt

         enddo ! aa
      enddo ! bb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
      return
      end
      
c  phase field part xkebe

!            xKebe(11,aa,bb) = xKebe(11,aa,bb) +
!     &           fact2*( kgTS*shhessgu(aa,1,1)
!     &                   *( du1dxidxj(1,1)+du1dxidxj(2,2) )
!     &                 )*DetJgwt
!
!            xKebe(12,aa,bb) = xKebe(12,aa,bb) +
!     &           fact2*( -kgTS*shhessgu(aa,1,2)
!     &                   *( du2dxidxj(1,1)+du2dxidxj(2,2) )
!     &                 )*DetJgwt
!
!            xKebe(16,aa,bb) = xKebe(16,aa,bb) +
!     &           fact2*( -kgTS*shhessgu(aa,1,2)
!     &                   *( du1dxidxj(1,1)+du1dxidxj(2,2) )
!     &                 )*DetJgwt
!
!            xKebe(17,aa,bb) = xKebe(17,aa,bb) +
!     &           fact2*( kgTS*shhessgu(aa,2,2)
!     &                   *( du2dxidxj(1,1)+du2dxidxj(2,2) )
!     &                 )*DetJgwt
