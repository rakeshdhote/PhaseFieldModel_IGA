      subroutine e3Rhs_3D (
     &     xi, ui, aci, rLi, chii,
     &     duidxi, duidxixj, gwt,
     &     shgu, shgradgu, shhessgu, 
     &     grade2energy,lapu1, lapu2,
     &     Rhsu)
      
      use aAdjKeep
      use common

      implicit none

      integer aa, idof

      real*8  shgu(NSHLu), 
     &        shgradgu(NSHLu,NSD), 
     &        shhessgu(NSHLu,NSD,NSD)

      real*8 fact1, fact2

      real*8
     &     ui(NDOF), 
     &     aci(NDOF), 
     &     rLi(NDOF), 
     &     chii(NSD),
     &     xi(NSD),
     &     timeAlpha

      real*8 duidxi(NDOF,NSD), duidxixj(NDOF,NSD,NSD)
      real*8 grade2energy

      real*8
     &     gwt,
     &     Rhsu(NDOF,NSHLu), 
     &     DetJgwt

      real*8 A1(NDOF,NDOF)
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
!    Constants

c      zero = 0d0   ! zero
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DetJgwt = DetJ*gwt

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
      du1dx = duidxi(1,1)! +const*timeAlpha
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

!     Laplacians in u1 and u2 directions
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
      S11 =  aa1*se1T+aa2*tau*se2T-aa4*se2T**3+aa6*se2T**5
      S12 =  aa3*se3T
      S21 =  aa3*se3T
      S22 =  aa1*se1T-aa2*tau*se2T+aa4*se2T**3-aa6*se2T**5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Time dependent matrix (initialization)
      A1   = 0d0
!     Populating the A1 matrix cells
      A1(1,1) = 1d0
      A1(2,2) = 1d0
      A1(3,3) = rhoT
      A1(4,4) = rhoT
      A1(5,5) = cvT

!     New variables and forcing functions (f2)
!     + coupling term (f3) added together (initialization)
      f2f3   = 0d0
!     Populating the f2f3 vector
      f2f3(1) = -1d0*v1
      f2f3(2) = -1d0*v2
      f2f3(3) = -1d0*fx
      f2f3(4) = -1d0*fy
      f2f3(5) = -1d0*gth-1d0/2d0*a2th*th*se2T*se2Tdot

!     Stress flux - x direction (initialization)
      f4x   = 0d0
!     Populating the f4 vector
!      f4x(1) = 0d0
!      f4x(2) = 0d0
      f4x(3) = S11+etaT*dv1dx
      f4x(4) = S21+etaT*dv2dx
      f4x(5) = kappaT*dthdx

!     Stress flux - y direction (initialization)
      f4y   = 0d0
!     Populating the f4y vector
!      f4y(1) = 0d0
!      f4y(2) = 0d0
      f4y(3) = S12+etaT*dv1dy
      f4y(4) = S22+etaT*dv2dy
      f4y(5) = kappaT*dthdy

!     Phase field part - 2nd derivative d^2/dx^2 terms (initialization)
      f5xx   = 0d0
!     Populating the f5xx vector
!      f5xx(1) = 0d0
!      f5xx(2) = 0d0
      f5xx(3) = kgTS*lapu1
!      f5xx(4) = 0d0
!      f5xx(5) = 0d0

!     Phase field part - 2nd derivative d^2/dxdy terms (initialization)
      f5xy   = 0d0
!     Populating the f5xy vector
!      f5xy(1) = 0d0
!      f5xy(2) = 0d0
      f5xy(3) = -1d0*kgTS*lapu2
      f5xy(4) = -1d0*kgTS*lapu1
!      f5xy(5) = 0d0

!     Phase field part - 2nd derivative d^2/dy^2 terms (initialization)
      f5yy   = 0d0
!     Populating the f5yy vector
!      f5yy(1) = 0d0
!      f5yy(2) = 0d0
!      f5yy(3) = 0d0
      f5yy(4) = kgTS*lapu2
!      f5yy(5) = 0d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Populate Residual Vector

      do aa = 1, NSHLu
         do idof = 1, NDOF
            
            Rhsu(idof,aa) = Rhsu(idof,aa)
     &         -(
     &              shgu(aa)*( A1(idof,1)*aci(1) + A1(idof,2)*aci(2) +
     &                         A1(idof,3)*aci(3) + A1(idof,4)*aci(4) +
     &                         A1(idof,5)*aci(5)
     &                       )
     &            +shgu(aa)*f2f3(idof)
     &            +shgradgu(aa,1)*f4x(idof)
     &            +shgradgu(aa,2)*f4y(idof)
     &            +shhessgu(aa,1,1)*f5xx(idof)
     &            +shhessgu(aa,1,2)*f5xy(idof)
     &            +shhessgu(aa,2,2)*f5yy(idof)
     &          )*DetJgwt

         enddo  ! idof
      enddo  ! aa
      
      return
      end

c_____________________________________________________________________________

      subroutine Energy(ui, duidxi, grade2energy, gwt)
      !This function is evaluated at each Gauss point

      use aAdjKeep
      use common

      implicit none     
      
      real*8 DetJgwt, ui(NDOF), gwt, DetJgwtc, tau
      real*8 du1dx, du2dx, dv1dx, dv2dx, dthdx, e1T, e2T, e3T
      real*8 du1dy, du2dy, dv1dy, dv2dy, dthdy
      real*8 grade2energy, duidxi(NDOF,NSD)
      real*8 u1, u2, v1, v2, th
      real*8 Sxx, Sxy, Syy, se1T, se2T, se3T
!      real*8 aa1T, aa2T, aa3T

!     Variables
      u1 = ui(1)
      u2 = ui(2)
      v1 = ui(3)
      v2 = ui(4)
      th = ui(5)

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

!     Strain definations
      e1T = (du1dx+du2dy)/sqrt(2d0)   ! Strains
      e2T = (du1dx-du2dy)/sqrt(2d0)
      e3T = (du1dy+du2dx)/2d0
      se1T = (du1dx+du2dy)            ! rescaled strains
      se2T = (du1dx-du2dy)
      se3T = (du1dy+du2dx)

!     Temperature Variable

      tau = (th-Tm)/Tm

!     Stress definations
      Sxx =  aa1*se1T+aa2*tau*se2T-aa4*se2T**3+aa6*se2T**5
      Sxy =  aa3*se3T
      Syy =  aa1*se1T-aa2*tau*se2T+aa4*se2T**3-aa6*se2T**5
      DetJgwt = DetJ*gwt
      DetJgwtc = DetJgwt*2d0/0.555555555555d0 ! Factor to account for
      ! linear element across z-direction

c      rho = ui(1)

c      pot = (8d0/27d0)*theta*rho*log(rho/(1d0-rho))-rho*rho

      Ener = Ener +
     &       (
     &          1d0/2d0*aa1T*e1T**2 + 1d0/2d0*aa3T*e3T**2
     &          + 1d0/2d0*aa2T*tau*e2T**2 - 1d0/4d0*e2T**4
     &          + 1d0/6d0*e2T**6 + 1d0/2d0*grade2energy
     &       )*DetJgwtc

      AvgTemp = AvgTemp + th*DetJgwtc ! define AvgTemp =0 (near Ener =0)

      AvgSxx = AvgSxx + Sxx*DetJgwtc
      AvgSxy = AvgSxy + Sxy*DetJgwtc
      AvgSyy = AvgSyy + Syy*DetJgwtc

      AvgExx = AvgExx + du1dx*DetJgwtc
      AvgExy = AvgExy + du1dy*DetJgwtc
      AvgEyy = AvgEyy + du2dy*DetJgwtc

      return
      end

c_____________________________________________________________________
!
!      subroutine Elliptic(ui, laprho, gwt)
!
!      use aAdjKeep
!      use common
!
!      implicit none
!
!      real*8 DetJgwt, ui(NDOF), gwt
!
!      real*8 laprho, h, rho
!
!      DetJgwt = DetJ*gwt
!
!      rho = ui(1)
!
!      h = (8d0*theta/3d0)*log(rho/(3d0-rho))
!     &  -  8d0*theta/(rho-3d0)-6d0*rho
!     &  - lambda*laprho
!
!      if (h.gt.hmax) hmax = h
!      if (h.lt.hmin) hmin = h
!
!
!      havg = havg + h*DetJgwt
!
!
!      return
!      end
!
!
