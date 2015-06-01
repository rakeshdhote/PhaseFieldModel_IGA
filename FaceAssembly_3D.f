c     This subroutine calculates the contributions to the left hand
c      side forcing vector due to boundary faces.
c
c     This is for the 3D case.
c     
c     
c     July 1, 2003
c     J. Austin Cottrell
c     CES Graduate Student
c     Texas Institute for Computational Engineering Science
c     University of Texas at Austin
c
c     Modified from 2D "BdryEdjAss.f"  
c     January 20, 2004
c     J. Austin Cottrell


      subroutine FaceAssembly_3D(ugAlpha, acgAlpha, pgAlpha,
     &  colu, rowu, colp, rowp)
      
      use aAdjKeep
      
      use common
      
      implicit none

      integer ifac, nshb, igaussb, jgaussb, i, g, ni, nj, 
     &  nk, niu, nju, nku, nip, njp, nkp, j,
     &  area_flag1, area_flag2, area_flag_m1, area_flag_m2,
     &  area_flag_m3, rowu(NNODZu*8*(Pu+1)*(Qu+1)*(Ru+1)),
     &  colu(NNODZu+1), rowp(NNODZp*8*(Pp+1)*(Qp+1)*(Rp+1)),
     &  colp(NNODZp+1)
      real*8  gp(NGAUSS), gw(NGAUSS), da, gwt
      real*8  Rhsu(NSD,NSHLu), Rhsp(NSHLp),
     &  xKebe(NSD*NSD,NSHLu,NSHLu), xGebe(NSD,NSHLu,NSHLp),
     &  xDebet(NSD,NSHLu,NSHLp),
     &  shbu(NSHLu), shbug(NSHLu,NSD), nor(NSD),
     &  Pres, norsurf,
     &  ugAlpha(NNODZu,NSD), acgAlpha(NNODZu,NSD), pgAlpha(NNODZp),
     &  ui(NSD), pri, fi(NSD), gi(NSD),
     &  duidxi(NSD,NSD), tauB, shbp(NSHLp)
      real*8 fl(NSHLu,NSD), ul(NSHLu,NSD), gl(NSHLu,NSD),
     &  acl(NSHLu,NSD), pl(NSHLp), dxidx(NSD,NSD)
      
c     shb will be the shape function array while shbg will hold the
c     gradients of the shape functions

      
      
c...  Loop over Faces

      norsurf = 0d+0
      
      do ifac = 1, NFACEu
        
        area_flag1 = 0          ! initialize flag
        area_flag2 = 0
        
c...    Check if face is of nonzero area
        
        niu = INNu(FACE_IENu(ifac,1),1)
        nju = INNu(FACE_IENu(ifac,1),2)
        nku = INNu(FACE_IENu(ifac,1),3)
        
        nip = INNp(FACE_IENp(FACE_CON(ifac),1),1)
        njp = INNp(FACE_IENp(FACE_CON(ifac),1),2)
        nkp = INNp(FACE_IENp(FACE_CON(ifac),1),3)
        
        if ((FACE_OR(ifac).eq.1).or.(FACE_OR(ifac).eq.6)) then
          da = (U_KNOTu(niu+1) - U_KNOTu(niu))* 
     &      (V_KNOTu(nju+1) - V_KNOTu(nju))/4d+0 ! change in gwt
                                ! due to mapping from
                                ! [-1,1] x [-1,1]
          if((U_KNOTu(niu).eq.(U_KNOTu(niu+1))).or.
     &      (V_KNOTu(nju).eq.(V_KNOTu(nju+1))) ) then
            area_flag1 = 1       ! indicates face has zero area
          endif
          
          
        else if ((FACE_OR(ifac).eq.2).or.(FACE_OR(ifac).eq.4)) then
          da = (U_KNOTu(niu+1) - U_KNOTu(niu))* 
     &      (W_KNOTu(nku+1) - W_KNOTu(nku))/4d+0 ! change in gwt
                                ! due to mapping from
                                ! [-1,1] x [-1,1]
          if((U_KNOTu(niu).eq.(U_KNOTu(niu+1))).or.
     &      (W_KNOTu(nku).eq.(W_KNOTu(nku+1))) ) then
            area_flag1 = 1       ! indicates face has zero area
          endif
          
        else if ((FACE_OR(ifac).eq.3).or.(FACE_OR(ifac).eq.5)) then
          da = (V_KNOTu(nju+1) - V_KNOTu(nju))* 
     &      (W_KNOTu(nku+1) - W_KNOTu(nku))/4d+0 ! change in gwt
                                ! due to mapping from
                                ! [-1,1] x [-1,1]
          if((V_KNOTu(nju).eq.(V_KNOTu(nju+1))).or.
     &      (W_KNOTu(nku).eq.(W_KNOTu(nku+1))) ) then
            area_flag1 = 1       ! indicates face has zero area
          endif
          
        endif

                                ! Perform assembly if Neumann or Weak Dir.

        area_flag_m1 = 0
        if((IBC_FACE(ifac,1).ne.2).and.(IBC_FACE(ifac,1).ne.3)) then
          area_flag_m1 = 1
        endif
        
        area_flag_m2 = 0
        if((IBC_FACE(ifac,2).ne.2).and.(IBC_FACE(ifac,2).ne.3)) then
          area_flag_m2 = 1
        endif

        area_flag_m3 = 0
        if((IBC_FACE(ifac,3).ne.2).and.(IBC_FACE(ifac,3).ne.3)) then
          area_flag_m3 = 1
        endif

        area_flag2 = area_flag_m1*area_flag_m2*area_flag_m3
        
c...    If face has zero area, skip entirely
        if ((area_flag1.eq.0).and.(area_flag2.eq.0)) then 
          
          
c...      get Gauss Point/Weight Arrays
          
          call genGPandGW(gp,gw,NGAUSS)
          
          shbu = 0d+0
          shbug = 0d+0
          shbp = 0d+0
          
                                ! Get local solution arrays
          
          do i = 1, NSHLu
            ul(i,:) = ugAlpha(FACE_IENu(ifac,i),:)
            gl(i,:) = DIR_BC(FACE_IENu(ifac,i),:)
          enddo
          
          do i = 1, NSHLp
            pl(i) = pgAlpha(FACE_IENp(FACE_CON(ifac),i))
          enddo
          
          mu = VISC(1)
          rho = DENS(1)
          
          do i = 1, NSHLu
            do j = 1, NSD
              fl(i,j) = LD_FACE(ifac,j)
            enddo
          enddo
          
          xKebe = 0d+0          ! initialize local stiffness matrix
          xGebe = 0d+0
          xDebet = 0d+0
          Rhsu = 0d+0           ! initialize local load vector
          Rhsp = 0d+0
          
!         Loop over integration points
          
          do igaussb = 1, NGAUSS
            do jgaussb = 1, NGAUSS
              
!             Get Boundary Shape functions and set Jacobian
              
              call eval_FACE_3D(ifac, gp(igaussb), gp(jgaussb),
     &          niu, nju, nku, nip, njp, nkp,
     &          shbu, shbug, shbp, dxidx, nor)
              
!             Get velocity, its gradient and f at int. point
              
!              call e3bint(ul, pl, fl, gl, shbu, shbug, shbp,
!     &          ui, pri, fi, gi, duidxi)
              
!              gwt = gw(igaussb)*gw(jgaussb)*da
              
!             Get boundary Stab. Parameter C*mu/h_n
              
!              call e3bSTAB_3D(mu, nor, dxidx, tauB)                            
              
c*              call e3bLHS_3D(ui, tauB, gwt, 
c*     &          shbu, shbug, shbp, xGebe, xDebet, xKebe,
c*     &          nor)
              
c*              if (IBC_FACE(ifac,1).eq.2) then
c*                xKebe(1:3,:,:) = 0d+0
c*                xKebe(4,:,:) = 0d+0
c*                xKebe(7,:,:) = 0d+0
c*                xGebe(1,:,:) = 0d+0
c*                xDebet(1,:,:) = 0d+0
c*              elseif (IBC_FACE(ifac,2).eq.2) then
c*                xKebe(4:6,:,:) = 0d+0
c*                xKebe(2,:,:) = 0d+0
c*                xKebe(8,:,:) = 0d+0
c*                xGebe(2,:,:) = 0d+0
c*                xDebet(2,:,:) = 0d+0
c*              elseif (IBC_FACE(ifac,3).eq.2) then   
c*                xKebe(7:9,:,:) = 0d+0
c*                xKebe(3,:,:) = 0d+0
c*                xKebe(6,:,:) = 0d+0
c*                xGebe(3,:,:) = 0d+0
c*                xDebet(3,:,:) = 0d+0
c*              endif
              
!              call e3bRHS_3D(ifac, nor, tauB, gwt, 
!     &          shbu, shbug, shbp, ui, pri, fi, gi, duidxi,
!     &          Rhsu, Rhsp)
              
              
                                    
            enddo
          enddo
            

!          call BCBLhs_3D(ifac, xKebe, xGebe, xDebet, Rhsu)
          
c...      assemble into the Sparse Global Stiffness Matrix and Rhs Vector
          
!          call FillSparseMatB_3D_u (ifac, colu, rowu, xKebe) ! Assemble stiffness
!          call LocaltoGlobalB_3D_u (ifac, Rhsu) ! Assemble load vector
          
!          call FillSparseMatB_3D_p (ifac, colp, rowp, xGebe) ! Assemble stiffness
!          call LocaltoGlobalB_3D_p (ifac, Rhsp) ! Assemble load vector
          
!          call FillSparseMatB_3D_p2 (ifac, colp, rowp, xDebet)
          
          
          
        endif
        
      enddo
      
ccc      write(*,*) "Surface Area=    ", norsurf
      
      return
      end
      



      
