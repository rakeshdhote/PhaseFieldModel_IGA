      subroutine VizGMV_3D(iproc,numprocs,istep,inistep)
      
      use aAdjKeep
      include "common.h"
      
c...  Local variables
      integer e, igauss, jgauss, i, icount, j, k, kgauss, aa,
     &  NEL_NZA_QUAD, NEL_LIN, NEL_LINm1, icount2, iln, iel, nnodz_lin
      integer ni,nj,nk,strss_file
      integer iproc,numprocs,istep,inistep
      real*8 shl(NSHL), shg(NSHL), shgradl(NSHL,NSD),
     &  shgradg(NSHL, NSD), QuantQuad(8,9), ytemp(NNODZ,NSD)
      real*8 gp(NGAUSS), gw(NGAUSS), skl(2,2,NSD), stl(NSHL,NSD)
      real*8 lambda, numod, Emod, mu, yl(NSHL,NSD), dist, distmin,
     &  summ, std, stdmax
      real*8 dcero
      integer*4 ipoin, idof, ipart 

      character*30 fname
      character*10 cname
      
      integer ifil, lstep
      
      integer, allocatable :: IEN_LIN(:,:)
      real*8, allocatable :: x_lin(:,:), y_lin(:,:), QuantLin(:,:,:)
      
c...  get Gaussian points and weights
      
      NGAUSS = 2
      gp(1) = -1d0 ! No subdivision, if subdivided gp(2)=0 and gp(3)=1
      gp(2) =  1d0
      gw(:) =  0d0
      
      NEL_NZA_QUAD = 0
      NEL_LIN      = 0
      NEL_LINm1    = 0

      NEL_LIN = NEL
      
      stdmax = 0d0

      allocate(QuantLin(NEL,8,9))

      do iel = 1, NEL
        
c...    Check to see if current element has nonzero area
        ni = INN(IEN(iel,1),1)  ! get NURB coordinates
        nj = INN(IEN(iel,1),2)
        nk = INN(IEN(iel,1),3)
        
        if ((U_KNOT(ni).ne.U_KNOT(ni+1)).and.
     &    (V_KNOT(nj).ne.V_KNOT(nj+1)).and.
     &    (W_KNOT(nk).ne.W_KNOT(nk+1))) then
          
          icount2 = 0
          
c...      Loop over integration points (NGAUSS in each direction)
          
          do igauss = 1, NGAUSS
            do jgauss = 1, NGAUSS
              do kgauss = 1, NGAUSS
                
                icount2 = icount2+1
                
c...            Get Element Shape functions and their gradients
                
                shl = 0d+0      ! initialize
                shg = 0d+0
                shgradl = 0d+0
                shgradg = 0d+0
                
                call eval_SHAPE_3D(iel, gp(igauss), gp(jgauss),
     &            gp(kgauss), shl, shgradl,shgradg)

                shg = shl                           
                
                do i = 1,NSHL
!                  yl(i, :) = yg(IEN(iel,i),1:NSD)
                  yl(i, 1) = yg(IEN(iel,i),1)  ! u1 
				  yl(i, 2) = yg(IEN(iel,i),2)  ! u2 
				  yl(i, 3) = yg(IEN(iel,i),5)  ! theta 				  
                enddo
                
                call eval_QUANT_3D(iel, yl, shg, shgradg, shgradl, 
     &            icount2, QuantQuad)
                
              enddo
            enddo
          enddo
          
      ! Collect Local Arrays for "linear" elements

          QuantLin(iel,1,:) = QuantQuad(1,:) ! Connectivity
          QuantLin(iel,2,:) = QuantQuad(2,:)
          QuantLin(iel,3,:) = QuantQuad(4,:)
          QuantLin(iel,4,:) = QuantQuad(3,:)
          QuantLin(iel,5,:) = QuantQuad(5,:)
          QuantLin(iel,6,:) = QuantQuad(6,:)
          QuantLin(iel,7,:) = QuantQuad(8,:)
          QuantLin(iel,8,:) = QuantQuad(7,:)
          
        endif
        
      enddo


      ! Bulid Linear IEN Based on Spatial Location
      
      allocate(IEN_LIN(NEL_LIN,8))
      allocate(x_lin(8*NEL_LIN,NSD))
      allocate(y_lin(8*NEL_LIN,6))

      
      icount = 1   ! Initialize Global Node Count
      x_lin(1,:) = QuantLin(1,1,7:9) ! Insert First Global Node ! x, y , z coordinates
      y_lin(1,:) = QuantLin(1,1,1:6) ! solution
      ! written for 6 stress of structral components

      IEN_LIN(1,1) = 1
      
      distmin = 1d-13
      
      do iel = 1, NEL_LIN
        
        do iln = 1, 8                 
          
          do k = 1, icount      ! Check for repeated nodes
            
      ! Distance to existing Global Nodes
            
            dist = (x_lin(k,1) - QuantLin(iel,iln,7))**2 +
     &             (x_lin(k,2) - QuantLin(iel,iln,8))**2 
            
            dist = sqrt(dist)
            
            if (dist.le.distmin) then
              IEN_LIN(iel,iln) = k ! Repeted Global Node
              goto 1000
            endif
            
          enddo

          icount = icount + 1   ! Add New Global Node
          IEN_LIN(iel,iln) = icount
          x_lin(icount,:) = QuantLin(iel,iln,7:9) ! x, y, z coordinates
          y_lin(icount,:) = QuantLin(iel,iln,1:6) ! solution vector

1000     continue

        enddo
      enddo

      nnodz_lin = icount ! total # nodes


      if (iproc.eq.1) then
        nnodz_glob_lin = 0
        NEL_GLOB_LIN  = 0
        allocate(x_glob_lin(nnodz_lin*2*numprocs,NSD))
        allocate(y_glob_lin(nnodz_lin*2*numprocs,6)) ! Structural Cauchy  ! sterss ! for subdivided mesh (finer mesh)
        allocate(IEN_GLOB_LIN(NEL_LIN*2*numprocs,8)) ! hex elem 8 nodes
        allocate(glob_mat_dat(NEL_LIN*2*numprocs))
      endif

      ! nnodz_lin - number of global linear nodes for this patch
      ! NEL_LIN - number of linear elements in this patch

      do i = 1, nnodz_lin
        x_glob_lin(nnodz_glob_lin+i,1:NSD) = x_lin(i,1:NSD)
        y_glob_lin(nnodz_glob_lin+i,1:6  ) = y_lin(i,1:6  )
      enddo

      do i = 1, NEL_LIN
        glob_mat_dat(NEL_GLOB_LIN+i) = real(iproc) +
     &    1d+2*sin(1.57*real(iproc))
        do j = 1, 8
          IEN_GLOB_LIN(NEL_GLOB_LIN+i,j) = IEN_LIN(i,j) + nnodz_glob_lin
        enddo
      enddo
      
      nnodz_glob_lin = nnodz_glob_lin + nnodz_lin
      NEL_GLOB_LIN = NEL_GLOB_LIN + NEL_LIN
      
      nnodz_lin = 0
      NEL_LIN = 0

!      if (iproc.eq.1) then

!        write(*,*) "Writing solution file"
        
        nnodz_lin = nnodz_glob_lin
        NEL_LIN = NEL_GLOB_LIN
        
        dcero = 0d0
        ifil = 99

      ! Writing ensight files
        if (istep.eq.inistep) then

        fname = trim ('SMA2D' // cname(istep)) // '.case' ! Writing ensight files
!        fname = 'SMA2D.case' ! Writing ensight files

           OPEN  (ifil,FILE=fname,STATUS='UNKNOWN',FORM='FORMATTED')

           write(ifil,200)
 200       format('#BOF: SMA2D.case'//
     &            'FORMAT'//
     &            'type: ensight'//
     &            'GEOMETRY'//
     &            'model: geometry.0.geo'//
     &            'VARIABLE'//
     &       'scalar per node: 1 Temperature temperature.*.res'//
     &       'vector per node: 2 Displacement displacment.*.res'//
     &       'vector per node: 3 Strain strains.*.res'//	 
     &            'TIME'//
     &            'time set:                    1'//
     &            'number of steps:             10'//
     &            'filename start number:       0'//
     &            'filename increment:         100'//
     &            'time values:')
           do ipoin = 1, 10000
              write(ifil,'(i8)') ipoin
           enddo  ! #nodes, x, y, z coordinates
           close(ifil)

!           write(*,*) "Writing geometry.*.geo file"
           fname = trim('geometry' // cname(istep)) // '.geo'
!           fname = 'geometry.geo'

           OPEN  (ifil,FILE=fname,STATUS='UNKNOWN',FORM='FORMATTED')
           
           write(ifil,'(a)')'SMA'
           write(ifil,'(a)')'node'
           write(ifil,'(a)')'node id given'
           write(ifil,'(a)')'element id given'
           write(ifil,'(a)')'coordinates'
           write(ifil,'(i8)') nnodz_lin
           
!           write(*,*) "Writing coordinates"

           do ipoin = 1, nnodz_lin
              write(ifil,'(i8,3e12.5)') ipoin, 
     &             x_glob_lin(ipoin,1), x_glob_lin(ipoin,2), dcero ! dcero -zcoordiante
           enddo  ! #nodes, x, y, z coordinates
           ipart = 1
           write(ifil,'(a,i8)') 'part', ipart  ! processor(or patch)
           write(ifil,'(a)') 'todo'
           write(ifil,'(a)') 'quad4'
           write(ifil,'(i8)') NEL_LIN
           
!           write(*,*) "Writing Connectivity"
           
           do i = 1, NEL_LIN
             write(ifil,'(5i8)') i, IEN_GLOB_LIN(i,1),IEN_GLOB_LIN(i,3),
     &                              IEN_GLOB_LIN(i,7),IEN_GLOB_LIN(i,5)
              
           enddo
           close(ifil)

        endif

!        write(*,*) "Writing solution - temperature.*.res"

        fname = trim ('temperature' // cname(istep)) // '.res'
        OPEN (ifil,FILE=fname,STATUS='UNKNOWN',FORM='FORMATTED')
        write(ifil,'(a)') 'temperature'
        write(ifil,'(6e12.5)') y_glob_lin(1:nnodz_lin,6) ! temperature
        
        close(ifil)

!        write(*,*) "Writing solution - displacement.*.res"

        fname = trim ('displacment' // cname(istep)) // '.res'
        OPEN (ifil,FILE=fname,STATUS='UNKNOWN',FORM='FORMATTED')
        write(ifil,'(a)') 'displacment'
        write(ifil,'(6e12.5)') ((y_glob_lin(ipoin,idof),idof=4,5),dcero,
     &                          ipoin=1,nnodz_lin) !u1, u2 disp
        close(ifil)

!        write(*,*) "Writing solution - strains.*.res"

        fname = trim ('strains' // cname(istep)) // '.res'
        OPEN (ifil,FILE=fname,STATUS='UNKNOWN',FORM='FORMATTED')
        write(ifil,'(a)') 'strains'
        write(ifil,'(6e12.5)') ((y_glob_lin(ipoin,idof),idof=1,2),dcero,
     &                          ipoin=1,nnodz_lin) !u1, u2 disp
        close(ifil)		
		
!      endif !filters out 4th, 5th column, add z-deformation dcero to it
      ! and extracts for all nodes,
      ! 6e12.5 - 6=# of entries in one line, 12=width, 5=# of decimal
      
      deallocate(IEN_LIN)
      deallocate(x_lin)
      deallocate(y_lin)
      deallocate(QuantLin)

      write(*,*)"Finished"
      
      
      return
      end
      

!#########################################################################


      subroutine eval_QUANT_3D(e, yl, shg, shgradg, shgradl,
     &  icount, Quant)
      
      use aAdjKeep
      
      include "common.h"
      
      real*8 yl(NSHL,NSD), shg(NSHL), Quant(NGAUSS**3,9),
     &  xl(NSD), disp(NSD), shgradg(NSHL,NSD), strain(6),
     &  du1dxi(3), du2dxi(3), du3dxi(3), dthdxi(3), shgradl(NSHL,NSD)	 
      
      integer aa, bb, icount, e
      
      xl = 0d+0
      disp = 0d+0
      du1dxi = 0d+0
      du2dxi = 0d+0
      du3dxi = 0d+0	  
      dthdxi = 0d+0		  
	  
c...  Calculate Symmetric Strain Tensor - 3D
      
      do aa = 1, NSHL
        
c...    Displacements and Stresses
        
        disp(1) = disp(1) + yl(aa,1)*shg(aa)
        disp(2) = disp(2) + yl(aa,2)*shg(aa)
        disp(3) = disp(3) + yl(aa,3)*shg(aa)
		
c...    Strains

! IN the following definations shgradg is changed to shgradl        
        du1dxi(1) = du1dxi(1) + yl(aa,1)*shgradl(aa,1)
        du1dxi(2) = du1dxi(2) + yl(aa,1)*shgradl(aa,2)
!        du1dxi(3) = du1dxi(3) + yl(aa,1)*shgradl(aa,3)	

        du2dxi(1) = du2dxi(1) + yl(aa,2)*shgradl(aa,1)
        du2dxi(2) = du2dxi(2) + yl(aa,2)*shgradl(aa,2)
!        du2dxi(3) = du2dxi(3) + yl(aa,2)*shgradl(aa,3)	

!        du3dxi(1) = du3dxi(1) + yl(aa,3)*shgradl(aa,1)
!        du3dxi(2) = du3dxi(2) + yl(aa,3)*shgradl(aa,2)
!        du3dxi(3) = du3dxi(3) + yl(aa,3)*shgradl(aa,3)		

!    here the y1(aa,3) represents temperature
        dthdxi(1) = dthdxi(1) + yl(aa,3)*shgradl(aa,1)
        dthdxi(2) = dthdxi(2) + yl(aa,3)*shgradl(aa,2)
!        dthdxi(3) = dthdxi(3) + yl(aa,3)*shgradl(aa,3)		

	
c...    Locations in reference configuration
        
        xl(1) = xl(1) + 
     &    shg(aa)*B_NET(INN(IEN(e,aa),1), INN(IEN(e,aa),2),
     &    INN(IEN(e,aa),3),1)
        xl(2) = xl(2) +
     &    shg(aa)*B_NET(INN(IEN(e,aa),1), INN(IEN(e,aa),2),
     &    INN(IEN(e,aa),3),2)
        xl(3) = xl(3) +
     &    shg(aa)*B_NET(INN(IEN(e,aa),1), INN(IEN(e,aa),2), 
     &    INN(IEN(e,aa),3),3)
        
      enddo
 
        strain(1) = du1dxi(1)  !  1d0/sqrt(2d0)*(du1dxi(1)+du2dxi(2))! e1=hydrostratic
        strain(2) = du2dxi(2)  ! 1d0/sqrt(2d0)*(du1dxi(1)-du2dxi(2))! e2=deviatoric
        strain(3) = du1dxi(2)  ! 1d0/2d0*(du1dxi(2)+du2dxi(1))  ! e3=shear 
      
c...  Format of the Elmstress is:
c...  for each integration point (indexed by icount) we have 
c...  {S11,S22,S33,S12,S23,S31,dx,dy,dz,x,y,z}
      
!      Quant(icount,1) = 0d0 ! Can use this
!      Quant(icount,2) = 0d0
!      Quant(icount,3) = 0d0

      Quant(icount,1) = strain(1)
      Quant(icount,2) = strain(2)
      Quant(icount,3) = strain(3)
      
      Quant(icount,4) = disp(1) ! u1
      Quant(icount,5) = disp(2) ! u2
      Quant(icount,6) = disp(3) ! temp
      
      Quant(icount,7) = xl(1)   ! x- coord
      Quant(icount,8) = xl(2)   ! y- coord
      Quant(icount,9) = xl(3)   ! z- coord
      
      return
      end


!#####################################################################

      function cname (i)
      
      logical      beg
      character*10 cname,cc
      integer      il(0:8)

      integer i, ic0, ii, k
      
      ic0 = ICHAR("0")
      cc = " "
      ii = i
      
      il(0) = mod(ii,10)
      do k = 1,8
        ii = (ii - il(k-1)) / 10
        il(k) = mod (ii,10)
      enddo
      
      beg = .false.
      
      do k = 8,1,-1
        if (il(k) .ne. 0 .or. beg) then
          beg = .true.
          cc  = TRIM(cc) // CHAR(ic0 + il(k))
        endif
      enddo
      
      cc = TRIM(cc)//CHAR(ic0 + il(0))
      cname = "." // cc
      
      return
      end
