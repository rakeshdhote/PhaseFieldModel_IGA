c------------------------------------------------------------------------c
c     Main routine to call all the subroutines
c------------------------------------------------------------------------c
      ! RPD Change
      program SMA2DSerial

      use maindef

      implicit none
      real*8 gami_f, tolt, avgdt, Dtold

      real*8 rtdum(2), rtdumg(2), rdum ! Norm for DUMMY u1 and u2
      real*8 rtmom(2), rtmomg(2), rmom ! Norm for v1 and v2
      real*8 rtTemp(2), rtTempg(2), rTemp ! Norm for Temperature
      real*8 rtot

      integer*4 ireject, Stmin, jstep

      character*10 cname

      time = 0d0
      avgdt= 0d0
      ireject = 0  ! # of iterations for error to go below tolerance
      myid = 0
      numnodes = 1

!=======================================================================

      call main_init            ! Initialize

!=======================================================================


      yg = ygold
      call mean
      call main_prn
      yg = 0d0

      print*, 'main program'
      print*, 'Stepcount',Stepcount

!      Delt  = 1d-1  ! Step increment - overridden by next initiliztion
      Stmin = 2
!      print*, 'Delt = ', Delt
      
!      Dtgl = 1d0/Delt
!      Dtgl = Dtgl*Dtgl
         gami   = 1d0

      gami_f = (gami-1d0)/gami

      do istep = 1, Nstep       ! LOOP OVER TIME STEPS

         StepCount = StepCount+1

 1001    continue

            write(*,'(a)') "###########################################"
            write(*,'(a)') "## Time Step Number and size, Time       ##"
            write(*,'(i7,2e16.6)') StepCount,Delt,time

c........Backward Euler iteration
!
         if (myid.eq.mpi_master) then
         write(*,'(a)')  "--Backward Euler --"
         endif		 
         method = 1

         gami   = 1d0
         alfi   = 1d0
         almi   = 1d0
         gami_f = 0d0

         yg  = ygold           ! Predict 
         acg = acgold * gami_f

         call solflowGMRES(
     &        rowu,  colu,
     &        icntu)

         ygBE  = yg
c........End Backward Euler iteration

         write(*,'(a)')  "-- Generalized Alpha --"

         method = 2

         almi   = (3d0-rhoinf)/(1d0+rhoinf)/2d0
         alfi   = 1d0/(1d0+rhoinf)
         gami   = 5d-1+almi-alfi
         gami_f = (gami-1d0)/gami

      ! Predict
         yg  = ygold
         acg = acgold * gami_f
 
 !!!!!!!!!!
!      solf = 15
!      fname = trim('yg_pre_GMRES')
!      open (solf, file=fname, status='unknown')
!      do i = 1, NNODZu
!         write(solf,*) yg(i,1),yg(i,2),yg(i,3),yg(i,4),yg(i,5)
!      enddo
!      close(solf)

!      fname = trim('acg_pre_GMRES')
!      open (solf, file=fname, status='unknown')
!      do i = 1, NNODZu
!         write(solf,*) acg(i,1),acg(i,2),acg(i,3),acg(i,4),acg(i,5)
!      enddo
!      close(solf)
!!!!!!!!!!


      ! Flow Solver For this time step
         call solflowGMRES(
     &        rowu,  colu,
     &        icntu)

 !!!!!!!!!!
!      solf = 15
!      fname = trim('yg_post_GMRES')
!      open (solf, file=fname, status='unknown')
!      do i = 1, NNODZu
!         write(solf,*) yg(i,1),yg(i,2),yg(i,3),yg(i,4),yg(i,5)
!      enddo
!      close(solf)

!      fname = trim('acg_post_GMRES')
!      open (solf, file=fname, status='unknown')
!      do i = 1, NNODZu
!         write(solf,*) acg(i,1),acg(i,2),acg(i,3),acg(i,4),acg(i,5)
!      enddo
!      close(solf)
!!!!!!!!!!
         Dtold = Delt  ! increment in timestep

!          if (Stepcount==1) Delt = 1d-1
!          if (Stepcount==100) Delt = 1d-1
!          if (Stepcount==200) Delt = 25d-2
!          if (Stepcount==250) Delt = 50d-2
!          if (Stepcount==300) Delt = 1d0
c........Time step correction
         if (mod(istep,1)==0) then  ! Backward Euler iteration
          tolt = 1d-3      ! Tolerance

      ! RPD Change
!          !  Structural variables DUMMY - u_1, u_2
          rtdum(1) =  sum((yg(:,1)-ygBE(:,1))*(yg(:,1)-ygBE(:,1)))
     &              + sum((yg(:,2)-ygBE(:,2))*(yg(:,2)-ygBE(:,2)))

          rtdum(2) = sum(yg(:,1)*yg(:,1)) + sum(yg(:,2)*yg(:,2))

!          rtdumg = rtdum
!          call MPI_ALLREDUCE (rtdumg, rtdum, 2,
!     &         MPI_DOUBLE_PRECISION,MPI_SUM, 
!     &         MPI_COMM_WORLD,mpi_err)		  
		  
          rdum = dsqrt(rtdum(1)/rtdum(2))

!      ! RPD Change
          !  Momentum equation - v_1, and v_2
          rtmom(1) =  sum((yg(:,3)-ygBE(:,3))*(yg(:,3)-ygBE(:,3)))
     &             + sum((yg(:,4)-ygBE(:,4))*(yg(:,4)-ygBE(:,4)))

          rtmom(2) = sum(yg(:,3)*yg(:,3))+ sum(yg(:,4)*yg(:,4))
!
!         rtmomg = rtmom
!          call MPI_ALLREDUCE (rtmomg, rtmom, 2,
!     &         MPI_DOUBLE_PRECISION,MPI_SUM, 
!     &         MPI_COMM_WORLD,mpi_err)		  
!		  
          rmom = dsqrt(rtmom(1)/rtmom(2))
!
!      ! RPD Change
!          ! Temperature
          rtTemp(1) = sum( (yg(:,5)-ygBE(:,5))*(yg(:,5)-ygBE(:,5)) )
          rtTemp(2) = sum( yg(:,5)*yg(:,5) )
!
!          rtTempg = rtTemp
!          call MPI_ALLREDUCE (rtTempg, rtTemp, 2,
!     &         MPI_DOUBLE_PRECISION,MPI_SUM, 
!     &         MPI_COMM_WORLD,mpi_err)			  
!		  
          rTemp = dsqrt(rtTemp(1)/rtTemp(2))
!
!          ! Maximum Error
          rtot = max(rTemp,rmom,rdum)

          write(*,*)  "Error", rtot

          ! First iteration
          if (Stepcount.lt.Stmin) then 
             rtot = 9d-1*tolt ! reduce tolerance by 0.9 factor
          endif

          ! Second iteration onwards
          if (Stepcount.ge.Stmin) Delt=Delt*dsqrt(7.d-01*tolt/rtot)

          if (rtot.gt.tolt) then  ! If error greater than tolerance
            ireject = ireject + 1
!           if (myid.eq.mpi_master) then			
            write(*,'(a)') "-- Time Step Rejected --"
            write(*,*) "Repeat Time Step: ", istep
!           endif			
            goto 1001
          else  ! If error less than tolerance (converged)
!      endif ! adaptive time step

                 time = time + Dtold  ! Step increment
                 ! write statistics

                 Vener(istep) = Ener   ! Dimensionless Energy
                 Vtime(istep) = time   ! Time
                 Vdelt(istep) = Delt   ! Step increment
                 Vtemp(istep) = AvgTemp   ! Avg Temp				 
                 VSxx(istep) = AvgSxx   ! Avg Sxx
                 VSxy(istep) = AvgSxy   ! Avg Sxy
                 VSyy(istep) = AvgSyy   ! Avg Syy

                 VExx(istep) = AvgExx   ! Avg Exx
                 VExy(istep) = AvgExy   ! Avg Exy
                 VEyy(istep) = AvgEyy   ! Avg Eyy

                 do jstep = 1, istep
                    fname = trim('AvgTemp.dat' // cname(Stepcount))
                    open(unit=29, file=fname,    status='unknown')
                    write(29,*) Vtime(jstep), Vtemp(jstep)
                 enddo

                 close(29)						 
				 
                 do jstep = 1, istep
                    fname = trim('tmstep.dat' // cname(Stepcount))
                    open(unit=24, file=fname,    status='unknown')
                    write(24,*) Vtime(jstep), Vdelt(jstep)
                 enddo

                 close(24)

                 do jstep = 1, istep
                    fname = trim('energy.dat' // cname(Stepcount))
                    open(unit=26, file=fname,    status='unknown')
                    write(26,*) Vtime(jstep), Vener(jstep)
                 enddo

                 close(26)
				 
                 do jstep = 1, istep
                    fname = trim('resultados.dat' // cname(Stepcount))
                    open(unit=37, file=fname,     status='unknown')
                    write(37,*) jstep, Vtime(jstep)
                 enddo

                 close(37)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Stress
      if ((StrStn == 1).or.(StrStn == 2).or.(StrStn == 3)) then ! save sxx, and exx
                 do jstep = 1, istep
                    fname = trim('AvgSxx.dat' // cname(Stepcount))
                    open(unit=41, file=fname,    status='unknown')
                    write(41,*) Vtime(jstep), VSxx(jstep)
                 enddo
                 close(41)

                 do jstep = 1, istep
                    fname = trim('AvgExx.dat' // cname(Stepcount))
                    open(unit=41, file=fname,    status='unknown')
                    write(41,*) Vtime(jstep), VExx(jstep)
                 enddo
                 close(41)
      endif ! save sxx, and exx


!      if ((StrStn == 2).or.(StrStn == 3)) then ! save syy, and eyy
!                 do jstep = 1, istep
!                    fname = trim('AvgSyy.dat' // cname(Stepcount))
!                    open(unit=34, file=fname,    status='unknown')
!                    write(34,*) Vtime(jstep), VSyy(jstep)
!                 enddo
!                 close(34)
!
!                 do jstep = 1, istep
!                    fname = trim('AvgEyy.dat' // cname(Stepcount))
!                    open(unit=44, file=fname,    status='unknown')
!                    write(44,*) Vtime(jstep), VEyy(jstep)
!                 enddo
!                 close(44)
!      endif ! save syy, and eyy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

               if (mod(istep,ifreq)==0) call main_prn ! Restart printout
               ! AFter every ifreq (=10) loops printout
               if (avg>0 )then        !sample simulation results
               if(mod(StepCount,avg)==0) call mean
               endif

             avgdt = avgdt + Delt  ! Increment

!      if (tstepadp == 1) then ! adaptive time step
            endif  ! End  - If error greater than tolerance
         endif      ! End Backward Euler iteration
!      endif ! adaptive time step

c........End time step correction

         ygold  = yg       ! Displacement and Temperature variables
         acgold = acg      ! Velocities - Variables

         print *,'Rejected time step =', ireject  ! rejected frequency
         print *,'Average time step =', avgdt/dfloat(istep)
         ! dfloat - convert integer to double precision
   
      enddo  ! Loop over time steps

          call main_prn
          istep = istep - 1

      if (mod(istep,ifreq)/=0) call main_prn
      ! If reminder is not equal to zero call main prn
      ! mod returns reminder  /= is logical not equal

      call main_end

      print*,'End of step istep            =', istep

      end ! end program SMA2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! main_init
 
      subroutine main_init

      use maindef

      implicit none

!      real*8 xran
      integer*4 iostat

      write(*,*) '--- Initialising ---'
      call processors_u  ! Read processors used in XYZ dir, IC and BCs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c...   Material Inputs

      ! Material Constants
      A11 = 140.0d9
      A33 = 280.0d9
      A22 = 212.0d9
      A24 = 1.70d13
      A26 = 3.00d16
!      TempIC = 270d0    ! Initial temperature
      Tm = 265d0
      Rho = 10000d0
      Eta = 0.025
      kg = 3.15d-8     ! Ginzburg Constant
      kappa = 50d0    ! Thermal conductivity
      Cv = 350d0  ! Sp. heat

      ! Loadings
      fx = 0d0
      fy = 0d0

      ! Rescaled Constants
      delta = sqrt(kg*A26/A24**2) ! scaling factor - dimension earlier
      H = 1d0   ! Scale factor for A22
      ec = sqrt(abs(A24)/(A26*H)) ! Strain scaling constant
      A0 = A26*ec**4  !  Scaling for elastic constants
      aa1T = A11/A0
      aa3T = A33/A0
      aa4T = 1d0
      aa6T = 1d0
      aa2T = A22*A26/A24**2 ! Structural
      Fz0 = A26*ec**6  !  Free Energy scale factor
      StressScale = A26*ec**5
      timec = sqrt(Rho*delta**2/A0)  !  Scaled time
      EtaM = Eta*timec/(Rho*delta**2)  !  Scaled Damping
      Tc = Tm  !  Temperature Constant
!      CvC = Rho*Cv*Tc/timec  !  CV constant
      CvC = Rho*Cv/timec  !  CV constant	  
!      kT = kappa*Tc/delta**2/CvC
      kT = kappa/delta**2/CvC	  
!      NLT = 1d0/sqrt(2d0)*A22*ec**2/Rho/Cv/Tm
      NLT = 1d0/sqrt(2d0)*A22*ec**2/timec/CvC/Tm	  
      gth = 0d0
      kgT = 1d0
!      scalev = ec*delta/timec

      ! Geometry
      length = lenx*delta   ! nanometer
      height = leny*delta

      ! Rescaled Geomtry
      lenT = length/delta ! scaled length
      htT = height/delta   ! scaled height
      scalex = lenT
      scaley = htT
      scalez = 1d0

      ! Loading
      tload = 1.0d-9
      tunload = 1.0d-9
      tloadT = tload/timec
      tunloadT = tunload/timec
      ttotT = tloadT+tunloadT
      tloadbyunloadT = tloadT/tunloadT
      StrainReq = 3d-2
      StrainT = StrainReq/ec
      uT = lenT*StrainT
      velT = uT/tloadT  !/2d0 - modification for one side pull

      ! Rescaled Constants - dimensionless Weak form
      rhoT = 1d0    ! Density (rescaled)
      cvT = 1d0     ! Sp. heat. (rescaled)
      a2th = NLT    ! thermal a2 term - rescaled
      etaT = EtaM   ! thermal a2 term - rescaled
      kappaT = kT  ! thermal a2 term - rescaled
      kgTS = kgT/2d0   ! gradient coefficient (rescaled)
      Ca2th = -1d0/2d0*a2th  !/Tm

      aa1        = aa1T/2d0
      aa2        = aa2T/2d0
      aa4        = aa4T/4d0
      aa6        = aa6T/8d0
      aa3        = aa3T/4d0

      write(*,*) '--- Printing Material Inputs ---'

        print *,'A11                  =', A11
        print *,'A33                  =', A33
        print *,'A22                  =', A22
        print *,'A24                  =', A24
        print *,'A26                  =', A26
        print *,'TempIC (TEMPERATURE) =', TempIC
        print *,'Tm                   =', Tm
        print *,'Rho                  =', Rho
        print *,'Eta                  =', Eta
        print *,'kg                   =', kg
        print *,'fx                   =', fx
        print *,'fy                   =', fy
        print *,'kappa                =', kappa
        print *,'Cv                   =', Cv

      write(*,*) '--- Printing Geometry Inputs ---'
        print *,'length               =', length
        print *,'height               =', height
        print *,'lenT                 =', lenT
        print *,'htT                  =', htT

      write(*,*) '--- Printing Rescaled Constants ---'
        print *,'delta                =', delta
        print *,'ec                   =', ec
        print *,'A0                   =', A0
        print *,'aa1T                 =', aa1T
        print *,'aa2T                 =', aa2T
        print *,'aa3T                 =', aa3T
        print *,'aa4T                 =', aa4T
        print *,'aa6T                 =', aa6T
        print *,'aa1                  =', aa1
        print *,'aa3                  =', aa3
        print *,'aa2                  =', aa2
        print *,'aa4                  =', aa4
        print *,'aa6                  =', aa6

        print *,'Fz0                  =', Fz0
        print *,'StressScale          =', StressScale
        print *,'timec                =', timec
        print *,'EtaM                 =', EtaM
        print *,'Tc (or Tm)           =', Tc
        print *,'CvC                  =', CvC
        print *,'CvT                  =', CvT
        print *,'kT                   =', kT
        print *,'NLT                  =', NLT

        print *,'gth                  =', gth
        print *,'rhoT                 =', rhoT
        print *,'cvT                  =', cvT
        print *,'a2th                 =', a2th
        print *,'kappaT               =', kappaT
        print *,'kgTS                 =', kgTS
        print *,'Ca2th                =', Ca2th

      write(*,*) '--- Printing Loading Inputs ---'
        print *,'tload                =', tload
        print *,'tloadT               =', tloadT
        print *,'tunloadT             =', tunloadT
        print *,'ttotT                =', ttotT
        print *,'tloadbyunloadT       =', tloadbyunloadT
        print *,'StrainReq            =', StrainReq
        print *,'StrainT              =', StrainT
        print *,'uT                   =', uT
        print *,'velT                 =', velT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Simulation Parameters
      Avgu1      = 0d0  ! Mean avg IC value of u1
      Avgu2      = 0d0  ! Mean avg IC value of u2
      distr      = 1d0 ! distribution

!      xnu        = 1d-1
      rpm        = 0d0 
c      theta      = 0.85d0
!      Delt       = 1d-2     ! Step increment- overrides earier init.
      rhoinf     = 5d-1
!      Utol       = 1d-3
!      NLtol      = 1d-1
      Nstep      = 1000000		! 100000
!      Nnewt      = 7
      NDOF       = 5
      ifreq      = samplesave   !50   ! Frequencies of restarts
      Kspaceu_mn = 15   ! Min # of iterations of GMRES
      Kspaceu    = 300	! Max # of iterations of GMRES
      Ksp_scale  = 300	 
      avg        = samplesave  !  50    ! FFT's frequency 0-Takes no samples
      iavg       = 0
      iavgfreq   = 1
      U_mx       = 1d-0

!      epsilon    = 0.013d0*xnu      ! Dummy variable assignment
!      lambda     = 0.001d0*xnu*xnu  ! Dummy variable assignment
!
!c      epsilon    = 15d-3                ! This is 1/Re! OK in 128 ^ 2
!c      lambda     = epsilon*epsilon/25d-1 !This is Ca ^ 2! OK in 128^2
!
!      epsilon    = 1d0/2d0/256d0       ! This is 1/Re! OK in 128 ^ 2
!      lambda     = 1d0/256d0/256d0     ! This is Ca ^ 2! OK in 128 ^ 2


      if (Kspaceu < Ksp_scale) Ksp_scale = Kspaceu
      if (Delt > Delt_mx) Delt_mx = Delt
      if (Delt < Delt_mn) Delt_mn = Delt

c      if (myid==mpi_master) then
         print *,'Delt             =', Delt
!         print *,'xnu              =', xnu
         print *,'rpm              =', rpm
!         print *,'theta            =', theta
!         print *,'lambda           =', lambda
!         print *,'epsilon          =', epsilon
         print *,'rhoinf           =', rhoinf
         print *,'Utol             =', Utol
         print *,'Nstep            =', Nstep
         print *,'Nnewt            =', Nnewt
         print *,'ifreq            =', ifreq
         print *,'Kspaceu_mn       =', Kspaceu_mn
         print *,'Ksp_scale        =', Ksp_scale
         print *,'Kspaceu          =', Kspaceu
         print *,'avg(0:no sample) =', avg
         print *,'iavgfreq         =', iavgfreq
c      endif

c...  get main input
      print*,'Reading Displ. Data'
      call input_3D_u              ! Read Displ. Data
!      print*,'Completed reading Displ. Data'

!!!!
!     restart file
!      restf = 99
!      fname = trim('IBCValues')
!      open (restf, file=fname, status='unknown')
!      do i = 1, NNODZu
!         write(restf,*) IBC(i,1),IBC(i,2),IBC(i,3),IBC(i,4),IBC(i,5)
!      enddo
!      close(restf)

!      restf = 99
!      fname = trim('DIRBCValues')
!      open (restf, file=fname, status='unknown')
!      do i = 1, NNODZu
!         write(restf,*) DIR_BC(i,1),DIR_BC(i,2),DIR_BC(i,3),
!     &                  DIR_BC(i,4),DIR_BC(i,5)
!      enddo
!      close(restf)
!!!!!
      
c      if (myid==mpi_master) print*,'Communication types'
c      call ctypes                  ! Communication types  (real)


	     dummyvar1 = MCPu - Pu
		    print *,'Mesh Size           =', dummyvar1
      NGAUSS  = max(Pu,Qu,Ru)+1
      print *,'NGAUSS           =', NGAUSS

      NSHLu = (Pu+1)*(Qu+1)*(Ru+1)
      print *,'Total DOFs           =', NSHLu
      print*,'Allocate Solution Vectors'

      allocate(rowu (NNODZu*8*NSHLu))
      allocate(colu (NNODZu+1))

      allocate(yg    (NNODZu,NDOF))
      allocate(ygBE  (NNODZu,NDOF))
      allocate(ygold (NNODZu,NDOF))
      allocate(acg   (NNODZu,NDOF))
      allocate(acgold(NNODZu,NDOF))

      allocate(Vener(Nstep))
      allocate(Vtime(Nstep))
      allocate(Vdelt(Nstep))
      allocate(Vtemp(Nstep))
      allocate(VSxx(Nstep))
      allocate(VExx(Nstep))
      allocate(VSxy(Nstep))
      allocate(VExy(Nstep))
      allocate(VSyy(Nstep))
      allocate(VEyy(Nstep))
	  
      !   CONTINUITY CONSTRAINT
      print*,'Continuity'
      call continuity
!      print*,'Completed reading - Continuity'

      print*,'Generating Data Structures'

c...  generate IEN and  Coordinates
      call genIEN_INN_3D_u
      print*,'Completed - genIEN_INN_3D_u'

c...  generate mesh face/interior loading attributes 
      call gen_FACE_IEN_3D_u
      print*,'Completed - gen_FACE_IEN_3D_u'

c...  generate Sparse Archit.
      call genSparStruc_3D_u  (colu,  rowu,  icntu)

c      if (myid==mpi_master)
      print*,'Allocate sparse matrices'

      allocate(RHSGu(NNODZu,NDOF)) ! global RHS
      allocate(LHSK(NDOF*NDOF,icntu))   ! global LHS

      mu  = VISC(1)
      allocate(xor(2))
      xor = 0d0

      !   GET INITIAL CONDITION
      print*,'Reading Initial Data'

      call initial_data
      if (initcond==3) then	  ! microstructure as an initial condition
         time = 0d0
         StepCount = 0
         Delt = 1d-2
      endif

!!!!!!!!!!
!      solf = 15
!      fname = trim('yg_post_init_data')
!      open (solf, file=fname, status='unknown')
!      do i = 1, NNODZu
!         write(solf,*) ygold(i,1),ygold(i,2),ygold(i,3),
!     &                 ygold(i,4),ygold(i,5)
!      enddo
!      close(solf)


!      fname = trim('acg_post_init_data')
!      open (solf, file=fname, status='unknown')
!      do i = 1, NNODZu
!         write(solf,*) acgold(i,1),acgold(i,2),acgold(i,3),
!     &                 acgold(i,4),acgold(i,5)
!      enddo
!      close(solf)
!!!!!!!!!!
	  
      print*,'Completed Reading Initial Data'

      write(*,*) '---------------------------------'
      write(*,*) '     Finished Initialising '
      write(*,*) '---------------------------------'

      return
      end

!=======================================================================

      subroutine main_end

      use maindef

      implicit none
      character*10 cname

c      if (myid == mpi_master) 
      print*,'Simulation complete'
c      call MPI_Finalize(mpi_err)

      return
      end

!=======================================================================
      subroutine continuity

      use maindef

      implicit none

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     U-direction
!     Continuity enforcement: Restriction operator

      allocate(TuMu(2*Pu,Pu))   ! Velocity
      allocate(TuSu(2*Pu,Pu))

      TuSu = 0.d0; TuMu = 0.d0

      call constraint_input(  ! Reads restriction operators
     *     TuMu,TuSu,           !         Master, Slave
     *     Pu,myid+1,lpu,"u","u")

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     V-direction
!     Continuity enforcement: Restriction operator

      allocate(TvMu(2*Qu,Qu))   ! Velocity
      allocate(TvSu(2*Qu,Qu))

      TvSu = 0.d0; TvMu = 0.d0

      call constraint_input(  ! Reads restriction operators
     *     TvMu,TvSu,           !         Master, Slave
     *     Qu,myid+1,lpv,"u","v")

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     W-direction
!     Continuity enforcement: Restriction operator

      allocate(TwMu(2*Ru,Ru))   ! Velocity
      allocate(TwSu(2*Ru,Ru))

      TwSu = 0.d0; TwMu = 0.d0

      ! Since 2D problem no function is called in z direction
       TwSu(1:2*Ru,1:Ru) = 1d0
       TwMu(1:2*Ru,1:Ru) = 1d0

c$$$      call constraint_input(  ! Reads restriction operators
c$$$     *     TwMu,TwSu,           !         Master, Slave
c$$$     *     Ru,myid+1,lpw,"u","w")

      return
      end

!=======================================================================
      subroutine initial_data    ! GET INITIAL CONDITION

      use maindef
      use common
      use aAdjKeep

      implicit none

      character*10 cname

      real*8 random, rru, rrul, da, funct

      real*8 c1x, c1y, c1z  ! centre of bubbles
!      real*8 d1, d2, d3  ! diameters of bubbles

      integer*4 idof, inod, iel, icount

      real*8 xl(NSHLu,NSD)

c------MEAN
      if (avg>0) call mean_read

c------INPUT

      ygold  = 0d0
      acgold = 0d0

      solf = 15

      fname = trim('restart.last' // cname(myid+1))
      open (solf, file=fname, status='old',iostat=ierr)
      ! ierr-error flag

      if (ierr == 0) then ! if ierr flag == 0

         print*, ' Reading File ', fname

         read(solf,*) StepCount,Delt,istep,time

         do i = 1, NNODZu
            read(solf,*)  ygold(i,1),ygold(i,2),ygold(i,3),ygold(i,4),
     &                    ygold(i,5)
         enddo

         do i = 1, NNODZu
          read(solf,*) acgold(i,1),acgold(i,2),acgold(i,3),acgold(i,4),
     &                 acgold(i,5)
        enddo

         close(solf)

      else

         restf = 99
         fname = trim('restart.fix' // cname(myid+1))
         open (restf, file=fname, status='old',iostat=ierr)

         if (ierr == 0) then

            print*, ' Reading File ', fname

            read(restf,*) StepCount,Delt,istep,time

            do i = 1, NNODZu
               read(restf,*) ygold(i,:)
            enddo

            close (restf)

         else


      do iel = 1, NELu
          
c...  Check to see if current element has nonzero area
         ni = INNu(IENu(iel,1),1) ! get NURB coordinates
         nj = INNu(IENu(iel,1),2)
         nk = INNu(IENu(iel,1),3)


          if ( (U_KNOTu(ni).ne.U_KNOTu(ni+1)).and.
     &         (V_KNOTu(nj).ne.V_KNOTu(nj+1)).and.
     &         (W_KNOTu(nk).ne.W_KNOTu(nk+1))) then
            
            da = (U_KNOTu(ni+1) - U_KNOTu(ni))*
     &           (V_KNOTu(nj+1) - V_KNOTu(nj))*
     &           (W_KNOTu(nk+1) - W_KNOTu(nk))/8d0 
                                
            icount = 0
        
            do k = 0, Ru
                do j = 0, Qu
                    do i = 0, Pu
                        icount = icount + 1
                        xl(icount,:) = B_NETu(ni-i,nj-j,nk-k,1:NSD)
                    enddo
                enddo
            enddo

		     if (initcond == 1) then ! initcond = 1 Random IC
               iseed = (3491+myid*17)*2+1
               call randomSeed(iseed)
               do i = 1, NNODZu
               ygold (i,1) =  Avgu1 + (random(iseed)-.5)*2d-2
               ygold (i,2) =  Avgu2 + (random(iseed)-.5)*2d-2
                 ygold (i,3:4) = 0d0
                 ygold (i,5) = TempIC
               enddo
		     endif ! initcond = 1 Random IC

		     if (initcond == 0) then ! initcond = 0 Parabolic IC
              do i = 1, NSHLu
               dd1 = (lenT-xl(i,1))*xl(i,1)*(htT-xl(i,2))*xl(i,2)*1d-5
               dd2 = (lenT-xl(i,1))*xl(i,1)*(htT-xl(i,2))*xl(i,2)*1d-5
               ygold (IENu(iel,i),1) = dd1
               ygold (IENu(iel,i),2) = dd2
               ygold (IENu(iel,i),3:4) = 0d0
               ygold (IENu(iel,i),5) = TempIC
              enddo
		     endif ! initcond = 0 Parabolic IC

            c1x = lenT/2d0
            c1y = htT/2d0
!            c1z = lenT/2d0
             if (initcond == 2) then ! initcond = 0 exponential
           do i = 1, NSHLu
            dd1 = exp(-(xl(i,1)-c1x)**2)*
     &            exp(-(xl(i,2)-c1y)**2)*5e-3
            dd2 = -exp(-(xl(i,1)-c1x)**2)*
     &            exp(-(xl(i,2)-c1y)**2)*5e-3
!            dd1 = exp(-(xl(i,1)-c1x)**2*5e4)*
!     &            exp(-(xl(i,2)-c1y)**2*5e5)*5e-4
!            dd2 = -exp(-(xl(i,1)-c1x)**2*5e4)*
!     &            exp(-(xl(i,2)-c1y)**2*5e5)*5e-4	
                 ygold (IENu(iel,i),1) = dd1
                 ygold (IENu(iel,i),2) = dd2
                 ygold (IENu(iel,i),3:4) = 0d0
                 ygold (IENu(iel,i),5) = TempIC
               enddo
             endif ! initcond = 0 Parabolic IC

             if (initcond == 4) then ! initcond = 0 zero
           do i = 1, NSHLu
                 ygold (IENu(iel,i),1) = 0d0
                 ygold (IENu(iel,i),2) = 0d0
                 ygold (IENu(iel,i),3:4) = 0d0
                 ygold (IENu(iel,i),5) = TempIC
               enddo
             endif ! initcond = 0 Parabolic IC
			 
         endif ! Non-zero area element If statement

        enddo ! iel

!!!!!!!!!!
!      solf = 15
!      fname = trim('ygold_initiated_pre_b4shuffle')
!      open (solf, file=fname, status='unknown')
!      do i = 1, NNODZu
!         write(solf,*) ygold(i,1),ygold(i,2),ygold(i,3),
!     &                 ygold(i,4),ygold(i,5)
!      enddo
!      close(solf)


!      fname = trim('acg_initiated_pre_b4shuffle')
!      open (solf, file=fname, status='unknown')
!      do i = 1, NNODZu
!         write(solf,*) acgold(i,1),acgold(i,2),acgold(i,3),
!     &                 acgold(i,4),acgold(i,5)
!      enddo
!      close(solf)
!!!!!!!!!!

            call shuffle_u(
     &           ygold,
     &           TuMu,TuSu,lpu, Pu,
     &           NNODZu, MCPu, NCPu, OCPu, NDOF,
     &           'out','u')

            call shuffle_v(
     &           ygold,
     &           TvMu,TvSu,lpv, Qu,
     &           NNODZu, MCPu, NCPu, OCPu, NDOF,
     &           'out','u')


            call shuffle_w(
     &           ygold,
     &           TwMu,TwSu,lpw, Ru,
     &           NNODZu, MCPu, NCPu, OCPu, NDOF,
     &           'out','u')


            StepCount = 0
            count = 0

         endif ! else statement restart.fix
      endif   ! if ierr flag == 0 restart.last

!!!!!!!!!!
!      solf = 15
!      fname = trim('yg_initiated_post_shuffle')
!      open (solf, file=fname, status='unknown')
!      do i = 1, NNODZu
!         write(solf,*) ygold(i,1),ygold(i,2),ygold(i,3),
!     &                 ygold(i,4),ygold(i,5)
!      enddo
!      close(solf)


!      fname = trim('acg_initiated_post_shuffle')
!      open (solf, file=fname, status='unknown')
!      do i = 1, NNODZu
!         write(solf,*) acgold(i,1),acgold(i,2),acgold(i,3),
!     &                 acgold(i,4),acgold(i,5)
!      enddo
!      close(solf)
!!!!!!!!!!

    ! Impose Dirichlet BC strongly
!      if((periodicx==0).or.(periodicy==0).or.(periodicz==0)) then !non-periodic
          where (IBC==1)
             ygold  = DIR_BC
             acgold = 0d0
          end where
!      endif

      ! ADJUST for Strong BCs on velocity
!            where (IBC==1)
!            ygold  = DIR_BC
!            acgold = 0d0
!            end where

!!!!!!!!!!	  
!      solf = 15
!      fname = trim('yg_initiated_post_IBC')
!      open (solf, file=fname, status='unknown')
!      do i = 1, NNODZu
!         write(solf,*) ygold(i,1),ygold(i,2),ygold(i,3),
!     &                 ygold(i,4),ygold(i,5)
!      enddo
!      close(solf)


!      fname = trim('acg_initiated_post_IBC')
!      open (solf, file=fname, status='unknown')
!      do i = 1, NNODZu
!         write(solf,*) acgold(i,1),acgold(i,2),acgold(i,3),
!     &                 acgold(i,4),acgold(i,5)
!      enddo
!      close(solf)  
!!!!!!!!!!



      return
      end  ! end subroutine initial_data

!=======================================================================
!                      Output routines
!=======================================================================

      subroutine main_prn

      use maindef

      implicit none
      character*10 cname

!     solution file
      solf = 15
      
      fname = trim('restart.last' // cname(myid+1))

      open (solf, file=fname, status='unknown')
      write(solf,*) StepCount,Delt,istep,time

      do i = 1, NNODZu
         write(solf,*) yg(i,1),yg(i,2),yg(i,3),yg(i,4),yg(i,5)
      enddo

      do i = 1, NNODZu
         write(solf,*) acg(i,1),acg(i,2),acg(i,3),acg(i,4),acg(i,5)
      enddo

      close(solf)

!     restart file
      restf = 99

      fname = trim('restart' // cname(StepCount)) // cname (myid+1)

      open (restf, file=fname, status='unknown')
      write(restf,*) StepCount,Delt,istep,time

      do i = 1, NNODZu
         write(restf,*) yg(i,1),yg(i,2),yg(i,3),yg(i,4),yg(i,5)
      enddo

      close (restf)

      return
      end


!=======================================================================
      subroutine mean_read    ! GET INITIAL STATISTICS


      use maindef

      implicit none

c...  Local variables
      integer igauss, jgauss, kgauss, icount, iel, fac
      real*8
     &     shl(NSHLu), shg(NSHLu), shgradl(NSHLu,NSD),
     &     shgradg(NSHLu, NSD),
     &     std1, std2, std3,
     &     da, divu, divv, divw, pi

      pi = acos(-1.)

                                ! Generate Evaluation points
      NPTSu = 2*(MCPu-Pu)
      NPTSv = 2*(NCPu-Qu)
      NPTSw = 2*(OCPu-Ru)

      allocate(PTSu(NPTSu))
      allocate(PTSv(NPTSv))
      allocate(PTSw(NPTSw))
      divu = 1d0/real(NPTSu)
      divv = 1d0/real(NPTSv)
      divw = 1d0/real(NPTSw)

      do i = 1, NPTSu
         PTSu(i) = U_KNOTu(1)
     &         + real(i-1)*(U_KNOTu(MCPu+Pu+1)-U_KNOTu(1))*divu
      enddo

      do i = 1, NPTSv
         PTSv(i) = V_KNOTu(1)
     &         + real(i-1)*(V_KNOTu(NCPu+Qu+1)-V_KNOTu(1))*divv
      enddo

      do i = 1, NPTSw
         PTSw(i) = W_KNOTu(1) +
     &        real(i-1)*(W_KNOTu(OCPu+Ru+1)-W_KNOTu(1))*divw
      enddo

c------MEAN
                                ! Initialize averages
      allocate(SOL_AR_AVG(NPTSu,NPTSv,NPTSw,NDOF))
      allocate(STRAIN_AVG(NPTSu,NPTSv,NPTSw,3))

      divfac = dble(NPTSw*NPTSv*NPTSu)*dble(numnodes)
      divfac = 1d0/divfac

      return
      end


c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine mean    ! compute mean

      use maindef

      implicit none

c...  Local variables
      integer  igauss, jgauss, kgauss,  iel
      real*8
     &     shg(NSHLu), shgradl(NSHLu,NSD), shhessg(NSHLu,NSD,NSD),
     &     dxidx(NSD,NSD), BXM, BXP, BYM, BYP, BZM, BZP,
     &     divu1, divv1, divw1, nurbcoordx, nurbcoordy, nurbcoordz
      real*8
     &     yl (NSHLu,NDOF),
     &     nurbcoord(3), std1, std2, std3,
     &     ui, vi, wi, pri, da, temp, ddu1dx, ddu1dy, ddu2dx, ddu2dy


      do igauss = 1, NPTSu
      do jgauss = 1, NPTSv
      do kgauss = 1, NPTSw

         nurbcoord(1) = PTSu(igauss)
         nurbcoord(2) = PTSv(jgauss)
         nurbcoord(3) = PTSw(kgauss)

         do iel = 1, NELu

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Scaling for correction

      NPTSu1 = 1*(MCPu-Pu)
      NPTSv1 = 1*(NCPu-Qu)
      NPTSw1 = 1*(OCPu-Ru)

      divu1 = 1d0/real(NPTSu1)
      divv1 = 1d0/real(NPTSv1)
      divw1 = 1d0/real(NPTSw1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! X
         if (mod(iel,(MCPu-Pu))==0) then
            BXM = U_KNOTu(1) + (U_KNOTu(MCPu+Pu+1)-U_KNOTu(1))*divu1
     &              *real((MCPu-Pu-1))
            BXP = U_KNOTu(1) + (U_KNOTu(MCPu+Pu+1)-U_KNOTu(1))*divu1
     &              *real((MCPu-Pu))
         end if
         if (mod(iel,(MCPu-Pu))> 0) then
         BXM = U_KNOTu(1) + (U_KNOTu(MCPu+Pu+1)-U_KNOTu(1))*divu1*
     &              real((mod(iel,(MCPu-Pu))-1))
         BXP = U_KNOTu(1) + (U_KNOTu(MCPu+Pu+1)-U_KNOTu(1))*divu1*
     &              real((mod(iel,(MCPu-Pu))))
         end if

!!!!!!!!! Y
         if (mod(iel,(MCPu-Pu)) > 0) then
             BYM = V_KNOTu(1) + (V_KNOTu(NCPu+Qu+1)-V_KNOTu(1))*divv1*
     &              real(((iel/(MCPu-Pu))-(iel/((MCPu-Pu)*(NCPu-Qu)))
     &              *(NCPu-Qu)))
             BYP = V_KNOTu(1) + (V_KNOTu(NCPu+Qu+1)-V_KNOTu(1))*divv1*
     &              real(((iel/(MCPu-Pu))+1-(iel/((MCPu-Pu)*(NCPu-Qu)))
     &              *(NCPu-Qu)))
          end if

         if ((mod(iel,(MCPu-Pu)) == 0).and.
     &              (mod(iel,(MCPu-Pu)*(NCPu-Qu))>0)) then
              BYM = V_KNOTu(1) + (V_KNOTu(NCPu+Qu+1)-V_KNOTu(1))*divv1*
     &              real(((iel/(MCPu-Pu))-(iel/((MCPu-Pu)*(NCPu-Qu)))
     &              *(NCPu-Qu)-1))
              BYP = V_KNOTu(1) + (V_KNOTu(NCPu+Qu+1)-V_KNOTu(1))*divv1*
     &              real(((iel/(MCPu-Pu))-(iel/((MCPu-Pu)*(NCPu-Qu)))
     &              *(NCPu-Qu)))
         end if

         if ((mod(iel,(MCPu-Pu))==0).and.
     &              (mod(iel,(MCPu-Pu)*(NCPu-Qu))==0)) then
               BYM = V_KNOTu(1)+(V_KNOTu(NCPu+Qu+1)-V_KNOTu(1))*divv1*
     &              real((NCPu-Qu-1))
               BYP = V_KNOTu(1)+(V_KNOTu(NCPu+Qu+1)-V_KNOTu(1))*divv1
     &              *real((NCPu-Qu))
         end if

!!!!!!!!! Z
         if ((mod(iel,((MCPu-Pu)*(NCPu-Qu))))==0) then
              BZM = W_KNOTu(1) + (W_KNOTu(OCPu+Ru+1)-W_KNOTu(1))*divw1*
     &              real(((iel/((MCPu-Pu)*(NCPu-Qu)))-1))
              BZP = W_KNOTu(1) + (W_KNOTu(OCPu+Ru+1)-W_KNOTu(1))*divw1*
     &              real(((iel/((MCPu-Pu)*(NCPu-Qu)))))
         end if
         if (mod(iel,(MCPu-Pu)*(NCPu-Qu))> 0) then
              BZM = W_KNOTu(1) + (W_KNOTu(OCPu+Ru+1)-W_KNOTu(1))*divw1*
     &              real((iel/((MCPu-Pu)*(NCPu-Qu))))
              BZP = W_KNOTu(1) + (W_KNOTu(OCPu+Ru+1)-W_KNOTu(1))*divw1*
     &              real((iel/((MCPu-Pu)*(NCPu-Qu))+1))
         end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         nurbcoordx = 0d0
         nurbcoordy = 0d0
         nurbcoordz = 0d0

         nurbcoordx = 2d0*( nurbcoord(1) - BXM)/( BXP - BXM ) - 1d0
         nurbcoordy = 2d0*( nurbcoord(2) - BYM)/( BYP - BYM ) - 1d0
         nurbcoordz = 2d0*( nurbcoord(3) - BZM)/( BZP - BZM ) - 1d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c...  Check to see if current element has nonzero area
            ni = INNu(IENu(iel,1),1) ! get NURB coordinates
            nj = INNu(IENu(iel,1),2)
            nk = INNu(IENu(iel,1),3)


            std1 = (nurbcoord(1)-U_KNOTu(ni))*
     &           (  nurbcoord(1)-U_KNOTu(ni+1))
            std2 = (nurbcoord(2)-V_KNOTu(nj))*
     &           (  nurbcoord(2)-V_KNOTu(nj+1))
            std3 = (nurbcoord(3)-W_KNOTu(nk))*
     &           (  nurbcoord(3)-W_KNOTu(nk+1))


            da = (U_KNOTu(ni+1) - U_KNOTu(ni))*
     &           (V_KNOTu(nj+1) - V_KNOTu(nj))*
     &           (W_KNOTu(nk+1) - W_KNOTu(nk))

                                ! I'm in the right element
            if ( (std1.le.0d0).and.
     &           (std2.le.0d0).and.
     &           (std3.le.0d0).and.
     &           (da  .ne.0d0)) then

c...  Get Element Shape functions and their gradients

               shg = 0d0        ! initialize

!               call eval_SHAPE_3D_int(
!     &              iel,
!     &              nurbcoord(1), nurbcoord(2), nurbcoord(3),
!     &              shg,shgradl)

               call eval_SHAPE_3D(
     &              iel,
     &              nurbcoordx, nurbcoordy, nurbcoordz,
     &              shg,shgradl,shhessg,dxidx)

               do i = 1,NSHLu
                  yl(i,:) = yg(IENu(iel,i),:)
               enddo

               ui  = sum(shg*yl(:,1))
               vi  = sum(shg*yl(:,2))	
               ddu1dx  = sum(shgradl(:,1)*yl(:,1))			
               ddu1dy  = sum(shgradl(:,2)*yl(:,1))	
               ddu2dx  = sum(shgradl(:,1)*yl(:,2))	
               ddu2dy  = sum(shgradl(:,2)*yl(:,2))				   
               temp= sum(shg*yl(:,5))

               goto 2000        ! Go to next point

            endif

         enddo

 2000        continue

        !  Variables to be printed
         SOL_AR_AVG(igauss,jgauss,kgauss,1) = ui
         SOL_AR_AVG(igauss,jgauss,kgauss,2) = vi
         SOL_AR_AVG(igauss,jgauss,kgauss,5) = temp
         STRAIN_AVG(igauss,jgauss,kgauss,1)=1d0/sqrt(2d0)
     &		                                *(ddu1dx+ddu2dy)
         STRAIN_AVG(igauss,jgauss,kgauss,2)=1d0/sqrt(2d0)
     &		                                *(ddu1dx-ddu2dy)
         STRAIN_AVG(igauss,jgauss,kgauss,3)=1d0/2d0
     &		                                *(ddu1dy+ddu2dx)		 
		 
! declare other array to store strain values

      enddo
      enddo
      enddo

      iavg = iavg + 1

      call prn4FFT
      open(unit=37, file='resultados.dat', status='unknown')
      write(37,*) iavg,time
      return
      end

!=======================================================================

      subroutine prn4FFT    ! compute mean

      use maindef

      implicit none

      character*10 cname

!!!!! udisp
      restf = 99
      fname = trim('udisp' // cname(iavg))
      open (restf, file=fname, status='unknown', form="unformatted")
      write(restf) SOL_AR_AVG(:,:,:,1)
      close(restf)

!!!!! vdisp
      restf = 97
      fname = trim('vdisp' // cname(iavg))
      open (restf, file=fname, status='unknown', form="unformatted")
      write(restf) SOL_AR_AVG(:,:,:,2)
      close(restf)
	  
!!!!! temperature
      restf = 98
      fname = trim('temperature' // cname(iavg))
      open (restf, file=fname, status='unknown', form="unformatted")
      write(restf) SOL_AR_AVG(:,:,:,5)
      close(restf)

!      restf = 86
!      fname = trim('e1strain' // cname(iavg))
!      open (restf, file=fname, status='unknown', form="unformatted")
!      write(restf) STRAIN_AVG(:,:,:,1)
!      close(restf)

!!!!! straindev = straindeviatoric
      restf = 83
      fname = trim('e2strain' // cname(iavg))
      open (restf, file=fname, status='unknown', form="unformatted")
      write(restf) STRAIN_AVG(:,:,:,2)
      close(restf)
	  
!      restf = 82
!      fname = trim('e3strain' // cname(iavg))
!      open (restf, file=fname, status='unknown', form="unformatted")
!      write(restf) STRAIN_AVG(:,:,:,3)
!      close(restf)
	  	  
	  ! define function to write strain files
      restf = 81
      fname = 'FFT.pts'
      open (restf, file=fname, status='unknown', form="formatted")
      write(restf,*) NPTSu, NPTSv, NPTSw
      close(restf)

      return
      end


