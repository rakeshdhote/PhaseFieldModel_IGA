c------------------------------------------------------------------------c
c                                                                        c
c        Main routine to call all the subroutines                        c
c                                                                        c
c                                                                        c
c------------------------------------------------------------------------c

      use aAdjKeep
      
      include "common.h"
      
      integer so, sf, count, icnt, i,ni,nj,nk,
     &  numproc, istep, inistep, fstep, inter,iunit,iunitg,iunitp
	 
c...  get main input
      iunit = 101
      iunitg = 112
      iunitp = 113	  

      open (iunit, file="step.dat", 
     &     form="formatted", status="old", iostat=ierr)
      read(iunit,*) numproc, inistep, fstep, inter ! 1  0   50000  10
      close(iunit)
	  
      open (iunitg, file="scalegeom.dat", 
     &     form="formatted", status="old", iostat=ierr)
      read(iunitg,*) scalex, scaley, scalez !  10   10  10	  
        print *,'scalex                  =', scalex		
        print *,'scaley                  =', scaley		
        print *,'scalez                  =', scalez			
      close(iunitg)
	  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	 		
      open (iunitp, file="meshupostpCASE.dat", 
     &     form="formatted", status="old", iostat=ierr)
      
      read(iunitp,*) Pu, Qu, Ru !  10   10  10	  	  
      read(iunitp,*) MCPu, NCPu, OCPu	  
        print *,'Pu                  =', Pu		
        print *,'Qu                  =', Qu		
        print *,'Ru                  =', Ru				  
        print *,'MCPu                  =', MCPu		
        print *,'NCPu                  =', NCPu		
        print *,'OCPu                  =', OCPu		

      allocate(U_KNOTu(MCPu+Pu+1))
      allocate(V_KNOTu(NCPu+Qu+1))
      allocate(W_KNOTu(OCPu+Ru+1))

      do i = 1,MCPu+Pu+1
         if (i < MCPu+Pu+1) then
            read (iunitp,"(F13.9)", ADVANCE = "NO") U_KNOTu(i) ! read U_KNOT
         else
            read (iunitp,"(F13.9)", ADVANCE = "YES") U_KNOTu(i)
         endif
      enddo

      do i = 1,NCPu+Qu+1
         if (i < NCPu+Qu+1) then
            read (iunitp,"(F13.9)", ADVANCE = "NO") V_KNOTu(i) ! read V_KNOT
         else
            read (iunitp,"(F13.9)", ADVANCE = "YES") V_KNOTu(i)
         endif
      enddo

      do i = 1,OCPu+Ru+1
         if (i < OCPu+Ru+1) then
            read (iunitp,"(F13.9)", ADVANCE = "NO") W_KNOTu(i) ! read W_KNOT
         else
            read (iunitp,"(F13.9)", ADVANCE = "YES") W_KNOTu(i)
         endif
      enddo	  
      close(iunitp)	  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	  
	  
      
      do istep = inistep, fstep, inter ! Loop every time step

        print *,'istep                  =', istep
!         write(*,*) istep
         
         do ni = 1, numproc  ! for all processors in parallel code

!        print *,'iproc                  =', ni				 
                                         
            call input_3D(ni,scalex,scaley,scalez) ! read mesh data
            
!            allocate(yg(NNODZ,NSD+1))
            allocate(yg(NNODZ,5))
            
            call results_3D(ni,istep) ! read result files
            
            NSHL = (P+1)*(Q+1)*(R+1) ! # of shape funcs
            NGAUSS = 3          ! Visualization
            
c...        Form projective control net

!            BW = B_NET
!            do i = 1,NSD
!               BW(:,:,:,i) = BW(:,:,:,i)*B_NET(:,:,:,NSD+1)
!            enddo ! Multiplied all columns by wieghing factor-4th colum
            
c...        generate IEN and  Coordinates
            
            call genIEN_INN_3D
!     IEN matrix, which relates element numbers and local node numbers
!     to the appropriate global node numbers
!     INN matrix, which relates global node
!     number to the "NURBS coordinates" of the node

            write(*,*) "Patch (Processor) No.   ", ni
            
            call VizGMV_3D(ni,numproc,istep,inistep)

            deallocate(yg)
            deallocate(U_KNOT)    
            deallocate(V_KNOT)
            deallocate(W_KNOT)
            deallocate(B_NET)
            deallocate(IPER)
            deallocate(IBC)
            deallocate(DIR_BC)
            deallocate(IBC_FACE)
            deallocate(LD_FACE)
            deallocate(LD_EL)
            deallocate(INN)
            deallocate(IEN)

         enddo  ! for all processors in parallel code

         deallocate(x_glob_lin)
         deallocate(y_glob_lin)
         deallocate(IEN_GLOB_LIN)
         deallocate(glob_mat_dat)
!             write(*,*) " Cleared all vectors for ni = ", ni
      enddo  ! ! Loop every time step
      
      end ! database_3D.f


