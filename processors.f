      ! # of processors in XYZ directions
      subroutine processors_u
      
      use common
!      use mpi

      implicit none
      
      integer :: meshf, ios
      character*30 fname	  

      meshf = 11
      
      fname = 'specs.dat'

      open (meshf, file=fname, status='old', IOSTAT=ios)
      IF (ios.NE.0) PRINT*, ios

      ! # of processors in XYZ directions
      read (meshf,*)  procx, procy, procz
      read (meshf,*)  periodicx, periodicy, periodicz

      read (meshf,*)  initcond	  
      read (meshf,*)  tstepadp

      read (meshf,"(F13.9)") lenx
      read (meshf,"(F13.9)") leny
      read (meshf,"(F13.9)") lenz

      read (meshf,"(F13.9)") xminusu1
      read (meshf,"(F13.9)") xminusu2
      read (meshf,"(F13.9)") xminusu3
      read (meshf,"(F13.9)") xminusu4
      read (meshf,"(F13.9)") xminusu5
      read (meshf,"(F13.9)") xminusu6
      read (meshf,"(F13.9)") xminusu7

      read (meshf,"(F13.9)") xplusu1
      read (meshf,"(F13.9)") xplusu2
      read (meshf,"(F13.9)") xplusu3
      read (meshf,"(F13.9)") xplusu4
      read (meshf,"(F13.9)") xplusu5
      read (meshf,"(F13.9)") xplusu6
      read (meshf,"(F13.9)") xplusu7

      read (meshf,"(F13.9)") yminusu1
      read (meshf,"(F13.9)") yminusu2
      read (meshf,"(F13.9)") yminusu3
      read (meshf,"(F13.9)") yminusu4
      read (meshf,"(F13.9)") yminusu5
      read (meshf,"(F13.9)") yminusu6
      read (meshf,"(F13.9)") yminusu7

      read (meshf,"(F13.9)") yplusu1
      read (meshf,"(F13.9)") yplusu2
      read (meshf,"(F13.9)") yplusu3
      read (meshf,"(F13.9)") yplusu4
      read (meshf,"(F13.9)") yplusu5
      read (meshf,"(F13.9)") yplusu6
      read (meshf,"(F13.9)") yplusu7

      read (meshf,"(F13.9)") zminusu1
      read (meshf,"(F13.9)") zminusu2
      read (meshf,"(F13.9)") zminusu3
      read (meshf,"(F13.9)") zminusu4
      read (meshf,"(F13.9)") zminusu5
      read (meshf,"(F13.9)") zminusu6
      read (meshf,"(F13.9)") zminusu7
	  
      read (meshf,"(F13.9)") zplusu1
      read (meshf,"(F13.9)") zplusu2
      read (meshf,"(F13.9)") zplusu3
      read (meshf,"(F13.9)") zplusu4
      read (meshf,"(F13.9)") zplusu5
      read (meshf,"(F13.9)") zplusu6
      read (meshf,"(F13.9)") zplusu7

      read (meshf,*)  StrStn

      read (meshf,*)  pullxminusu1
      read (meshf,*)  pullxminusu2
      read (meshf,*)  pullxminusu3	  
      read (meshf,*)  pullxminusu4
      read (meshf,*)  pullxminusu5
      read (meshf,*)  pullxminusu6	 
      read (meshf,*)  pullxminusu7

      read (meshf,*)  pullxplusu1
      read (meshf,*)  pullxplusu2
      read (meshf,*)  pullxplusu3
      read (meshf,*)  pullxplusu4
      read (meshf,*)  pullxplusu5
      read (meshf,*)  pullxplusu6
      read (meshf,*)  pullxplusu7

      read (meshf,*)  pullyminusu1
      read (meshf,*)  pullyminusu2
      read (meshf,*)  pullyminusu3	  
      read (meshf,*)  pullyminusu4
      read (meshf,*)  pullyminusu5
      read (meshf,*)  pullyminusu6	
      read (meshf,*)  pullyminusu7

      read (meshf,*)  pullyplusu1
      read (meshf,*)  pullyplusu2
      read (meshf,*)  pullyplusu3	  
      read (meshf,*)  pullyplusu4
      read (meshf,*)  pullyplusu5
      read (meshf,*)  pullyplusu6	
      read (meshf,*)  pullyplusu7

      read (meshf,*)  pullzminusu1
      read (meshf,*)  pullzminusu2
      read (meshf,*)  pullzminusu3	  
      read (meshf,*)  pullzminusu4
      read (meshf,*)  pullzminusu5
      read (meshf,*)  pullzminusu6	
      read (meshf,*)  pullzminusu7

      read (meshf,*)  pullzplusu1
      read (meshf,*)  pullzplusu2
      read (meshf,*)  pullzplusu3	  
      read (meshf,*)  pullzplusu4
      read (meshf,*)  pullzplusu5
      read (meshf,*)  pullzplusu6	
      read (meshf,*)  pullzplusu7
	  
      read (meshf,"(F13.9)") TempIC
      read (meshf,"(F13.9)") Delt
      read (meshf,"(F13.9)") Utol	
      read (meshf,"(F13.9)") NLtol
      read (meshf,*)  Nnewt	
      read (meshf,*)  samplesave

!         if (myid.eq.mpi_master) then	       print *,'=============================================='
      print *,'Processors in X-direction           =', procx
      print *,'Processors in Y-direction           =', procy
      print *,'Processors in Z-direction           =', procz

      IF (periodicx == 1) PRINT*, 'Periodicity in X-direction: YES'
      IF (periodicx == 0) PRINT*, 'Periodicity in X-direction: NO'

      IF (periodicy == 1) PRINT*, 'Periodicity in Y-direction: YES'
      IF (periodicy == 0) PRINT*, 'Periodicity in Y-direction: NO'

      IF (periodicz == 1) PRINT*, 'Periodicity in Z-direction: YES'
      IF (periodicz == 0) PRINT*, 'Periodicity in Z-direction: NO'

      IF (initcond == 1) PRINT*, 'Initial Condition:     Random'
      IF (initcond == 0) PRINT*, 'Initial Condition:     Parabolic'
      IF (initcond == 3) PRINT*, 'Initial Condition:    Microstructure'
      IF (initcond == 2) PRINT*, 'Initial Condition:    Exponential'

      IF (tstepadp == 1) PRINT*, 'Time Stepping:     Adaptive'
      IF (tstepadp == 0) PRINT*, 'Time Stepping:     Fixed'

      print *,'X-minus - u1          =', xminusu1
      print *,'X-minus - u2          =', xminusu2
      print *,'X-minus - u3          =', xminusu3
      print *,'X-minus - u4          =', xminusu4
      print *,'X-minus - u5          =', xminusu5
      print *,'X-minus - u6          =', xminusu6
      print *,'X-minus - u7          =', xminusu7
	  
      print *,'X-plus - u1           =', xplusu1
      print *,'X-plus - u2           =', xplusu2
      print *,'X-plus - u3           =', xplusu3
      print *,'X-plus - u4           =', xplusu4
      print *,'X-plus - u5           =', xplusu5
      print *,'X-plus - u6           =', xplusu6
      print *,'X-plus - u7           =', xplusu7
	  
      print *,'Y-minus - u1          =', yminusu1
      print *,'Y-minus - u2          =', yminusu2
      print *,'Y-minus - u3          =', yminusu3
      print *,'Y-minus - u4          =', yminusu4
      print *,'Y-minus - u5          =', yminusu5
      print *,'Y-minus - u6          =', yminusu6
      print *,'Y-minus - u7          =', yminusu7
	  
      print *,'Y-plus - u1           =', yplusu1
      print *,'Y-plus - u2           =', yplusu2
      print *,'Y-plus - u3           =', yplusu3
      print *,'Y-plus - u4           =', yplusu4
      print *,'Y-plus - u5           =', yplusu5
      print *,'Y-plus - u6           =', yplusu6
      print *,'Y-plus - u7           =', yplusu7
	  
      print *,'Z-minus - u1          =', zminusu1
      print *,'Z-minus - u2          =', zminusu2
      print *,'Z-minus - u3          =', zminusu3
      print *,'Z-minus - u4          =', zminusu4
      print *,'Z-minus - u5          =', zminusu5
      print *,'Z-minus - u6          =', zminusu6	  
      print *,'Z-minus - u7          =', zminusu7
	  
      print *,'Z-plus - u1           =', zplusu1
      print *,'Z-plus - u2           =', zplusu2
      print *,'Z-plus - u3           =', zplusu3
      print *,'Z-plus - u4           =', zplusu4
      print *,'Z-plus - u5           =', zplusu5
      print *,'Z-plus - u6           =', zplusu6
      print *,'Z-plus - u7           =', zplusu7

      print *,'pullxminusu1           =', pullxminusu1
      print *,'pullxminusu2           =', pullxminusu2
      print *,'pullxminusu3           =', pullxminusu3
      print *,'pullxminusu4           =', pullxminusu4
      print *,'pullxminusu5           =', pullxminusu5
      print *,'pullxminusu6           =', pullxminusu6
      print *,'pullxminusu7           =', pullxminusu7

      print *,'pullxplusu1           =', pullxplusu1
      print *,'pullxplusu2           =', pullxplusu2
      print *,'pullxplusu3           =', pullxplusu3
      print *,'pullxplusu4           =', pullxplusu4
      print *,'pullxplusu5           =', pullxplusu5
      print *,'pullxplusu6           =', pullxplusu6
      print *,'pullxplusu7           =', pullxplusu7

      print *,'pullyminusu1           =', pullyminusu1
      print *,'pullyminusu2           =', pullyminusu2
      print *,'pullyminusu3           =', pullyminusu3
      print *,'pullyminusu4           =', pullyminusu4
      print *,'pullyminusu5           =', pullyminusu5
      print *,'pullyminusu6           =', pullyminusu6
      print *,'pullyminusu7           =', pullyminusu7

      print *,'pullyplusu1           =', pullyplusu1
      print *,'pullyplusu2           =', pullyplusu2
      print *,'pullyplusu3           =', pullyplusu3
      print *,'pullyplusu4           =', pullyplusu4
      print *,'pullyplusu5           =', pullyplusu5
      print *,'pullyplusu6           =', pullyplusu6
      print *,'pullyplusu7           =', pullyplusu7

      print *,'pullzminusu1           =', pullzminusu1
      print *,'pullzminusu2           =', pullzminusu2
      print *,'pullzminusu3           =', pullzminusu3
      print *,'pullzminusu4           =', pullzminusu4
      print *,'pullzminusu5           =', pullzminusu5
      print *,'pullzminusu6           =', pullzminusu6
      print *,'pullzminusu7           =', pullzminusu7

      print *,'pullzplusu1           =', pullzplusu1
      print *,'pullzplusu2           =', pullzplusu2
      print *,'pullzplusu3           =', pullzplusu3
      print *,'pullzplusu4           =', pullzplusu4
      print *,'pullzplusu5           =', pullzplusu5
      print *,'pullzplusu6           =', pullzplusu6
      print *,'pullzplusu7           =', pullzplusu7

      IF (pullxminusu1 == 1) PRINT*, 'Time dependent Dirichlet loading
     &  X-minus boundary u1 direction: YES'
      IF (pullxminusu1 == 0) PRINT*, 'Time dependent Dirichlet loading
     &  X-minus boundary u1 direction: No'	 

      IF (pullxminusu2 == 1) PRINT*, 'Time dependent Dirichlet loading
     &  X-minus boundary u2 direction: YES'
      IF (pullxminusu2 == 0) PRINT*, 'Time dependent Dirichlet loading
     &  X-minus boundary u2 direction: No'	 
	 
      IF (pullxminusu3 == 1) PRINT*, 'Time dependent Dirichlet loading
     &  X-minus boundary u3 direction: YES'
      IF (pullxminusu3 == 0) PRINT*, 'Time dependent Dirichlet loading
     &  X-minus boundary u3 direction: No'	 	 
	 
      IF (pullxminusu4 == 1) PRINT*, 'Time dependent Dirichlet loading
     &  X-minus boundary u4 direction: YES'
      IF (pullxminusu4 == 0) PRINT*, 'Time dependent Dirichlet loading
     &  X-minus boundary u4 direction: No'	

      IF (pullxminusu5 == 1) PRINT*, 'Time dependent Dirichlet loading
     &  X-minus boundary u5 direction: YES'
      IF (pullxminusu5 == 0) PRINT*, 'Time dependent Dirichlet loading
     &  X-minus boundary u5 direction: No'	

      IF (pullxminusu6 == 1) PRINT*, 'Time dependent Dirichlet loading
     &  X-minus boundary u6 direction: YES'
      IF (pullxminusu6 == 0) PRINT*, 'Time dependent Dirichlet loading
     &  X-minus boundary u6 direction: No'	

      IF (pullxminusu7 == 1) PRINT*, 'Time dependent Dirichlet loading
     &  X-minus boundary u7 direction: YES'
      IF (pullxminusu7 == 0) PRINT*, 'Time dependent Dirichlet loading
     &  X-minus boundary u7 direction: No'	


	 
      IF (pullxplusu1 == 1) PRINT*, 'Time dependent Dirichlet loading
     &  X-plus boundary u1 direction: YES'
      IF (pullxplusu1 == 0) PRINT*, 'Time dependent Dirichlet loading
     &  X-plus boundary u1 direction: NO'
      IF (pullxplusu1 == 2) PRINT*, 'Time dependent Dirichlet loading
     &  X-plus boundary u1 direction: YES & Pull'

      IF (pullxplusu2 == 1) PRINT*, 'Time dependent Dirichlet loading
     &  X-plus boundary u2 direction: YES'
      IF (pullxplusu2 == 0) PRINT*, 'Time dependent Dirichlet loading
     &  X-plus boundary u2 direction: NO'

      IF (pullxplusu3 == 1) PRINT*, 'Time dependent Dirichlet loading
     &  X-plus boundary u3 direction: YES'
      IF (pullxplusu3 == 0) PRINT*, 'Time dependent Dirichlet loading
     &  X-plus boundary u3 direction: NO'	 
	 
      IF (pullxplusu4 == 1) PRINT*, 'Time dependent Dirichlet loading
     &  X-plus boundary u4 direction: YES'
      IF (pullxplusu4 == 0) PRINT*, 'Time dependent Dirichlet loading
     &  X-plus boundary u4 direction: NO'

      IF (pullxplusu5 == 1) PRINT*, 'Time dependent Dirichlet loading
     &  X-plus boundary u5 direction: YES'
      IF (pullxplusu5 == 0) PRINT*, 'Time dependent Dirichlet loading
     &  X-plus boundary u5 direction: NO'	

      IF (pullxplusu6 == 1) PRINT*, 'Time dependent Dirichlet loading
     &  X-plus boundary u6 direction: YES'
      IF (pullxplusu6 == 0) PRINT*, 'Time dependent Dirichlet loading
     &  X-plus boundary u6 direction: NO'

      IF (pullxplusu7 == 1) PRINT*, 'Time dependent Dirichlet loading
     &  X-plus boundary u7 direction: YES'
      IF (pullxplusu7 == 0) PRINT*, 'Time dependent Dirichlet loading
     &  X-plus boundary u7 direction: NO'	
	 
	 
	 
      IF (pullyminusu1 == 1) PRINT*, 'Time dependent Dirichlet loading
     &  Y-minus boundary u1 direction: YES'
      IF (pullyminusu1 == 0) PRINT*, 'Time dependent Dirichlet loading
     &  Y-minus boundary u1 direction: NO'

      IF (pullyminusu2 == 1) PRINT*, 'Time dependent Dirichlet loading
     &  Y-minus boundary u2 direction: YES'
      IF (pullyminusu2 == 0) PRINT*, 'Time dependent Dirichlet loading
     &  Y-minus boundary u2 direction: NO'

      IF (pullyminusu3 == 1) PRINT*, 'Time dependent Dirichlet loading
     &  Y-minus boundary u3 direction: YES'
      IF (pullyminusu3 == 0) PRINT*, 'Time dependent Dirichlet loading
     &  Y-minus boundary u3 direction: NO'	 

      IF (pullyminusu4 == 1) PRINT*, 'Time dependent Dirichlet loading
     &  Y-minus boundary u4 direction: YES'
      IF (pullyminusu4 == 0) PRINT*, 'Time dependent Dirichlet loading
     &  Y-minus boundary u4 direction: NO'

      IF (pullyminusu5 == 1) PRINT*, 'Time dependent Dirichlet loading
     &  Y-minus boundary u5 direction: YES'
      IF (pullyminusu5 == 0) PRINT*, 'Time dependent Dirichlet loading
     &  Y-minus boundary u5 direction: NO'	 

      IF (pullyminusu6 == 1) PRINT*, 'Time dependent Dirichlet loading
     &  Y-minus boundary u6 direction: YES'
      IF (pullyminusu6 == 0) PRINT*, 'Time dependent Dirichlet loading
     &  Y-minus boundary u6 direction: NO'

      IF (pullyminusu7 == 1) PRINT*, 'Time dependent Dirichlet loading
     &  Y-minus boundary u7 direction: YES'
      IF (pullyminusu7 == 0) PRINT*, 'Time dependent Dirichlet loading
     &  Y-minus boundary u7 direction: NO'	 

	 
	 
      IF (pullyplusu1 == 1) PRINT*, 'Time dependent Dirichlet loading
     &  Y-plus boundary u1 direction: YES'
      IF (pullyplusu1 == 0) PRINT*, 'Time dependent Dirichlet loading
     &  Y-plus boundary u1 direction: NO'

      IF (pullyplusu2 == 1) PRINT*, 'Time dependent Dirichlet loading
     &  Y-plus boundary u2 direction: YES'
      IF (pullyplusu2 == 0) PRINT*, 'Time dependent Dirichlet loading
     &  Y-plus boundary u2 direction: NO'

      IF (pullyplusu3 == 1) PRINT*, 'Time dependent Dirichlet loading
     &  Y-plus boundary u3 direction: YES'
      IF (pullyplusu3 == 0) PRINT*, 'Time dependent Dirichlet loading
     &  Y-plus boundary u3 direction: NO'

      IF (pullyplusu4 == 1) PRINT*, 'Time dependent Dirichlet loading
     &  Y-plus boundary u4 direction: YES'
      IF (pullyplusu4 == 0) PRINT*, 'Time dependent Dirichlet loading
     &  Y-plus boundary u4 direction: NO'

      IF (pullyplusu5 == 1) PRINT*, 'Time dependent Dirichlet loading
     &  Y-plus boundary u5 direction: YES'
      IF (pullyplusu5 == 0) PRINT*, 'Time dependent Dirichlet loading
     &  Y-plus boundary u5 direction: NO'

      IF (pullyplusu6 == 1) PRINT*, 'Time dependent Dirichlet loading
     &  Y-plus boundary u6 direction: YES'
      IF (pullyplusu6 == 0) PRINT*, 'Time dependent Dirichlet loading
     &  Y-plus boundary u6 direction: NO'

      IF (pullyplusu7 == 1) PRINT*, 'Time dependent Dirichlet loading
     &  Y-plus boundary u7 direction: YES'
      IF (pullyplusu7 == 0) PRINT*, 'Time dependent Dirichlet loading
     &  Y-plus boundary u7 direction: NO'	 
	 

      IF (pullzminusu1 == 1) PRINT*, 'Time dependent Dirichlet loading
     &  Z-minus boundary u1 direction: YES'
      IF (pullzminusu1 == 0) PRINT*, 'Time dependent Dirichlet loading
     &  Z-minus boundary u1 direction: NO'

      IF (pullzminusu2 == 1) PRINT*, 'Time dependent Dirichlet loading
     &  Z-minus boundary u2 direction: YES'
      IF (pullzminusu2 == 0) PRINT*, 'Time dependent Dirichlet loading
     &  Z-minus boundary u2 direction: NO'

      IF (pullzminusu3 == 1) PRINT*, 'Time dependent Dirichlet loading
     &  Z-minus boundary u3 direction: YES'
      IF (pullzminusu3 == 0) PRINT*, 'Time dependent Dirichlet loading
     &  Z-minus boundary u3 direction: NO'	 

      IF (pullzminusu4 == 1) PRINT*, 'Time dependent Dirichlet loading
     &  Z-minus boundary u4 direction: YES'
      IF (pullzminusu4 == 0) PRINT*, 'Time dependent Dirichlet loading
     &  Z-minus boundary u4 direction: NO'

      IF (pullzminusu5 == 1) PRINT*, 'Time dependent Dirichlet loading
     &  Z-minus boundary u5 direction: YES'
      IF (pullzminusu5 == 0) PRINT*, 'Time dependent Dirichlet loading
     &  Z-minus boundary u5 direction: NO'	 

      IF (pullzminusu6 == 1) PRINT*, 'Time dependent Dirichlet loading
     &  Z-minus boundary u6 direction: YES'
      IF (pullzminusu6 == 0) PRINT*, 'Time dependent Dirichlet loading
     &  Z-minus boundary u6 direction: NO'

      IF (pullzminusu7 == 1) PRINT*, 'Time dependent Dirichlet loading
     &  Z-minus boundary u7 direction: YES'
      IF (pullzminusu7 == 0) PRINT*, 'Time dependent Dirichlet loading
     &  Z-minus boundary u7 direction: NO'	 	 

	 
      IF (pullzplusu1 == 1) PRINT*, 'Time dependent Dirichlet loading
     &  Z-plus boundary u1 direction: YES'
      IF (pullzplusu1 == 0) PRINT*, 'Time dependent Dirichlet loading
     &  Z-plus boundary u1 direction: NO'

      IF (pullzplusu2 == 1) PRINT*, 'Time dependent Dirichlet loading
     &  Z-plus boundary u2 direction: YES'
      IF (pullzplusu2 == 0) PRINT*, 'Time dependent Dirichlet loading
     &  Z-plus boundary u2 direction: NO'

      IF (pullzplusu3 == 1) PRINT*, 'Time dependent Dirichlet loading
     &  Z-plus boundary u3 direction: YES'
      IF (pullzplusu3 == 0) PRINT*, 'Time dependent Dirichlet loading
     &  Z-plus boundary u3 direction: NO'	 
	 
      IF (pullzplusu4 == 1) PRINT*, 'Time dependent Dirichlet loading
     &  Z-plus boundary u4 direction: YES'
      IF (pullzplusu4 == 0) PRINT*, 'Time dependent Dirichlet loading
     &  Z-plus boundary u4 direction: NO'

      IF (pullzplusu5 == 1) PRINT*, 'Time dependent Dirichlet loading
     &  Z-plus boundary u5 direction: YES'
      IF (pullzplusu5 == 0) PRINT*, 'Time dependent Dirichlet loading
     &  Z-plus boundary u5 direction: NO'	 

      IF (pullzplusu6 == 1) PRINT*, 'Time dependent Dirichlet loading
     &  Z-plus boundary u6 direction: YES'
      IF (pullzplusu6 == 0) PRINT*, 'Time dependent Dirichlet loading
     &  Z-plus boundary u6 direction: NO'

      IF (pullzplusu7 == 1) PRINT*, 'Time dependent Dirichlet loading
     &  Z-plus boundary u7 direction: YES'
      IF (pullzplusu7 == 0) PRINT*, 'Time dependent Dirichlet loading
     &  Z-plus boundary u7 direction: NO'	

      print *,'=============================================='
!         endif	
	  
      return

      end
