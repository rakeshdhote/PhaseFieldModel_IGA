
      subroutine results_3D (iproc,istep)
      
      use aAdjKeep
      
      include "common.h"
      
      integer i, j, k,kk, g, sol_file
      integer iproc, istep
      character*30 fname
      character*10 cname

      integer*4 i1, i2
      real*8    x1, x2
      
      sol_file = 12+iproc
      fname = trim('restart' // cname(istep)) // cname (iproc)
      open (sol_file, file=fname, status='old')
      read(sol_file,*) i1, x1, i2, x2 !  StepCount,Delt,istep,time


      g = 0
      do k = 1, OCP             ! Read Velocity Data
        do j = 1, NCP
          do i = 1, MCP
            g = g+1
            
c*            do kk = 1,NSD
c*              if (kk < NSD) then
c*                read (sol_file,"(F13.9)", ADVANCE = "NO") 
c*     &            yg(g,kk)
c*              else
c*                read (sol_file,"(F13.9)", ADVANCE = "YES") 
c*     &            yg(g,kk)
c*              endif
c*              
c*            enddo
            
           read(sol_file,*) yg(g,1), yg(g,2), yg(g,3), yg(g,4), yg(g,5)
            
          enddo         
        enddo
      enddo
      
c$$$      g = 0
c$$$      do k = 1, OCP             ! Read Pressure Data
c$$$        do j = 1, NCP
c$$$          do i = 1, MCP
c$$$            g = g+1
c$$$            
c$$$c*            do kk = 1,NSD
c$$$c*              if (kk < NSD) then
c$$$c*                read (sol_file,"(F13.9)", ADVANCE = "NO") 
c$$$c*     &            yg(g,NSD+1)
c$$$c*              else
c$$$c*                read (sol_file,"(F13.9)", ADVANCE = "YES") 
c$$$c*     &            yg(g,NSD+1)
c$$$c*              endif
c$$$c*              
c$$$c*            enddo
c$$$
c$$$            read(sol_file,*) yg(g,4)
c$$$            
c$$$          enddo         
c$$$        enddo
c$$$      enddo

      close(sol_file)
      
      return
      end
      


