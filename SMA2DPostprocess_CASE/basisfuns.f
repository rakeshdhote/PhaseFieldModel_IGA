c
c  Subroutine basisfuns.f consumes a knot index, parameter value, and a knot
c     vector and returns an vector containing all nonzero 1D b-spline shape
c     functions evaluated at that parameter. The last entry in the vector is
c     the shape function associated with the given knot index (i.e. the first
c     local node).
c
c
c Algorithm from Piegl, Les. "The NURBS Book". Springer-Verlag: 
c    Berlin 1995; p. 70.
c
c
c  modified June 30, 2003
c
c
c  J. Austin Cottrell
c  CAM Graduate Student
c  Institute for Computational Engineering Science
c  The University of Texas at Austin


      subroutine basisfuns(i,pl,mcpl,u,u_knotl,N)


c --------------VARIABLE DECLARATIONS--------------------------------
c...  knot span index, degree of curve, number of control points
      integer i,pl,mcpl

c...  parameter value, knot vector, vector of evaluated basis functions
      real*8 u, u_knotl(mcpl+pl+1),N(pl+1),temp

c...  local variables
c
c     counters for loops
      integer j, r

c     temporary storage
      real*8 saved,leftl(pl+1),rightl(pl+1)
c -------------------------------------------------------------------

c     form vector
      
      N(1) = 1
      do j = 1,pl
         leftl(j+1) = u - u_knotl(i+1-j)
         rightl(j+1) = u_knotl(i+j) - u
         saved = 0
         do r = 0,(j-1)
            temp = N(r+1)/(rightl(r+2)+leftl(j-r+1))
            N(r+1) = saved + rightl(r+2)*temp
            saved = leftl(j-r+1)*temp
         enddo
         N(j+1) = saved
      enddo

      return 
      end
  
