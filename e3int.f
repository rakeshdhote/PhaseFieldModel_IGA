      subroutine e3int (
     &     ul, acl, xl, chil,                    
     &     shgu, shgradgu, shhessgu,
     &     ui, aci, xi, chii, grade2energy,
     &     duidxi, duidxixj, rLi)
      
      use aAdjKeep
      
      use common

      implicit none

      integer i, j, k 
!     &        ik, jk, kk
      real*8 
     &     shgu(NSHLu), 
     &     shgradgu(NSHLu,NSD), 
     &     shhessgu(NSHLu,NSD,NSD)
      real*8
     &     xl(NSHLu,NSD),
     &     ul(NSHLu,NDOF),
     &     acl(NSHLu,NDOF),
     &     chil(NSHLu,NSD)
      real*8 
     &     xi(NSD),
     &     ui(NDOF),
     &     aci(NDOF),
     &     chii(NSD),
     &     duidxi(NDOF,NSD),
     &     duidxixj(NDOF,NSD,NSD),
     &     rLi(NDOF)

      real*8 grade2energy, de2dx, de2dy

!      real*8 rho, u, v, w

      real*8 u1, u2, v1, v2, th

      ui       = 0d0
      aci      = 0d0
      xi       = 0d0
      chii     = 0d0
      duidxi   = 0d0
      duidxixj = 0d0

!.....Get quantities and their grads and hessians at Gauss points
      ! Get quantities and their grads and hessians

      do i = 1, NSHLu
         ui   = ui   + ul  (i,:)*shgu(i)
         aci  = aci  + acl (i,:)*shgu(i)
         xi   = xi   + xl  (i,:)*shgu(i)
         chii = chii + chil(i,:)*shgu(i)
         do j = 1, NSD
            duidxi(:,j) = duidxi(:,j)+ul(i,:)*shgradgu(i,j)
            do k = 1, NSD
               duidxixj(:,j,k) = duidxixj(:,j,k)+ul(i,:)*shhessgu(i,j,k)
            enddo ! k
         enddo ! j
      enddo !  i

!     Variables
      u1 = ui(1)
      u2 = ui(2)
      v1 = ui(3)
      v2 = ui(4)
      th = ui(5)

c      laprho  = duidxixj(1,1,1) + duidxixj(1,2,2) + duidxixj(1,3,3)
c      gradnor = duidxi(1,1)**2  + duidxi(1,2)**2  + duidxi(1,3)**2
c      divu    = duidxi(2,1)     + duidxi(3,2)     + duidxi(4,3)

!      laprho = 0d0
!      divu = 0d0

      de2dx = (duidxixj(1,1,1)-duidxixj(2,1,2))/sqrt(2d0)
      de2dy = (duidxixj(1,1,2)-duidxixj(2,2,2))/sqrt(2d0)

      grade2energy = de2dx*de2dx+de2dy*de2dy

      rLi = 0d0                 ! PDE Residual

      return
      end

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

c$$$      subroutine e3bint (ul, fl, pl, gl, shgu, shgradgu, shgp, ui,
c$$$     &     pri, fi, gi, duidxi)
c$$$      
c$$$      
c$$$      use aAdjKeep
c$$$      
c$$$      use common
c$$$      implicit none
c$$$
c$$$      integer i,j,k
c$$$      real*8 shgu(NSHLu), shgp(NSHLp),
c$$$     &     shgradgu(NSHLu,NSD)
c$$$      real*8 ul(NSHLu,NSD), pl(NSHLp),
c$$$     &     fl(NSHLu,NSD), gl(NSHLu,NSD)      ! Y la aceleracion?
c$$$      real*8 ui(NSD), pri,
c$$$     &     duidxi(NSD,NSD), fi(NSD), gi(NSD)
c$$$      
c$$$      
c$$$      ui  = 0d0
c$$$      pri = 0d0
c$$$      fi  = 0d0
c$$$      gi  = 0d0
c$$$      duidxi = 0d0
c$$$                             ! Get quantities and their grads and hessians
c$$$      
c$$$      do i = 1, NSHLu
c$$$         ui(:) = ui(:) + ul(i,:)*shgu(i)
c$$$         fi(:) = fi(:) + fl(i,:)*shgu(i)
c$$$         gi(:) = gi(:) + gl(i,:)*shgu(i)
c$$$         do j = 1, NSD
c$$$            duidxi(:,j) = duidxi(:,j) + ul(i,:)*shgradgu(i,j)
c$$$         enddo
c$$$      enddo
c$$$      
c$$$      do i = 1, NSHLp
c$$$         pri = pri + shgp(i)*pl(i)
c$$$      enddo
c$$$      
c$$$      return
c$$$      end
