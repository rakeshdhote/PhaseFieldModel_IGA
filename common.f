      module common

      implicit none

      integer 
     &     NSD, NDOF,
     &     Pu, Qu, Ru, MCPu, NCPu, OCPu, NNODZu, NELu,
     &     Pp, Qp, Rp, MCPp, NCPp, OCPp, NNODZp, NELp,	 
     &     Kspaceu, Kspaceu_mn,
     &     Ksp_usd, Ksp_scale,
     &     CLOSED_U_flag, CLOSED_V_flag, CLOSED_W_flag,
     &     NGAUSS, NGAUSSf, NEL_NZA,
     &     NSHLu, NSHLp,
     &     NFACEu, NFACEp,
     &     isolver, Nstep, Nnewt, ifreq, iper, istep,
     &     lpu,lpv,lpw, dummyvar1,procx,procy,procz,
     &     periodicx, periodicy, periodicz, initcond, tstepadp, StrStn

      real*8
     &     DetJ, DetJb, ENRGY, ERROR,
     &     rhoinf, almi, alfi, gami, method,
     &     Dtgl, Delt, Ener, Elli, AvgMass, AvgTemp,
     &     hmax, hmin, havg, rpm,
     &     mu, E_kf, lambda, epsilon, xnu, theta, 
     &     Delt_mn, Delt_mx, Delt_fact, Delt_facti, time,
     &     Utol, NLtol, U_mx, pint_save, re_lambda_save

      integer k_f,k_f2
      real*8  E_in, absomega, q_in

      integer iseed, myid, mpi_master, numnodes

      ! Material constants
      real*8 A11, A33, A22, A24, A26, TempIC, Tm
      real*8 Rho, Eta, kg, length, height, fx, fy, kappa
      real*8 Cv, delta, lenT, htT, H, ec, A0
      real*8 Fz0, StressScale, timec, EtaM, Tc, CvC, kT, NLT, gth
      real*8 tload, tunload, tloadT, tunloadT, ttotT, tloadbyunloadT
      real*8 StrainReq, StrainT, uT, velT
      real*8 aa1, aa2, aa3, aa4, aa6
      real*8 rhoT, cvT, a2th, etaT, kappaT, kgTS, Ca2th
      real*8 aa1T, aa2T, aa3T, aa4T, aa6T, kgT
      real*8 Avgu1, Avgu2, scalev, distr ! distribution
      real*8 scalex, scaley, scalez, dd1, dd2
      real*8 xminusu1, xminusu2, xminusu3, xminusu4, xminusu5,xminusu6,
     &        xminusu7, yminusu1, yminusu2,yminusu3,yminusu4,yminusu5,
     &        yminusu6, yminusu7 
      real*8 zminusu1, zminusu2, zminusu3, zminusu4, zminusu5,zminusu6,
     &        zminusu7, xplusu1, xplusu2, xplusu3, xplusu4, xplusu5,  
     &         xplusu6, xplusu7
      real*8 yplusu1, yplusu2, yplusu3, yplusu4, yplusu5, yplusu6, 
     &         yplusu7, zplusu1, zplusu2, zplusu3, zplusu4, zplusu5,
     &         zplusu6, zplusu7
	  
      real*8 lenx, leny, lenz, VXn, VXnp1, VXnpAlpha

      integer pullxminusu1, pullxminusu2, pullxminusu3, pullxminusu4, 
     &        pullxminusu5, pullxminusu6, pullxminusu7, 
     &        pullxplusu1, pullxplusu2, pullxplusu3,  pullxplusu4,
     &        pullxplusu5, pullxplusu6,  pullxplusu7,
     &        pullyminusu1, pullyminusu2, pullyminusu3, pullyminusu4, 
     &        pullyminusu5, pullyminusu6, pullyminusu7, 
     &        pullyplusu1, pullyplusu2, pullyplusu3, pullyplusu4,
     &        pullyplusu5, pullyplusu6, pullyplusu7,
     &        pullzminusu1, pullzminusu2, pullzminusu3, pullzminusu4, 
     &        pullzminusu5, pullzminusu6, pullzminusu7, 
     &        pullzplusu1, pullzplusu2, pullzplusu3, pullzplusu4,
     &        pullzplusu5, pullzplusu6, pullzplusu7
	 
      real*8 AvgSxx, AvgSxy, AvgSyy
      real*8 AvgExx, AvgExy, AvgEyy

      integer samplesave	 

      end module common
