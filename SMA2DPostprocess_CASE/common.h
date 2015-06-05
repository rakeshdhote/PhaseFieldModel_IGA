
       integer NSD, P, Q, R, MCP, NCP, OCP, NNODZ, NEL, 
     &             CLOSED_U_flag, CLOSED_V_flag, NGAUSS, NSHL,
     &             NEL_NZA, CLOSED_W_flag, NFACE, isolver,
     &             nnodz_glob_lin, NEL_GLOB_LIN
	 
	 real*8 DetJ, DetJb, ENRGY, ERROR, scalex, scaley, scalez
       integer Pu, Qu, Ru, MCPu, NCPu, OCPu
	 
       common /comparint/ NSD, P, Q, R, MCP, NCP, OCP, NNODZ, NEL,
     &             CLOSED_U_flag, CLOSED_V_flag, NGAUSS, NSHL,
     &             NEL_NZA, CLOSED_W_flag, NFACE, isolver,
     &             nnodz_glob_lin, NEL_GLOB_LIN

	 common /comparrel/ DetJ, DetJb, ENRGY, ERROR
       
