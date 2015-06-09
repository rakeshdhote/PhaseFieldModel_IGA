clear
clc
%% %%%%%%%%%%%%%%%%%%%%%%% 
% Details about the cluster - Just used in folder name
% Orca - MPI jobs should use multiples of 24 cores
% Saw - MPI jobs should use multiples of 8 cores
dummy00 = 'Orca_';
%dummy00 = 'Saw_';
% dummy00 = 'Scinet_';
%dummy00 = 'Goblin_';
% dummy00 = 'NAME_';
%% %%%%%%%%%%%%%%%%%%%%%%% 
% Periodic or Non-Periodic Boundary condition apply in X, Y, and Z directioni
% Periodic? = 1 (Periodic BCs); = 0 (Non-Periodic BCs)
PeriodicX = 0;
PeriodicY = 0;
PeriodicZ = 0; % For 2D the periodic z is 1
%% %%%%%%%%%%%%%%%%%%%%%%% 
% Dimensionality of the problem - 2D or 3D
dimension = 3; % 2 = 2D ;  3 = 3D
typerun = 1; % 1 =  evolution; 2 = tensile test
%% %%%%%%%%%%%%%%%%%%%%%%% 
NSD = 3 ;% DO NOT TOUCH input('Input Number of space dimension: ');
%% %%%%%%%%%%%%%%%%%%%%%%% 
% Degree of the curve
% 2D: Pu = qu = 2; Ru = 1   %% 3D: Pu = qu = Ru = 2
% DO NOT TOUCH
Pu = 2 ;% DO NOT TOUCH % input('Input degree of curve in the u direction: ');
Qu = 2 ;% DO NOT TOUCH % input('Input degree of curve in the v direction: ');
Ru = 2 ;% DO NOT TOUCH % input('Input degree of curve in the w direction: ');
if (dimension == 2) % 2D
    Ru = 1; %m
end
%% %%%%%%%%%%%%%%%%%%%%%%% 
% # of cores in X, Y, and Z direction
% 2D: CoresZ = 1   %% 3D: CoresZ != 1
% Core refers to rank of the processor in N size cluster
CoresX = 8;% input('Input number of cores in the X direction: ');
CoresY = 2;% input('Input number of cores in the Y direction: ');
CoresZ = 2;% 
if (dimension == 2) % 2D
    CoresZ = 1; %m
end
TotCores = CoresX * CoresY * CoresZ;
%% %%%%%%%%%%%%%%%%%%%%%%% 
spacex = 1; % Scaled space in x-dimension
spacey = 1; % Scaled space in y-dimension
spacez = 1; % DEFAULT =1 % Scaled space in z-dimension
%% %%%%%%%%%%%%%%%%%%%%%%% 
% length scaling constant delta 
% calculated from the model rescaling in spatial scale
delta = 0;
if (dimension == 2) % 2D
    delta = 1.8083e-9; %m
end

kGG = 3.15e-8 ; % kg value
Azero = 1.97e10; 

if (dimension == 3) % 3D
    delta = sqrt(kGG/Azero);
%     delta = 1.26451e-9; %m
end
%% %%%%%%%%%%%%%%%%%%%%%%% 
% Physical length
 physicalx = 200 * 10^-9; %m
 physicaly = (physicalx/CoresX)*CoresY; %m   
 physicalz = (physicalx/CoresX)*CoresZ; %m   
%physicalx = 100*delta; %m
%physicaly = 20*delta; %m
%physicalz = 1*delta; %m

if (dimension == 2) % 2D
    physicalz = 1; %m
end
%% %%%%%%%%%%%%%%%%%%%%%%% 
% scaled geometry
scalex = physicalx/delta;
scaley = physicaly/delta;
scalez = physicalz/delta;
if (dimension == 2) % 2D
    scalez = 1; %m
end
%% %%%%%%%%%%%%%%%%%%%%%%% 
% elements in one core in X, Y, and Z direction

% APPROXIMATELY 1 MESH IN 1 nm OR 1 MESH IN 1 SCALED GEOMETRY LENGTH
% FOR MAXIMUM ERROR TO BE LESS THAN 0.5 %
unitelemsPerCore = 20; % 16 for 3D; 32 for 2D
elemx = unitelemsPerCore*1; % Elements in the X direction in each processor
elemy = unitelemsPerCore*1; % Elements in the Y direction
elemz = unitelemsPerCore*1; % Elements in the Z direction
if (dimension == 2) %(CoresZ == 1)
    elemz = 1; % 2D
end
%% %%%%%%%%%%%%%%%%%%%%%%% 
%unitdistancePerCore = 56; % 8 for 3D; 16 for 2D
% length in terms of delta
lenx = scalex ;% unitdistancePerCore * CoresX;
leny = scaley; % unitdistancePerCore * CoresY;
lenz = scalez; % unitdistancePerCore * CoresZ;
if (dimension == 2)
    lenz = 1; % 2D
end
%% %%%%%%%%%%%%%%%%%%%%%%% 
%Initial condition: 
% Random = 1; Parabolic = 0; Exponential = 2; 
%3 - initial conditions; 4 = zero ic
% 5 - restart runs evolution
initcond = 1;

%% %%%%%%%%%%%%%%%%%%%%%%% 
% Time stepping: Adaptive = 1, fixed = 0
tstepadp = 0; % adaptive

%% %%%%%%%%%%%%%%%%%%%%%%% 
% Initial temperature 
TempIC = 240; 

%% %%%%%%%%%%%%%%%%%%%%%%%  
% Solver parameters
Delt = 1; % Delta t %%%0.5 1e-3;
Utol = 1e-3; % 1e-3 for larger time steps; 
NLtol = 1e-1; 
Nnewt = 5;  % 5 for faster
% sample frequency to save data
if(typerun == 1) 
    samplesave = 10; % evolution
end
if(typerun == 2) 
    samplesave = 10;  % loading
end

%% %%%%%%%%%%%%%%%%%%%%%%% 
% IBC  Values setup
% Following IBC values are used for 
    % Dynamic BC (called Pull): IBC == 2 (if IBC == 2, DIRBC ==1 at the
    % corresponding DOF)
    % Static BC (called Pull) : IBC == 1   
% Pullx  = 2 (pull); 1 = constraint ; 0 (no pull) 
% IBC Values::::::::: ! IBC = 2 pull; = 1 constraint
% VVVIMP ::::: ! if IBC == 2; DIRBC==1
% Setting up predefined DIR_BC values: set IBC = 5, DIR_BC = value

% pullxminus refers to the left edge/surface in the X direction
pullxminusu1 = 1; % DOF:   2D = u1      ; 3D = u1
pullxminusu2 = 1; %        2D = u2      ; 3D = u2
pullxminusu3 = 1; %        2D = v1      ; 3D = u3
pullxminusu4 = 1; %        2D = v2      ; 3D = v1
pullxminusu5 = 1; %        2D = theta   ; 3D = v2
pullxminusu6 = 1; %        2D = --      ; 3D = v3
pullxminusu7 = 0; %        2D = --      ; 3D = theta

% pullxminus refers to the right edge/surface in the X direction
pullxplusu1  = 0;  %2
pullxplusu2  = 0;  % 3
pullxplusu3  = 4; %  4 2D IBC = 1
pullxplusu4  = 0; % 12 3D IBC = 1
pullxplusu5  = 0; % 13
pullxplusu6  = 14; % 14
pullxplusu7  = 0;

pullyminusu1 = 0;
pullyminusu2 = 0;
pullyminusu3 = 0;
pullyminusu4 = 0;
pullyminusu5 = 0;
pullyminusu6 = 0;
pullyminusu7 = 0;

pullyplusu1  = 0;
pullyplusu2  = 0;
pullyplusu3  = 0;
pullyplusu4  = 0;
pullyplusu5  = 0;
pullyplusu6  = 0;
pullyplusu7  = 0;

pullzminusu1 = 0;
pullzminusu2 = 0;
pullzminusu3 = 0;
pullzminusu4 = 0;
pullzminusu5 = 0;
pullzminusu6 = 0;
pullzminusu7 = 0;

pullzplusu1  = 0;
pullzplusu2  = 0;
pullzplusu3  = 0;
pullzplusu4  = 0;
pullzplusu5  = 0;
pullzplusu6  = 0;
pullzplusu7  = 0;

%% %%%%%%%%%%%%%%%%%%%%%%% 
% DIRBC - dirichlet values
% X- side   ! put DIRBC = 1 for pulling
% % VVVIMP ::::: ! if IBC == 2; DIRBC==1

% Following DIRBC values are used for 
    % Dynamic BC (called Pull): IBC == 2 (if IBC == 2, DIRBC ==1 at the
    % corresponding DOF)
    % Static BC (called Pull) : IBC == 1   
% Pullx  = 2 (pull); 1 = constraint ; 0 (no pull) 
% IBC Values::::::::: ! IBC = 2 pull; = 1 constraint
% VVVIMP ::::: ! if IBC == 2; DIRBC==1

% Setting up predefined DIR_BC values: set IBC = 5, DIR_BC = value

dummy = 0;

xminusu1 = 0;   
xminusu2 = 0;     
xminusu3 = 0;
xminusu4 = 0;   
xminusu5 = 0;     
xminusu6 = 0;
xminusu7 = 0;   

% X+ side
xplusu1 = 0;        % =1  if ibc = 2
xplusu2 = 0;        % =1  if ibc = 3
xplusu3 = 1;        % =1 if ibc = 4      2D if IBC == 2; DIRBC==1
xplusu4 = 0;        % =1 if ibc = 12      3D if IBC == 2; DIRBC==1
xplusu5 = 0;        % =1 if ibc = 13     
xplusu6 = 1;        % =1 if ibc = 14     
xplusu7 = 0;    

% Y- side
yminusu1 = 0;   
yminusu2 = 0;   
yminusu3 = 0;
yminusu4 = 0;   
yminusu5 = 0;   
yminusu6 = 0;
yminusu7 = 0;   

% Y+ side
yplusu1 = 0;    
yplusu2 = 0;    
yplusu3 = 0;
yplusu4 = 0;    
yplusu5 = 0;    
yplusu6 = 0;
yplusu7 = 0;    

% Z- side
zminusu1 = 0;   
zminusu2 = 0;   
zminusu3 = 0;
zminusu4 = 0;   
zminusu5 = 0;   
zminusu6 = 0;
zminusu7 = 0;   

% Z+ side
zplusu1 = 0;    
zplusu2 = 0;    
zplusu3 = 0;
zplusu4 = 0;    
zplusu5 = 0;    
zplusu6 = 0;
zplusu7 = 0;    

%% %%%%%%%%%%%%%%%%%%%%%%% 
% Results saving Stress and strains
% StrStn = 0 NONE; 
%        = 1 Sxx and Exx; 
%        = 2 Sxx, Exx, Syy, Eyy, Szz, Ezz;
%        = 3 ALL
StrStn = 1; % none
%% %%%%%%%%%%%%%%%%%%%%%%% 
% Loading time fractions - Multiaxial loading

% X-direction
tlxfrac = 1; % DEFAULT = 1; end of loading as a fraction of max loadT
tsxfrac = 0; % DEFAULT = 0; start of loading as a fraction of max loadT

% Y-direction
tlyfrac = 1; % DEFAULT = 1; end of loading as a fraction of max loadT
tsyfrac = 0; % DEFAULT = 0; start of loading as a fraction of max loadT

% Z-direction
tlzfrac = 1; % DEFAULT = 1; end of loading as a fraction of max loadT
tszfrac = 0; % DEFAULT = 0; start of loading as a fraction of max loadT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
StrainReqX = 3e-2; % 3e-2
StrainReqY = 0; %3e-2; % 3e-2
StrainReqZ = 0; %3e-2; % 3e-2
%% %%%%%%%%%%%%%%%%%%%%%%% %% %%%%%%%%%%%%%%%%%%%%%%% 
tload = 1e-9;    % 1ns
tunload = 1e-9;  % 1ns
%% %%%%%%%%%%%%%%%%%%%%%%% %% %%%%%%%%%%%%%%%%%%%%%%% 
% Z-direction
twXlfrac = 1; % DEFAULT = 1; end of loading as a fraction of max loadT
twXsfrac = 0; % DEFAULT = 0; start of loading as a fraction of max loadT
twistReqX = 0; %5*pi/180; % ! 5 degrees
%% %%%%%%%%%%%%%%%%%%%%%%% %% %%%%%%%%%%%%%%%%%%%%%%% 
fx = 0; %0.001;  % positive for pulling in x+ direction
fy = 0;
fz = 0;
gth = 0;
%% %%%%%%%%%%%%%%%%%%%%%%% %% %%%%%%%%%%%%%%%%%%%%%%% 
Abulk = 19.23e10;
A0 = 1.97e10;
Ashear = 28e10; 
E0 = 0.0;
Ab = Abulk/A0;
As = Ashear/A0;
%%%%%%%%%%%%%%%%%%%%%%%%% %% %%%%%%%%%%%%%%%%%%%%%%% 
wt = 1; % weight

dummy0 = date; %'MeshFiles';
if (elemz == 1) 
    dummy1 = '_2D';
end
if (elemz > 1) 
    dummy1 = '_3D';
end

if ((PeriodicX & PeriodicY & PeriodicZ) == 1) 
    dummy2 = '_Periodic';
else
    dummy2 = '_NonPeriodic';
end

currentdir = pwd;
directory=strcat(dummy00,dummy0,dummy1,dummy2,'_MeshT_',num2str(elemx*CoresX),'x',num2str(elemy*CoresY),'x',num2str(elemz*CoresZ), ...
    '_GeomT_',num2str(lenx),'x',num2str(leny),'x',num2str(lenz), ...
    '_Cores_',num2str(CoresX),'x',num2str(CoresY),'x',num2str(CoresZ),'Temp_',num2str(TempIC)); % Name of the directory,where files to be saved
mkdir(directory);
cd(directory);
cdir = pwd;

%% %%%%%%%%%%%%%%%%%%%%%%% 

TotalElemX = elemx*CoresX;
TotalElemY = elemy*CoresY;
TotalElemZ = elemz*CoresZ;

MCPu = elemx+Pu ;% input('Input control points in the X direction: ');
NCPu = elemy+Qu ;% input('Input control points in the Y direction: ');
OCPu = elemz+Ru ;% input('Input control points in the Z direction: ');
TotCP = MCPu*NCPu*OCPu;

sxmax = spacex/CoresX; % space increment in each cores
symax = spacey/CoresY; 
szmax = spacez/CoresZ; 

neighbor = 2; % communication with the neighbouring two processors in x,y,z direction
nlworkuu_u = 5+4+2*2*NCPu*OCPu; % 5 lines for send intro + ...
% 4 lines for receive intro + NCpu*OCPu for 2 rows (C1 continuity) and 2 sides (East-West, 
% North-South, Top-Bottom)
nlworkuu_v = 5+4+2*2*OCPu;
nlworkuu_w = 5+4+2*2;
itag = 100;
iacc_send = 0; % send
iacc_receieve = 1; % receive

length_u = 2;  %length segment in X-direction; 2 for C1 continuity
length_v = 2*MCPu;
length_w = 2*MCPu*NCPu;

numsegment_u = NCPu*OCPu;
numsegment_v = 1*OCPu;
numsegment_w = 1*1;

NonPeriodic_Send_Receive = -2;
%%  Generating mesh_Files

meshu = '\mesh_File_Generator.m';
run(strcat(currentdir,meshu));

%%  Generating lwork_Files

lwork = '\lwork_File_Generator.m';
run(strcat(currentdir,lwork));

%%  Generating mesh_Files

uknot = '\uknot_File_Generator.m';
run(strcat(currentdir,uknot));

%%  Generating postp_CASE.dat file

uknot = '\mesh_File_Generator_postpCASE.m';
run(strcat(currentdir,uknot));
%% Writing specs data in the dat file
file = fopen('specs.dat','a+');
fprintf(file,'%5d%5d%5d\n',CoresX, CoresY, CoresZ);
fprintf(file,'%5d%5d%5d\n',PeriodicX, PeriodicY, PeriodicZ);
fprintf(file,'%5d\n',initcond);
fprintf(file,'%5d\n',tstepadp);
fprintf(file,'%13.9f\n',lenx);
fprintf(file,'%13.9f\n',leny);
fprintf(file,'%13.9f\n',lenz);

fprintf(file,'%13.9f\n',xminusu1);
fprintf(file,'%13.9f\n',xminusu2);
fprintf(file,'%13.9f\n',xminusu3);
fprintf(file,'%13.9f\n',xminusu4);
fprintf(file,'%13.9f\n',xminusu5);
fprintf(file,'%13.9f\n',xminusu6);
fprintf(file,'%13.9f\n',xminusu7);

fprintf(file,'%13.9f\n',xplusu1);
fprintf(file,'%13.9f\n',xplusu2);
fprintf(file,'%13.9f\n',xplusu3);
fprintf(file,'%13.9f\n',xplusu4);
fprintf(file,'%13.9f\n',xplusu5);
fprintf(file,'%13.9f\n',xplusu6);
fprintf(file,'%13.9f\n',xplusu7);

fprintf(file,'%13.9f\n',yminusu1);
fprintf(file,'%13.9f\n',yminusu2);
fprintf(file,'%13.9f\n',yminusu3);
fprintf(file,'%13.9f\n',yminusu4);
fprintf(file,'%13.9f\n',yminusu5);
fprintf(file,'%13.9f\n',yminusu6);
fprintf(file,'%13.9f\n',yminusu7);

fprintf(file,'%13.9f\n',yplusu1);
fprintf(file,'%13.9f\n',yplusu2);
fprintf(file,'%13.9f\n',yplusu3);
fprintf(file,'%13.9f\n',yplusu4);
fprintf(file,'%13.9f\n',yplusu5);
fprintf(file,'%13.9f\n',yplusu6);
fprintf(file,'%13.9f\n',yplusu7);

fprintf(file,'%13.9f\n',zminusu1);
fprintf(file,'%13.9f\n',zminusu2);
fprintf(file,'%13.9f\n',zminusu3);
fprintf(file,'%13.9f\n',zminusu4);
fprintf(file,'%13.9f\n',zminusu5);
fprintf(file,'%13.9f\n',zminusu6);
fprintf(file,'%13.9f\n',zminusu7);

fprintf(file,'%13.9f\n',zplusu1);
fprintf(file,'%13.9f\n',zplusu2);
fprintf(file,'%13.9f\n',zplusu3);
fprintf(file,'%13.9f\n',zplusu4);
fprintf(file,'%13.9f\n',zplusu5);
fprintf(file,'%13.9f\n',zplusu6);
fprintf(file,'%13.9f\n',zplusu7);

fprintf(file,'%5d\n',StrStn); 

fprintf(file,'%5d\n',pullxminusu1); 
fprintf(file,'%5d\n',pullxminusu2); 
fprintf(file,'%5d\n',pullxminusu3); 
fprintf(file,'%5d\n',pullxminusu4); 
fprintf(file,'%5d\n',pullxminusu5); 
fprintf(file,'%5d\n',pullxminusu6); 
fprintf(file,'%5d\n',pullxminusu7); 

fprintf(file,'%5d\n',pullxplusu1); 
fprintf(file,'%5d\n',pullxplusu2); 
fprintf(file,'%5d\n',pullxplusu3); 
fprintf(file,'%5d\n',pullxplusu4); 
fprintf(file,'%5d\n',pullxplusu5); 
fprintf(file,'%5d\n',pullxplusu6); 
fprintf(file,'%5d\n',pullxplusu7); 

fprintf(file,'%5d\n',pullyminusu1); 
fprintf(file,'%5d\n',pullyminusu2); 
fprintf(file,'%5d\n',pullyminusu3); 
fprintf(file,'%5d\n',pullyminusu4); 
fprintf(file,'%5d\n',pullyminusu5); 
fprintf(file,'%5d\n',pullyminusu6); 
fprintf(file,'%5d\n',pullyminusu7); 

fprintf(file,'%5d\n',pullyplusu1); 
fprintf(file,'%5d\n',pullyplusu2); 
fprintf(file,'%5d\n',pullyplusu3); 
fprintf(file,'%5d\n',pullyplusu4); 
fprintf(file,'%5d\n',pullyplusu5); 
fprintf(file,'%5d\n',pullyplusu6); 
fprintf(file,'%5d\n',pullyplusu7); 

fprintf(file,'%5d\n',pullzminusu1); 
fprintf(file,'%5d\n',pullzminusu2); 
fprintf(file,'%5d\n',pullzminusu3); 
fprintf(file,'%5d\n',pullzminusu4); 
fprintf(file,'%5d\n',pullzminusu5); 
fprintf(file,'%5d\n',pullzminusu6); 
fprintf(file,'%5d\n',pullzminusu7); 

fprintf(file,'%5d\n',pullzplusu1); 
fprintf(file,'%5d\n',pullzplusu2); 
fprintf(file,'%5d\n',pullzplusu3); 
fprintf(file,'%5d\n',pullzplusu4); 
fprintf(file,'%5d\n',pullzplusu5); 
fprintf(file,'%5d\n',pullzplusu6); 
fprintf(file,'%5d\n',pullzplusu7); 

fprintf(file,'%13.9f\n',TempIC);
fprintf(file,'%13.9f\n',Delt);
fprintf(file,'%13.9f\n',Utol);
fprintf(file,'%13.9f\n',NLtol);
fprintf(file,'%5d\n',Nnewt); 
fprintf(file,'%5d\n',samplesave); 
fprintf(file,'%13.9f\n',tlxfrac); 
fprintf(file,'%13.9f\n',tsxfrac); 
fprintf(file,'%13.9f\n',tlyfrac); 
fprintf(file,'%13.9f\n',tsyfrac); 
fprintf(file,'%13.9f\n',tlzfrac); 
fprintf(file,'%13.9f\n',tszfrac); 
fprintf(file,'%13.9f\n',StrainReqX); 
fprintf(file,'%13.9f\n',StrainReqY); 
fprintf(file,'%13.9f\n',StrainReqZ); 
fprintf(file,'%13.9f\n',tload); 
fprintf(file,'%13.9f\n',tunload); 
fprintf(file,'%13.9f\n',twXlfrac); 
fprintf(file,'%13.9f\n',twXsfrac); 
fprintf(file,'%13.9f\n',twistReqX);
fprintf(file,'%13.9f\n',fx); 
fprintf(file,'%13.9f\n',fy); 
fprintf(file,'%13.9f\n',fz); 
fprintf(file,'%13.9f\n',gth);
fclose('all');

%% Writing specs details for quick understanding of settings in the dat file
file = fopen('SimulationDetails.dat','a+');
fprintf(file,'dimension = %5d\n',dimension);
fprintf(file,'CoresX, CoresY, CoresZ = %5d%5d%5d\n',CoresX, CoresY, CoresZ);
fprintf(file,'Total Cores = %5d\n',TotCores);
fprintf(file,'PeriodicX, PeriodicY, PeriodicZ = %5d%5d%5d\n',PeriodicX, PeriodicY, PeriodicZ);
fprintf(file,'initcondition (1 = random; 0 = parabolic; 2 = others) = %5d\n',initcond);
fprintf(file,'time stepping (1 = adaptive; 0 = fixed) = %5d\n',tstepadp);
if (dimension == 2) % 2D
fprintf(file,'physicalx, physicaly, physicalz (nm) = %13.9f %13.9f %13.9f\n',physicalx/1e-9,physicaly/1e-9,1); 
end
if (dimension == 3) % 2D
fprintf(file,'physicalx, physicaly, physicalz (nm) = %13.9f %13.9f %13.9f\n',physicalx/1e-9,physicaly/1e-9,physicalz/1e-9); 
end
fprintf(file,'lengthX = %13.9f\n',lenx);
fprintf(file,'lengthY = %13.9f\n',leny);
fprintf(file,'lengthZ = %13.9f\n',lenz);
fprintf(file,'xminusu1 = %13.9f\n',xminusu1);
fprintf(file,'xminusu2 = %13.9f\n',xminusu2);
fprintf(file,'xminusu3 = %13.9f\n',xminusu3);
fprintf(file,'xminusu4 = %13.9f\n',xminusu4);
fprintf(file,'xminusu5 = %13.9f\n',xminusu5);
fprintf(file,'xminusu6 = %13.9f\n',xminusu6);
fprintf(file,'xminusu7 = %13.9f\n',xminusu7);

fprintf(file,'xplusu1 = %13.9f\n',xplusu1);
fprintf(file,'xplusu2 = %13.9f\n',xplusu2);
fprintf(file,'xplusu3 = %13.9f\n',xplusu3);
fprintf(file,'xplusu4 = %13.9f\n',xplusu4);
fprintf(file,'xplusu5 = %13.9f\n',xplusu5);
fprintf(file,'xplusu6 = %13.9f\n',xplusu6);
fprintf(file,'xplusu7 = %13.9f\n',xplusu7);

fprintf(file,'yminusu1 = %13.9f\n',yminusu1);
fprintf(file,'yminusu2 = %13.9f\n',yminusu2);
fprintf(file,'yminusu3 = %13.9f\n',yminusu3);
fprintf(file,'yminusu4 = %13.9f\n',yminusu4);
fprintf(file,'yminusu5 = %13.9f\n',yminusu5);
fprintf(file,'yminusu6 = %13.9f\n',yminusu6);
fprintf(file,'yminusu7 = %13.9f\n',yminusu7);

fprintf(file,'yplusu1 = %13.9f\n',yplusu1);
fprintf(file,'yplusu2 = %13.9f\n',yplusu2);
fprintf(file,'yplusu3 = %13.9f\n',yplusu3);
fprintf(file,'yplusu4 = %13.9f\n',yplusu4);
fprintf(file,'yplusu5 = %13.9f\n',yplusu5);
fprintf(file,'yplusu6 = %13.9f\n',yplusu6);
fprintf(file,'yplusu7 = %13.9f\n',yplusu7);

fprintf(file,'zminusu1 = %13.9f\n',zminusu1);
fprintf(file,'zminusu2 = %13.9f\n',zminusu2);
fprintf(file,'zminusu3 = %13.9f\n',zminusu3);
fprintf(file,'zminusu4 = %13.9f\n',zminusu4);
fprintf(file,'zminusu5 = %13.9f\n',zminusu5);
fprintf(file,'zminusu6 = %13.9f\n',zminusu6);
fprintf(file,'zminusu7 = %13.9f\n',zminusu7);

fprintf(file,'zplusu1 = %13.9f\n',zplusu1);
fprintf(file,'zplusu2 = %13.9f\n',zplusu2);
fprintf(file,'zplusu3 = %13.9f\n',zplusu3);
fprintf(file,'zplusu4 = %13.9f\n',zplusu4);
fprintf(file,'zplusu5 = %13.9f\n',zplusu5);
fprintf(file,'zplusu6 = %13.9f\n',zplusu6);
fprintf(file,'zplusu7 = %13.9f\n',zplusu7);

fprintf(file,'StrStn = %5d\n',StrStn); 
fprintf(file,'pullxminusu1 = %5d\n',pullxminusu1); 
fprintf(file,'pullxminusu2 = %5d\n',pullxminusu2); 
fprintf(file,'pullxminusu3 = %5d\n',pullxminusu3); 
fprintf(file,'pullxminusu4 = %5d\n',pullxminusu4); 
fprintf(file,'pullxminusu5 = %5d\n',pullxminusu5); 
fprintf(file,'pullxminusu6 = %5d\n',pullxminusu6); 
fprintf(file,'pullxminusu7 = %5d\n',pullxminusu7); 

fprintf(file,'pullxplusu1 = %5d\n',pullxplusu1); 
fprintf(file,'pullxplusu2 = %5d\n',pullxplusu2); 
fprintf(file,'pullxplusu3 = %5d\n',pullxplusu3); 
fprintf(file,'pullxplusu4 = %5d\n',pullxplusu4); 
fprintf(file,'pullxplusu5 = %5d\n',pullxplusu5); 
fprintf(file,'pullxplusu6 = %5d\n',pullxplusu6); 
fprintf(file,'pullxplusu7 = %5d\n',pullxplusu7); 

fprintf(file,'pullyminusu1 = %5d\n',pullyminusu1); 
fprintf(file,'pullyminusu2 = %5d\n',pullyminusu2); 
fprintf(file,'pullyminusu3 = %5d\n',pullyminusu3); 
fprintf(file,'pullyminusu4 = %5d\n',pullyminusu4); 
fprintf(file,'pullyminusu5 = %5d\n',pullyminusu5); 
fprintf(file,'pullyminusu6 = %5d\n',pullyminusu6); 
fprintf(file,'pullyminusu7 = %5d\n',pullyminusu7); 

fprintf(file,'pullyplusu1 = %5d\n',pullyplusu1); 
fprintf(file,'pullyplusu2 = %5d\n',pullyplusu2); 
fprintf(file,'pullyplusu3 = %5d\n',pullyplusu3); 
fprintf(file,'pullyplusu4 = %5d\n',pullyplusu4); 
fprintf(file,'pullyplusu5 = %5d\n',pullyplusu5); 
fprintf(file,'pullyplusu6 = %5d\n',pullyplusu6); 
fprintf(file,'pullyplusu7 = %5d\n',pullyplusu7); 

fprintf(file,'pullzminusu1 = %5d\n',pullzminusu1); 
fprintf(file,'pullzminusu2 = %5d\n',pullzminusu2); 
fprintf(file,'pullzminusu3 = %5d\n',pullzminusu3); 
fprintf(file,'pullzminusu4 = %5d\n',pullzminusu4); 
fprintf(file,'pullzminusu5 = %5d\n',pullzminusu5); 
fprintf(file,'pullzminusu6 = %5d\n',pullzminusu6); 
fprintf(file,'pullzminusu7 = %5d\n',pullzminusu7); 

fprintf(file,'pullzplusu1 = %5d\n',pullzplusu1); 
fprintf(file,'pullzplusu2 = %5d\n',pullzplusu2); 
fprintf(file,'pullzplusu3 = %5d\n',pullzplusu3); 
fprintf(file,'pullzplusu4 = %5d\n',pullzplusu4); 
fprintf(file,'pullzplusu5 = %5d\n',pullzplusu5); 
fprintf(file,'pullzplusu6 = %5d\n',pullzplusu6); 
fprintf(file,'pullzplusu7 = %5d\n',pullzplusu7); 

fprintf(file,'TempIC = %13.9f\n',TempIC);
fprintf(file,'Delt = %13.9f\n',Delt);
fprintf(file,'Utol = %13.9f\n',Utol);
fprintf(file,'NLtol = %13.9f\n',NLtol);
fprintf(file,'Nnewt = %5d\n',Nnewt); 
fprintf(file,'samplesave = %5d\n',samplesave);  
fprintf(file,'tlxfrac = %13.9f\n',tlxfrac);
fprintf(file,'tsxfrac = %13.9f\n',tsxfrac);
fprintf(file,'tlyfrac = %13.9f\n',tlyfrac);
fprintf(file,'tsyfrac = %13.9f\n',tsyfrac);
fprintf(file,'tlzfrac = %13.9f\n',tlzfrac);
fprintf(file,'tszfrac = %13.9f\n',tszfrac);
fprintf(file,'StrainReqX = %13.9f\n',StrainReqX);
fprintf(file,'StrainReqY = %13.9f\n',StrainReqY);
fprintf(file,'StrainReqZ = %13.9f\n',StrainReqZ);
fprintf(file,'tload      = %13.9f\n',tload);
fprintf(file,'tunload    = %13.9f\n',tunload);
fprintf(file,' twXlfrac  = %13.9f\n',twXlfrac);
fprintf(file,' twXsfrac  = %13.9f\n',twXsfrac);
fprintf(file,' twistReqX = %13.9f\n',twistReqX);
fprintf(file,'fx         = %13.9f\n',fx);
fprintf(file,' fy        = %13.9f\n',fy);
fprintf(file,' fz        = %13.9f\n',fz);
fprintf(file,' gth       = %13.9f\n',gth);
fclose('all');
%%
tinit = 0;
tfinal = 10000;
tinc = 100;

file = fopen('step.dat','a+');
fprintf(file,'%5d %5d %5d %5d\n',TotCores, tinit, tfinal, tinc);
fclose('all');

%%
file = fopen('postp.dat','a+');
fprintf(file,'%5d %5d %5d %5d\n',Pu, Qu, Ru);
fprintf(file,'%5d %5d %5d %5d\n',MCPu, NCPu, OCPu);
fclose('all');
%%

file = fopen('scalegeom.dat','a+');
fprintf(file,'%13.9f %13.9f %13.9f\n',scalex, scaley, scalez);
%fprintf(file,'%13.9f %13.9f %13.9f\n',Ab, As, E0);
fclose('all');
%%
copyfile('..\*.m','.');
%%
message = 'Check files in the Folder:   ';
currentdir = pwd;
disp(strcat(message,currentdir))
cd ..
winopen(currentdir);
% 
% message = sprintf('Total Elements in X = %d\nTotal Elements in Y = %d\nTotal Elements in Z = %d\nLength in X = %d\nLength in Y = %d\nLength in Z = %d', ...
%     TotalElemX,TotalElemY,TotalElemZ,lenx,leny,lenz)
% uiwait(msgbox(message));
