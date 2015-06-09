% mesh?u.?.dat file Generator

cd(cdir);

for k = 1:CoresZ  % Current Cores No in the Z direction
    for j = 1:CoresY % Current Cores No in the Y direction
        for i = 1:CoresX % Current Cores No in the X direction

cx = sxmax*(i-1); % coordinate x at the bottom left corner
cy = symax*(j-1); % coordinate y at the bottom left corner
cz = szmax*(k-1); % coordinate z at the bottom left corner

currentcore = i + CoresX*(j-1)+CoresX*CoresY*(k-1);

file = fopen(strcat('meshu.',int2str(currentcore),'.dat'),'a+');

fprintf(file,'%12d\n',NSD);
fprintf(file,'%5d%5d%5d\n',Pu, Qu, Ru);
fprintf(file,'%5d%5d%5d\n',MCPu, NCPu, OCPu);

% Knot vectors
% Uknotu = [i-1 i-1 linspace(i-1,i,MCPu-Pu+1) i i];
% Vknotu = [j-1 j-1 linspace(j-1,j,NCPu-Qu+1) j j];

Uknotu = [repmat(i-1,1,Pu) linspace(i-1,i,MCPu-Pu+1) repmat(i,1,Pu)];
Vknotu = [repmat(j-1,1,Qu) linspace(j-1,j,NCPu-Qu+1) repmat(j,1,Qu)];

% 2D code - Z direction knot vector
 if (dimension == 2) % 2D
        Wknotu = [0 0 1 1]; 
 end
    
 % 3D code - Z direction knot vector
 if OCPu > 2  % 3D code
        Wknotu = [repmat(k-1,1,Ru) linspace(k-1,k,OCPu-Ru+1) repmat(k,1,Ru)]; 
 end

fprintf(file,'%13.9f',Uknotu);
fprintf(file,'\n');
fprintf(file,'%13.9f',Vknotu);
fprintf(file,'\n');
fprintf(file,'%13.9f',Wknotu);
fprintf(file,'\n');

% Control Net
for r = 1:1:OCPu % z direction
    for q = 1:1:NCPu % y direction
        for p = 1:1:MCPu % x direction
            
     xterm1 = 1/((MCPu-Pu)*2)*(p>=2);
     xterm2 = 1/(MCPu-Pu)*(p-2)*(p>2)*(p<=(MCPu-1));
     xterm3 = (1/((MCPu-Pu)*2)+(MCPu-3)*1/(MCPu-Pu))*(p>(MCPu-1));            
            
     yterm1 = 1/((NCPu-Qu)*2)*(q>=2);
     yterm2 = 1/(NCPu-Qu)*(q-2)*(q>2)*(q<=(MCPu-1));
     yterm3 = (1/((NCPu-Qu)*2)+(NCPu-3)*1/(NCPu-Qu))*(q>(NCPu-1));         
     
     if OCPu > 2 %3D code
         zterm1 = 1/((OCPu-Ru)*2)*(r>=2);
         zterm2 = 1/(OCPu-Ru)*(r-2)*(r>2)*(r<=(OCPu-1));
         zterm3 = (1/((OCPu-Ru)*2)+(OCPu-3)*1/(OCPu-Ru))*(r>(OCPu-1));         
     end     
    
     if OCPu == 2  % 2D code
         zterm1 = r-1;
         zterm2 = 0;
         zterm3 = 0;    
     end
            
            array = [cx+(xterm1+xterm2+xterm3)*sxmax  ...
                cy+(yterm1+yterm2+yterm3)*symax ...
                cz+(zterm1+zterm2+zterm3)*szmax  wt];
            fprintf(file,'%13.9f',array);
            fprintf(file,'\n');
        end
    end
end

for IPERu = 1:1:TotCP
            fprintf(file,'%12d',IPERu);
            fprintf(file,'\n');
end

CLOSED_U_flag = 0;
CLOSED_V_flag = 0;
CLOSED_W_flag = 0;

fprintf(file,'%2d%2d%2d\n',CLOSED_U_flag, CLOSED_V_flag, CLOSED_W_flag);

% Number of global nodes
NNODZu = MCPu*NCPu*OCPu;

% Number of elements
NELu = (MCPu-Pu)*(NCPu-Qu)*(OCPu-Ru);

% Number of faces
NFACEu = 2*(MCPu-Pu)*(NCPu-Qu)*(1-CLOSED_W_flag) ...
    + 2*(NCPu-Qu)*(OCPu-Ru)*(1-CLOSED_U_flag) ...
    + 2*(MCPu-Pu)*(OCPu-Ru)*(1-CLOSED_V_flag);

IBC = 0;
for ii = 1:1:NNODZu
    fprintf(file,'%2d%2d%2d\n',IBC, IBC, IBC);
end

DIR_BC = 0.0;
for ii = 1:1:NNODZu
    fprintf(file,'%13.9f%13.9f%13.9f\n',DIR_BC, DIR_BC, DIR_BC);
end

IBC_FACE = 0;
for ii = 1:1:NFACEu
    fprintf(file,'%2d%2d%2d\n',IBC_FACE, IBC_FACE, IBC_FACE);
end

LD_FACE = 0.0;
for ii = 1:1:NFACEu
    fprintf(file,'%13.9f%13.9f%13.9f\n',LD_FACE, LD_FACE, LD_FACE);
end

LD_EL = 0.0;
for ii = 1:1:NELu
    fprintf(file,'%13.9f%13.9f%13.9f\n',LD_EL, LD_EL, LD_EL);
end

DENS = 1.0;
for ii = 1:1:NELu
    fprintf(file,'%19.14f     \n',DENS);
end

VISC = 0.006666667;
for ii = 1:1:NELu
    fprintf(file,'%13.9f\n',VISC);
end

fclose('all');
disp(strcat('Generated File:    mesh.',int2str(currentcore),'.dat'));

      end % X direction cores
    end % Y direction cores
end % Z direction cores
