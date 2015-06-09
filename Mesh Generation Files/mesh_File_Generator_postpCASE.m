% mesh?u.?.dat file Generator

cd(cdir);

%for k = 1:CoresZ  % Current Cores No in the Z direction
%    for j = 1:CoresY % Current Cores No in the Y direction
%        for i = 1:CoresX % Current Cores No in the X direction

cx = 0; % sxmax*(i-1); % coordinate x at the bottom left corner
cy = 0; % symax*(j-1); % coordinate y at the bottom left corner
cz = 0; % szmax*(k-1); % coordinate z at the bottom left corner

file = fopen('meshupostpCASE.dat','a+');

fprintf(file,'%5d%5d%5d\n',Pu, Qu, Ru);
fprintf(file,'%5d%5d%5d\n',MCPu, NCPu, OCPu);

MCPu = CoresX*(MCPu-Pu)+2;
NCPu = CoresY*(NCPu-Qu)+2;
OCPu = CoresZ*(OCPu-Ru)+2;

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

fclose('all');
disp(strcat('Generated File: meshupostpCASE.dat'));

%      end % X direction cores
%    end % Y direction cores
%end % Z direction cores
