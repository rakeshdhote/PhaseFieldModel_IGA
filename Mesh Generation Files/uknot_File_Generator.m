% uknot?u.?.dat file Generator

cd(cdir);

ri = 4;
rj = 2;
lprt1 = 1;
lprt0 = 0;
lprt2 = 2;
tw1 = 1.0;
twp5 = 0.5;
tw0 = 0;

for k = 1:CoresZ  % Current Cores No in the Z direction
    for j = 1:CoresY % Current Cores No in the Y direction
        for i = 1:CoresX % Current Cores No in the X direction
            
currentcore = i + CoresX*(j-1)+CoresX*CoresY*(k-1);

%%           
fileuknot = fopen(strcat('uknotfu.',int2str(currentcore),'.dat'),'a+');

if ((i == 1) && (PeriodicX == 1) && (TotCores > 1)) % periodic and MPI
    fprintf(fileuknot,'%5d%5d%5d\n',ri, rj, lprt1);
end

if ((i == 1) && (PeriodicX == 1) && (TotCores == 1)) % periodic and serial
    fprintf(fileuknot,'%5d%5d%5d\n',ri, rj, lprt1);
end

if ((i == 1) && (PeriodicX == 0) && (TotCores > 1)) % non-periodic and MPI
    fprintf(fileuknot,'%5d%5d%5d\n',ri, rj, lprt0);
end

if ((i == 1) && (PeriodicX == 0) && (TotCores == 1)) % non-periodic and serial
    fprintf(fileuknot,'%5d%5d%5d\n',ri, rj, lprt0);
end

if ((i > 1) && (i < CoresX) && (TotCores > 1) ) % periodic & non-periodic MPI
    fprintf(fileuknot,'%5d%5d%5d\n',ri, rj, lprt1);
end

if ((i == CoresX) && (PeriodicX == 1) && (TotCores > 1)) % periodic and MPI
    fprintf(fileuknot,'%5d%5d%5d\n',ri, rj, lprt1);
end

if ((i == CoresX) && (PeriodicX == 0) && (TotCores > 1)) % non-periodic and MPI
    fprintf(fileuknot,'%5d%5d%5d\n',ri, rj, lprt2);
end

fprintf(fileuknot,'%13.9f%13.9f\n',tw1, tw0);
fprintf(fileuknot,'%13.9f%13.9f\n',twp5, twp5);
fprintf(fileuknot,'%13.9f%13.9f\n',twp5, twp5);
fprintf(fileuknot,'%13.9f%13.9f\n',tw0, tw1);

if ((PeriodicX == 1) || ((i>1)&&(i<CoresX)))
    fprintf(fileuknot,'%13.9f%13.9f\n',tw1, tw0);
    fprintf(fileuknot,'%13.9f%13.9f\n',twp5, twp5);
    fprintf(fileuknot,'%13.9f%13.9f\n',twp5, twp5);
    fprintf(fileuknot,'%13.9f%13.9f\n',tw0, tw1);
end
fclose('all');
disp(strcat('Generated File:    uknotfu.',int2str(currentcore),'.dat'));

%%
fileuknot = fopen(strcat('vknotfu.',int2str(currentcore),'.dat'),'a+');

if ((j == 1) && (PeriodicY == 1) && (TotCores > 1)) % periodic and MPI
    fprintf(fileuknot,'%5d%5d%5d\n',ri, rj, lprt1);
end

if ((j == 1) && (PeriodicY == 1) && (TotCores == 1)) % periodic and serial
    fprintf(fileuknot,'%5d%5d%5d\n',ri, rj, lprt1);
end

if ((j == 1) && (PeriodicY == 0) && (TotCores > 1)) % non-periodic and MPI
    fprintf(fileuknot,'%5d%5d%5d\n',ri, rj, lprt0);
end

if ((j == 1) && (PeriodicY == 0) && (TotCores == 1)) % non-periodic and serial
    fprintf(fileuknot,'%5d%5d%5d\n',ri, rj, lprt0);
end

if ((j > 1) && (j < CoresY) && (TotCores > 1)) % periodic & non-periodic MPI
    fprintf(fileuknot,'%5d%5d%5d\n',ri, rj, lprt1);
end

if ((j == CoresY) && (PeriodicY == 1) && (TotCores > 1)) % periodic and MPI
    fprintf(fileuknot,'%5d%5d%5d\n',ri, rj, lprt1);
end

if ((j == CoresY) && (PeriodicY == 0) && (TotCores > 1)) % non-periodic and MPI
    fprintf(fileuknot,'%5d%5d%5d\n',ri, rj, lprt2);
end

fprintf(fileuknot,'%13.9f%13.9f\n',tw1, tw0);
fprintf(fileuknot,'%13.9f%13.9f\n',twp5, twp5);
fprintf(fileuknot,'%13.9f%13.9f\n',twp5, twp5);
fprintf(fileuknot,'%13.9f%13.9f\n',tw0, tw1);
if ((PeriodicY == 1) || ((j>1)&&(j<CoresY)))
    fprintf(fileuknot,'%13.9f%13.9f\n',tw1, tw0);
    fprintf(fileuknot,'%13.9f%13.9f\n',twp5, twp5);
    fprintf(fileuknot,'%13.9f%13.9f\n',twp5, twp5);
    fprintf(fileuknot,'%13.9f%13.9f\n',tw0, tw1);
end
fclose('all');
disp(strcat('Generated File:    vknotfu.',int2str(currentcore),'.dat'));

%%

if ((OCPu==2)&&(CoresX>1)) % 2D case parallel
    continue
end

fileuknot = fopen(strcat('wknotfu.',int2str(currentcore),'.dat'),'a+');

if (dimension == 3) % 3D
    
    if ((k == 1) && (PeriodicZ == 1) && (TotCores > 1)) % periodic and MPI
        fprintf(fileuknot,'%5d%5d%5d\n',ri, rj, lprt1);
    end

    if ((k == 1) && (PeriodicZ == 0) && (TotCores > 1)) % non-periodic and MPI
        fprintf(fileuknot,'%5d%5d%5d\n',ri, rj, lprt0);
    end

    if ((k > 1) && (k < CoresZ) && (TotCores > 1)) % periodic & non-periodic MPI
        fprintf(fileuknot,'%5d%5d%5d\n',ri, rj, lprt1);
    end

    if ((k == CoresZ) && (PeriodicZ == 1) && (TotCores > 1)) % periodic and MPI
        fprintf(fileuknot,'%5d%5d%5d\n',ri, rj, lprt1);
    end

    if ((k == CoresZ) && (PeriodicZ == 0) && (TotCores > 1)) % non-periodic and MPI
        fprintf(fileuknot,'%5d%5d%5d\n',ri, rj, lprt2);
    end

  fprintf(fileuknot,'%13.9f%13.9f\n',tw1, tw0);
  fprintf(fileuknot,'%13.9f%13.9f\n',twp5, twp5);
  fprintf(fileuknot,'%13.9f%13.9f\n',twp5, twp5);
  fprintf(fileuknot,'%13.9f%13.9f\n',tw0, tw1);
  if ((PeriodicZ == 1) || ((k>1)&&(k<CoresZ)))
    fprintf(fileuknot,'%13.9f%13.9f\n',tw1, tw0);
    fprintf(fileuknot,'%13.9f%13.9f\n',twp5, twp5);
    fprintf(fileuknot,'%13.9f%13.9f\n',twp5, twp5);
    fprintf(fileuknot,'%13.9f%13.9f\n',tw0, tw1);
  end
    fclose('all');
    disp(strcat('Generated File:    wknotfu.',int2str(currentcore),'.dat'));
end

if (OCPu==2) %2D
    
    fprintf(fileuknot,'%5d%5d%5d\n',rj, lprt1, lprt1);
    fprintf(fileuknot,'%13.9f\n',tw1);
    fprintf(fileuknot,'%13.9f\n',tw1);
if (PeriodicZ == 1)    
    fprintf(fileuknot,'%13.9f\n',tw1);
    fprintf(fileuknot,'%13.9f\n',tw1);
end    
    fclose('all');
    disp(strcat('Generated File:    wknotfu.',int2str(currentcore),'.dat'));
end
      end % X direction cores
    end % Y direction cores
end % Z direction cores