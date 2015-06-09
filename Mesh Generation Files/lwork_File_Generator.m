% lknot?u.?.dat file Generator
% http://www.mpitutorial.com/mpi-send-and-receive/
% MPI_Send(void* data, int count, MPI_Datatype datatype, int destination, 
% int tag, MPI_Comm communicator)
% MPI_Recv(void* data, int count, MPI_Datatype datatype, int destination, 
% int tag, MPI_Comm communicator, MPI_Status* status)

cd(cdir);

for k = 1:CoresZ  % Current Cores No in the Z direction
    for j = 1:CoresY % Current Cores No in the Y direction
        for i = 1:CoresX % Current Cores No in the X direction

currentcore = i + CoresX*(j-1)+CoresX*CoresY*(k-1);
Start_CoreX = 1 + CoresX*(j-1)+ CoresX*CoresY*(k-1);
Start_CoreY = i + CoresX*CoresY*(k-1);
Start_CoreZ = i + CoresX*(j-1);  
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Core Communication in X direction

% Communicaiton - Send
if (i<CoresX) 
    neighbor_sendto_u = currentcore + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PeriodicX 
if (PeriodicX==1)
    if (i==CoresX)                         
        neighbor_sendto_u = Start_CoreX;
    end
end
% Communicaiton - Receive
if (PeriodicX==1)
    if (i==1) 
        neighbor_recvfrom_u = currentcore + (CoresX-1);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Non-PeriodicX

if (PeriodicX==0)
    if (i==CoresX)                         
        neighbor_sendto_u = NonPeriodic_Send_Receive;
    end
end
% Communicaiton - Receive
if (PeriodicX==0)
    if (i==1) 
        neighbor_recvfrom_u = NonPeriodic_Send_Receive;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (i>1) 
    neighbor_recvfrom_u = currentcore - 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file = fopen(strcat('lworkuu.',int2str(currentcore),'.dat'),'a+');
fprintf(file,'%12d%12d\n',neighbor,nlworkuu_u);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% odd core  communication
if mod(i,2)==1  % mod(number,2) = 0 (even); 1(odd)
        
    fprintf(file,'%12d\n',itag);
    fprintf(file,'%12d\n',iacc_send);
    
    if ((OCPu==2)&&(CoresX==1)) % 2D case serial
        if (PeriodicX == 1)
            fprintf(file,'%12d\n',2);  
        end
        if (PeriodicX == 0)
            fprintf(file,'%12d\n',neighbor_sendto_u);  
        end        
    end    
    
    if ((OCPu>2)||(CoresX>1)) % 2D parallel and 3D Case
        fprintf(file,'%12d\n',neighbor_sendto_u);
    end
    
    fprintf(file,'%12d\n',numsegment_u);
    
    for ii = 1:1:numsegment_u
        fprintf(file,'%12d\n',(MCPu-1)+(ii-1)*MCPu);
        fprintf(file,'%12d\n',length_u);
    end

    fprintf(file,'%12d\n',itag);
    fprintf(file,'%12d\n',iacc_receieve);
    fprintf(file,'%12d\n',neighbor_recvfrom_u);
    fprintf(file,'%12d\n',numsegment_u);
    
    for ii = 1:1:numsegment_u
        fprintf(file,'%12d\n',1+(ii-1)*MCPu);
        fprintf(file,'%12d\n',length_u);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% even core communication
if mod(i,2)==0
    
    fprintf(file,'%12d\n',itag);
    fprintf(file,'%12d\n',iacc_receieve);
    fprintf(file,'%12d\n',neighbor_recvfrom_u);
    fprintf(file,'%12d\n',numsegment_u);
    
    for ii = 1:1:numsegment_u
        fprintf(file,'%12d\n',1+(ii-1)*MCPu);
        fprintf(file,'%12d\n',length_u);
    end
    
    fprintf(file,'%12d\n',itag);
    if ( (PeriodicX == 0) && (i == CoresX)) % to account non-periodic bc
        fprintf(file,'%12d\n',1);
    else
        fprintf(file,'%12d\n',iacc_send);
    end
    fprintf(file,'%12d\n',neighbor_sendto_u);
    fprintf(file,'%12d\n',numsegment_u);
    
    for ii = 1:1:numsegment_u
        fprintf(file,'%12d\n',(MCPu-1)+(ii-1)*MCPu);
        fprintf(file,'%12d\n',length_u);
    end
end    

fclose('all');
disp(strcat('Generated File:    lworkuu.',int2str(currentcore),'.dat'));
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Core Communication in Y direction
% Communicaiton - Send
if (j<CoresY) 
    neighbor_sendto_v = currentcore + CoresX;
end

if (PeriodicY==1)
    if (j==CoresY)                         
        neighbor_sendto_v = Start_CoreY;
    end
end
% Communicaiton - Receive
if (PeriodicY==1)
    if (j==1) 
        neighbor_recvfrom_v = currentcore + (CoresY-1)*CoresX;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Non-Periodic

if (PeriodicY==0)
    if (j==CoresY)                         
        neighbor_sendto_v = NonPeriodic_Send_Receive;
    end
end
% Communicaiton - Receive
if (PeriodicY==0)
    if (j==1) 
        neighbor_recvfrom_v = NonPeriodic_Send_Receive;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (j>1) 
    neighbor_recvfrom_v = currentcore - CoresX;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file = fopen(strcat('lworkvu.',int2str(currentcore),'.dat'),'a+');
fprintf(file,'%12d%12d\n',neighbor,nlworkuu_v);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% odd core  communication
if mod(j,2)==1  % mod(number,2) = 0 (even); 1(odd)
        
    fprintf(file,'%12d\n',itag);
    fprintf(file,'%12d\n',iacc_send);
    
    if ((OCPu==2)&&(CoresX==1)) % 2D case serial
        if (PeriodicY == 1)
            fprintf(file,'%12d\n',2);  
        end
        if (PeriodicY == 0)
            fprintf(file,'%12d\n',neighbor_sendto_v);  
        end   
    end    
    
    if ((OCPu>2)||(CoresX>1)) % 2D parallel and 3D Case
        fprintf(file,'%12d\n',neighbor_sendto_v);
    end
    
    fprintf(file,'%12d\n',numsegment_v);
    
    for ii = 1:1:numsegment_v
        fprintf(file,'%12d\n',(MCPu*(NCPu-2)+1)+(ii-1)*MCPu*NCPu);
        fprintf(file,'%12d\n',length_v);
    end

    fprintf(file,'%12d\n',itag);
    fprintf(file,'%12d\n',iacc_receieve);
    fprintf(file,'%12d\n',neighbor_recvfrom_v);
    fprintf(file,'%12d\n',numsegment_v);
    
    for ii = 1:1:numsegment_v
        fprintf(file,'%12d\n',1+(ii-1)*MCPu*NCPu);
        fprintf(file,'%12d\n',length_v);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% even core communication
if mod(j,2)==0
    
    fprintf(file,'%12d\n',itag);
    fprintf(file,'%12d\n',iacc_receieve);
    fprintf(file,'%12d\n',neighbor_recvfrom_v);
    fprintf(file,'%12d\n',numsegment_v);
    
    for ii = 1:1:numsegment_v
        fprintf(file,'%12d\n',1+(ii-1)*MCPu*NCPu);
        fprintf(file,'%12d\n',length_v);
    end
    
    fprintf(file,'%12d\n',itag);
    if ( (PeriodicY ==0) && (j == CoresY)) 
        fprintf(file,'%12d\n',1);
    else
        fprintf(file,'%12d\n',iacc_send);
    end
    fprintf(file,'%12d\n',neighbor_sendto_v);
    fprintf(file,'%12d\n',numsegment_v);
    
    for ii = 1:1:numsegment_v
        fprintf(file,'%12d\n',(MCPu*(NCPu-2)+1)+(ii-1)*MCPu*NCPu);
        fprintf(file,'%12d\n',length_v);
    end
end    
fclose('all');
disp(strcat('Generated File:    lworkvu.',int2str(currentcore),'.dat'));

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Core Communication in Z direction
% Communicaiton - Send
if (k<CoresZ) 
    neighbor_sendto_w = currentcore + CoresX*CoresY;
end

if (PeriodicZ==1)
    if (k==CoresZ) 
        neighbor_sendto_w = Start_CoreZ;
    end
end

% Communicaiton - Receive
if (PeriodicZ==1)
    if (k==1) 
        neighbor_recvfrom_w = currentcore + (CoresZ-1)*CoresX*CoresY;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Non-PeriodicZ

if (PeriodicZ==0)
    if (k==CoresZ)                         
        neighbor_sendto_w = NonPeriodic_Send_Receive;
    end
end
% Communicaiton - Receive
if (PeriodicZ==0)
    if (k==1) 
        neighbor_recvfrom_w = NonPeriodic_Send_Receive;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (k>1) 
    neighbor_recvfrom_w = currentcore - CoresX*CoresY;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ((OCPu==2)&&(CoresX>1)) % 2D case parallel
    continue
end

file = fopen(strcat('lworkwu.',int2str(currentcore),'.dat'),'a+');
fprintf(file,'%12d%12d\n',neighbor,nlworkuu_w);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% odd core  communication
if mod(k,2)==1  % mod(number,2) = 0 (even); 1(odd)
        
    fprintf(file,'%12d\n',itag);
    fprintf(file,'%12d\n',iacc_send);
    
    if ((OCPu==2)&&(CoresX==1)) % 2D case serial
        if (PeriodicZ == 1)
            fprintf(file,'%12d\n',2);  
        end
        if (PeriodicZ == 0)
            fprintf(file,'%12d\n',neighbor_sendto_w);  
        end   
    end    
    
    if ((OCPu>2)||(CoresX>1)) % 2D parallel and 3D Case
        fprintf(file,'%12d\n',neighbor_sendto_w);
    end

    fprintf(file,'%12d\n',numsegment_w);

    if (OCPu==2) % 2D Case
        for ii = 1:1:numsegment_w
            fprintf(file,'%12d\n',MCPu*NCPu+1);
            fprintf(file,'%12d\n',MCPu*NCPu);
        end
    end
    
    if (OCPu>2) % 3D Case    
        for ii = 1:1:numsegment_w
            fprintf(file,'%12d\n',(MCPu*NCPu*OCPu-2*MCPu*NCPu+1)+(ii-1)*MCPu*NCPu*2);
            fprintf(file,'%12d\n',length_w);
        end
    end

    fprintf(file,'%12d\n',itag);
    fprintf(file,'%12d\n',iacc_receieve);
    fprintf(file,'%12d\n',neighbor_recvfrom_w);
    fprintf(file,'%12d\n',numsegment_w);

    if (OCPu==2) % 2D Case      
        for ii = 1:1:numsegment_w
            fprintf(file,'%12d\n',1);
            fprintf(file,'%12d\n',MCPu*NCPu);
        end
    end    
    
    if (OCPu>2) % 3D Case      
        for ii = 1:1:numsegment_w
            fprintf(file,'%12d\n',1+(ii-1)*MCPu*NCPu);
            fprintf(file,'%12d\n',length_w);
        end
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% even core communication
if mod(k,2)==0
    
    fprintf(file,'%12d\n',itag);
    fprintf(file,'%12d\n',iacc_receieve);
    fprintf(file,'%12d\n',neighbor_recvfrom_w);
    fprintf(file,'%12d\n',numsegment_w);
    
    for ii = 1:1:numsegment_w
        fprintf(file,'%12d\n',1+(ii-1)*MCPu*NCPu);
        fprintf(file,'%12d\n',length_w);
    end
    
    fprintf(file,'%12d\n',itag);
    if ( (PeriodicZ ==0) && (k == CoresZ)) 
        fprintf(file,'%12d\n',1);
    else
        fprintf(file,'%12d\n',iacc_send);
    end
    fprintf(file,'%12d\n',neighbor_sendto_w);
    fprintf(file,'%12d\n',numsegment_w);
    
    for ii = 1:1:numsegment_w
        fprintf(file,'%12d\n',(MCPu*NCPu*OCPu-2*MCPu*NCPu+1)+(ii-1)*MCPu*NCPu*2);
        fprintf(file,'%12d\n',length_w);
    end
end    
fclose('all');
disp(strcat('Generated File:    lworkwu.',int2str(currentcore),'.dat'));

      end % X direction cores
    end % Y direction cores
end % Z direction cores