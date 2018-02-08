function [Vx,Vy] = Right_Edge_sites(U,Nr,Ns)

Vx(:,1)=U(:,Nr);

if Ns==3
    NV=floor(Ns/3)+1;
elseif Ns>3
    NV=floor(Ns/3);
end

Vy(:,1) = U(:,Nr*NV+Nr);
   
Hstep = NV*2*Nr;

Vx(:,2) = U(:,Nr+Hstep); Vy(:,2) = U(:,Nr*NV+Nr+Hstep);
Vx(:,3) = U(:,Nr+Hstep*2); Vy(:,3) = U(:,Nr*NV+Nr+Hstep*2);

for i=4:Ns
if mod(i,3)==1
Vx(:,i) = U(:,Nr+floor(i/3)*Nr);        
Vy(:,i) = U(:,Nr*NV+Nr+floor(i/3)*Nr);
elseif mod(i,3)==2
Vx(:,i) = U(:,Nr+Hstep+floor(i/3)*Nr);
Vy(:,i) = U(:,Nr*NV+Nr+Hstep+floor(i/3)*Nr);
elseif mod(i,3)==0
Vx(:,i) = U(:,Nr+Hstep*2+(floor(i/3)-1)*Nr);
Vy(:,i) = U(:,Nr*NV+Nr+Hstep*2+(floor(i/3)-1)*Nr);    
end
end
end