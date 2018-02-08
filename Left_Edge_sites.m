function [Vx,Vy] = Left_Edge_sites(U,Nr,Ns)

Vx(:,1)=U(:,1);

NV=floor(Ns/3);

Vy(:,1) = U(:,Nr*NV+1);
   
Hstep = NV*2*Nr;

Vx(:,2) = U(:,1+Hstep); Vy(:,2) = U(:,Nr*NV+1+Hstep);
Vx(:,3) = U(:,1+Hstep*2); Vy(:,3) = U(:,Nr*NV+1+Hstep*2);

for i=4:Ns
if mod(i,3)==1
Vx(:,i) = U(:,1+floor(i/3)*Nr);        
Vy(:,i) = U(:,Nr*NV+1+floor(i/3)*Nr);
elseif mod(i,3)==2
Vx(:,i) = U(:,1+Hstep+floor(i/3)*Nr);
Vy(:,i) = U(:,Nr*NV+1+Hstep+floor(i/3)*Nr);
elseif mod(i,3)==0
Vx(:,i) = U(:,1+Hstep*2+(floor(i/3)-1)*Nr);
Vy(:,i) = U(:,Nr*NV+1+Hstep*2+(floor(i/3)-1)*Nr);    
end
end
end