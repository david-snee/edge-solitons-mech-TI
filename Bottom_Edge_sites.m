function [Vx,Vy] = Bottom_Edge_sites(U,Nr,Ns)

Vx(:,1:Nr)=U(:,1:Nr);

if Ns==3
    NV=floor(Ns/3)+1;
elseif Ns>3
    NV=floor(Ns/3);
end

Vy(:,1:Nr) = U(:,Nr*NV+1:Nr*NV+Nr);
end