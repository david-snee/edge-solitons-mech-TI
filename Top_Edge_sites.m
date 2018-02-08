function [Vx,Vy] = Top_Edge_sites(U,Nr,Ns)

if Ns==3
    NV=floor(Ns/3)+1;
elseif Ns>3
    NV=floor(Ns/3);
end

Vx(:,1:Nr)=U(:,5*Nr*NV-Nr+1:5*Nr*NV);
Vy(:,1:Nr) = U(:,6*Nr*NV-Nr+1:6*Nr*NV);
end