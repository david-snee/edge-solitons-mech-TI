
function EDGE_dispersion_field_mechanical_TI(N,par,k_step)


% N is the number of lattice sites in the x direction.  This can be
% arbitrarily large.

%par is an array holding values for F,A,omega0.

%k_step defines how smooth/detailed the plot will be.


    NT = 6*N;    %6 variables.

    f=par(3);    %free parameters
    A=par(2);
    omega0=par(1);
    coef=-(omega0^2+A*f);

function MAT=mat_build(S,D1,D2,D3) %size of required matrix
   M1=spdiags(D1*ones(S,1),0,S,S); %main diagonal
   M2=spdiags(D2*ones(S,1),1,S,S);   %above main diagonal
   M3=spdiags(D3*ones(S,1),-1,S,S);   %below main diagonal

   MAT=M1+M2+M3;
end

k_vec= 0:(2*pi)/(k_step):2*pi;

    vals = zeros(length(k_vec),NT);
    vecs = zeros(length(k_vec),NT,NT);
    
    % Compute dispersion relationship from monodromy matrix for all frequency values omega.
    
    time = clock;

    tic
    
    mono_mat=zeros(NT);
    
for j=0:1
    mono_mat(j*N+1:(j+1)*N,j*N+1:(j+1)*N)=mat_build(N,coef,f,f); 
    %x0-x0,y0-y0
end
    
for i=2:5
    mono_mat(i*N+1:(i+1)*N,i*N+1:(i+1)*N)=mat_build(N,coef,-f/2,-f/2);
    %x1-x1,y1-y1,x2-x2,y2-y2
end
    
for ii=[0 4]
    mono_mat(ii*N+1:(2+ii)*N,2*N+1:4*N)=mat_build(2*N,f,0,0);   
    %x0-x1,y0-y1,x2-x1,y2-y1
end

for iii=[0 4]
    mono_mat(2*N+1:4*N,iii*N+1:(2+iii)*N)=mat_build(2*N,f,0,0);
    %x1-x0,y1-y0,x1-x2,y1-y2
end
 
%x1-y1,y1-x1,x2-y2,y2-x2
mono_mat(2*N+1:3*N,3*N+1:4*N)=mat_build(N,0,sqrt(3)*f/2,-sqrt(3)*f/2);
mono_mat(3*N+1:4*N,2*N+1:3*N)=mat_build(N,0,-sqrt(3)*f/2,sqrt(3)*f/2);
mono_mat(4*N+1:5*N,5*N+1:6*N)=mat_build(N,0,-sqrt(3)*f/2,sqrt(3)*f/2);
mono_mat(5*N+1:6*N,4*N+1:5*N)=mat_build(N,0,sqrt(3)*f/2,-sqrt(3)*f/2);
 
    
        
    for kk=1:length(k_vec)
        
        % Compute monodromy matrix for given k value. 
     mono_mat(1:N,4*N+1:5*N)=mat_build(N,f*exp(-1i*k_vec(kk)),0,0);
     mono_mat(4*N+1:5*N,1:N)=mat_build(N,f*exp(1i*k_vec(kk)),0,0);
     mono_mat(N+1:2*N,5*N+1:6*N)=mat_build(N,f*exp(-1i*k_vec(kk)),0,0);
     mono_mat(5*N+1:6*N,N+1:2*N)=mat_build(N,f*exp(1i*k_vec(kk)),0,0);
     
        
    %Eigenvalues and eigenfunctions are computed.            
    [V,D] = eig(mono_mat);
    D=diag(D);
    [B,IX] = sort(abs(real(D)));
    D = D(IX); V = V(:,IX);
    vals(kk,:) = D; vecs(kk,:,:) = V;
    
    %time related
        eta = datenum([0 0 0 0 0 toc/kk*(length(k_vec)-kk)]); 
    days_left = floor(eta); 
    time_left = datevec(eta-days_left); 
    fprintf([num2str(kk) ' out of ' num2str(length(k_vec))...
        ' - Est time left: ' num2str(days_left) ' d ' ... 
        num2str(time_left(4)) ' h ' ...
        num2str(time_left(5)) ' m ' ...
        num2str(time_left(6)) ' s \n'])
    
        
    end
    
    if N>60
        val(:,1:10) = vals(:,2*N-4:2*N+5);
        vec(:,:,1:10)=vecs(:,:,2*N-4:2*N+5);
        vals=val; vecs=vec;
    else 
        vals=val; vecs=vec;
    end
save(['EDGE_N' int2str(N) '_f' num2str(f) '_omega0' num2str(omega0) '_A' num2str(A) '.mat'],'k_vec','vals','vecs');

    
end