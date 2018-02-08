%% Performs the time evolution. This code relies on Initial condition
%  coming from the EDGE_nls_ampfac.m, Dispersion relation and profiles
%  of 1D edge states numerically read in from linear semi-infinite
%  problem, and the nonliner network in NETW_nonlin.m and on a specified
%  domain from either NETW_generate_graph or NETW_generate_graph_holhex.

function mech_NETW_nonlin_solver(shape,Nr,Ns,tf,epsilon,par,ttheta,amp,pos,dt,dispersion,nu,option)

omega0=par(1);
A=par(2);
f=par(3);
coef=-(omega0^2+A*f);

sigma= omega0^2/6; %coefficient due to sin expansion

if option==1 
% FULL DOMAIN
rsx0=load(['NETW_' shape int2str(Nr) 'x' int2str(Ns) '_rsx0.dat']);
rsx1=load(['NETW_' shape int2str(Nr) 'x' int2str(Ns) '_rsx1.dat']);
rsx2=load(['NETW_' shape int2str(Nr) 'x' int2str(Ns) '_rsx2.dat']);
rsy0=load(['NETW_' shape int2str(Nr) 'x' int2str(Ns) '_rsy0.dat']);
rsy1=load(['NETW_' shape int2str(Nr) 'x' int2str(Ns) '_rsy1.dat']);
rsy2=load(['NETW_' shape int2str(Nr) 'x' int2str(Ns) '_rsy2.dat']);
bonds=load(['NETW_' shape int2str(Nr) 'x' int2str(Ns) '_bonds.dat']);
elseif option==2 || option==4 || option==5 || option==6 || option==7 || option==8
% PERIODIC DOMAIN
rsx0=load(['NETW_periodic_' shape int2str(Nr) 'x' int2str(Ns) '_rsx0.dat']);
rsx1=load(['NETW_periodic_' shape int2str(Nr) 'x' int2str(Ns) '_rsx1.dat']);
rsx2=load(['NETW_periodic_' shape int2str(Nr) 'x' int2str(Ns) '_rsx2.dat']);
rsy0=load(['NETW_periodic_' shape int2str(Nr) 'x' int2str(Ns) '_rsy0.dat']);
rsy1=load(['NETW_periodic_' shape int2str(Nr) 'x' int2str(Ns) '_rsy1.dat']);
rsy2=load(['NETW_periodic_' shape int2str(Nr) 'x' int2str(Ns) '_rsy2.dat']);
bonds=load(['NETW_periodic_' shape int2str(Nr) 'x' int2str(Ns) '_bonds.dat']);
elseif option==3
% SITES REMOVED
rsx0=load(['NETW_sites_removed' int2str(Nr) 'x' int2str(Ns) '_rsx0.dat']);
rsx1=load(['NETW_sites_removed' int2str(Nr) 'x' int2str(Ns) '_rsx1.dat']);
rsx2=load(['NETW_sites_removed' int2str(Nr) 'x' int2str(Ns) '_rsx2.dat']);
rsy0=load(['NETW_sites_removed' int2str(Nr) 'x' int2str(Ns) '_rsy0.dat']);
rsy1=load(['NETW_sites_removed' int2str(Nr) 'x' int2str(Ns) '_rsy1.dat']);
rsy2=load(['NETW_sites_removed' int2str(Nr) 'x' int2str(Ns) '_rsy2.dat']);
bonds=load(['NETW_sites_removed' int2str(Nr) 'x' int2str(Ns) '_bonds.dat']);
end


if isempty(amp)
    amp = ttheta;
    for i = 1:length(ttheta)
        amp(i) = epsilon(i)*EDGE_mech_nls_ampfac(par,ttheta(i),Nr,dispersion);
    end
end

N1 = size(rsx0,1);
N2 = size(rsy0,1);
N3 = size(rsx1,1);
N4 = size(rsy1,1);
N5 = size(rsx2,1);
N6 = size(rsy2,1);
N=[N1 N2 N3 N4 N5 N6]; TN=sum(N);
asym = load(['EDGE_N' int2str(Nr) '_f' num2str(f) '_omega0' num2str(omega0) '_A' num2str(A) '.mat']);
asym.vals=(sqrt(-real(asym.vals)));
[~,index] = min(abs(asym.k_vec-ttheta));
vec_matrix=asym.vecs(index,:,dispersion)';
vo = csapi(asym.k_vec,asym.vals(:,dispersion)); %value chosen from edge problem.
ndi = 3;
vdi = zeros(ndi+1,1);
for i = 1:ndi+1
    tmp = fnder(vo,i-1);
    vdi(i) = fnval(tmp,ttheta);
end
alsec = vdi(3); alfirst= vdi(2); alpha0=vdi(1);
ndi = 3;
initpos = zeros(TN,1);
initvel = zeros(TN,1);
%% Initialises condition on all cells (sites).
for i = 1:TN
   
      
   if i<=N1 
       rr = rsx0(i,1);
       ss = rsx0(i,2);
       nullind=rr+1; 
   elseif i>N1 && i<=sum(N(1:2))
       rr = rsy0((i-N1),1);
       ss = rsy0((i-N1),2);
       nullind=Nr+rr+1;
   elseif i>sum(N(1:2)) && i<=sum(N(1:3))
       rr = rsx1((i-sum(N(1:2))),1);
       ss = rsx1((i-sum(N(1:2))),2); 
       nullind=2*Nr+rr+1;
   elseif i>sum(N(1:3)) && i<=sum(N(1:4))
       rr = rsy1((i-sum(N(1:3))),1);
       ss = rsy1((i-sum(N(1:3))),2); 
       nullind=3*Nr+rr+1;
   elseif i>sum(N(1:4)) && i<=sum(N(1:5))
       rr = rsx2((i-sum(N(1:4))),1);
       ss = rsx2((i-sum(N(1:4))),2); 
       nullind=4*Nr+rr+1;
   elseif i>sum(N(1:5)) && i<=TN
       rr = rsy2((i-sum(N(1:5))),1);
       ss = rsy2((i-sum(N(1:5))),2); 
       nullind=5*Nr+rr+1;
   end    
   
    scell=floor(ss/3); 
            
        if option==1 || option==3
            % Bright Soliton solution
           initpos(i) = initpos(i)+amp.*sech(epsilon.*(scell-pos)).*exp(1i.*ttheta.*scell).*vec_matrix(nullind) + ... %term.
            amp.*sech(epsilon.*(scell-pos)).*exp(-1i.*ttheta.*scell).*conj(vec_matrix(nullind)); %conjugate term.
        initvel(i) = initvel(i)-amp.*alpha0.*1i.*sech(epsilon.*(scell-pos)).*exp(1i.*ttheta.*scell).*vec_matrix(nullind) + ... %term.
            amp.*alpha0.*1i.*sech(epsilon.*(scell-pos)).*exp(-1i.*ttheta.*scell).*conj(vec_matrix(nullind)); %conjugate term
        elseif option==2
        % Peregrine soliton solution
        inittime=-300*epsilon^2*alsec;
        
        num=4*(1+2*1i*inittime);
        den=1+4*epsilon^2*(scell-pos)^2+4*inittime^2;
        
       initpos(i) = initpos(i)+amp*(1-num./den).*exp(1i*inittime*(1-alpha0)).*exp(1i.*ttheta.*scell).*vec_matrix(nullind) + ... %term.
           conj(amp*(1-num./den).*exp(1i*inittime*(1-alpha0)).*exp(1i.*ttheta.*scell).*vec_matrix(nullind)); %conjugate term.
       initvel(i) = initvel(i)-amp*alpha0.*1i.*(1-num./den).*exp(1i*inittime*(1-alpha0)).*exp(1i.*ttheta.*scell).*vec_matrix(nullind) + ... %term.
          amp*alpha0.*1i.*(1-conj(num)./den).*exp(-1i*inittime*(1-alpha0)).*exp(-1i.*ttheta.*scell).*conj(vec_matrix(nullind)); %conjugate term            
      
        elseif option==4
       % K-M soliton solution
       inittime=-0*alsec/2; phi=0.5;
           p=2*sinh(phi); OMEGA=2*sinh(2*phi);

              num=cos(OMEGA*inittime-2*1i*phi)-cosh(phi)*cosh(p*epsilon*(scell-pos));
              den=cos(OMEGA*inittime)-cosh(phi)*cosh(p*epsilon*(scell-pos));

              initpos(i) = initpos(i)+amp.*num./den.*exp(1i*inittime*(2-alpha0)).*exp(1i.*ttheta.*scell).*vec_matrix(nullind) + ... %term.
              amp.*conj(num./den).*exp(-1i*inittime*(2-alpha0)).*exp(-1i.*ttheta.*scell).*conj(vec_matrix(nullind)); %conjugate term.
              initvel(i) = initvel(i)-amp.*alpha0.*1i.*num./den.*exp(1i*inittime*(2-alpha0)).*exp(1i.*ttheta.*scell).*vec_matrix(nullind) + ... %term.
              amp.*alpha0.*1i.*conj(num./den).*exp(-1i*inittime*(2-alpha0)).*exp(-1i.*ttheta.*scell).*conj(vec_matrix(nullind));                 
        elseif option==5
            % A- breather solution
            inittime=-0*alsec/2; phi=2;
            p=2*sin(phi); OMEGA=2*sin(2*phi);

            num=cosh(OMEGA*inittime-2*1i*phi)-cos(phi)*cos(p*epsilon*(scell-pos));
            den=cosh(OMEGA*inittime)-cos(phi)*cos(p*epsilon*(scell-pos));

            initpos(i) = initpos(i)+amp.*num./den.*exp(2*1i*inittime).*exp(1i.*ttheta.*scell).*vec_matrix(nullind) + ... %term.
            amp.*conj(num)./den.*exp(2*-1i*inittime).*exp(-1i.*ttheta.*scell).*conj(vec_matrix(nullind)); %conjugate term.
            initvel(i) = initvel(i)-amp.*alpha0.*1i.*num./den.*exp(2*1i*inittime).*exp(1i.*ttheta.*scell).*vec_matrix(nullind) + ... %term.
            amp.*alpha0.*1i.*conj(num)./den.*exp(2*-1i*inittime).*exp(-1i.*ttheta.*scell).*conj(vec_matrix(nullind));
        elseif option==6 
         % Dark soliton solution
         inittime=0*alsec*epsilon^2;
         
         initpos(i) = initpos(i)+amp.*tanh(epsilon.*(scell-pos)).*exp(1i.*ttheta.*scell).*exp(1i*(epsilon^2-alpha0)*inittime).*vec_matrix(nullind) + ... %term.
            conj(amp.*tanh(epsilon.*(scell-pos)).*exp(1i.*ttheta.*scell).*exp(1i*(epsilon^2-alpha0)*inittime).*vec_matrix(nullind)); %conjugate term.
        initvel(i) = initvel(i)-amp.*alpha0.*1i.*tanh(epsilon.*(scell-pos)).*exp(1i.*ttheta.*scell).*exp(1i*(epsilon^2-alpha0)*inittime).*vec_matrix(nullind) - ... %term.
            conj(amp.*alpha0.*1i.*tanh(epsilon.*(scell-pos)).*exp(1i.*ttheta.*scell).*exp(1i*(epsilon^2-alpha0)*inittime).*vec_matrix(nullind)); %conjugate term
        elseif option==7 
            % Two soliton collision
        if scell<(Ns/3)/2    
              initpos(i) = initpos(i)+amp.*sech(epsilon.*(scell-pos+35)).*exp(1i.*ttheta.*scell).*vec_matrix(nullind) + ... %term.
                amp.*sech(epsilon.*(scell-pos+35)).*exp(-1i.*ttheta.*scell).*conj(vec_matrix(nullind)); %conjugate term.
              initvel(i) = initvel(i)-amp.*alpha0.*1i.*sech(epsilon.*(scell-pos+35)).*exp(1i.*ttheta.*scell).*vec_matrix(nullind) + ... %term.
                amp.*alpha0.*1i.*sech(epsilon.*(scell-pos+35)).*exp(-1i.*ttheta.*scell).*conj(vec_matrix(nullind)); %conjugate term
        elseif scell>=(Ns/3)/2
              initpos(i) = initpos(i)+amp.*sech(epsilon.*(scell-pos-35)).*exp(1i.*ttheta.*scell).*vec_matrix(nullind) + ... %term.
               amp.*sech(epsilon.*(scell-pos-35)).*exp(-1i.*ttheta.*scell).*conj(vec_matrix(nullind)); %conjugate term.
              initvel(i) = initvel(i)+amp.*alpha0.*1i.*sech(epsilon.*(scell-pos-35)).*exp(1i.*ttheta.*scell).*vec_matrix(nullind) - ... %term.
               amp.*alpha0.*1i.*sech(epsilon.*(scell-pos-35)).*exp(-1i.*ttheta.*scell).*conj(vec_matrix(nullind)); %conjugate term
        end
        else 
           % Gaussian envelope
            inittime=-0*alsec/2*epsilon^2;
        
       initpos(i) = initpos(i)+amp*exp(-epsilon^2*(scell-pos)^2/(4*nu^2)).*exp(1i.*(ttheta.*scell-alpha0*inittime)).*vec_matrix(nullind) + ... %term.
           conj(amp*exp(-epsilon^2*(scell-pos)^2/(4*nu^2)).*exp(1i.*(ttheta.*scell-alpha0*inittime)).*vec_matrix(nullind)); %conjugate term.
       initvel(i) = initvel(i)-amp*alpha0.*1i.*exp(-epsilon^2*(scell-pos)^2/(4*nu^2)).*exp(1i.*(ttheta.*scell-alpha0*inittime)).*vec_matrix(nullind) - ... %term.
          conj(amp*alpha0.*1i.*exp(-epsilon^2*(scell-pos)^2/(4*nu^2)).*exp(1i.*(ttheta.*scell-alpha0*inittime)).*vec_matrix(nullind)); %conjugate term            

        end  
end


if option==1 || option==3 
M= mech_NETW_matrix(rsx0,rsy0,rsx1,rsy1,rsx2,rsy2,bonds,par);
else
M= mech_NETW_matrix_periodic_coupling(rsx0,rsy0,rsx1,rsy1,rsx2,rsy2,bonds,par,Ns);    
end



sigmamat=spdiags(sigma*ones(TN,1),0,TN,TN);

MS=[zeros(size(M)) eye(size(M));M zeros(size(M))];
S=[zeros(size(M)) zeros(size(M)); sigmamat zeros(size(M))];
init=[initpos ; initvel];

% Solve nonlinear system
[T,U] = ode45(@(t,lhs) MS*lhs+S*(lhs.^3), 0:dt:tf, init);


epsstr = [];
tthetastr = [];
ampstr = [];
posstr = [];
for i = 1:length(ttheta)
    epsstr = [epsstr num2str(epsilon(i)) 't'];
    tthetastr = [tthetastr num2str(ttheta(i)) 't'];
    ampstr = [ampstr num2str(amp(i)) 't'];
    posstr = [posstr num2str(pos(i)) 't'];
end
epsstr = epsstr(1:end-1);
tthetastr = tthetastr(1:end-1);
ampstr = ampstr(1:end-1);
posstr = posstr(1:end-1);

if option==1
%% SOLITON
U1=U(1:(tf/dt)/8,:);
U2=U((tf/dt)/8+1:2*(tf/dt)/8,:);
U3=U(2*(tf/dt)/8+1:3*(tf/dt)/8,:);
U4=U(3*(tf/dt)/8+1:4*(tf/dt)/8,:);
U5=U(4*(tf/dt)/8+1:5*(tf/dt)/8,:);
U6=U(5*(tf/dt)/8+1:6*(tf/dt)/8,:);
U7=U(6*(tf/dt)/8+1:7*(tf/dt)/8,:);
U8=U(7*(tf/dt)/8+1:(tf/dt)+1,:);
fn1 = ['mech_NETW_nlin_' shape int2str(Nr) 'x' int2str(Ns) '_dispersion' num2str(dispersion) '_eps' epsstr '_ttheta' tthetastr '_amp' ampstr 'U1' int2str(tf) '_tf'];
fn2 = ['mech_NETW_nlin_' shape int2str(Nr) 'x' int2str(Ns) '_dispersion' num2str(dispersion) '_eps' epsstr '_ttheta' tthetastr '_amp' ampstr 'U2' int2str(tf) '_tf'];
fn3 = ['mech_NETW_nlin_' shape int2str(Nr) 'x' int2str(Ns) '_dispersion' num2str(dispersion) '_eps' epsstr '_ttheta' tthetastr '_amp' ampstr 'U3' int2str(tf) '_tf'];
fn4 = ['mech_NETW_nlin_' shape int2str(Nr) 'x' int2str(Ns) '_dispersion' num2str(dispersion) '_eps' epsstr '_ttheta' tthetastr '_amp' ampstr 'U4' int2str(tf) '_tf'];
fn5 = ['mech_NETW_nlin_' shape int2str(Nr) 'x' int2str(Ns) '_dispersion' num2str(dispersion) '_eps' epsstr '_ttheta' tthetastr '_amp' ampstr 'U5' int2str(tf) '_tf'];
fn6 = ['mech_NETW_nlin_' shape int2str(Nr) 'x' int2str(Ns) '_dispersion' num2str(dispersion) '_eps' epsstr '_ttheta' tthetastr '_amp' ampstr 'U6' int2str(tf) '_tf'];
fn7 = ['mech_NETW_nlin_' shape int2str(Nr) 'x' int2str(Ns) '_dispersion' num2str(dispersion) '_eps' epsstr '_ttheta' tthetastr '_amp' ampstr 'U7' int2str(tf) '_tf'];
fn8 = ['mech_NETW_nlin_' shape int2str(Nr) 'x' int2str(Ns) '_dispersion' num2str(dispersion) '_eps' epsstr '_ttheta' tthetastr '_amp' ampstr 'U8' int2str(tf) '_tf'];
save([fn1 '.mat'],'T','U1');
save([fn2 '.mat'],'U2');
save([fn3 '.mat'],'U3');
save([fn4 '.mat'],'U4');
save([fn5 '.mat'],'U5');
save([fn6 '.mat'],'U6');
save([fn7 '.mat'],'U7');
save([fn8 '.mat'],'U8');
elseif option==2
%% PEREGRINE Soliton
fn1 = ['peregrine_' shape int2str(Nr) 'x' int2str(Ns) '_dispersion' num2str(dispersion) '_eps' epsstr '_ttheta' tthetastr '_amp' ampstr 'U' int2str(tf) '_tf'];
save([fn1 '.mat'],'T','U');
elseif option==3
%% SITES REMOVED 
U1=U(1:(tf/dt)/8,:);
U2=U((tf/dt)/8+1:2*(tf/dt)/8,:);
U3=U(2*(tf/dt)/8+1:3*(tf/dt)/8,:);
U4=U(3*(tf/dt)/8+1:4*(tf/dt)/8,:);
U5=U(4*(tf/dt)/8+1:5*(tf/dt)/8,:);
U6=U(5*(tf/dt)/8+1:6*(tf/dt)/8,:);
U7=U(6*(tf/dt)/8+1:7*(tf/dt)/8,:);
U8=U(7*(tf/dt)/8+1:(tf/dt)+1,:);
fn1 = ['sites_removed_' int2str(Nr) 'x' int2str(Ns) '_dispersion' num2str(dispersion) '_eps' epsstr '_ttheta' tthetastr '_amp' ampstr 'U1' int2str(tf) '_tf'];
fn2 = ['sites_removed_' int2str(Nr) 'x' int2str(Ns) '_dispersion' num2str(dispersion) '_eps' epsstr '_ttheta' tthetastr '_amp' ampstr 'U2' int2str(tf) '_tf'];
fn3 = ['sites_removed_'  int2str(Nr) 'x' int2str(Ns) '_dispersion' num2str(dispersion) '_eps' epsstr '_ttheta' tthetastr '_amp' ampstr 'U3' int2str(tf) '_tf'];
fn4 = ['sites_removed_'  int2str(Nr) 'x' int2str(Ns) '_dispersion' num2str(dispersion) '_eps' epsstr '_ttheta' tthetastr '_amp' ampstr 'U4' int2str(tf) '_tf'];
fn5 = ['sites_removed_' int2str(Nr) 'x' int2str(Ns) '_dispersion' num2str(dispersion) '_eps' epsstr '_ttheta' tthetastr '_amp' ampstr 'U5' int2str(tf) '_tf'];
fn6 = ['sites_removed_' int2str(Nr) 'x' int2str(Ns) '_dispersion' num2str(dispersion) '_eps' epsstr '_ttheta' tthetastr '_amp' ampstr 'U6' int2str(tf) '_tf'];
fn7 = ['sites_removed_' int2str(Nr) 'x' int2str(Ns) '_dispersion' num2str(dispersion) '_eps' epsstr '_ttheta' tthetastr '_amp' ampstr 'U7' int2str(tf) '_tf'];
fn8 = ['sites_removed_'  int2str(Nr) 'x' int2str(Ns) '_dispersion' num2str(dispersion) '_eps' epsstr '_ttheta' tthetastr '_amp' ampstr 'U8' int2str(tf) '_tf'];
save([fn1 '.mat'],'T','U1');
save([fn2 '.mat'],'U2');
save([fn3 '.mat'],'U3');
save([fn4 '.mat'],'U4');
save([fn5 '.mat'],'U5');
save([fn6 '.mat'],'U6');
save([fn7 '.mat'],'U7');
save([fn8 '.mat'],'U8');    
elseif option==4
%% K-M Soliton
fn1 = ['K_M_' shape int2str(Nr) 'x' int2str(Ns) '_dispersion' num2str(dispersion) '_eps' epsstr '_ttheta' tthetastr '_amp' ampstr 'U' int2str(tf) '_tf'];
save([fn1 '.mat'],'T','U');
elseif option==5
%% A- breather
fn1 = ['A_breather_' shape int2str(Nr) 'x' int2str(Ns) '_dispersion' num2str(dispersion) '_eps' epsstr '_ttheta' tthetastr '_amp' ampstr 'U' int2str(tf) '_tf'];
save([fn1 '.mat'],'T','U');
elseif option==6 
%% DARK soliton
U1=U(1:(tf/dt)/8,:);
U2=U((tf/dt)/8+1:2*(tf/dt)/8,:);
U3=U(2*(tf/dt)/8+1:3*(tf/dt)/8,:);
U4=U(3*(tf/dt)/8+1:4*(tf/dt)/8,:);
U5=U(4*(tf/dt)/8+1:5*(tf/dt)/8,:);
U6=U(5*(tf/dt)/8+1:6*(tf/dt)/8,:);
U7=U(6*(tf/dt)/8+1:7*(tf/dt)/8,:);
U8=U(7*(tf/dt)/8+1:(tf/dt)+1,:);
fn1 = ['dark_' shape int2str(Nr) 'x' int2str(Ns) '_dispersion' num2str(dispersion) '_eps' epsstr '_ttheta' tthetastr '_amp' ampstr 'U1' int2str(tf) '_tf'];
fn2 = ['dark_' shape int2str(Nr) 'x' int2str(Ns) '_dispersion' num2str(dispersion) '_eps' epsstr '_ttheta' tthetastr '_amp' ampstr 'U2' int2str(tf) '_tf'];
fn3 = ['dark_' shape int2str(Nr) 'x' int2str(Ns) '_dispersion' num2str(dispersion) '_eps' epsstr '_ttheta' tthetastr '_amp' ampstr 'U3' int2str(tf) '_tf'];
fn4 = ['dark_' shape int2str(Nr) 'x' int2str(Ns) '_dispersion' num2str(dispersion) '_eps' epsstr '_ttheta' tthetastr '_amp' ampstr 'U4' int2str(tf) '_tf'];
fn5 = ['dark_' shape int2str(Nr) 'x' int2str(Ns) '_dispersion' num2str(dispersion) '_eps' epsstr '_ttheta' tthetastr '_amp' ampstr 'U5' int2str(tf) '_tf'];
fn6 = ['dark_' shape int2str(Nr) 'x' int2str(Ns) '_dispersion' num2str(dispersion) '_eps' epsstr '_ttheta' tthetastr '_amp' ampstr 'U6' int2str(tf) '_tf'];
fn7 = ['dark_' shape int2str(Nr) 'x' int2str(Ns) '_dispersion' num2str(dispersion) '_eps' epsstr '_ttheta' tthetastr '_amp' ampstr 'U7' int2str(tf) '_tf'];
fn8 = ['dark_' shape int2str(Nr) 'x' int2str(Ns) '_dispersion' num2str(dispersion) '_eps' epsstr '_ttheta' tthetastr '_amp' ampstr 'U8' int2str(tf) '_tf'];
save([fn1 '.mat'],'T','U1');
save([fn2 '.mat'],'U2');
save([fn3 '.mat'],'U3');
save([fn4 '.mat'],'U4');
save([fn5 '.mat'],'U5');
save([fn6 '.mat'],'U6');
save([fn7 '.mat'],'U7');
save([fn8 '.mat'],'U8');
elseif option==7 
%% Two soliton collision
U1=U(1:(tf/dt)/8,:);
U2=U((tf/dt)/8+1:2*(tf/dt)/8,:);
U3=U(2*(tf/dt)/8+1:3*(tf/dt)/8,:);
U4=U(3*(tf/dt)/8+1:4*(tf/dt)/8,:);
U5=U(4*(tf/dt)/8+1:5*(tf/dt)/8,:);
U6=U(5*(tf/dt)/8+1:6*(tf/dt)/8,:);
U7=U(6*(tf/dt)/8+1:7*(tf/dt)/8,:);
U8=U(7*(tf/dt)/8+1:(tf/dt)+1,:);
fn1 = ['Two_Soliton_' int2str(Nr) 'x' int2str(Ns) '_dispersion' num2str(dispersion) '_eps' epsstr '_ttheta' tthetastr '_amp' ampstr 'U1' int2str(tf) '_tf'];
fn2 = ['Two_Soliton_' int2str(Nr) 'x' int2str(Ns) '_dispersion' num2str(dispersion) '_eps' epsstr '_ttheta' tthetastr '_amp' ampstr 'U2' int2str(tf) '_tf'];
fn3 = ['Two_Soliton_' int2str(Nr) 'x' int2str(Ns) '_dispersion' num2str(dispersion) '_eps' epsstr '_ttheta' tthetastr '_amp' ampstr 'U3' int2str(tf) '_tf'];
fn4 = ['Two_Soliton_' int2str(Nr) 'x' int2str(Ns) '_dispersion' num2str(dispersion) '_eps' epsstr '_ttheta' tthetastr '_amp' ampstr 'U4' int2str(tf) '_tf'];
fn5 = ['Two_Soliton_' int2str(Nr) 'x' int2str(Ns) '_dispersion' num2str(dispersion) '_eps' epsstr '_ttheta' tthetastr '_amp' ampstr 'U5' int2str(tf) '_tf'];
fn6 = ['Two_Soliton_' int2str(Nr) 'x' int2str(Ns) '_dispersion' num2str(dispersion) '_eps' epsstr '_ttheta' tthetastr '_amp' ampstr 'U6' int2str(tf) '_tf'];
fn7 = ['Two_Soliton_' int2str(Nr) 'x' int2str(Ns) '_dispersion' num2str(dispersion) '_eps' epsstr '_ttheta' tthetastr '_amp' ampstr 'U7' int2str(tf) '_tf'];
fn8 = ['Two_Soliton_' int2str(Nr) 'x' int2str(Ns) '_dispersion' num2str(dispersion) '_eps' epsstr '_ttheta' tthetastr '_amp' ampstr 'U8' int2str(tf) '_tf'];
save([fn1 '.mat'],'T','U1');
save([fn2 '.mat'],'U2');
save([fn3 '.mat'],'U3');
save([fn4 '.mat'],'U4');
save([fn5 '.mat'],'U5');
save([fn6 '.mat'],'U6');
save([fn7 '.mat'],'U7');
save([fn8 '.mat'],'U8');
else
%% Guassian initial data  
U1=U(1:(tf/dt)/4,:);
U2=U((tf/dt)/4+1:2*(tf/dt)/4,:);
U3=U(2*(tf/dt)/4+1:3*(tf/dt)/4,:);
U4=U(3*(tf/dt)/4+1:(tf/dt)+1,:);
fn1 = ['gaussian_' shape int2str(Nr) 'x' int2str(Ns) '_dispersion' num2str(dispersion) '_eps' epsstr '_ttheta' tthetastr '_amp' ampstr 'U1' int2str(tf) '_tf_' num2str(nu) 'nu'];
fn2 = ['gaussian_' shape int2str(Nr) 'x' int2str(Ns) '_dispersion' num2str(dispersion) '_eps' epsstr '_ttheta' tthetastr '_amp' ampstr 'U2' int2str(tf) '_tf_' num2str(nu) 'nu'];
fn3 = ['gaussian_' shape int2str(Nr) 'x' int2str(Ns) '_dispersion' num2str(dispersion) '_eps' epsstr '_ttheta' tthetastr '_amp' ampstr 'U3' int2str(tf) '_tf_' num2str(nu) 'nu'];
fn4 = ['gaussian_' shape int2str(Nr) 'x' int2str(Ns) '_dispersion' num2str(dispersion) '_eps' epsstr '_ttheta' tthetastr '_amp' ampstr 'U4' int2str(tf) '_tf_' num2str(nu) 'nu'];
save([fn1 '.mat'],'T','U1');
save([fn2 '.mat'],'U2');
save([fn3 '.mat'],'U3');
save([fn4 '.mat'],'U4');
end
disp([fn1 ' complete;']);

end
