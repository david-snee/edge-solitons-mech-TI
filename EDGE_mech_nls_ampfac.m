%% Fixes the amplitude to match the exact solution of 1D NLS equation.

function ampfac = EDGE_mech_nls_ampfac(par,ttheta,Nr,dispersion)

omega0=par(1);
A=par(2);
f=par(3);
coef=-(omega0^2+A*f);

sigma= omega0^2/6; %coefficient due to sin expansion

% tf = time to run simulation out to.

% sigma = magnitude of nonlinearity

% epsilon = width of wavpacket.
asym = load(['EDGE_N' int2str(Nr) '_f' num2str(f) '_omega0' num2str(omega0) '_A' num2str(A) '.mat']);
asym.vals=(sqrt(-real(asym.vals)));
ndi = 3;
vo = csapi(asym.k_vec,asym.vals(:,dispersion));  %value chosen from edge problem.
vdi = zeros(ndi+1,1);
for i = 1:ndi+1
    tmp = fnder(vo,i-1);
    vdi(i) = fnval(tmp,ttheta);
end
alsec = vdi(3);


% Find nonlinear component in Fourier space
[~,index] = min(abs(asym.k_vec-ttheta));
vec_matrix=asym.vecs(index,:,dispersion)';

sigmaeff = sigma*(norm(vec_matrix,4)^4)*3/(2*vdi(1));

ampfac = sqrt(alsec/sigmaeff);

end
