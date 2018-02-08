%% Performs the numerical simulation for differing parameters.

function mech_NETW_nlin_data()
shape = 'rectangle';
Nr = 60; Ns=180;   %Choose Ns such that mod(Ns,3)=0.
tf =2e3;
dt = 1;
epsilon = 0.1;
par = [(0.75*2*pi)  (3+sqrt(3)) (1.02*2*pi)^2];
ttheta = 9*pi/10; 
amp = [];
pos = 30;
dispersion =119;  %  59,Nr=30  119,Nr=60  (DECREASED VECS SIZE)=> 4,Nr=120
option=8;  %1=FULL DOMAIN  2=PERIODIC DOMAIN (peregrine)  3=REMOVAL OF SITES  4=K-M  5=A breather. 6=Dark  7=Two-soliton  8=Gaussian
nu=0.1;
mech_NETW_nonlin_solver(shape,Nr,Ns,tf,epsilon,par,ttheta,amp,pos,dt,dispersion,nu,option);
end
  
