%% Plots a movie of multiple frames relying on the code NETW_plot_mnz.m.

function mech_NETW_nlin_movie()
shape = 'rectangle';
Nr = 60; Ns=180;   %Choose Ns such that mod(Ns,3)=0.
epsilon = 0.1;
par = [(0.75*2*pi)  (3+sqrt(3)) (1.02*2*pi)^2];
ttheta = 9*pi/10;
tf =5e2;
amp = [];
pos = 30;
dispersion = 119;  
option=3;

if option==1
% FULL DOMAIN
rsx0=load(['NETW_' shape int2str(Nr) 'x' int2str(Ns) '_rsx0.dat']);
rsx1=load(['NETW_' shape int2str(Nr) 'x' int2str(Ns) '_rsx1.dat']);
rsx2=load(['NETW_' shape int2str(Nr) 'x' int2str(Ns) '_rsx2.dat']);
bonds=load(['NETW_' shape int2str(Nr) 'x' int2str(Ns) '_bonds.dat']);
elseif option==2
% PERIODIC DOMAIN
rsx0=load(['NETW_periodic_' shape int2str(Nr) 'x' int2str(Ns) '_rsx0.dat']);
rsx1=load(['NETW_periodic_' shape int2str(Nr) 'x' int2str(Ns) '_rsx1.dat']);
rsx2=load(['NETW_periodic_' shape int2str(Nr) 'x' int2str(Ns) '_rsx2.dat']);
bonds=load(['NETW_periodic_' shape int2str(Nr) 'x' int2str(Ns) '_bonds.dat']);
elseif option==3
% SITES REMOVED
rsx0=load(['NETW_sites_removed' int2str(Nr) 'x' int2str(Ns) '_rsx0.dat']);
rsx1=load(['NETW_sites_removed' int2str(Nr) 'x' int2str(Ns) '_rsx1.dat']);
rsx2=load(['NETW_sites_removed' int2str(Nr) 'x' int2str(Ns) '_rsx2.dat']);
bonds=load(['NETW_sites_removed' int2str(Nr) 'x' int2str(Ns) '_bonds.dat']);
end

if isempty(amp)
    amp = ttheta;
    for p = 1:length(ttheta)
        amp(p) = epsilon(p)*EDGE_mech_nls_ampfac(par,ttheta(p),Nr,dispersion);
    end
end
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
    fn1 = ['mech_NETW_nlin_' shape int2str(Nr) 'x' int2str(Ns) '_dispersion' num2str(dispersion) '_eps' epsstr '_ttheta' tthetastr '_amp' ampstr 'U1' int2str(tf) '_tf'];
fn2 = ['mech_NETW_nlin_' shape int2str(Nr) 'x' int2str(Ns) '_dispersion' num2str(dispersion) '_eps' epsstr '_ttheta' tthetastr '_amp' ampstr 'U2' int2str(tf) '_tf'];
fn3 = ['mech_NETW_nlin_' shape int2str(Nr) 'x' int2str(Ns) '_dispersion' num2str(dispersion) '_eps' epsstr '_ttheta' tthetastr '_amp' ampstr 'U3' int2str(tf) '_tf'];
fn4 = ['mech_NETW_nlin_' shape int2str(Nr) 'x' int2str(Ns) '_dispersion' num2str(dispersion) '_eps' epsstr '_ttheta' tthetastr '_amp' ampstr 'U4' int2str(tf) '_tf'];
fn5 = ['mech_NETW_nlin_' shape int2str(Nr) 'x' int2str(Ns) '_dispersion' num2str(dispersion) '_eps' epsstr '_ttheta' tthetastr '_amp' ampstr 'U5' int2str(tf) '_tf'];
fn6 = ['mech_NETW_nlin_' shape int2str(Nr) 'x' int2str(Ns) '_dispersion' num2str(dispersion) '_eps' epsstr '_ttheta' tthetastr '_amp' ampstr 'U6' int2str(tf) '_tf'];
fn7 = ['mech_NETW_nlin_' shape int2str(Nr) 'x' int2str(Ns) '_dispersion' num2str(dispersion) '_eps' epsstr '_ttheta' tthetastr '_amp' ampstr 'U7' int2str(tf) '_tf'];
fn8 = ['mech_NETW_nlin_' shape int2str(Nr) 'x' int2str(Ns) '_dispersion' num2str(dispersion) '_eps' epsstr '_ttheta' tthetastr '_amp' ampstr 'U8' int2str(tf) '_tf'];
load([fn1 '.mat']);
load([fn2 '.mat']);
load([fn3 '.mat']);
load([fn4 '.mat']);
load([fn5 '.mat']);
load([fn6 '.mat']);
load([fn7 '.mat']);
load([fn8 '.mat']);
U = comb_U(U1,U2,U3,U4,U5,U6,U7,U8,8);
elseif option==2
fn1 = ['peregrine_' shape int2str(Nr) 'x' int2str(Ns) '_dispersion' num2str(dispersion) '_eps' epsstr '_ttheta' tthetastr '_amp' ampstr 'U' int2str(tf) '_tf'];
elseif option==3
    fn1 = ['sites_removed_' int2str(Nr) 'x' int2str(Ns) '_dispersion' num2str(dispersion) '_eps' epsstr '_ttheta' tthetastr '_amp' ampstr 'U1' int2str(tf) '_tf'];
fn2 = ['sites_removed_' int2str(Nr) 'x' int2str(Ns) '_dispersion' num2str(dispersion) '_eps' epsstr '_ttheta' tthetastr '_amp' ampstr 'U2' int2str(tf) '_tf'];
fn3 = ['sites_removed_'  int2str(Nr) 'x' int2str(Ns) '_dispersion' num2str(dispersion) '_eps' epsstr '_ttheta' tthetastr '_amp' ampstr 'U3' int2str(tf) '_tf'];
fn4 = ['sites_removed_'  int2str(Nr) 'x' int2str(Ns) '_dispersion' num2str(dispersion) '_eps' epsstr '_ttheta' tthetastr '_amp' ampstr 'U4' int2str(tf) '_tf'];
fn5 = ['sites_removed_' int2str(Nr) 'x' int2str(Ns) '_dispersion' num2str(dispersion) '_eps' epsstr '_ttheta' tthetastr '_amp' ampstr 'U5' int2str(tf) '_tf'];
fn6 = ['sites_removed_' int2str(Nr) 'x' int2str(Ns) '_dispersion' num2str(dispersion) '_eps' epsstr '_ttheta' tthetastr '_amp' ampstr 'U6' int2str(tf) '_tf'];
fn7 = ['sites_removed_' int2str(Nr) 'x' int2str(Ns) '_dispersion' num2str(dispersion) '_eps' epsstr '_ttheta' tthetastr '_amp' ampstr 'U7' int2str(tf) '_tf'];
fn8 = ['sites_removed_'  int2str(Nr) 'x' int2str(Ns) '_dispersion' num2str(dispersion) '_eps' epsstr '_ttheta' tthetastr '_amp' ampstr 'U8' int2str(tf) '_tf'];
load([fn1 '.mat']);
load([fn2 '.mat']);
load([fn3 '.mat']);
load([fn4 '.mat']);
load([fn5 '.mat']);
load([fn6 '.mat']);
load([fn7 '.mat']);
load([fn8 '.mat']);
U = comb_U(U1,U2,U3,U4,U5,U6,U7,U8,8);
end

U=U(:,1:size(U,2)/2);

itv = 10;
frac = 1;
tstart=101;
figure(1)

wobj = VideoWriter([fn1 '.avi']);
wobj.FrameRate = 5;
open(wobj);

for tin = tstart+(0:itv:frac*(size(U,1)-tstart))
textcolor = [1,1,0];

mech_NETW_plot_mnz(rsx0,rsx1,rsx2,U,Nr,Ns,tin,textcolor);
text(2,2,['t = ' num2str(T(tin))], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 20, 'FontName', 'Times', 'Color', 'g');
    
frame = getframe(gca);
writeVideo(wobj,frame);
end
close(wobj);


end
