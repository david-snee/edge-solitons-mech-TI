%% Plots a single frame of the problem.

function mech_NETW_plot_mnz(rsx0,rsx1,rsx2,U,Nr,Ns,tin,textcolor)

L=Nr;
clb=0.1;

% the horizontal layout
x1 = 0;
x2 = 600;
x3 = 0;
Lx = x1+x2+x3;

% the vertical layout
y1 = 0;
y2 = 700;
y3 = 0;
Ly = y1+y2+y3;

% NOTE: gcf=grab current figure. gca=grab current axes.

figure(1)
clf
set(gcf, 'Units', 'points')  % set units to be points.
set(gcf, 'Position', [23 48 Lx Ly])  % specify position which window will pop up.
set(gcf, 'PaperPositionMode', 'auto') % auto sizes the print/save to be size of plot. 
set(0, 'DefaultTextFontAngle', 'italic', 'DefaultTextFontName', 'Times', 'DefaultTextFontSize', 16)  % text related.

axes
axis([0 1 0 1])  %defines the axis to be 1x1 plot.
hold on
box off  %removes the box defined on the axes. (top and right outlines)
axis off  %removes the axis if off
set(gca, 'Position', [0, 0, 1, 1])  %defines the plot to be the whole window.
rectangle('Position', [0.001 0.001 0.998 0.998], 'FaceColor', [0 0 0], 'EdgeColor', [.99 .99 .99])

axes;
eps = 2;
axis([0-eps Nr+eps 0-eps Ns+eps])
hold on
box off
axis off
set(gca, 'FontSize', 12, 'FontName', 'Times');
set(gca, 'Position', [(x1)/Lx, (y1)/Ly, (x2)/Lx, (y2)/Ly])

N1 = size(rsx0,1);
N3 = size(rsx1,1);
N5 = size(rsx2,1);
N=[N1 N1 N3 N3 N5 N5]; TN=sum(N);

ux0=U(tin,1:N1); uy0=U(tin,N1+1:2*N1);
U0=sqrt((ux0.^2+uy0.^2));
for i1 = 1:N1
    rr = rsx0(i1,1);
    ss = rsx0(i1,2);
    A = [rr ss];
    col = U0(i1);
    plot(A(1),A(2),'o','MarkerEdgeColor',min(col/clb,1)*textcolor,'MarkerFaceColor',min(col/clb,1)*textcolor,'MarkerSize',3);
end
ux1=U(tin,sum(N(1:2))+1:sum(N(1:3))); uy1=U(tin,sum(N(1:3))+1:sum(N(1:4)));
U1=sqrt((ux1.^2+uy1.^2));
for i3 = 1:N3
    rr = rsx1(i3,1);
    ss = rsx1(i3,2);
    B = [rr ss];
    col = U1(i3);
    plot(B(1),B(2),'o','MarkerEdgeColor',min(col/clb,1)*textcolor,'MarkerFaceColor',min(col/clb,1)*textcolor,'MarkerSize',3);
end
ux2=U(tin,sum(N(1:4))+1:sum(N(1:5))); uy2=U(tin,sum(N(1:5))+1:TN);
U2=sqrt((ux2.^2+uy2.^2));
for i5 = 1:N5
    rr = rsx2(i5,1);
    ss = rsx2(i5,2);
    C = [rr ss];
    col = U2(i5);
    plot(C(1),C(2),'o','MarkerEdgeColor',min(col/clb,1)*textcolor,'MarkerFaceColor',min(col/clb,1)*textcolor,'MarkerSize',3);
end
end
