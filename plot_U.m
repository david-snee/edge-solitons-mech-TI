function plot_U(T,U,Nr,Ns,tin)
Upos=U(:,1:2*Nr*Ns);
Uvel=U(:,2*Nr*Ns+1:4*Nr*Ns);

[Px1,Py1] = Left_Edge_sites(Upos,Nr,Ns);
[Vx1,Vy1] = Left_Edge_sites(Uvel,Nr,Ns);


Ap1=sqrt((Px1.^2+Py1.^2));
Av1=sqrt((Vx1.^2+Vy1.^2));

[Px2,Py2] = Top_Edge_sites(Upos,Nr,Ns);
[Vx2,Vy2] = Top_Edge_sites(Uvel,Nr,Ns);


Ap2=sqrt((Px2.^2+Py2.^2));
Av2=sqrt((Vx2.^2+Vy2.^2));

[Px3,Py3] = Right_Edge_sites(Upos,Nr,Ns);
[Vx3,Vy3] = Right_Edge_sites(Uvel,Nr,Ns);


Ap3=sqrt((Px3.^2+Py3.^2));
Av3=sqrt((Vx3.^2+Vy3.^2));

Px3=fliplr(Px3); Py3=fliplr(Py3); Vx3=fliplr(Vx3); Vy3=fliplr(Vy3);
Ap3=fliplr(Ap3); Av3=fliplr(Av3);


[Px4,Py4] = Bottom_Edge_sites(Upos,Nr,Ns);
[Vx4,Vy4] = Bottom_Edge_sites(Uvel,Nr,Ns);


Ap4=sqrt((Px4.^2+Py4.^2));
Av4=sqrt((Vx4.^2+Vy4.^2));

Px4=fliplr(Px4); Py4=fliplr(Py4); Vx4=fliplr(Vx4); Vy4=fliplr(Vy4);
Ap4=fliplr(Ap4); Av4=fliplr(Av4);

Acomb(tin,:)=[Ap1(tin,:) Ap2(tin,:) Ap3(tin,:) Ap4(tin,:)];
X=0:Ns-1+Nr+Ns+Nr;

figure(3)
clf
f5=pcolor(X,T(tin),Acomb(tin,:));
set(f5, 'EdgeColor', 'none');
g1=line([Ns-1 Ns-1],[T(tin(1)) T(tin(end))],[max(max(Acomb)) max(max(Acomb))]);
g2=line([Ns-1+Nr Ns-1+Nr],[T(tin(1)) T(tin(end))],[max(max(Acomb)) max(max(Acomb))]);
g3=line([Ns-1+Nr+Ns Ns-1+Nr+Ns],[T(tin(1)) T(tin(end))],[max(max(Acomb)) max(max(Acomb))]);
set(g1,'color','black','LineWidth',2);
set(g2,'color','black','LineWidth',2);
set(g3,'color','black','LineWidth',2);
set(gcf, 'Units', 'points') 
ylabel('t','FontSize',16)
set(get(gca,'ylabel'),'rotation',0)
set(get(gca,'ylabel'), 'Position', [-50, 700]);
set(gca,'xTick',[])
set(gca,'yTick',[0 1000])
set(gca,'yTickLabel',{'0','1000'},'FontSize',16)
x1 = 90;
y1 = -50;
txt1 = 'L';
text(x1,y1,txt1,'FontWeight','bold','FontSize',16);
x1 = 210;
y1 = -50;
txt1 = 'T';
text(x1,y1,txt1,'FontWeight','bold','FontSize',16);
x1 = 330;
y1 = -50;
txt1 = 'R';
text(x1,y1,txt1,'FontWeight','bold','FontSize',16);
x1 = 450;
y1 = -50;
txt1 = 'B';
text(x1,y1,txt1,'FontWeight','bold','FontSize',14);
colorbar

figure(4)
clf
f5=surf(X,T(tin),Acomb(tin,:));
title('Position of the average A')
set(f5, 'EdgeColor', 'none');

end

