%% Generates a domain on the lattice sites. Note that removal of a site
%  removes the x and y pendula from that site.


function mechanical_TI_domain(Nr,Ns)
str = ['NETW_rectangle' int2str(Nr) 'x' int2str(Ns)];
%% Create lattice sites for each variable.
function adm = admit_x0y0(r,s)
        adm = mod(s,3)==0;
end
rsx0 = [];
rsy0 = [];
for s = 0:Ns-1
    for r = 0:Nr-1
        if admit_x0y0(r,s)
            rsx0 = [rsx0;r,s];
            rsy0 = [rsy0;r,s];
        end
    end
end
save([str '_rsx0.dat'],'rsx0','-ASCII');
save([str '_rsy0.dat'],'rsy0','-ASCII');

function adm = admit_x1y1(r,s)
        adm = mod(s,3)==1;
end

rsx1 = [];
rsy1 = [];
for s = 0:Ns-1
    for r = 0:Nr-1
        if admit_x1y1(r,s)
            rsx1 = [rsx1;r,s];
            rsy1 = [rsy1;r,s];
        end
    end
end
save([str '_rsx1.dat'],'rsx1','-ASCII');
save([str '_rsy1.dat'],'rsy1','-ASCII');

function adm = admit_x2y2(r,s)
        adm = mod(s,3)==2;
end

rsx2 = [];
rsy2 = [];
for s = 0:Ns-1
    for r = 0:Nr-1
        if admit_x2y2(r,s)
            rsx2 = [rsx2;r,s];
            rsy2 = [rsy2;r,s];
        end
    end
end
save([str '_rsx2.dat'],'rsx2','-ASCII');
save([str '_rsy2.dat'],'rsy2','-ASCII');

bonds=[];
N1 = size(rsx0,1);
N2 = size(rsy0,1);
N3 = size(rsx1,1);
N4 = size(rsy1,1);
N5 = size(rsx2,1);
N6 = size(rsy2,1);

%% Define x0 bonds.
function adm = bond_x0x0(x0L,x0R)
    rx0 = rsx0(x0L,1);
    sx0 = rsx0(x0L,2);
    rx0i = rsx0(x0R,1);
    sx0i = rsx0(x0R,2);
    x1L=0; x2L=0; y0L=0; y1L=0; y2L=0; x2R=0; x1R=0; y0R=0; y1R=0; y2R=0;
    adm = (rx0-rx0i==0 && sx0-sx0i==0) || (rx0-rx0i==-1 && sx0-sx0i==0) || (rx0-rx0i==1 && sx0-sx0i==0);
end

for x0L = 1:N1
    for x0R = 1:N1
        if bond_x0x0(x0L,x0R)
            bonds = [bonds;x0L,y0L,x1L,y1L,x2L,y2L,x0R,y0R,x1R,y1R,x2R,y2R];
        end
    end    
end

function adm = bond_x0x1(x0L,x1R)
    rx0 = rsx0(x0L,1);
    sx0 = rsx0(x0L,2);
    rx1 = rsx1(x1R,1);
    sx1 = rsx1(x1R,2);
    x1L=0; x2L=0; y0L=0; y1L=0; y2L=0; x0R=0; x2R=0; y0R=0; y1R=0; y2R=0;
    adm = (rx0-rx1==0 && sx0-sx1==-1);
end

for x0L = 1:N1
    for x1R = 1:N3
        if bond_x0x1(x0L,x1R)
            bonds = [bonds;x0L,y0L,x1L,y1L,x2L,y2L,x0R,y0R,x1R,y1R,x2R,y2R];
        end
    end    
end

function adm = bond_x0x2(x0L,x2R)
    rx0 = rsx0(x0L,1);
    sx0 = rsx0(x0L,2);
    rx2 = rsx2(x2R,1);
    sx2 = rsx2(x2R,2);
    x1L=0; x2L=0; y0L=0; y1L=0; y2L=0; x0R=0; x1R=0; y0R=0; y1R=0; y2R=0;
    adm = (rx0-rx2==0 && sx0-sx2==1);
end

for x0L = 1:N1
    for x2R = 1:N5
        if bond_x0x2(x0L,x2R)
            bonds = [bonds;x0L,y0L,x1L,y1L,x2L,y2L,x0R,y0R,x1R,y1R,x2R,y2R];
        end
    end    
end
% %% Define y0 bonds.
function adm = bond_y0y0(y0L,y0R)
    ry0 = rsy0(y0L,1);
    sy0 = rsy0(y0L,2);
    ry0i = rsy0(y0R,1);
    sy0i = rsy0(y0R,2);
    x1L=0; x2L=0; x0L=0; y1L=0; y2L=0; x2R=0; x1R=0; x0R=0; y1R=0; y2R=0;
    adm = (ry0-ry0i==0 && sy0-sy0i==0) || (ry0-ry0i==-1 && sy0-sy0i==0) || (ry0-ry0i==1 && sy0-sy0i==0);
end

for y0L = 1:N2
    for y0R = 1:N2
        if bond_y0y0(y0L,y0R)
            bonds = [bonds;x0L,y0L,x1L,y1L,x2L,y2L,x0R,y0R,x1R,y1R,x2R,y2R];
        end
    end    
end

function adm = bond_y0y1(y0L,y1R)
    ry0 = rsy0(y0L,1);
    sy0 = rsy0(y0L,2);
    ry1 = rsy1(y1R,1);
    sy1 = rsy1(y1R,2);
    x1L=0; x2L=0; x0L=0; y1L=0; y2L=0; x0R=0; x2R=0; y0R=0; x1R=0; y2R=0;
    adm = (ry0-ry1==0 && sy0-sy1==-1);
end

for y0L = 1:N2
    for y1R = 1:N4
        if bond_y0y1(y0L,y1R)
            bonds = [bonds;x0L,y0L,x1L,y1L,x2L,y2L,x0R,y0R,x1R,y1R,x2R,y2R];
        end
    end    
end

function adm = bond_y0y2(y0L,y2R)
    ry0 = rsy0(y0L,1);
    sy0 = rsy0(y0L,2);
    ry2 = rsy2(y2R,1);
    sy2 = rsy2(y2R,2);
    x1L=0; x2L=0; x0L=0; y1L=0; y2L=0; x0R=0; x2R=0; y0R=0; x1R=0; y1R=0;
    adm = (ry0-ry2==0 && sy0-sy2==1);
end

for y0L = 1:N2
    for y2R = 1:N6
        if bond_y0y2(y0L,y2R)
            bonds = [bonds;x0L,y0L,x1L,y1L,x2L,y2L,x0R,y0R,x1R,y1R,x2R,y2R];
        end
    end    
end


%% Define x1 bonds.
function adm = bond_x1x0(x1L,x0R)
    rx1 = rsx1(x1L,1);
    sx1 = rsx1(x1L,2);
    rx0 = rsx0(x0R,1);
    sx0 = rsx0(x0R,2);
    x0L=0; x2L=0; y0L=0; y1L=0; y2L=0; x2R=0; x1R=0; y0R=0; y1R=0; y2R=0;
    adm = (rx1-rx0==0 && sx1-sx0==1);
end

for x1L = 1:N3
    for x0R = 1:N1
        if bond_x1x0(x1L,x0R)
            bonds = [bonds;x0L,y0L,x1L,y1L,x2L,y2L,x0R,y0R,x1R,y1R,x2R,y2R];
        end
    end    
end

function adm = bond_x1x1(x1L,x1R)
    rx1 = rsx1(x1L,1);
    sx1 = rsx1(x1L,2);
    rx1i = rsx1(x1R,1);
    sx1i = rsx1(x1R,2);
    x0L=0; x2L=0; y0L=0; y1L=0; y2L=0; x2R=0; x0R=0; y0R=0; y1R=0; y2R=0;
    adm = (rx1-rx1i==0 && sx1-sx1i==0) || (rx1-rx1i==-1 && sx1-sx1i==0) || (rx1-rx1i==1 && sx1-sx1i==0);
end

for x1L = 1:N3
    for x1R = 1:N3
        if bond_x1x1(x1L,x1R)
            bonds = [bonds;x0L,y0L,x1L,y1L,x2L,y2L,x0R,y0R,x1R,y1R,x2R,y2R];
        end
    end    
end

function adm = bond_x1y1(x1L,y1R)
    rx1 = rsx1(x1L,1);
    sx1 = rsx1(x1L,2);
    ry1 = rsy1(y1R,1);
    sy1 = rsy1(y1R,2);
    x0L=0; x2L=0; y0L=0; y1L=0; y2L=0; x2R=0; x1R=0; y0R=0; x2R=0; y2R=0;
    adm = (rx1-ry1==-1 && sx1-sy1==0) || (ry1-rx1==-1 && sy1-sx1==0);
end

for x1L = 1:N3
    for y1R = 1:N4
        if bond_x1y1(x1L,y1R)
            bonds = [bonds;x0L,y0L,x1L,y1L,x2L,y2L,x0R,y0R,x1R,y1R,x2R,y2R];
        end
    end    
end

function adm = bond_x1x2(x1L,x2R)
    rx1 = rsx1(x1L,1);
    sx1 = rsx1(x1L,2);
    rx2 = rsx2(x2R,1);
    sx2 = rsx2(x2R,2);
    x0L=0; x2L=0; y0L=0; y1L=0; y2L=0; x0R=0; x1R=0; y0R=0; y1R=0; y2R=0;
    adm = (rx1-rx2==0 && sx1-sx2==-1);
end

for x1L = 1:N3
    for x2R = 1:N5
        if bond_x1x2(x1L,x2R)
            bonds = [bonds;x0L,y0L,x1L,y1L,x2L,y2L,x0R,y0R,x1R,y1R,x2R,y2R];
        end
    end    
end



% %% Define y1 bonds.
function adm = bond_y1y0(y1L,y0R)
    ry1 = rsy1(y1L,1);
    sy1 = rsy1(y1L,2);
    ry0 = rsy0(y0R,1);
    sy0 = rsy0(y0R,2);
    x0L=0; x2L=0; y0L=0; x1L=0; y2L=0; x2R=0; x1R=0; x0R=0; y1R=0; y2R=0;
    adm = (ry1-ry0==0 && sy1-sy0==1);
end

for y1L = 1:N4
    for y0R = 1:N2
        if bond_y1y0(y1L,y0R)
            bonds = [bonds;x0L,y0L,x1L,y1L,x2L,y2L,x0R,y0R,x1R,y1R,x2R,y2R];
        end
    end    
end

function adm = bond_y1x1(y1L,x1R)
    rx1 = rsx1(x1R,1);
    sx1 = rsx1(x1R,2);
    ry1 = rsy1(y1L,1);
    sy1 = rsy1(y1L,2);
    x0L=0; x2L=0; y0L=0; y1R=0; y2L=0; x2R=0; x1L=0; y0R=0; x2R=0; y2R=0;
    adm = (ry1-rx1==-1 && sy1-sx1==0) || (rx1-ry1==-1 && sx1-sy1==0);
end

for x1R = 1:N3
    for y1L = 1:N4
        if bond_y1x1(y1L,x1R)
            bonds = [bonds;x0L,y0L,x1L,y1L,x2L,y2L,x0R,y0R,x1R,y1R,x2R,y2R];
        end
    end    
end


function adm = bond_y1y1(y1L,y1R)
    ry1 = rsy1(y1L,1);
    sy1 = rsy1(y1L,2);
    ry1i = rsy1(y1R,1);
    sy1i = rsy1(y1R,2);
    x0L=0; x2L=0; y0L=0; x1L=0; y2L=0; x2R=0; x0R=0; y0R=0; x1R=0; y2R=0;
    adm = (ry1-ry1i==0 && sy1-sy1i==0) || (ry1-ry1i==-1 && sy1-sy1i==0) || (ry1-ry1i==1 && sy1-sy1i==0);
end

for y1L = 1:N4
    for y1R = 1:N4
        if bond_y1y1(y1L,y1R)
            bonds = [bonds;x0L,y0L,x1L,y1L,x2L,y2L,x0R,y0R,x1R,y1R,x2R,y2R];
        end
    end    
end

function adm = bond_y1y2(y1L,y2R)
    ry1 = rsy1(y1L,1);
    sy1 = rsy1(y1L,2);
    ry2 = rsy2(y2R,1);
    sy2 = rsy2(y2R,2);
    x0L=0; x2L=0; y0L=0; x1L=0; y2L=0; x0R=0; x1R=0; y0R=0; y1R=0; x2R=0;
    adm = (ry1-ry2==0 && sy1-sy2==-1);
end

for y1L = 1:N4
    for y2R = 1:N6
        if bond_y1y2(y1L,y2R)
            bonds = [bonds;x0L,y0L,x1L,y1L,x2L,y2L,x0R,y0R,x1R,y1R,x2R,y2R];
        end
    end    
end

% %% Define x2 bonds.

function adm = bond_x2x0(x2L,x0R)
    rx2 = rsx2(x2L,1);
    sx2 = rsx2(x2L,2);
    rx0 = rsx0(x0R,1);
    sx0 = rsx0(x0R,2);
    x0L=0; x1L=0; y0L=0; y1L=0; y2L=0; x2R=0; x1R=0; y0R=0; y1R=0; y2R=0;
    adm = (rx2-rx0==0 && sx2-sx0==-1);
end

for x2L = 1:N5
    for x0R = 1:N1
        if bond_x2x0(x2L,x0R)
            bonds = [bonds;x0L,y0L,x1L,y1L,x2L,y2L,x0R,y0R,x1R,y1R,x2R,y2R];
        end
    end    
end

function adm = bond_x2x1(x2L,x1R)
    rx2 = rsx2(x2L,1);
    sx2 = rsx2(x2L,2);
    rx1 = rsx1(x1R,1);
    sx1 = rsx1(x1R,2);
    x0L=0; x1L=0; y0L=0; y1L=0; y2L=0; x2R=0; x0R=0; y0R=0; y1R=0; y2R=0;
    adm = (rx2-rx1==0 && sx2-sx1==1);
end

for x2L = 1:N5
    for x1R = 1:N3
        if bond_x2x1(x2L,x1R)
            bonds = [bonds;x0L,y0L,x1L,y1L,x2L,y2L,x0R,y0R,x1R,y1R,x2R,y2R];
        end
    end    
end

function adm = bond_x2x2(x2L,x2R)
    rx2 = rsx2(x2L,1);
    sx2 = rsx2(x2L,2);
    rx2i = rsx2(x2R,1);
    sx2i = rsx2(x2R,2);
    x0L=0; x1L=0; y0L=0; y1L=0; y2L=0; x0R=0; x1R=0; y0R=0; y1R=0; y2R=0;
    adm = (rx2-rx2i==0 && sx2-sx2i==0) || (rx2-rx2i==-1 && sx2-sx2i==0) || (rx2-rx2i==1 && sx2-sx2i==0);
end

for x2L = 1:N5
    for x2R = 1:N5
        if bond_x2x2(x2L,x2R)
            bonds = [bonds;x0L,y0L,x1L,y1L,x2L,y2L,x0R,y0R,x1R,y1R,x2R,y2R];
        end
    end    
end


function adm = bond_x2y2(x2L,y2R)
    rx2 = rsx2(x2L,1);
    sx2 = rsx2(x2L,2);
    ry2 = rsy2(y2R,1);
    sy2 = rsy2(y2R,2);
    x0L=0; x1L=0; y0L=0; y1L=0; y2L=0; x2R=0; x1R=0; y0R=0; x0R=0; y1R=0;
    adm = (rx2-ry2==-1 && sx2-sy2==0) || (ry2-rx2==-1 && sy2-sx2==0);
end

for x2L = 1:N5
    for y2R = 1:N6
        if bond_x2y2(x2L,y2R)
            bonds = [bonds;x0L,y0L,x1L,y1L,x2L,y2L,x0R,y0R,x1R,y1R,x2R,y2R];
        end
    end    
end

%% Define y2 bonds
function adm = bond_y2y0(y2L,y0R)
    ry2 = rsy2(y2L,1);
    sy2 = rsy2(y2L,2);
    ry0 = rsy0(y0R,1);
    sy0 = rsy0(y0R,2);
    x0L=0; x1L=0; y0L=0; y1L=0; x2L=0; x2R=0; x1R=0; x0R=0; y1R=0; y2R=0;
    adm = (ry2-ry0==0 && sy2-sy0==-1);
end

for y2L = 1:N6
    for y0R = 1:N2
        if bond_y2y0(y2L,y0R)
            bonds = [bonds;x0L,y0L,x1L,y1L,x2L,y2L,x0R,y0R,x1R,y1R,x2R,y2R];
        end
    end    
end

function adm = bond_y2y1(y2L,y1R)
    ry2 = rsy2(y2L,1);
    sy2 = rsy2(y2L,2);
    ry1 = rsy1(y1R,1);
    sy1 = rsy1(y1R,2);
    x0L=0; x1L=0; y0L=0; y1L=0; x2L=0; x2R=0; x0R=0; y0R=0; x1R=0; y2R=0;
    adm = (ry2-ry1==0 && sy2-sy1==1);
end

for y2L = 1:N6
    for y1R = 1:N4
        if bond_y2y1(y2L,y1R)
            bonds = [bonds;x0L,y0L,x1L,y1L,x2L,y2L,x0R,y0R,x1R,y1R,x2R,y2R];
        end
    end    
end

function adm = bond_y2x2(y2L,x2R)
    rx2 = rsx2(x2R,1);
    sx2 = rsx2(x2R,2);
    ry2 = rsy2(y2L,1);
    sy2 = rsy2(y2L,2);
    x0L=0; x1L=0; y0L=0; y1L=0; x2L=0; x0R=0; x1R=0; y0R=0; y2R=0; y1R=0;
    adm = (ry2-rx2==-1 && sy2-sx2==0) || (rx2-ry2==-1 && sx2-sy2==0);
end

for x2R = 1:N5
    for y2L = 1:N6
        if bond_y2x2(y2L,x2R)
            bonds = [bonds;x0L,y0L,x1L,y1L,x2L,y2L,x0R,y0R,x1R,y1R,x2R,y2R];
        end
    end    
end

function adm = bond_y2y2(y2L,y2R)
    ry2 = rsy2(y2L,1);
    sy2 = rsy2(y2L,2);
    ry2i = rsy2(y2R,1);
    sy2i = rsy2(y2R,2);
    x0L=0; x1L=0; y0L=0; y1L=0; x2L=0; x0R=0; x1R=0; y0R=0; y1R=0; x2R=0;
    adm = (ry2-ry2i==0 && sy2-sy2i==0) || (ry2-ry2i==-1 && sy2-sy2i==0) || (ry2-ry2i==1 && sy2-sy2i==0);
end

for y2L = 1:N6
    for y2R = 1:N6
        if bond_y2y2(y2L,y2R)
            bonds = [bonds;x0L,y0L,x1L,y1L,x2L,y2L,x0R,y0R,x1R,y1R,x2R,y2R];
        end
    end    
end

save([str '_bonds.dat'],'bonds','-ASCII');
end