%% Generates the matrix defining the linear RHS of the equations.

function rhs = mech_NETW_matrix(rsx0,rsy0,rsx1,rsy1,rsx2,rsy2,bonds,par)
N1 = size(rsx0,1);
N2 = size(rsy0,1);
N3 = size(rsx1,1);
N4 = size(rsy1,1);
N5 = size(rsx2,1);
N6 = size(rsy2,1);
N = [N1 N2 N3 N4 N5 N6]; TN=sum(N);
Nb = size(bonds,1);
MAT = sparse(TN,TN);
for i = 1:Nb
    M1=[]; M2=[];
    k=find(bonds(i,:));
    i1=k(1);
    i2=k(2);
    b1=bonds(i,i1);
    b2=bonds(i,i2);
    if i1==1
        rs1=rsx0;
        M1=0;
    elseif i1==2
        rs1=rsy0;
        M1=N1;
    elseif i1==3
        rs1=rsx1;
        M1=sum(N(1:2));
    elseif i1==4
        rs1=rsy1;
        M1=sum(N(1:3));
    elseif i1==5
        rs1=rsx2;
        M1=sum(N(1:4));
    elseif i1==6
        rs1=rsy2;
        M1=sum(N(1:5));
    end
    if i2==7
        rs2=rsx0;
        M2=0;
    elseif i2==8
        rs2=rsy0;
        M2=N1;
    elseif i2==9
        rs2=rsx1;
        M2=sum(N(1:2));
    elseif i2==10
        rs2=rsy1;
        M2=sum(N(1:3));
    elseif i2==11
        rs2=rsx2;
        M2=sum(N(1:4));
    elseif i2==12
        rs2=rsy2;
        M2=sum(N(1:5));
    end
    RS1=rs1(b1,:); RS2=rs2(b2,:);
    r1 = RS1(1); s1 = RS1(2); r2 = RS2(1); s2 = RS2(2);
    
    
    MAT(M1+b1,M2+b2) = mech_NETW_coupling(r1,r2,s1,s2,par,i1,i2);
end

rhs = MAT; 
end
