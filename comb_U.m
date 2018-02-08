function U = comb_U(U1,U2,U3,U4,U5,U6,U7,U8,num)
U=[];
for i=1:num
    if i==1
        V=U1;
    elseif i==2
        V=U2;
    elseif i==3
        V=U3;
    elseif i==4
        V=U4;
     elseif i==5
        V=U5;
     elseif i==6
        V=U6;
     elseif i==7
        V=U7;
     elseif i==8
        V=U8;
    end
    U=[U ; V];
end

