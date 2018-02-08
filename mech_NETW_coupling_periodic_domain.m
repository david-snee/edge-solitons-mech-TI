%% Contains all of the coupling coefficients for the linear ODE. 

function el = mech_NETW_coupling_periodic_domain(r1,r2,s1,s2,par,i1,i2,Ns)
omega0=par(1);
A=par(2);
f=par(3);
coef=-(omega0^2+A*f);

if i1==1 || i1==2
    if (r1-r2==-1 && s1-s2==0) || (r1-r2==1 && s1-s2==0) || (r1-r2==0 && s1-s2==-1) || (r1-r2==0 && s1-s2==1)
      el = f;
    elseif r1-r2==0 && s1-s2==0
      el = coef;
    elseif r1-r2==0 && s1-s2==-(Ns-1)
      el = f;  
    else
      el = 0;
    end  
elseif i1==3
    if (r1-r2==0 && s1-s2==1) || (r1-r2==0 && s1-s2==-1) 
      el = f;
    elseif (r1-r2==-1 && s1-s2==0 && mod(i2,2)==1) || (r1-r2==1 && s1-s2==0 && mod(i2,2)==1)
      el = -f/2; 
    elseif r1-r2==-1 && s1-s2==0 && mod(i2,2)==0
      el = sqrt(3)*f/2;
    elseif r2-r1==-1 && s2-s1==0
      el = -sqrt(3)*f/2;
    elseif r1-r2==0 && s1-s2==0
      el = coef;
    else
      el = 0;
    end
elseif i1==4
    if (r1-r2==0 && s1-s2==1) || (r1-r2==0 && s1-s2==-1) 
      el = f;
    elseif (r1-r2==-1 && s1-s2==0 && mod(i2,2)==0) || (r1-r2==1 && s1-s2==0 && mod(i2,2)==0)
      el = -f/2; 
    elseif r2-r1==-1 && s2-s1==0
      el = sqrt(3)*f/2;
    elseif r1-r2==-1 && s1-s2==0 && mod(i2,2)==1
      el = -sqrt(3)*f/2;
    elseif r1-r2==0 && s1-s2==0
      el = coef;
    else
      el = 0;
    end
elseif i1==5
    if (r1-r2==0 && s1-s2==1) || (r1-r2==0 && s1-s2==-1) 
      el = f;
    elseif (r1-r2==-1 && s1-s2==0 && mod(i2,2)==1) || (r1-r2==1 && s1-s2==0 && mod(i2,2)==1)
      el = -f/2; 
    elseif r2-r1==-1 && s2-s1==0
      el = sqrt(3)*f/2;
    elseif r1-r2==-1 && s1-s2==0 && mod(i2,2)==0
      el = -sqrt(3)*f/2;
    elseif r1-r2==0 && s1-s2==0
      el = coef;
    elseif r1-r2==0 && s1-s2==Ns-1
      el = f;
    else
      el = 0;
    end
elseif i1==6
    if (r1-r2==0 && s1-s2==1) || (r1-r2==0 && s1-s2==-1) 
      el = f;
    elseif (r1-r2==-1 && s1-s2==0 && mod(i2,2)==0) || (r1-r2==1 && s1-s2==0 && mod(i2,2)==0)
      el = -f/2; 
    elseif r1-r2==-1 && s1-s2==0 && mod(i2,2)==1
      el = sqrt(3)*f/2;
    elseif r2-r1==-1 && s2-s1==0
      el = -sqrt(3)*f/2;
    elseif r1-r2==0 && s1-s2==0
      el = coef;
    elseif r1-r2==0 && s1-s2==Ns-1
      el = f;
    else
      el = 0;
    end
end
end
