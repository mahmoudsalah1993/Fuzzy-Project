function [xdot] = ATFC( w )
x1 = w(1);
x2 = w(2);
xh(1) = w(3);
xh(2) = w(4);
mh(1) = w(5);
mh(2) = w(6);
hh(1) = w(7);
hh(2) = w(8);
for i = 9:13
    theta11(i-8)=w(i);
end
for i = 14:18
    theta21(i-13)=w(i);
end
for i = 19:23
    theta12(i-18)=w(i);
end
for i = 24:28
    theta22(i-23)=w(i);
end
t = w(29);
%generate k
if t>=0 && t<5
    k=1;
elseif t>=5 && t<15
    k=2;
elseif t>=15 && t<25
    k=1;
elseif t>=25 && t<35
    k=2;
elseif t>=35 && t<45
    k=1;
elseif t>=45 && t<55
    k=2;
elseif t>=55 && t<65
    k=1;
elseif t>=65 && t<75
    k=2;
elseif t>=75 && t<85
    k=1;
elseif t>=85 && t<95
    k=2;
elseif t>=95 && t<105
    k=1;
elseif t>=105 && t<115
    k=2;
elseif t>=115 && t<125
    k=1;
elseif t>=125 && t<135
    k=2;
elseif t>=135 && t<150
    k=1;
end;
for L=1:5
    mu1(L)= exp(-(xh(1) - 3 + L)^2/2);
end
%phi1 , phi2
s=0;
for L=1:5
    s=s+ mu1(L);
end
for L=1:5
    phi1(L)=mu1(L)/s;
end
for L=1:5
    mu2(L)=exp(-(xh(1) - 3 + L)^2/2)*exp(-(xh(2) - 3 + L)^2/2);
end
d= 0;
for n=1:5
    d = d + mu2(n);
end
for L=1:5
    phi2(L) = mu2(L)/d;
end

%y and yr
y = x1;
yr=sin(0.5*t);
yrd=0.5 * cos(0.5 * t);

%z
z(1)= y-yr;
%alphas
if(k==1)
    alpha1= -z(1) - 2*z(1)- theta11 * transpose(phi1) + 0.5 * cos(0.5*t);
elseif(k==2)
    alpha1= -z(1) - 2*z(1)- theta12 * transpose(phi1)  + 0.5 * cos(0.5*t);
end
%z2 equation 17
z(2)= xh(2)-alpha1;

%dalpha
dalpha1dx= -3;
dalpha1dtheta=-transpose(phi1);
dalpha1dyr=3;

syms x ncount;

%dphi
for L=1:5
    dphix = diff(exp(-(x - 3 + L)^2/2)/symsum(exp(-(x - 3 + ncount)^2/2), ncount, 1, 5));
    dphi1(L)= subs(dphix,xh(1));
end

%N = exp(-(xh(1) - 2 )^2/2) + exp(-(xh(1) - 1 )^2/2) + exp(-(xh(1))^2/2) + exp(-(xh(1) + 1 )^2/2) + exp(-(xh(1) + 2 )^2/2);
%ND = -(xh(1) - 2 ) * exp(-(xh(1) - 2 )^2/2) -(xh(1) - 1 ) * exp(-(xh(1) - 1 )^2/2) -(xh(1)) * exp(-(xh(1))^2/2) -(xh(1) + 1 ) * exp(-(xh(1) + 1 )^2/2) -(xh(1) + 2 ) * exp(-(xh(1) + 2 )^2/2);
%for L=1:5
%    D= exp(-(xh(1) - 3 + L)^2/2);
%    DD = -(xh(1) - 3 + L) * D;
%    dphi1(L) = (N * DD - D * ND)/(N^2);
%end
%disp(dphi1);

%dalphadhx
dalpha1dxh= - theta11 *  transpose(dphi1);

%dtheta
if k==1
    dtheta11 = 0.001* phi1 * z(1) - 0.05 * theta11; %switch with k?
    dtheta12 = [0,0,0,0,0];
    dtheta21 = 0.001 * z(2) * phi2 - 0.05 * theta21;
    dtheta22 = [0,0,0,0,0];
elseif k==2
    dtheta11 = [0,0,0,0,0];
    dtheta12 = 0.001* phi1 * z(1) - 0.05 * theta12;
    dtheta21 = [0,0,0,0,0];
    dtheta22 = 0.001 * z(2) * phi2 - 0.05 * theta22;
end

%xhdot(1) for H
if k==1
    fh = theta11* transpose(phi1);
    xhdot(1)= xh(2) + fh + 4 * (x1 - xh(1));
elseif k==2
    fh = theta12 * transpose(phi1);
    xhdot(1)= xh(2) + fh + 8 * (x1 - xh(1));
end

%H
if k==1
    H2= 12 * (x1 - xh(1)) + theta21* transpose(phi2) - dalpha1dxh *xhdot(1) - dtheta11 * dalpha1dtheta - dalpha1dyr * yrd - dalpha1dx *(xh(2)+theta11*transpose(phi1));
    disp(H2);
    v=(1/mh(1)) * (-20*z(2) - z(2) - 2 * (dalpha1dx ^2) * z(2)- sign(z(2)) * hh(1) - H2);
    hhdot(1)= 0.1 * abs(z(2)) - 0.5 * hh(1);
    hhdot(2)= 0;
    %mhdot
    mhdot(2)=0;
    if abs(mh(1)) < 4 || (abs(mh(1)) == 4 && mh(1) * z(2) * v <=0 )
        mhdot(1)= 0.01*z(2)-0.5*mh(1);
    elseif abs(mh(1))==4 && mh(1)*z(2)*v>0
        mhdot(1)= 0.01*z(2)-0.5*mh(1)- (0.1*mh(1)*z(2)*v*mh(1))/(abs(mh(1))^2);
    end
    u(1) = mh(1) * v + hh(1);
    u(2)=0;
elseif k==2
    H2= 12 * (x1 - xh(1)) + theta22 *transpose(phi2) - dalpha1dxh * xhdot(1) -  dtheta12 * dalpha1dtheta - dalpha1dyr * yrd - dalpha1dx *(xh(2)+theta12*transpose(phi1));    
    v=(1/mh(2)) * (-20*z(2) - z(2) - 2 * (dalpha1dx) ^2 * z(2)- sign(z(2)) * hh(2) - H2);
    hhdot(1)= 0;
    hhdot(2)= 0.1 * abs(z(2)) - 0.5 * hh(2);
    mhdot(1)=0;
    if abs(mh(2)) < 4 || (abs(mh(2)) == 4 && mh(2) * z(2) * v <=0 )
        mhdot(2)= 0.01*z(2)-0.5*mh(2);
    elseif abs(mh(2))==4 && mh(2)*z(2)*v(2)>0
        mhdot(2)= 0.01*z(2)-0.5*mh(2)- (0.1*mh(2)*z(2)*v*mh(2))/(abs(mh(2))^2);
    end
    u(2) = mh(2) * v + hh(2);
    u(1)=0;
end

%xdot and xhdot2
if k==1
    fh = theta21 * transpose(phi2);
    x1dot = x2 - x1 * exp(-0.5 * x1) + cos(t);
    x2dot = u(1) + x1 * sin(x2^2) + cos(2*t) ;
    xhdot(2)= u(1)+ fh + 12 * (x1 - xh(1));
elseif k==2
    fh = theta22 * transpose(phi2);
    x1dot = x2 - x1^2 + cos(t^2);
    x2dot = u(2) + 0.1 * x1 * x2^2 + sin(t);
    xhdot(2)= u(2)+ fh + 12 * (x1 - xh(1));
end
xdot = [x1dot, x2dot, xhdot, mhdot,hhdot, dtheta11, dtheta21, dtheta12, dtheta22, u, yr, k];
xdot = double(xdot);
end

