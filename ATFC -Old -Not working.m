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
k=1;
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
disp(k);
for L=1:5
    mu11(L)= exp(-(xh(1) - 3 + L)^2/2);
    mu12(L)= exp(-(xh(1) - 3 + L)^2/2);
end
%phi1 , phi2
s1=0;
s2=0;
for L=1:5
    s1=s1+ mu11(L);
    s2=s2+ mu12(L);
end
for L=1:5
    phi11(L)=mu11(L)/s1;
    phi12(L)=mu12(L)/s2;
end
for L=1:5
    mu21(L)=exp(-(xh(1) - 3 + L)^2/2)*exp(-(xh(2) - 3 + L)^2/2);
    mu22(L)=exp(-(xh(1) - 3 + L)^2/2)*exp(-(xh(2) - 3 + L)^2/2);
end
d1= 0;
d2= 0;
for n=1:5
    d1 = d1+ mu21(n);
    d2 = d2 + mu22(n);
end
for L=1:5
    phi21(L) = mu21(L)/d1;
    phi22(L) = mu22(L)/d2;
end

%y and yr
y = x1;
yr=sin(0.5*t);
yrd=0.5 * cos(0.5 * t);
yrddot=-0.25 * sin(0.5 * t);

%z
z(1)= y-yr;
%alphas
alpha1(1)= -z(1) - 2*z(1)- theta11 * transpose(phi11) + 0.5 * cos(0.5*t);
alpha1(2)= -z(1) - 2*z(1)- theta12 * transpose(phi12)  + 0.5 * cos(0.5*t);

%z2 equation 17
z(2)= xh(2)-alpha1(k);

%dalpha
dalpha1dx(1)= -3;
dalpha1dtheta(1)=phi11(1);
dalpha1dyr(1)=3;

dalpha1dx(2)= -3;
dalpha1dtheta(2)=phi12(1);
dalpha1dyr(2)=3;

syms x ncount;

%dphi
for L=1:5
    dphix = diff(exp(-(x - 3 + L)^2/2)/symsum(exp(-(x - 3 + ncount)^2/2), ncount, 1, 5));
    dphi1(L)= double(subs(dphix,xh(1)));
end
for L=1:5
    dphix = diff((exp(-(x - 3 + L)^2/2)*exp(-(xh(2) - 3 + L)^2/2))/symsum((exp(-(x - 3 + ncount)^2/2)*exp(-(xh(2) - 3 + ncount)^2/2)),ncount,1,5));
    dphi2(L)= double(subs(dphix, xh(1)));
end

%dalphadhx
dalpha1dxh(1)= - theta11 *  transpose(dphi1);
dalpha1dxh(2)= -theta12 * transpose(dphi2);

%dtheta
dtheta11 = 0.001* phi11 * z(1) - 0.05 * theta11;
dtheta12 = 0.001* phi12 * z(1) - 0.05 * theta12;
dtheta21 = 0.001 * z(2) * phi21 - 0.05 * theta21;
dtheta22 = 0.001 * z(2) * phi22 - 0.05 * theta22;

%xhdot(1) for H
fh11 = theta11* transpose(phi11);
fh12 = theta12 * transpose(phi12);

if k==1
    xhdot(1)= xh(2) + fh11 + 4;
elseif k==2
    xhdot(1)= xh(2) + fh12 + 8;
end

%H
H2(1)= 12 + theta21* transpose(phi21) - dalpha1dxh(1)*xhdot(1) - dalpha1dtheta(1)*dtheta11(1) - dalpha1dyr(1)*yrd - dalpha1dxh(1)*(xh(2)+theta11*transpose(phi11));
H2(2)= 12 + theta22 *transpose(phi22) - dalpha1dxh(2)*xhdot(1) - dalpha1dtheta(2)*dtheta12(2) - dalpha1dyr(2)*yrd - dalpha1dxh(2)*(xh(2)+theta12*transpose(phi12));

%v
v(1)=(1/mh(1)) * (-20*z(2) - z(2) - 2 * (dalpha1dx(1)) ^2 * z(2)- sign(z(2)) * hh(1) - H2(1));
v(2)=(1/mh(2)) * (-20*z(2) - z(2) - 2 * (dalpha1dx(2)) ^2 * z(2)- sign(z(2)) * hh(2) - H2(2));

%hhdot
hhdot(1)= 0.1 * abs(z(2)) - 0.5 * hh(1);
hhdot(2)= 0.1 * abs(z(2)) - 0.5 * hh(2);

mhdot(1)=mh(1);
mhdot(2)=mh(2);

%mhdot
if abs(mh(1)) < 4 || (abs(mh(1)) == 4 && mh(1) * z(2) * v(1) <=0 )
    mhdot(1)= 0.01*z(2)-0.5*mh(1);
elseif abs(mh(1))==4 && mh(1)*z(2)*v(1)>0
    mhdot(1)= 0.01*z(2)-0.5*mh(1)- (0.1*mh(1)*z(2)*mh(1))/(abs(mh(1))^2);
end

if abs(mh(2)) < 4 || (abs(mh(2)) == 4 && mh(2) * z(2) * v(2) <=0 )
    mhdot(2)= 0.01*z(2)-0.5*mh(2);
elseif abs(mh(2))==4 && mh(2)*z(2)*v(2)>0
    mhdot(2)= 0.01*z(2)-0.5*mh(2)- (0.1*mh(2)*z(2)*mh(2))/(abs(mh(2))^2);
end

%u
u(1) = mh(1) * v(1) + hh(1);
u(2) = mh(2) * v(2) + hh(2);

%fh
fh11 = theta11 * transpose(phi11);
fh21 = theta21 * transpose(phi21);
fh12 = theta12 * transpose(phi12);
fh22 = theta22* transpose(phi22);

%xdot and xhdot
if k==1
    x1dot = x2 - x1 * exp(-0.5 * x1) + cos(t);
    x2dot = u(1) + x1 * sin(x2^2) + cos(2*t) ;
    xhdot(1)= xh(2) + fh11 + 4;
    xhdot(2)= u(1)+ fh21 + 12;
elseif k==2
    x1dot = x2 - x1^2 + cos(t^2);
    x2dot = u(2) + 0.1 * x1 * x2^2 + sin(t);
    xhdot(1)= xh(2) + fh12 + 8;
    xhdot(2)= u(2)+ fh22 + 12 ;
end
xdot = [x1dot, x2dot, xhdot, mhdot,hhdot, dtheta11, dtheta21, dtheta12, dtheta22, u, yr];
xdot= double(xdot);
disp(xdot);
end

