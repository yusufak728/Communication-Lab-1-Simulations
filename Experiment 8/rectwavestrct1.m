function rectwavestrct1(ao,bo,d,H0,f,t)
%矩形波???构 所有?算?位?米 ?入?毫米
%f l0   operator frequency/wavelength
%lg    waveguide wavelength   %lc    TE10 mode cut frequency
%a b   waveguide dimention   %c     Propagation direction 這裡取為波長 
%d     采?精度   %t     t?刻的??构?
a=ao/1000;
b=bo/1000;
lc=2*a;    %TE10截止?率
l0=(3*10^8)/f;
u=4*pi*10^(-7);

if(l0>lc)
        return;
else
        clf;
lg=l0/((1-(l0/lc)^2)^0.5);
c=lg;
B=2*pi/lg;
w=B/(3*10^8);
x=0:a/d:a;
y=0:b/d:b;
z=0:c/d:c;
[x1,y1,z1]=meshgrid(x,y,z);
%mesh(x1,y1,z1);
hx=-B.*a.*H0.*sin(pi./a.*x1).*sin(w*t-B.*z1)./pi;
hz=H0.*cos(pi./a.*x1).*cos(w*t-z1.*B);
hy=zeros(size(y1));

quiver3(z1,x1,y1,hz,hx,hy,'b');
hold on;
x2=x1-0.001;
y2=y1-0.001;
z2=z1-0.001;
ex=zeros(size(x2));
ey=w.*u.*a.*H0.*sin(pi./a.*x2).*sin(w*t-B.*z2)./pi;
ez=zeros(size(z2));
quiver3(z2,x2,y2,ez,ex,ey,'r');
xlabel('??方向');
ylabel('波???a');
zlabel('波?窄?b');
hold off;
end
a=22.86*1e-3;
b=10.16*1e-3;
f=9.84*1e9;
m=1;
n=0;
miu=4*pi*1e-7;
eps=8.854*1e-12;
kc=((m*pi/a)^2+(n*pi/b)^2)^0.5;
w=2*pi*f;
beta=(miu*eps*w^2-kc^2)^0.5;

t=0;
ngrid=20;
x=[0:a/ngrid:a];y=[0:b/2:b];z=[0:0.04/ngrid:0.04];

ex=randn(3,ngrid+1,ngrid+1);
ey=randn(3,ngrid+1,ngrid+1);
ez=randn(3,ngrid+1,ngrid+1);
hx=randn(3,ngrid+1,ngrid+1);
hy=randn(3,ngrid+1,ngrid+1);
hz=randn(3,ngrid+1,ngrid+1);

for q=1:3
    for p=1:ngrid+1
        for r=1:ngrid+1
    ex(q,p,r)=1i*(w*miu/kc^2)*(n*pi/b)*cos((m*pi/a)*x(p))*sin((n*pi/b)*y(q))*exp(1i*((w*t)-(beta*z(r))));
    ey(q,p,r)=1i*(w*miu/kc^2)*(m*pi/a)*cos((m*pi/a)*x(p))*cos((n*pi/b)*y(q))*exp(1i*((w*t)-(beta*z(r))));
    ez(q,p,r)=0;
    hx(q,p,r)=1i*(beta/kc^2)*(m*pi/a)*sin((m*pi/a)*x(p))*cos((n*pi/b)*y(q))*exp(1i*((w*t)-(beta*z(r))));
    hy(q,p,r)=1i*(beta/kc^2)*(m*pi/a)*cos((m*pi/a)*x(p))*sin((n*pi/b)*y(q))*exp(1i*((w*t)-(beta*z(r))));
    hz(q,p,r)=cos((m*pi/a)*x(p))*cos((n*pi/b)*y(q))*exp(1i*((w*t)-(beta*z(r))));
        end
    end
end
[X,Y,Z]=meshgrid(0:a/ngrid:a,0:b/2:b,0:0.04/ngrid:0.04);
quiver3(X,Y,Z,abs(hx),abs(hz),abs(hy));
hold on;
quiver3(X,Y,Z,abs(ex),abs(ez),abs(ey),'r');
