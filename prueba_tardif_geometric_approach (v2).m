function prueba_tardif_geometric_approach()
%
% prueba del método de Tardif et. al.: Calibration of Cameras with Radially
% Symmetric Distortion, IEEE PAMI 2009
% simulación partiendo de puntos distribuidos alrededor del centro de
% distorsión en circuferencias de diferentes radios (rd)

% genera puntos en la imagen con distorsión radial
f = 800;  %longitud focal
cc = [800, 800];  %centro de distorsión
r = (0.2:0.2:1)';  %radios normalizados de los puntos sin distorsión
theta = (0:10:350)';  %distribución de los puntos a cada 10 grados 
kc = [-0.2; 0.05];
[u,v] = image_points(theta,r,f,cc,kc);

% agrega ruido a las coordenadas de los puntos
rng(1);
sigma = 5;
u = u + sigma*randn(size(u));
v = v + sigma*randn(size(v));

alpha = 20;  %ángulo de inclinación de la cámara
R = [cosd(alpha) 0 -sind(alpha); 0 1 0; sind(alpha) 0 cosd(alpha)];
C = [3; 0; 5];  %centro de proyección C = (xc, yc=0, zc)'
[X,Y] = plane_points(u,v,R,C,f,cc);  %puntos en el plano de calibración

plot_points(u,v);
plot_points(X,Y);

% encuentra los parámetros de las elipses
% epsilon_d = [ad 0 -ad*kd; 0 bd 0; -ad*kd 0 ad*kd^2-1]  ec.(16)
% [X Y 1] * epsilon_d * [X Y 1]' = 0
% ad*(X-kd)^2 + b*Y^2 = 1
[ad,bd,kd] = calibration_conics(X,Y);

% calcula el centro de proyección como el punto más cercano a la
% intersección de las hipérbolas
[xc1,zc1] = closest_point(ad,bd,kd);
plot_viewpoint_conics(ad,bd,kd);
plot(xc1,zc1,'r.','markersize',30);

% calcula el centro de proyección C = (xc,0,zc)'
% alpha_est - ángulo entre el eje óptico y el plano de calibración
% theta_est - semi-ángulo de apertura del cono
[xc,zc,alpha_est,theta_est] = camera_parameters(ad,bd,kd,xc1,zc1);
[a,b,k] = ellipse_parameters(xc,zc,alpha_est,theta_est);
plot_viewpoint_conics(a,b,k);
plot(xc1,zc1,'r.',xc,zc,'b.','markersize',30);

rd = tand(theta_est);  %radios de los puntos normalizados con distorsión
figure(); plot(r,rd,'bo'); grid on;
set(gca,'fontsize',12);
xlabel('$r$','interpreter','latex','fontsize',18);
ylabel('$r_d$','interpreter','latex','fontsize',18);


%=====================================================
function [u,v] = image_points(theta,r,f,cc,kc)
% puntos en la imagen distribuidos en círculos de radios rd con centro cc y
% en ángulos theta
nr = length(r);  np = length(theta);
u = zeros(np,nr);  v = zeros(np,nr);
for i=1:nr
    % aplica distorsión radial
    rd = r(i).*(1 + kc(1)*r(i).^2 + kc(2)*r(i).^4);
    x = rd*cosd(theta);  y = rd*sind(theta);
    u(:,i) = f*x + cc(1);  v(:,i) = f*y + cc(2);
end

%======================================================
function [X,Y] = plane_points(u,v,R,C,f,cc)
% puntos proyectados de la imagen al plano de calibración
% centro de proyección C = [Xc, Yc=0, Zc]
[np,nr] = size(u);
X = zeros(np,nr);  Y = zeros(np,nr);
for i=1:nr
    x = (u(:,i)-cc(1))/f;  y = (v(:,i)-cc(2))/f;
    P = R*[x, y, ones(np,1)]';
    xr = P(1,:)';  yr = P(2,:)';  zr = P(3,:)';
    X(:,i) = C(3)*xr./zr + C(1);
    Y(:,i) = C(3)*yr./zr + C(2);
end

%===========================================================
function [ad,bd,kd] = calibration_conics(X,Y)
nr = size(X,2);
ad = zeros(nr,1);  bd = zeros(nr,1);  kd = zeros(nr,1);
options = optimset('display','off');
for i=1:nr
    p0 = [1 1 0];
    %ajusta una elipse a los puntos (X(:,i), Y(:,i))
    % ad*(X-kd)^2 + bd*Yd^2 - 1 = 0
    fun_ellipse = @(p) (p(1)*(X(:,i)-p(3)).^2 + p(2)*Y(:,i).^2 - 1);
    p = lsqnonlin(fun_ellipse,p0,[],[],options);
    ad(i) = p(1);  bd(i) = p(2);  kd(i) = p(3);
end

%===========================================================
function [xc,zc] = closest_point(ad,bd,kd)
xc = 5;  zc = 10;
p0 = [xc, zc];
options = optimset('display','off');
p = fminsearch(@fun_closest_point,p0,options,ad,bd,kd);
xc = p(1);  zc = p(2);

%==================================================
function res = fun_closest_point(p,ad,bd,kd)
xc = p(1);  zc = p(2);
n = length(ad);
res = 0;
xp = (xc-3:0.001:xc+3)';
for i=1:n
    Psi = [ad(i)*bd(i), 0, -ad(i)*bd(i)*kd(i); 0, bd(i)*(ad(i)-bd(i)), 0; ...
        -ad(i)*bd(i)*kd(i), 0, ad(i)*bd(i)*kd(i)^2 + ad(i) - bd(i)];
    zp = sqrt((-Psi(3,3) - Psi(1,1)*xp.^2 - (Psi(1,3)+Psi(3,1))*xp) / Psi(2,2));
    ind = (imag(zp)==0);
    dist2 = min((xc-xp(ind)).^2 + (zc-zp(ind)).^2);
    res = res + dist2;
end

%===================================================================
function [xc,zc,alpha,theta] = camera_parameters(a,b,k,xc1,zc1)
n = length(a);
alpha = 80;  theta = 10*ones(n,1);
p0 = [xc1; zc1; alpha; theta];
pmin = [0; 0; 0; zeros(n,1)];
pmax = [inf; inf; 90; 90*ones(n,1)];
options = optimset('display','off');
p = lsqnonlin(@fun_parameters,p0,pmin,pmax,options,a,b,k);
xc = p(1);  zc = p(2);
alpha = p(3);
theta = p(4:end);

%===================================================
function res = fun_parameters(p,a,b,k)
xc = p(1);  zc = p(2);
alpha = p(3);
theta = p(4:end);
[ap,bp,kp] = ellipse_parameters(xc,zc,alpha,theta);
res = [(ap-a)./a; (bp-b)./b; (kp-k)./k];

%===========================================================
function [a,b,k] = ellipse_parameters(xc,zc,alpha,theta)
x0 = xc - zc/tand(alpha);
y0 = 0;
ym = sqrt(zc^2 + (xc-x0)^2) * tand(theta) + y0;
lambda = sind(theta).^2 / sind(alpha)^2;
am = 1 - lambda;
bm = 2*lambda*xc - 2*x0;
cm = x0^2 - lambda*xc^2 - lambda*zc^2;
rc = sqrt(bm.^2 - 4*am.*cm);
xm1 = (-bm - rc)./(2*am);
xm2 = (-bm + rc)./(2*am);
k = 0.5*(xm1 + xm2);
a = 1./(xm1 - k).^2;
b = (1-a.*(x0-k).^2)./ym.^2;

%================================================
function plot_points(x,y)
figure(); hold on; grid on; axis('equal');
rng(1);
for i=1:size(x,2)
    plot(x(:,i),y(:,i),'o','color',rand(1,3));
end

%================================================
function plot_viewpoint_conics(ad,bd,kd)
nr = length(ad);
Z = (-10:0.1:10);
figure(); hold on; grid on;
rng(1);
for i=1:nr
    Psi = [ad(i)*bd(i), 0, -ad(i)*bd(i)*kd(i); 0, bd(i)*(ad(i)-bd(i)), 0; ...
        -ad(i)*bd(i)*kd(i), 0, ad(i)*bd(i)*kd(i)^2 + ad(i) - bd(i)];
    a = Psi(1,1);
    b = Psi(1,3) + Psi(3,1);
    c = Psi(2,2)*Z.^2 + Psi(3,3);
    X1 = (-b + sqrt(b^2 - 4*a*c))/(2*a);
    X2 = (-b - sqrt(b^2 - 4*a*c))/(2*a);
    ind = (b^2 - 4*a*c)>=0;
    plot(X1(ind),Z(ind),'-',X2(ind),Z(ind),'-','color',rand(1,3));
end