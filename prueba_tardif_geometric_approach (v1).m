function prueba_tardif_geometric_approach()
%
% prueba del m�todo de Tardif et. al.: Calibration of Cameras with Radially
% Symmetric Distortion, IEEE PAMI 2009
% simulaci�n partiendo de puntos distribuidos alrededor del centro de
% distorsi�n en circuferencias de diferentes radios (rd)

% genera puntos en la imagen con distorsi�n radial
f = 800;  %longitud focal
cc = [800, 800];  %centro de distorsi�n
r = (0.2:0.2:1)';  %radios normalizados de los puntos sin distorsi�n
theta = (0:10:350)';  %distribuci�n de los puntos a cada 10 grados
kc = [-0.2; 0.05];
[u,v] = image_points(theta,r,f,cc,kc);
% agrega ruido a las coordenadas de los puntos
rng(1);
sigma = 1;
u = u + sigma*randn(size(u));
v = v + sigma*randn(size(v));

alpha = 20;  %�ngulo de inclinaci�n de la c�mara
R = [cosd(alpha) 0 -sind(alpha); 0 1 0; sind(alpha) 0 cosd(alpha)];
C = [3; 0; 5];
[X,Y] = plane_points(u,v,R,C,f,cc);  %puntos en el plano de calibraci�n

plot_points(u,v);
plot_points(X,Y);

% encuentra los par�metros de las elipses
% epsilon_d = [ad 0 -ad*kd; 0 bd 0; -ad*kd 0 ad*kd^2-1]  ec.(16)
% [X Y 1] * epsilon_d * [X Y 1]' = 0
% ad*(X-kd)^2 + b*Y^2 = 1
[ad,bd,kd] = calibration_conics(X,Y);

% viewpoint conics
% Psi = [ad*bd 0 -ad*bd*kd; 0 bd*(ad-bd) 0; -ad*bd*kd 0 ad*bd*kd^2+ad-bd]  ec (17)
plot_viewpoint_conics(ad,bd,kd);

% par�metros de las elipses  ec. (18)
% epsilon_d_hat = [gamma^2*rho 0 epsilon*rho*gamma+gamma; 0 delta^2*(2*rho+1) 0; ...
%                  epsilon*rho*gamma+gamma 0 rho*epsilon^2+2*epsilon-2]
% *solo rho var�a con el radio de los puntos en la imagen
[gamma,delta,epsilon,rho] = calibration_conics2(X,Y);
plot_viewpoint_conics2(gamma,epsilon,delta,rho);
% *hay un error en los par�metros calculados


%============================================================
function [u,v] = image_points(theta,r,f,cc,kc)
% puntos en la imagen distribuidos en c�rculos de radios rd con centro cc y
% en �ngulos theta
nr = length(r);  np = length(theta);
u = zeros(np,nr);  v = zeros(np,nr);
for i=1:nr
    % aplica distorsi�n radial
    rd = r(i).*(1 + kc(1)*r(i).^2 + kc(2)*r(i).^4);
    x = rd*cosd(theta);  y = rd*sind(theta);
    u(:,i) = f*x + cc(1);  v(:,i) = f*y + cc(2);
end

%======================================================
function [X,Y] = plane_points(u,v,R,C,f,cc)
% puntos proyectados de la imagen al plano de calibraci�n
% centro de proyecci�n C = [Xc, Yc=0, Zc]
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

%====================================================================
function [gamma,delta,epsilon,rho] = calibration_conics2(X,Y)
% ajusta elipses a los puntos (X,Y) usando los par�metros de la ec. (18)
nr = size(X,2);
options = optimset('display','off');
p0 = [1; 1; 1; ones(nr,1)];
p = lsqnonlin(@fun_conics,p0,[],[],options,X,Y);
gamma = p(1);  delta = p(2);  epsilon = p(3);
rho = p(4:end,1);

%=========================================================
function res = fun_conics(p,X,Y)
nr = size(X,2);
gamma = p(1);  delta = p(2);  epsilon = p(3);
rho = p(4:end,1);
res = [];
for i=1:nr
    % ec. (18) para epsilon_d_hat
    e = [gamma^2*rho(i), 0, gamma*(epsilon*rho(i)+1); ...
        0, delta^2*(2*rho(i)+1), 0; ...
        gamma*(epsilon*rho(i)+1), 0, rho(i)*epsilon^2+2*epsilon-2];
    res = [res; e(1,1)*X(:,i).^2 + e(2,2)*Y(:,i).^2 + (e(1,3)*e(3,1))*X(:,i) + e(3,3)];
end

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

%================================================================
function plot_viewpoint_conics2(gamma,epsilon,delta,rho)
nr = length(rho);
Z = (-10:0.1:10);
figure(); hold on; grid on;
rng(1);
for i=1:nr
    Psi = [gamma, 0, epsilon + 1/rho(i); 0, (gamma^2-2*delta^2-delta^2/rho(i))/gamma, 0; ...
        epsilon + 1/rho(i), 0, (gamma^2+2*delta^2*(epsilon-1)+delta^2*epsilon^2*rho(i))/...
        (gamma*delta^2*rho(i))];
    a = Psi(1,1);
    b = Psi(1,3) + Psi(3,1);
    c = Psi(2,2)*Z.^2 + Psi(3,3);
    X1 = (-b + sqrt(b^2 - 4*a*c))/(2*a);
    X2 = (-b - sqrt(b^2 - 4*a*c))/(2*a);
    ind = (b^2 - 4*a*c)>=0;
    plot(X1(ind),Z(ind),'-',X2(ind),Z(ind),'-','color',rand(1,3));
end