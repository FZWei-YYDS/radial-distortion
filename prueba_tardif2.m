function prueba_tardif2()
%
% prueba del método de Tardif et. al.: Calibration of Cameras with Radially
% Symmetric Distortion, IEEE PAMI 2009
% simulación partiendo de puntos distribuidos alrededor del centro de
% distorsión en circuferencias de diferentes radios (rd)

r = (0.2:0.2:1)';  %radios de los puntos sin distorsión
nrad = length(r);
f = 1;  %longitud focal de los puntos sin distorsión
kc = [-0.2, 0.05];
rd = r.*(1 + kc(1)*r.^2 + kc(2)*r.^4);  %aplica la distorsión radial
fd = f*rd./r;  %longitud focal para los puntos con distorsión

C = [0; 0; 10];  %vértice de los conos C = (X, Y=0, Z)'
theta = (0:10:350)';
npoints = length(theta);
alpha = 20;  %ángulo de inclinación de la cámara con respecto al plano de referencia
R = [cosd(alpha) 0 -sind(alpha); 0 1 0; sind(alpha) 0 cosd(alpha)];
%parámetros de las elipses
a = zeros(nrad,1);  b = zeros(nrad,1);  X0 = zeros(nrad,1);
for i=1:nrad
    %puntos en la imagen
    x = rd(i)*cosd(theta);  y = rd(i)*sind(theta);  z = fd(i)*ones(npoints,1);
    P = R*[x y z]';
    xr = P(1,:)';  yr = P(2,:)';  zr = P(3,:)';
    %coordenadas de los puntos proyectados al plano de referencia Z=0
    X = C(3)*xr./zr;  Y = C(3)*yr./zr;
    [a(i),b(i),X0(i)] = fit_ellipse(X,Y);
    q{i} = [x/f, y/f, ones(npoints,1)];
    Q{i} = [X, Y, ones(npoints,1)];
end
plot_ellipses(a,b,X0,0.01);
Zmin = -5;  Zmax = 15;
plot_hyperbolas(a,b,X0,Zmin,Zmax,0.01);

[fd,mu] = homography_calibration(q,Q);
r_est = f*rd./fd;
figure(3); plot(r,rd,'bo-',r_est,rd,'ro-','linewidth',2); grid on;

%==========================================================
function [a,b,X0] = fit_ellipse(X,Y)
% ajusta los puntos (X,Y)' a una elipse dada por a*(X-X0)^2 + b*Y^2 = 1
% vector de parámetros p = (a, b, X0)'
p0 = [1, 2, 0];
fun_ellipse = @(p) (p(1)*(X-p(3)).^2 + p(2)*Y.^2 - 1);
options = optimset('display','off');
p = lsqnonlin(fun_ellipse,p0,[],[],options);
a = p(1);  b = p(2);  X0 = p(3);

%=============================================================
function [fd,mu] = homography_calibration(q,Q)
n = length(q);
fd = zeros(n,1);  mu = zeros(n,1);
for i=1:n
    q1 = q{i}(:,1);  q2 = q{i}(:,2);  q3 = q{i}(:,3);
    Q1 = Q{i}(:,1);  Q2 = Q{i}(:,2);  Q3 = Q{i}(:,3);
    A = [q2.*Q1, q2.*Q2, q2.*Q3, -q1.*Q1, -q1.*Q2, -q1.*Q3];
    [U,S,V] = svd(A);
    M = reshape(V(:,6),3,2)';
    alpha = acos(M(1,1)/M(2,2));
    R = [cos(alpha) 0 sin(alpha); 0 1 0; -sin(alpha) 0 cos(alpha)];
    t0 = null(M);
    options = optimset('display','off');
    t0 = lsqnonlin(@minfunt0,t0,[],[],options,q1,q2,Q1,Q2,Q3,R);
    S = R * [[1 0; 0 1; 0 0], -t0] * [Q1, Q2, Q3]';
    A = [q3.*S(2,:)', q2.*Q3; q3.*S(1,:)', q1.*Q3];
    b = [q2.*S(3,:)'; q1.*S(3,:)'];
    r = (A'*A)\(A'*b);
    fd(i) = r(1);  mu(i) = r(2);
end

%=========================================================
function res = minfunt0(t0,q1,q2,Q1,Q2,Q3,R)
S = R * [[1 0; 0 1; 0 0], -t0] * [Q1, Q2, Q3]';
res = q1.*S(2,:)' - q2.*S(1,:)';

%=============================================================
function plot_ellipses(a,b,X0,delta)
figure(1); hold on; grid on;
rng(1);
for i=1:length(a)
    Xp = (X0(i)-sqrt(1/a(i)) : delta : X0(i)+sqrt(1/a(i)))';
    Yp = sqrt((1-a(i)*(Xp-X0(i)).^2)/b(i));
    ind = find(imag(Yp)==0);
    Xp = [Xp(ind); Xp(ind(end:-1:1))];
    Yp = [Yp(ind); -Yp(ind(end:-1:1))];
    plot(Xp,Yp,'-','color',rand(1,3));
end
axis('equal');  set(gca,'fontsize',12);
xlabel('X');  ylabel('Y');

%==================================================
function plot_hyperbolas(a,b,X0,Zmin,Zmax,delta)
Z = (Zmin : delta : Zmax)';
figure(2); hold on; grid on;
rng(1);
for i=1:length(a)
    X = sqrt((b(i)-a(i))*(b(i)*Z.^2 + 1)/(a(i)*b(i)));
    plot(X+X0(i),Z,'-',-X+X0(i),Z,'-','color',rand(1,3));
end
set(gca,'fontsize',12);
xlabel('X');  ylabel('Z');