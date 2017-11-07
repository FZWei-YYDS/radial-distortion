function prueba_tardif()
%
% prueba del metodo de Tardif et. al.: Calibration of Cameras with Radially
% Symmetric Distortion, IEEE PAMI 2009
% simulacion partiendo de puntos distribuidos alrededor del centro de
% distorsion en circuferencias de diferentes radios (rd)

r = (0.2:0.2:1)';  %radios de los puntos sin distorsion
nrad = length(r);
f = 1;  %longitud focal de los puntos sin distorsion
kc = [-0.2, 0.05];
rd = r.*(1 + kc(1)*r.^2 + kc(2)*r.^4);  %aplica la distorsion radial
fd = f*rd./r;  %longitud focal para los puntos con distorsion

C = [0; 0; 10];  %vertice de los conos C = (X, Y=0, Z)'
theta = (0:10:350)';
npoints = length(theta);
alpha = 30;  %angulo de inclinacion de la camara con respecto al plano de referencia
R = [cosd(alpha) 0 -sind(alpha); 0 1 0; sind(alpha) 0 cosd(alpha)];
%parametros de las elipses
a = zeros(nrad,1);  b = zeros(nrad,1);  X0 = zeros(nrad,1);
rng(1);
figure(1); hold on;  grid on;
figure(2); hold on;  grid on;
for i=1:nrad
    %puntos en la imagen
    x = rd(i)*cosd(theta);  y = rd(i)*sind(theta);  z = fd(i)*ones(npoints,1);
    P = R*[x y z]';
    xr = P(1,:)';  yr = P(2,:)';  zr = P(3,:)';
    %coordenadas de los puntos proyectados al plano de referencia Z=0
    X = C(3)*xr./zr;  Y = C(3)*yr./zr;
    cx = rand(1,3);
    [a(i),b(i),X0(i),Xp,Yp] = fit_ellipse(X,Y,0.01);
    figure(1); plot(X,Y,'o',Xp,Yp,'-','color',cx);
    %coordenadas de los puntos de las hiperbolas (Xc, Yc=0, Zc)' donde se
    %encuentra C
    Zc = (-5:0.1:15)';
    Xc = sqrt((b(i)-a(i))*(b(i)*Zc.^2 + 1)/(a(i)*b(i)));
    figure(2); plot(Xc+X0(i),Zc,'-',-Xc+X0(i),Zc,'linewidth',2,'color',cx);
end
figure(1); axis('equal'); set(gca,'fontsize',12); xlabel('X');  ylabel('Y');
figure(2); set(gca,'fontsize',12); xlabel('X');  ylabel('Z');

% encuentra el vertice de los conos en la interseccion de las hiperbolas
C = find_intersection(a,b,X0);
figure(2); plot(C(1),C(3),'ko','linewidth',2);

% encuentra el radio de los puntos sin distorsion a partir de
% (fd/rd)^2 = -(b * (a*b*(X-X0)^2 + a-b) / (a * (a-b)) (ec.(7) en el articulo)
% y r/f = rd/fd
r = f * sqrt(a.*(b-a)./(b.*(a.*b.*(C(1)-X0).^2 + a - b)));
figure(3); plot(r,rd,'bo-','linewidth',2); grid on;
axis([0 1.2 0 1.2]);
xlabel('r');  ylabel('rd');
set(gca,'fontsize',12);

%==========================================================
function [a,b,X0,Xp,Yp] = fit_ellipse(X,Y,delta)
% ajusta los puntos (X,Y)' a una elipse dada por a*(X-X0)^2 + b*Y^2 = 1
% vector de parametros p = (a, b, X0)'
p0 = [1, 2, 0];
fun_ellipse = @(p) (p(1)*(X-p(3)).^2 + p(2)*Y.^2 - 1);
options = optimset('display','off');
p = lsqnonlin(fun_ellipse,p0,[],[],options);
a = p(1);  b = p(2);  X0 = p(3);
%encuentra puntos (Xp,Yp)' en la elipse
Xp = (min(X):delta:max(X))';
Yp = sqrt((1-a*(Xp-X0).^2)/b);
ind = find(imag(Yp)==0);
Xp = [Xp(ind); Xp(ind(end:-1:1))];
Yp = [Yp(ind); -Yp(ind(end:-1:1))];

%===============================================================
function C = find_intersection(a,b,X0)
% encuentra la interseccion C = (X,Y=0,Z)' de las hiperbolas dadas por
% a*b*(X-X0)^2 + b*(a-b)*Z^2 + (a-b) = 0
p0 = [10; 10];
fun_intersection = @(p) (a.*b.*(p(1)-X0).^2 + b.*(a-b).*p(2)^2 + (a-b));
options = optimset('display','off');
p = lsqnonlin(fun_intersection,p0,[],[],options);
C = [p(1); 0; p(2)];