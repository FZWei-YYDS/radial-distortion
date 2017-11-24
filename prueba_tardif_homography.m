function [r,rd] = prueba_tardif_homography()
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
[u,v] = ideal_points(theta,r,f,cc);

% agrega ruido a las coordenadas de los puntos
rng(1);
sigma = 3;
u = u + sigma*randn(size(u));
v = v + sigma*randn(size(v));
[ud,vd] = distorted_points(u,v,f,kc,cc);

alpha = 20;  %ángulo de inclinación de la cámara
R = [cosd(alpha) 0 sind(alpha); 0 1 0; -sind(alpha) 0 cosd(alpha)];
C = [3; 0; 5];  %centro de proyección C = (xc, yc=0, zc)'
[X,Y] = plane_points(u,v,R,C,f,cc);  %puntos en el plano de calibración

%usa el método de la homografía para calcular la longitud focal aparente
%(fd) y el radio de los puntos con distorsión (rd)
for i=1:size(u,2)
    q{i} = [ud(:,i)-cc(1), vd(:,i)-cc(2), ones(size(u(:,i)))];
    Q{i} = [X(:,i), Y(:,i), ones(size(X(:,i)))];
end
[fd,mu] = homography_calibration(q,Q);

rd = fd.*r/f;
figure(); plot(r,rd,'bo','linewidth',2); grid on;
set(gca,'fontsize',12);
xlabel('$r_u$','interpreter','latex','fontsize',18);
ylabel('$r_d$','interpreter','latex','fontsize',18);


%==============================================================
function [fd,mu] = homography_calibration(q,Q)
A = [];
for i=1:length(Q)
    q1 = q{i}(:,1);  q2 = q{i}(:,2);  q3 = q{i}(:,3);
    Q1 = Q{i}(:,1);  Q2 = Q{i}(:,2);  Q3 = Q{i}(:,3);
    A = [q2.*Q1, q2.*Q2, q2.*Q3, -q1.*Q1, -q1.*Q2, -q1.*Q3];
end
[U,S,V] = svd(A);
M = reshape(V(:,6),3,2)';
alpha = acosd(M(1,1)/M(2,2));
R = [cosd(alpha) 0 -sind(alpha); 0 1 0; sind(alpha) 0 cosd(alpha)];
% calcula t0: el artículo menciona que t0 se obtiene del espacio nulo de M,
% pero no obtuve resultados correctos de esta manera, así que usé el
% término q1*S2 - q2*S1 = 0 donde S = R*T0*Q, mencionados también en el
% artículo y calculé t0 usando mínimos cuadrados para minimizar el error en
% este término
% t0 = null(M);  
A = [];  B = [];
for i=1:length(Q)
    q1 = q{i}(:,1);  q2 = q{i}(:,2);  q3 = q{i}(:,3);
    Q1 = Q{i}(:,1);  Q2 = Q{i}(:,2);  Q3 = Q{i}(:,3);
    A = [A; R(1,1)*q2.*Q3 - R(2,1)*q1.*Q3, R(1,2)*q2.*Q3 - R(2,2)*q1.*Q3, ...
        R(1,3)*q2.*Q3 - R(2,3)*q1.*Q3];
    B = [B; R(1,1)*q2.*Q1 + R(1,2)*q2.*Q2 - R(2,1)*q1.*Q1 - R(2,2)*q1.*Q2];
end
t0 = (A'*A)\(A'*B);
T0 = [[1 0; 0 1; 0 0],-t0];
% calcula fd para cada radio rd y una mu para todos los valores de rd
% usando un único sistema de ecuaciones para una vista modificando la
% ecuación (14c) en la forma en que indica el artículo
n = length(Q);
m = size(Q{1},1);
B = [];
for i=1:n
    q1 = q{i}(:,1);  q2 = q{i}(:,2);  q3 = q{i}(:,3);
    Q1 = Q{i}(:,1);  Q2 = Q{i}(:,2);  Q3 = Q{i}(:,3);
    S = R * T0 * [Q1, Q2, Q3]';
    if( i==1 )
        A = [q3.*S(2,:)', zeros(m,n-1), q2.*Q3; q3.*S(1,:)', zeros(m,n-1), q1.*Q3];
    elseif( i>1 && i<n )
        A = [A; zeros(m,i-1), q3.*S(2,:)', zeros(m,n-i), q2.*Q3; ...
            zeros(m,i-1), q3.*S(1,:)', zeros(m,n-i), q1.*Q3];
    else
        A = [A; zeros(m,n-1), q3.*S(2,:)', q2.*Q3; zeros(m,n-1), q3.*S(1,:)', q1.*Q3];
    end
    B = [B; q2.*S(3,:)'; q1.*S(3,:)'];
end
s = (A'*A)\(A'*B);
fd = s(1:end-1);
mu = s(end);

%=====================================================
function [u,v] = ideal_points(theta,r,f,cc)
%calcula las coordenadas de puntos ideales (sin distorsión) distribuidos en
%circunferencias de radio r con centro en cc
nr = length(r);  np = length(theta);
u = zeros(np,nr);  v = zeros(np,nr);
for i=1:nr
    x = r(i)*cosd(theta);  y = r(i)*sind(theta);
    u(:,i) = f*x + cc(1);  v(:,i) = f*y + cc(2);
end

%======================================================
function [ud,vd] = distorted_points(u,v,f,kc,cc)
%calcula las coordenadas de los puntos con distorsión de acuerdo al modelo
%de Brown usando dos coeficientes de distorsión radial kc
x = (u-cc(1))/f;  y = (v-cc(2))/f;
r = sqrt(x.^2 + y.^2);
xd = x.*(1 + kc(1)*r.^2 + kc(2)*r.^4);
yd = y.*(1 + kc(1)*r.^2 + kc(2)*r.^4);
ud = f*xd + cc(1);  vd = f*yd + cc(2);

%======================================================
function [X,Y] = plane_points(u,v,R,C,f,cc)
%calcula las coordenadas de los puntos en el plano de referencia a partir
%de la proyección de los puntos ideales en la imagen
[np,nr] = size(u);
X = zeros(np,nr);  Y = zeros(np,nr);
for i=1:nr
    x = (u(:,i)-cc(1))/f;  y = (v(:,i)-cc(2))/f;
    P = R*[x, y, ones(np,1)]';
    xr = P(1,:)';  yr = P(2,:)';  zr = P(3,:)';
    X(:,i) = C(3)*xr./zr + C(1);
    Y(:,i) = C(3)*yr./zr + C(2);
end
