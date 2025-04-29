%coefficient SOS
% Paso 1: Inicializar datos
clc;
clear;
%%
load CShape_UQ.mat
data1 = demoUQ{1};
disp(data1);
x1=data1.quat;
x2=x1';
vts=data1.tsVel;
load CShape.mat
t=demos{1}.t;

origin=[1,0,0]';
vq=Exp(vts,1,1); % add a 0 in the firs term
% vq2=zeros(4, 1000);
% vq2(1:3, :)=vq;
x2(:,1)=[];
 
%% 
%Configurar variables de polinomio (usaremos un grado 2 como ejemplo)
sdpvar c1 c2
binvar b1 b2

vqt=cell(1,999);
for i=1:999
    x=Log(origin,x2(i,:)');
    a=c1.*x + c2.*x;
    vqt{i} =a;  % Polinomio de grado 2
end

objective = 0;

% convertir A (3x1000 double) en celda de sdpvar escalares
vqt = reshape(vqt, size(vqt));   % celda 3x1000 de sdpvar, mismos valores que A
for i = 1:999
    bi = vqt{1,i}';  % 3x1 sdpvar
    for j = 1:3
        diff = vts(j,i) - bi(j);  
        objective = objective + diff^2;  % Acumulamos en sdpvar
    end
end
e=0.2;
% constrains=[c1>=e, c2>=e];
Constrains = [];

Constrains = [Constrains, implies(b1 == 0, c1 <= -e)];
Constrains = [Constrains, implies(b1 == 1, c1 >= e)];

Constrains = [Constrains, implies(b2 == 0, c2 <= -e)];
Constrains = [Constrains, implies(b2 == 1, c2 >= e)];
% Paso 6: Resolver con YALMIP + MOSEK (u otro solver SDP)
opts = sdpsettings('solver', 'mosek', 'verbose', 1);  % Configura solver
dignostics=optimize(Constrains, objective, opts);  % Resuelve (sin función objetivo, solo restricciones)

% Paso 7: Extraer resultados
c1_val = value(c1);
c2_val = value(c2);

disp('Constantes del polinomio:');
fprintf('c1 = %.4f\n', c1_val);
fprintf('c2 = %.4f\n', c2_val);
%% functions
function v = Log(origin, quaternion)
    u = origin;
    q = quaternion;

    di = dist(q,u);
    v = proj(u, q);
    % If the two points are "far apart", correct the norm.
    if di > 1e-6
        nv = norm(v);
        v = v * (di / nv);

    end
    % v=u.^(1/2).*log(u.^(-1/2).*q.*u.^(-1/2))*u.^(1/2)
    
end
function b = Exp(pos,w,t)
    if nargin == 2
            % t = 1
        td = w;
    else
        td = t*w;
    end

        nrm_td = norm(td);

        % Former versions of Manopt avoided the computation of sin(a)/a for
        % small a, but further investigations suggest this computation is
        % well-behaved numerically.
    if nrm_td > 0
        b = pos*cos(nrm_td) + td*sin(nrm_td)/nrm_td;
    else
        b = pos;
    end

end
function val= proj(x,d)
    val = d - x*(x(:)'*d(:));
end
function d = dist(x, y) 
    chordal_distance = norm(x - y, 'fro');
    d = real(2*asin(.5*chordal_distance));
end
