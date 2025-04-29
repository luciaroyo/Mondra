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

origin=[1,0,0,0]';
vq=Exp(vts,1,1); % add a 0 in the firs term
vq2=zeros(4, 1000);
vq2(1:3, :)=vq;


%% 
% Paso 2: Configurar variables de polinomio (usaremos un grado 2 como ejemplo)
vqar=cell(1,1000);
for i=1:1000
    x=Log(origin,x2(i,:)');
    vqt = 0.67+0.8042.*x - 0.74.*x.*x;
    vqt1= 0.6561+0.8042.*x + 0.1.*x.*x;
    vqar{i}=vqt;
    vqar1{i}=vqt1;
end

%% Display
% Primer componente
figure;
plot(cellfun(@(x) x(1), vqar),'r');
hold on;
plot(cellfun(@(x) x(1), vqar1),'b');
plot(vq2(1, :));
hold off;
title('Component 1');

% Segundo componente
figure;
plot(cellfun(@(x) x(2), vqar)),'r';
hold on;
plot(cellfun(@(x) x(2), vqar1),'b');
plot(vq2(2, :));
hold off;
title('Component 2');

% Tercer componente
figure;
plot(cellfun(@(x) x(3), vqar),'r');
hold on;
plot(cellfun(@(x) x(3), vqar1),'b');
plot(vq2(3, :));
hold off;
title('Component 3');

% Cuarto componente
figure;
plot(cellfun(@(x) x(4), vqar),'r');
hold on;
plot(cellfun(@(x) x(4), vqar1),'b');
plot(vq2(4, :));
hold off;
title('Component 4');

h1=cellfun(@(x) x(1), vqar);
h2=vq2(1, :);
%% 

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
function sum= sumsquare(v1,v2)
    sum=0;
    for i=1:4
        for j=1:1000
            a=v1{j}(i);
            b=v2(i,j);
            sum=sum+(a-b)^2;
        end
    end
end