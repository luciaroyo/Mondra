%%Lucia Royo Numerical aproach 4
clear 
close all
%% Polynoliam creation
n = 1:4; % polinomial order
num_const = length(n)+1;%amount of constant
const = randn(num_const, 1); % random coefficient
const = const ./ sqrt(sum(const.^2)); % Normalización

%% numerical method
n_point = 100;
% function
origin=[1,0,0,0]';
for i = 1:n_point
    pos = rand(4,1);
    pos = pos/norm(pos);
    fR3 = polyR3(pos, origin, const, n, true);
    integral_result = Exp(origin, fR3);
    n_(i) = norm(integral_result);
    %origin=integral_result;
end

plot(n_)
% try
%     chol(A); % Si no lanza un error, la matriz es definida positiva
%     isPositiveDefinite = true;
% catch
%     isPositiveDefinite = false;
% end
% isPositiveDefinite


%% funtions:
function f = polyR3(pos, origin, const, n, use_zero_coef)
    if use_zero_coef
        f = const(end) * [0,1,1,1]';
    else
        f = zeros(length(pos),1); % Inicialización
    end

    % Cálculo de los ángulos y logaritmos
    Log_ = Log(origin, pos);
    % Construcción del polinomio
    for i = n
        f = f + const(i)*Log_.^i;
    end
end		

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