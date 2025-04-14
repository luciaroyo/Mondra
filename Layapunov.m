clc; 
clear; 
close all;

%% Cargar datos 
load('datos_trayectoria.mat'); % Asume variables: posicion, velocidad, aceleracion

% Definir variables de decisión para la función de Lyapunov polinómica
sdpvar a0 a1 a2 
x1 = posicion; % Vector de posiciones desde el archivo
v1 = velocidad; % Vector de velocidades desde el archivo
sis = [x1, v1]; % Sistema con datos discretos
n = length(x1); % Número de puntos en la trayectoria

% Definir la función de Lyapunov como una cuadrática general
V = a0+ a1*x1.^2 + a2*x1.^2;

% Derivada de la función de Lyapunov (estimada usando datos)
dVdt = (a1 + 2*a1*x1).*v1;

% Establecer restricciones de positividad y decrecimiento
epsilon = 0.1;
sos_constraint1 = V - epsilon*(x1.^2 + v1.^2); % V(x) > 0
sos_constraint2 = -dVdt - epsilon*(x1.^2 + v1.^2); % -dVdt < 0

% Configurar la optimización usando MOSEK
ops = sdpsettings('solver','mosek','verbose',1);

% Definir restricciones para todos los puntos de datos
% F = [];
% for i = 1:n
%     F = [F, sos(sos_constraint1(i)), sos(sos_constraint2(i))];
% end
F = [F, sos(sos_constraint1), sos(sos_constraint2)];

% Resolver el problema
sol = optimize(F, [], ops);

% Mostrar resultados
if sol.problem == 0
    a0_val = value(a0)
    a1_val = value(a1);
    a2_val = value(a2);
    fprintf('Coeficientes obtenidos: a1 = %.4f, a2 = %.4f, a3 = %.4f, a4 = %.4f, a5 = %.4f, a6 = %.4f\n', ...
        a1_val, a2_val, a3_val, a4_val, a5_val, a6_val);
else
    disp('Error en la optimización.');
    sol.info
end

% Graficar la función de Lyapunov en el espacio de estados
figure;
[X1, X2] = meshgrid(linspace(min(x1),max(x1),50), linspace(min(x2),max(x2),50));
V_plot = a1_val*X1.^2 + a2_val*X2.^2 + a3_val*X1.*X2 + a4_val*X1 + a5_val*X2 + a6_val;
surf(X1, X2, V_plot);
xlabel('Posición');
ylabel('Velocidad');
zlabel('V(x)');
title('Función de Lyapunov basada en datos reales');
