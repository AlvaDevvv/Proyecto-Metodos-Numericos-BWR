% =========================================================================
% Archivo: Main_Simulacion.m
% Descripción: Simulación transitoria 1D de la barra de combustible (BWR)
% =========================================================================
clear; clc; close all;

% 1. Cargar Parámetros
Parametros_y_Propiedades;

% 2. Parámetros de Discretización
Nr = 50; % Número total de nodos espaciales
r = linspace(0, R3, Nr)';
dr = r(2) - r(1);

dt = 0.5;        % Paso de tiempo [s]
t_final = 100;   % Tiempo total de simulación [s]
t_vec = 0:dt:t_final;
Nt = length(t_vec);

% Identificación de las regiones (Índices de la malla)
idx_c = r <= R1;                       % Nodos en pastilla
idx_h = (r > R1) & (r <= R2);          % Nodos en huelgo
idx_v = r > R2;                        % Nodos en vaina

% 3. Inicialización
T_ant = T_ini * ones(Nr, 1); % Perfil en el instante t
T_hist = zeros(Nr, Nt);      % Historial de temperaturas
T_hist(:,1) = T_ant;

tol_picard = 1e-4; % Tolerancia para la iteración no lineal

% 4. Bucle Temporal
for n = 2:Nt
    % Variables para el método de Picard
    T_iter = T_ant; 
    error_picard = 1;
    iteraciones = 0;
    
    while error_picard > tol_picard && iteraciones < 50
        iteraciones = iteraciones + 1;
        
        % Pre-asignar matrices del sistema tridiagonal
        a_diag = zeros(Nr, 1); % Subdiagonal
        b_diag = zeros(Nr, 1); % Diagonal principal
        c_diag = zeros(Nr, 1); % Superdiagonal
        d_vec  = zeros(Nr, 1); % Vector independiente
        
        % Evaluar propiedades con la temperatura de la iteración actual
        k_eval = zeros(Nr, 1);
        rhoC_eval = zeros(Nr, 1);
        q_eval = zeros(Nr, 1);
        
        for i = 1:Nr
            if idx_c(i)
                k_eval(i) = k_c(T_iter(i));
                rhoC_eval(i) = rhoC_c(T_iter(i));
                q_eval(i) = q_gen; % Solo hay generación en la pastilla
            elseif idx_h(i)
                k_eval(i) = k_h(T_iter(i));
                rhoC_eval(i) = rhoC_h(T_iter(i));
                q_eval(i) = 0;
            else
                k_eval(i) = k_v(T_iter(i));
                rhoC_eval(i) = rhoC_v(T_iter(i));
                q_eval(i) = 0;
            end
        end
        
        % Ensamblaje del sistema:  -a*T_{i-1} + b*T_i - c*T_{i+1} = d
        for i = 1:Nr
            if i == 1 % Condición de frontera r=0 (Simetría dt/dr = 0)
                % Media armónica para k en la frontera
                k_int = (2 * k_eval(1) * k_eval(2)) / (k_eval(1) + k_eval(2));
                b_diag(1) = rhoC_eval(1)/dt + 4*k_int/(dr^2);
                c_diag(1) = 4*k_int/(dr^2);
                d_vec(1)  = rhoC_eval(1)/dt * T_ant(1) + q_eval(1);
                
            elseif i == Nr % Condición de frontera convectiva r=R3
                % Media armónica para k en la frontera
                k_int = (2 * k_eval(Nr-1) * k_eval(Nr)) / (k_eval(Nr-1) + k_eval(Nr));
                A_cond = k_int * (r(Nr) - dr/2) / (r(Nr) * dr);
                A_conv = h_inf / dr;
                
                a_diag(Nr) = A_cond;
                b_diag(Nr) = rhoC_eval(Nr)/(2*dt) + A_cond + A_conv;
                d_vec(Nr)  = rhoC_eval(Nr)/(2*dt) * T_ant(Nr) + A_conv * T_inf;
                
            else % Nodos internos (Aproximación central con Media Armónica)
                r_mas  = r(i) + dr/2;
                r_menos = r(i) - dr/2;
                
                % ¡CORRECCIÓN CRÍTICA: Media Armónica para conservar flujo!
                k_mas   = (2 * k_eval(i) * k_eval(i+1)) / (k_eval(i) + k_eval(i+1));
                k_menos = (2 * k_eval(i-1) * k_eval(i)) / (k_eval(i-1) + k_eval(i));
                
                coef_E = (k_mas * r_mas) / (r(i) * dr^2);
                coef_W = (k_menos * r_menos) / (r(i) * dr^2);
                
                a_diag(i) = coef_W;
                c_diag(i) = coef_E;
                b_diag(i) = rhoC_eval(i)/dt + coef_W + coef_E;
                d_vec(i)  = rhoC_eval(i)/dt * T_ant(i) + q_eval(i);
            end
        end
        
        % Resolver sistema lineal con TDMA (Se envían las subdiagonales en negativo)
        T_new = Solver_TDMA(-a_diag, b_diag, -c_diag, d_vec);
        
        % Calcular error de Picard
        error_picard = norm(T_new - T_iter, inf);
        T_iter = T_new; % Actualizar para siguiente iteración
    end
    
    % Guardar resultados del paso de tiempo
    T_ant = T_new;
    T_hist(:,n) = T_new;
end

% 5. Gráficos de Resultados
figure;
plot(r*1000, T_hist(:,1), '--k', 'DisplayName', 't = 0 s'); hold on;
plot(r*1000, T_hist(:, floor(Nt/4)), 'b', 'DisplayName', ['t = ', num2str(t_vec(floor(Nt/4))), ' s']);
plot(r*1000, T_hist(:, floor(Nt/2)), 'g', 'DisplayName', ['t = ', num2str(t_vec(floor(Nt/2))), ' s']);
plot(r*1000, T_hist(:, end), 'r', 'LineWidth', 1.5, 'DisplayName', ['t = ', num2str(t_final), ' s']);

% Líneas divisorias de regiones
xline(R1*1000, '--', 'Color', [0.5 0.5 0.5], 'HandleVisibility','off');
xline(R2*1000, '--', 'Color', [0.5 0.5 0.5], 'HandleVisibility','off');

% ¡CORRECCIÓN: Etiquetas dinámicas basadas en la temperatura máxima!
max_T = max(T_hist(:));
text(R1*1000/2, max_T*0.75, 'Pastilla UO_2', 'HorizontalAlignment','center');
text((R1+R2)*1000/2, max_T*0.85, 'Huelgo', 'HorizontalAlignment','center');
text((R2+R3)*1000/2, max_T*0.75, 'Vaina', 'HorizontalAlignment','center');

xlabel('Radio (mm)');
ylabel('Temperatura (K)');
title('Evolución Transitoria de la Temperatura en la Barra BWR');
legend('Location', 'best');
grid on;